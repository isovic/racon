// Implementation file for CUDA POA kernels.

#include "cudapoa_kernels.cuh"
#include "cudapoa_nw.cu"
#include "cudapoa_topsort.cu"
#include "cudapoa_add_alignment.cu"

#include <stdio.h>

namespace nvidia {

namespace cudapoa {

// Kernel for running POA.
__global__
void generatePOAKernel(uint8_t* consensus_d,
                       size_t consensus_pitch,
                       uint8_t* sequences_d,
                       size_t sequences_pitch,
                       uint32_t max_sequence_size,
                       uint16_t * num_sequences_per_window_d,
                       uint16_t * sequence_lengths_d,
                       uint32_t max_depth_per_window,
                       uint32_t total_windows,
                       int16_t* scores_d, int16_t* traceback_i_d, int16_t* traceback_j_d,
                       uint8_t* nodes_d,  uint16_t* incoming_edges_d, uint16_t* incoming_edge_count_d,
                       uint16_t* outgoing_edges_d, uint16_t* outgoing_edge_count_d,
                       uint16_t* incoming_edge_w_d, uint16_t* outgoing_edge_w_d,
                       uint16_t* sorted_poa_d, uint16_t* node_id_to_pos_d,
                       uint16_t* node_alignments_d, uint16_t* node_alignment_count_d)
{

    uint32_t block_idx = blockIdx.x;
    uint32_t thread_idx = threadIdx.x;

    long long int back_time = 0;
    long long int nw_time = 0;
    long long int add_time = 0;
    long long int top_time = 0;

    if (block_idx > total_windows)
        return;

    // Loop over all the windows that are assigned to a particular thread.
    for(uint32_t window_idx = block_idx; window_idx < total_windows; window_idx += gridDim.x)
    {
        // Find the buffer offsets for each thread within the global memory buffers.
        uint8_t* nodes = &nodes_d[CUDAPOA_MAX_NODES_PER_WINDOW * block_idx];
        uint16_t* incoming_edges = &incoming_edges_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
        uint16_t* incoming_edge_count = &incoming_edge_count_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW];
        uint16_t* outoing_edges = &outgoing_edges_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
        uint16_t* outgoing_edge_count = &outgoing_edge_count_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW];
        uint16_t* incoming_edge_weights = &incoming_edge_w_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
        uint16_t* outgoing_edge_weights = &outgoing_edge_w_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
        uint16_t* sorted_poa = &sorted_poa_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW];
        uint16_t* node_id_to_pos = &node_id_to_pos_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW];
        uint16_t* node_alignments = &node_alignments_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_ALIGNMENTS];
        uint16_t* node_alignment_count = &node_alignment_count_d[block_idx * CUDAPOA_MAX_NODES_PER_WINDOW];

        int16_t* scores = &scores_d[CUDAPOA_MAX_MATRIX_DIMENSION * CUDAPOA_MAX_MATRIX_DIMENSION * block_idx];
        int16_t* traceback_i = &traceback_i_d[CUDAPOA_MAX_MATRIX_DIMENSION * block_idx];
        int16_t* traceback_j = &traceback_j_d[CUDAPOA_MAX_MATRIX_DIMENSION * block_idx];

        // Fetch the sequence data and sequence length sub-arrays for specific window ID.
        uint32_t input_row_idx = window_idx * max_depth_per_window;
        uint8_t* window_data = &sequences_d[input_row_idx * sequences_pitch];

        uint16_t num_sequences_in_window = num_sequences_per_window_d[window_idx];
        uint16_t* sequence_length_data = &sequence_lengths_d[window_idx * max_depth_per_window];

        long long int t0 = clock64();

        if (thread_idx == 0)
        {

            if (num_sequences_in_window == 0)
            {
                printf("num sequences in window %d\n", num_sequences_in_window);
                continue;
            }

            // Create backbone for window based on first sequence in window.
            uint16_t sequence_0_length = sequence_length_data[0];
            nodes[0] = window_data[0];
            sorted_poa[0] = 0;
            incoming_edge_count[0] = 0;
            node_alignment_count[0] = 0;
            node_id_to_pos[0] = 0;
            outgoing_edge_count[sequence_0_length - 1] = 0;
            //Build the rest of the graphs
            for (int nucleotide_idx=1; nucleotide_idx<sequence_0_length; nucleotide_idx++){
                nodes[nucleotide_idx] = window_data[nucleotide_idx];
                sorted_poa[nucleotide_idx] = nucleotide_idx;
                outoing_edges[(nucleotide_idx-1) * CUDAPOA_MAX_NODE_EDGES] = nucleotide_idx;
                outgoing_edge_count[nucleotide_idx-1] = 1;
                incoming_edges[nucleotide_idx * CUDAPOA_MAX_NODE_EDGES] = nucleotide_idx - 1;
                incoming_edge_count[nucleotide_idx] = 1;
                node_alignment_count[nucleotide_idx] = 0;
                node_id_to_pos[nucleotide_idx] = nucleotide_idx;
            }

        }

        __syncthreads();

        back_time += (clock64() - t0);

        //for(uint16_t i = 0; i < sequence_0_length; i++)
        //{
        //    printf("%c ", nodes[i]);
        //}

        //printf("window id %d, sequences %d\n", window_idx, num_sequences_in_window - 1);

        // Align each subsequent read, add alignment to graph, run topoligical sort.
        for(uint16_t s = 1; s < num_sequences_in_window; s++)
        {
            //printf("running window %d seq %d / %d\n", window_idx, s, num_sequences_in_window);
            uint8_t* seq = &window_data[s * max_sequence_size];
            uint16_t seq_len = sequence_length_data[s];

            //for(uint16_t i = 0; i < seq_len; i++)
            //{
            //    printf("%c ", seq[i]);
            //}

            //return;
            // Run DP step and fetch traceback.
            //bool found_node = false;
            //for(uint16_t i = 0; i < sequence_length_data[0]; i++)
            //{
            //    if (outgoing_edge_count[i] == 0)
            //    {
            //        printf("node %d has 0 oe\n", i);
            //        found_node = true;
            //    }
            //}
            //if (!found_node)
            //{
            //    printf("DID NOT FIND A NODE WITH NO OUTGOING EDGE before alignment!!!!\n");
            //    return;
            //}

            // print sorted graph
            //for(uint16_t i = 0; i < sequence_length_data[0]; i++)
            //{
            //    printf("%d ", sorted_poa[i]);
            //}
            //printf("\n");

            if (thread_idx == 0)
            {

                if (sequence_length_data[0] >= CUDAPOA_MAX_NODES_PER_WINDOW)
                {
                    printf("Node count %d is greater than max matrix size %d\n", sequence_length_data[0], CUDAPOA_MAX_NODES_PER_WINDOW);
                    return;
                }
                if (seq_len >= CUDAPOA_MAX_NODES_PER_WINDOW)
                {
                    printf("Sequence len %d is greater than max matrix size %d\n", seq_len, CUDAPOA_MAX_NODES_PER_WINDOW);
                    return;
                }

            }
            long long int start = clock64();

            // Run Needleman-Wunsch alignment between graph and new sequence.
            //printf("running nw\n");
            uint16_t alignment_length = runNeedlemanWunsch(nodes,
                               sorted_poa,
                               node_id_to_pos,
                               sequence_length_data[0],
                               incoming_edge_count,
                               incoming_edges,
                               outgoing_edge_count,
                               outoing_edges,
                               seq,
                               seq_len,
                               scores,
                               traceback_i,
                               traceback_j);

            long long int nw_end = clock64();
            nw_time += (nw_end - start);
            __syncthreads();

            //found_node = false;
            //for(uint16_t i = 0; i < sequence_length_data[0]; i++)
            //{
            //    if (outgoing_edge_count[i] == 0)
            //    {
            //        printf("node %d has 0 oe\n", i);
            //        found_node = true;
            //    }
            //}
            //if (!found_node)
            //{
            //    printf("DID NOT FIND A NODE WITH NO OUTGOING EDGE before addition!!!!\n");
            //    return;
            //}

            start = clock64();

            if (thread_idx == 0)
            {

                // Add alignment to graph.
                //printf("running add\n");
                sequence_length_data[0] = addAlignmentToGraph(nodes, sequence_length_data[0],
                        node_alignments, node_alignment_count,
                        incoming_edges, incoming_edge_count,
                        outoing_edges, outgoing_edge_count,
                        incoming_edge_weights, outgoing_edge_weights,
                        alignment_length,
                        sorted_poa, traceback_i, 
                        seq, traceback_j);

                long long int add_end = clock64();
                add_time += (add_end - start);

                // Verify that each graph has at least one node with no outgoing edges.
                //bool found_node = false;
                //for(uint16_t i = 0; i < sequence_length_data[0]; i++)
                //{
                //    //printf("node id %d ie %d oe %d\n ", i, incoming_edge_count[i], outgoing_edge_count[i]);
                //    if (outgoing_edge_count[i] == 0)
                //        found_node = true;
                //}
                //if (!found_node)
                //{
                //    printf("DID NOT FIND A NODE WITH NO OUTGOING EDGE after addition!!!!\n");
                //    return;
                //}


                // Run a topsort on the graph. Not strictly necessary at this point
                //printf("running topsort\n");
                topologicalSortDeviceUtil(sorted_poa,
                        node_id_to_pos,
                        sequence_length_data[0],
                        incoming_edge_count,
                        outoing_edges, outgoing_edge_count);

                long long int top_end = clock64();
                top_time += (top_end - add_end);
                //printf("done loop\n");
            }

            __syncthreads();
        }

        // Dummy kernel code to copy first sequence as output.
        //uint8_t *input_row = &sequences_d[input_row_idx * sequences_pitch];
        //uint8_t *output_row = &consensus_d[window_idx * consensus_pitch];
        //for(uint32_t c = 0; c < sequence_lengths_d[window_idx * max_depth_per_window]; c++)
        //{
        //    output_row[c] = input_row[c];
        //}
    }

    //if (thread_idx == 0)
    //{
    //long long int total = back_time + nw_time + add_time + top_time;
    //printf("Total time of backbone generation is %lf %\n", ((double)back_time / total) * 100.f);
    //printf("Total time of nw is %lf %\n", ((double)nw_time / total) * 100.f);
    //printf("Total time of addition is %lf %\n", ((double)add_time / total) * 100.f);
    //printf("Total time of topsort is %lf %\n", ((double)top_time / total) * 100.f);
    //}

}

// Host function call for POA kernel.
void generatePOA(uint8_t* consensus_d,
                 size_t consensus_pitch,
                 uint8_t* sequences_d,
                 size_t sequences_pitch,
                 uint32_t max_sequence_size,
                 uint16_t* num_sequences_per_window_d,
                 uint16_t * sequence_lengths_d,
                 uint32_t max_depth_per_window,
                 uint32_t total_windows,
                 uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream,
                 int16_t* scores, int16_t* traceback_i, int16_t* traceback_j,
                 uint8_t* nodes,  uint16_t* incoming_edges, uint16_t* incoming_edge_count,
                 uint16_t* outgoing_edges, uint16_t* outgoing_edge_count,
                 uint16_t* incoming_edge_w, uint16_t* outgoing_edge_w,
                 uint16_t* sorted_poa, uint16_t* node_id_to_pos,
                 uint16_t* node_alignments, uint16_t* node_alignment_count)
{
    generatePOAKernel<<<num_blocks, num_threads, 0, stream>>>(consensus_d,
                                                              consensus_pitch,
                                                              sequences_d,
                                                              sequences_pitch,
                                                              max_sequence_size,
                                                              num_sequences_per_window_d,
                                                              sequence_lengths_d,
                                                              max_depth_per_window,
                                                              total_windows,
                                                              scores, traceback_i, traceback_j,
                                                              nodes, incoming_edges, incoming_edge_count,
                                                              outgoing_edges, outgoing_edge_count,
                                                              incoming_edge_w, outgoing_edge_w,
                                                              sorted_poa, node_id_to_pos,
                                                              node_alignments, node_alignment_count);
}

} // namespace cudapoa

} // namespace nvidia
