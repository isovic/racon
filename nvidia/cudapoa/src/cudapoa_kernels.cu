// Implementation file for CUDA POA kernels.

#include "cudapoa_kernels.cuh"
#include <stdio.h>

namespace nvidia {

namespace cudapoa {

const uint32_t MAX_DIMENSION = CUDAPOA_MAX_NODES_PER_WINDOW + 1;

// Device function for running topoligical sort on graph.
__device__
void topologicalSortDeviceUtil(uint16_t* sorted_poa,
                               uint16_t node_count,
                               uint16_t* incoming_edge_count,
                               uint16_t* outgoing_edges,
                               uint16_t* outgoing_edge_count)
{
    // Clear the incoming edge count for each node.
    uint16_t local_incoming_edge_count[CUDAPOA_MAX_NODES_PER_WINDOW];
    memset(local_incoming_edge_count, -1, 1000);
    uint16_t sorted_poa_position = 0;

    // Iterate through node IDs (since nodes are from 0
    // through node_count -1, a simple loop works) and fill 
    // out the incoming edge count.
    for(uint16_t n = 0; n < node_count; n++)
    {
        local_incoming_edge_count[n] = incoming_edge_count[n];
        // If we find a node ID has 0 incoming edges, add it to sorted nodes list.
        if (local_incoming_edge_count[n] == 0)
        {
            sorted_poa[sorted_poa_position++] = n;
        }
    }

    // Loop through set of node IDs with no incoming edges,
    // then iterate through their children. For each child decrement their 
    // incoming edge count. If incoming edge count of child == 0, 
    // add its node ID to the sorted order list.
    for(uint16_t n = 0; n < sorted_poa_position; n++)
    {
        uint16_t node = sorted_poa[n];
        for(uint16_t edge = 0; edge < outgoing_edge_count[node]; edge++)
        {
            uint16_t out_node = outgoing_edges[node * CUDAPOA_MAX_NODE_EDGES + edge];
            //printf("%d\n", out_node);
            local_incoming_edge_count[out_node]--;
            if (local_incoming_edge_count[out_node] == 0)
            {
                sorted_poa[sorted_poa_position++] = out_node;
            }
        }
    }

    // sorted_poa will have final ordering of node IDs.
}

// Device function for running Needleman-Wunsch dynamic programming loop.
__device__
void runNeedlemanWunsch(uint8_t* nodes,
                        uint16_t* graph,
                        uint16_t graph_count,
                        uint16_t* incoming_edge_count,
                        uint16_t* incoming_edges,
                        uint16_t* outgoing_edge_count,
                        uint16_t* outgoing_edges,
                        uint8_t* read,
                        uint16_t read_count,
                        int32_t* scores,
                        int16_t* traceback_i,
                        int16_t* traceback_j)
{
    // Assuming gap/mismatch penalty of 1, match rewards of -1.

    // Init boundary conditions.
    for(uint16_t i = 1; i < graph_count + 1; i++)
    {
        scores[i * MAX_DIMENSION + 0] = i;
    }
    for(uint16_t j = 1; j < read_count + 1; j++)
    {
        scores[0 * MAX_DIMENSION + j] = j;
    }

    // Run DP loop.
    for(uint16_t i = 1; i < graph_count + 1; i++)
    {
        for(uint16_t j = 1; j < read_count + 1; j++)
        {
            int32_t hor_val = scores[i * MAX_DIMENSION + (j-1)] + 1;
            int32_t ver_val = scores[(i-1) * MAX_DIMENSION + j] + 1;
            int32_t cell_val = (nodes[graph[i]] == read[j] ? -1 : 1);
            int32_t diag_val = scores[(i-1) * MAX_DIMENSION + (j-1)] + cell_val;

            int32_t final_val = min(hor_val, min(diag_val, ver_val));
            scores[i * MAX_DIMENSION + j] = final_val;
            if (hor_val <= final_val)
            {
                // Insert horizontal.
                traceback_i[i * MAX_DIMENSION + j] = i;
                traceback_j[i * MAX_DIMENSION + j] = j - 1;
            }
            else if (ver_val <= final_val)
            {
                // Insert vertical.
                traceback_i[i * MAX_DIMENSION + j] = i - 1;
                traceback_j[i * MAX_DIMENSION + j] = j;
            }
            else
            {
                // Insert diagonal.
                traceback_i[i * MAX_DIMENSION + j] = i - 1;
                traceback_j[i * MAX_DIMENSION + j] = j - 1;
            }
        }
    }
    
}

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
                       int32_t* scores_d, int16_t* traceback_i_d, int16_t* traceback_j_d,
                       uint8_t* nodes_d,  uint16_t* incoming_edges_d, uint16_t* incoming_edge_count_d,
                       uint16_t* outgoing_edges_d, uint16_t* outgoing_edge_count_d,
                       uint16_t* incoming_edge_w_d, uint16_t* outgoing_edge_w_d,
                       uint16_t* sorted_poa_d)
{
    uint32_t thread_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (thread_idx > total_windows)
        return;

    for(uint32_t window_idx = thread_idx; window_idx < total_windows; window_idx += blockDim.x)
    {
        // Memory layout for graph in adjacency list format.
        uint8_t* nodes = &nodes_d[CUDAPOA_MAX_NODES_PER_WINDOW * thread_idx];
        uint16_t* incoming_edges = &incoming_edges_d[thread_idx * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
        uint16_t* incoming_edge_count = &incoming_edge_count_d[thread_idx * CUDAPOA_MAX_NODES_PER_WINDOW];
        uint16_t* outoing_edges = &outgoing_edges_d[thread_idx * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
        uint16_t* outgoing_edge_count = &outgoing_edge_count_d[thread_idx * CUDAPOA_MAX_NODES_PER_WINDOW];
        uint16_t* incoming_edges_weights = &incoming_edge_w_d[thread_idx * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
        uint16_t* outoing_edges_weights = &outgoing_edge_w_d[thread_idx * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
        uint16_t* sorted_poa = &sorted_poa_d[thread_idx * CUDAPOA_MAX_NODES_PER_WINDOW];

        int32_t* scores = &scores_d[MAX_DIMENSION * MAX_DIMENSION * thread_idx];
        int16_t* traceback_i = &traceback_i_d[MAX_DIMENSION * MAX_DIMENSION * thread_idx];
        int16_t* traceback_j = &traceback_j_d[MAX_DIMENSION * MAX_DIMENSION * thread_idx];

        // Fetch the sequence data and sequence length sub-arrays for specific window ID.
        uint32_t input_row_idx = window_idx * max_depth_per_window;
        uint8_t* window_data = &sequences_d[input_row_idx * sequences_pitch];

        uint8_t num_sequences_in_window = num_sequences_per_window_d[window_idx];
        uint16_t* sequence_length_data = &sequence_lengths_d[window_idx * max_depth_per_window];

        // Create backbone for window based on first sequence in window.
        uint16_t node_count = 0;
        uint16_t sequence_0_length = sequence_length_data[0];
        nodes[0] = window_data[0];
        node_count++;
        sorted_poa[0] = 0;
        //Build the rest of the graphs
        for (int nucleotide_idx=1; nucleotide_idx<sequence_0_length; nucleotide_idx++){
            nodes[nucleotide_idx] = window_data[nucleotide_idx];
            node_count++;
            sorted_poa[nucleotide_idx] = nucleotide_idx;
            outoing_edges[(nucleotide_idx-1) * CUDAPOA_MAX_NODE_EDGES] = nucleotide_idx;
            outgoing_edge_count[nucleotide_idx-1] = 1;
            incoming_edges[nucleotide_idx * CUDAPOA_MAX_NODE_EDGES] = nucleotide_idx - 1;
            incoming_edge_count[nucleotide_idx] = 1;
        }

        //printf("node count %d\n", node_count);

        //// Align each subsequent read, add alignment to graph, run topoligical sort.
        for(uint16_t s = 1; s < num_sequences_in_window; s++)
        {
            uint8_t* seq = &window_data[s * max_sequence_size];
            uint16_t seq_len = sequence_length_data[s];

            // Run DP step and fetch traceback.
            runNeedlemanWunsch(nodes,
                               sorted_poa,
                               node_count,
                               incoming_edge_count,
                               incoming_edges,
                               outgoing_edge_count,
                               outoing_edges,
                               seq,
                               seq_len,
                               scores,
                               traceback_i,
                               traceback_j);

            // Fetch trackback alignment.

            // Add alignment to graph.

            // Run a topsort on the graph. Not strictly necessary at this point
            topologicalSortDeviceUtil(sorted_poa,
                                      node_count,
                                      incoming_edge_count,
                                      outoing_edges, outgoing_edge_count);
        }

        // Dummy kernel code to copy first sequence as output.
        uint8_t *input_row = &sequences_d[input_row_idx * sequences_pitch];
        uint8_t *output_row = &consensus_d[window_idx * consensus_pitch];
        for(uint32_t c = 0; c < sequence_lengths_d[window_idx * max_depth_per_window]; c++)
        {
            output_row[c] = input_row[c];
        }
    }

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
                 int32_t* scores, int16_t* traceback_i, int16_t* traceback_j,
                 uint8_t* nodes,  uint16_t* incoming_edges, uint16_t* incoming_edge_count,
                 uint16_t* outgoing_edges, uint16_t* outgoing_edge_count,
                 uint16_t* incoming_edge_w, uint16_t* outgoing_edge_w,
                 uint16_t* sorted_poa)
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
                                                              sorted_poa);
}


// Kernel for running topological independently.
__global__
void topologicalSortKernel(uint16_t* sorted_poa_d,
                           uint16_t node_count,
                           uint16_t* incoming_edge_count_d,
                           uint16_t* outgoing_edges_d,
                           uint16_t* outgoing_edge_count_d)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
    {
        topologicalSortDeviceUtil(sorted_poa_d,
                                  node_count,
                                  incoming_edge_count_d,
                                  outgoing_edges_d,
                                  outgoing_edge_count_d);
    }
}

// Host function for running topological sort kernel.
void topologicalSort(uint16_t* sorted_poa_d,
                     uint16_t node_count,
                     uint16_t* incoming_edge_count_d,
                     uint16_t* outgoing_edges_d,
                     uint16_t* outgoing_edge_count_d,
                     uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream)
{
    topologicalSortKernel<<<num_blocks, num_threads, 0, stream>>>(sorted_poa_d,
                                                                  node_count,
                                                                  incoming_edge_count_d,
                                                                  outgoing_edges_d,
                                                                  outgoing_edge_count_d);
}

} // namespace cudapoa

} // namespace nvidia
