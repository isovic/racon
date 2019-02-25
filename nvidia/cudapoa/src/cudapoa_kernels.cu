// Implementation file for CUDA POA kernels.

#include "cudapoa_kernels.cuh"

namespace nvidia {

namespace cudapoa {

// Device function for running topoligical sort on graph.
__device__
void topologicalSortDeviceUtil(uint8_t* sorted_poa,
                               uint16_t node_count,
                               uint16_t* incoming_edge_count,
                               uint16_t* outgoing_edges,
                               uint16_t* outgoing_edge_count)
{
    // Clear the incoming edge count for each node.
    __shared__ uint16_t local_incoming_edge_count[CUDAPOA_MAX_NODES_PER_WINDOW];
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
        uint8_t node = sorted_poa[n];
        for(uint16_t edge = 0; edge < outgoing_edge_count[node]; edge++)
        {
            uint8_t out_node = outgoing_edges[node * CUDAPOA_MAX_NODE_EDGES + edge];
            local_incoming_edge_count[out_node]--;
            if (local_incoming_edge_count[out_node] == 0)
            {
                sorted_poa[sorted_poa_position++] = out_node;
            }
        }
    }

    // sorted_poa will have final ordering of node IDs.
}

// Kernel for running POA.
__global__
void generatePOAKernel(uint8_t* consensus_d,
                       size_t consensus_pitch,
                       uint8_t* sequences_d,
                       size_t sequences_pitch,
                       uint32_t max_sequence_size,
                       uint8_t * num_sequences_per_window_d,
                       uint16_t * sequence_lengths_d,
                       uint32_t max_depth_per_window,
                       uint32_t total_windows)
{
    // Memory layout for graph in adjacency list format.
    __shared__ uint8_t nodes[CUDAPOA_MAX_NODES_PER_WINDOW];
    __shared__ uint16_t incoming_edges[CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
    __shared__ uint16_t incoming_edge_count[CUDAPOA_MAX_NODES_PER_WINDOW];
    __shared__ uint16_t outoing_edges[CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
    __shared__ uint16_t outgoing_edge_count[CUDAPOA_MAX_NODES_PER_WINDOW];
    __shared__ uint16_t incoming_edges_weights[CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
    __shared__ uint16_t outoing_edges_weights[CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
    __shared__ uint8_t sorted_poa[CUDAPOA_MAX_NODES_PER_WINDOW];

    uint16_t band_idx = blockIdx.x * blockDim.x + threadIdx.x; // The ID of the thread within the band

    for (uint16_t window_idx = 0; window_idx<total_windows; window_idx++){
        uint16_t node_count = 0;
        if (threadIdx.x == 0){
            uint16_t sequence_0_length = sequence_lengths_d[0];
            uint32_t input_row_idx = window_idx * max_depth_per_window;
            nodes[0] = sequences_d[input_row_idx];
            node_count++;
            sorted_poa[0] = 0;
            //Build the rest of the graphs
            for (int nucleotide_idx=1; nucleotide_idx<sequence_0_length; nucleotide_idx++){
                    nodes[nucleotide_idx] = sequences_d[input_row_idx + nucleotide_idx];
                    node_count++;
                    sorted_poa[nucleotide_idx] = nucleotide_idx;
                    outoing_edges[nucleotide_idx-1] = nucleotide_idx;
                    outgoing_edge_count[nucleotide_idx-1] = 1;
                    incoming_edges[nucleotide_idx] = nucleotide_idx - 1;
                    incoming_edge_count[nucleotide_idx] = 1;
            }

            //Run a topsort on the graph. Not strictly necessary at this point
            topologicalSortDeviceUtil(sorted_poa,
                                      node_count,
                                      incoming_edge_count,
                                      outoing_edges, outgoing_edge_count);
        }
    }


    // Dummy kernel code to copy first sequence as output.
    uint32_t window_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (window_id >= total_windows)
        return;
    uint32_t input_row_idx = window_id * max_depth_per_window;
    uint8_t *input_row = &sequences_d[input_row_idx * sequences_pitch];
    uint8_t *output_row = &consensus_d[window_id * consensus_pitch];
    for(uint32_t c = 0; c < max_sequence_size; c++)
    {
        output_row[c] = input_row[c];
    }
}

// Host function call for POA kernel.
void generatePOA(uint8_t* consensus_d,
                 size_t consensus_pitch,
                 uint8_t* sequences_d,
                 size_t sequences_pitch,
                 uint32_t max_sequence_size,
                 uint8_t* num_sequences_per_window_d,
                 uint16_t * sequence_lengths_d,
                 uint32_t max_depth_per_window,
                 uint32_t total_windows,
                 uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream)
{
    generatePOAKernel<<<num_blocks, num_threads, 0, stream>>>(consensus_d,
                                                              consensus_pitch,
                                                              sequences_d,
                                                              sequences_pitch,
                                                              max_sequence_size,
                                                              num_sequences_per_window_d,
                                                              sequence_lengths_d,
                                                              max_depth_per_window,
                                                              total_windows);
}


// Kernel for running topological independently.
__global__
void topologicalSortKernel(uint8_t* sorted_poa_d,
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
void topologicalSort(uint8_t* sorted_poa_d,
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
