// Implementation file for CUDA POA kernels.

#include "cudapoa_kernels.cuh"
#include <stdio.h>

namespace nvidia {

namespace cudapoa {

__global__
void generatePOAKernel(uint8_t* consensus_d,
                 uint8_t* sequences_d,
                 uint32_t max_sequence_size,
                 uint32_t max_depth_per_window,
                 uint32_t total_windows)
{
    uint32_t window_id = blockIdx.x * blockDim.x + threadIdx.x;

    if (window_id >= total_windows)
        return;

    uint32_t input_row_idx = window_id * max_depth_per_window;
    uint8_t *input_row = &sequences_d[input_row_idx * max_sequence_size];
    uint8_t *output_row = &consensus_d[window_id * max_sequence_size];

    int num_nodes = 0;

    // represent the graph as an adjacency list
    __shared__ uint16_t from_edges[MAX_NODES_PER_WINDOW * MAX_NODE_EDGES * 32];
    __shared__ uint16_t num_from_edges[MAX_NODES_PER_WINDOW * 32];
    __shared__ uint16_t to_edges[MAX_NODES_PER_WINDOW * MAX_NODE_EDGES * 32];
    __shared__ uint16_t num_to_edges[MAX_NODES_PER_WINDOW * 32];
    __shared__ uint16_t nodes[MAX_NODES_PER_WINDOW * 32];

    nodes[0] = input_row[0];
    num_nodes++;

    for(uint16_t node_id = 1; node_id < 100; node_id++) //TODO don't hardcode as 100
    {
        nodes[node_id] = input_row[node_id];
        num_nodes++;
        from_edges[node_id * MAX_NODE_EDGES] = node_id -1; // set from edge
        num_from_edges[node_id * MAX_NODE_EDGES]++;
        to_edges[(node_id - 1) * MAX_NODE_EDGES] = node_id; // set from edge
        num_to_edges[(node_id - 1) * MAX_NODE_EDGES]++;
    }

    for(uint32_t c = 0; c < max_sequence_size; c++)
    {
        output_row[c] = input_row[c];
    }
}

void generatePOA(uint8_t* consensus_d,
                 uint8_t* sequences_d,
                 uint32_t max_sequence_size,
                 uint32_t max_depth_per_window,
                 uint32_t total_windows,
                 uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream)
{
    generatePOAKernel<<<num_blocks, num_threads, 0, stream>>>(consensus_d,
                                                              sequences_d,
                                                              max_sequence_size,
                                                              max_depth_per_window,
                                                              total_windows);
}

}

}
