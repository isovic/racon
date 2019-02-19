// Implementation file for CUDA POA kernels.

#include "cudapoa_kernels.cuh"
#include <stdio.h>

namespace nvidia {

namespace cudapoa {

__global__
void generatePOAKernel(uint8_t* consensus_d,
                 uint8_t* sequences_d,
                 uint32_t max_sequence_size,
                 uint8_t * num_sequences_per_window,
                 uint16_t * sequence_lengths,
                 uint32_t max_depth_per_window,
                 uint32_t total_windows)
{
    uint32_t window_id = blockIdx.x * blockDim.x + threadIdx.x;

    if (window_id >= total_windows)
        return;

    uint32_t input_row_idx = window_id * max_depth_per_window;
    uint8_t *input_row = &sequences_d[input_row_idx * max_sequence_size];
    uint8_t *output_row = &consensus_d[window_id * max_sequence_size];

    for(uint32_t c = 0; c < max_sequence_size; c++)
    {
        output_row[c] = input_row[c];
    }
}

void generatePOA(uint8_t* consensus_d,
                 uint8_t* sequences_d,
                 uint32_t max_sequence_size,
                 uint8_t* num_sequences_per_window_d,
                 uint16_t * sequence_lengths_d,
                 uint32_t max_depth_per_window,
                 uint32_t total_windows,
                 uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream)
{
    generatePOAKernel<<<num_blocks, num_threads, 0, stream>>>(consensus_d,
                                                              sequences_d,
                                                              max_sequence_size,
                                                              num_sequences_per_window_d,
                                                              sequence_lengths_d,
                                                              max_depth_per_window,
                                                              total_windows);
}

}

}
