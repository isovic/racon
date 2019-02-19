// Header for for CUDA POA host kernel wrappers.

#pragma once

#include <stdint.h>
#include <cuda_runtime_api.h>

#define MAX_NODE_EDGES 5
#define MAX_NODES_PER_WINDOW 1000

namespace nvidia {

namespace cudapoa {

void generatePOA(uint8_t* consensus_d,
                 uint8_t* sequences_d,
                 uint32_t max_sequence_size,
                 uint8_t * num_sequences_per_window_d,
                 uint16_t * sequence_lengths_d,
                 uint32_t max_depth_per_window,
                 uint32_t total_windows,
                 uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream);

}

}
