// Header for for CUDA POA host kernel wrappers.

#pragma once

#include <stdint.h>
#include <cuda_runtime_api.h>

namespace nvidia {

namespace cudapoa {

void generatePOA(uint8_t* consensus_d,
                 uint8_t* sequences_d,
                 uint32_t max_sequence_size,
                 uint32_t max_depth_per_window,
                 uint32_t total_windows,
                 uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream);

}

}
