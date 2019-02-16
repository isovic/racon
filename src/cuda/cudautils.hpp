// Implementation file for CUDA POA utilities.

#pragma once

#include <stdlib.h>
#include <cuda_runtime_api.h>

namespace racon {

void cudaCheckError()
{
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)
    {
        fprintf(stderr, "Encountered CUDA error ; %s\n", cudaGetErrorString(error));
        exit(-1);
    }
}

} // namespace racon
