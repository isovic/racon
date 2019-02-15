/*!
 * @file cudabatch.cpp
 *
 * @brief CUDABatch class source file
 */

#include "cudabatch.hpp"

#include <iostream>

namespace racon {

std::unique_ptr<CUDABatch> createCUDABatch()
{
    return std::unique_ptr<CUDABatch>(new CUDABatch());
}

CUDABatch::CUDABatch()
    : windows_()
    , sequence_count_()
    , stream_()
{
    // Create new CUDA stream.
    cudaStreamCreate(&stream_);

    // Allocate host memory and CUDA memory based on max sequence and target counts.
}

CUDABatch::~CUDABatch()
{
    // Destroy CUDA stream.
    cudaStreamDestroy(stream_);

    // Free all the host and CUDA memory.
}

bool CUDABatch::doesWindowFit(std::shared_ptr<Window> window) const
{
    // Checks if adding new window will go over either the MAX_SEQUENCES
    // or MAX_WINDOWS count of the batch.
    return (window->sequences_.size() + sequence_count_ < MAX_SEQUENCES) &&
        (windows_.size() + 1 < MAX_WINDOWS);
}

bool CUDABatch::addWindow(std::shared_ptr<Window> window)
{
    if (doesWindowFit(window))
    {
        windows_.push_back(window);
        sequence_count_ += window->sequences_.size();
        return true;
    }
    else{
        return false;
    }
}

bool CUDABatch::hasWindows() const
{
    return (windows_.size() != 0);
}

void CUDABatch::generateMemoryMap()
{
    // Fill host/cuda memory with sequence information.
}

void CUDABatch::generatePOA()
{
    // Launch kernel to run 1 POA per thread in thread block.
}

bool CUDABatch::generateConsensus()
{
    // Generate consensus for all windows in the batch
    generateMemoryMap();
    generatePOA();

    return true;
}

void CUDABatch::reset()
{
    windows_.clear();
    sequence_count_ = 0;
}

} // namespace racon
