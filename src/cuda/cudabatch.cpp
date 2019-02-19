/*!
 * @file cudabatch.cpp
 *
 * @brief CUDABatch class source file
 */

#include <iostream>
#include <cstring>

#include "cudapoa_kernels.cuh"

#include "cudabatch.hpp"
#include "cudautils.hpp"

namespace racon {

std::unique_ptr<CUDABatchProcessor> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth)
{
    return std::unique_ptr<CUDABatchProcessor>(new CUDABatchProcessor(max_windows, max_window_depth));
}

CUDABatchProcessor::CUDABatchProcessor(uint32_t max_windows, uint32_t max_window_depth)
    : max_windows_(max_windows)
    , max_depth_per_window_(max_window_depth)
    , windows_()
    , sequence_count_()
    , stream_()
    , consensus_pitch_()
{
    // Create new CUDA stream.
    cudaStreamCreate(&stream_);

    // Allocate host memory and CUDA memory based on max sequence and target counts.

    // Input buffers.
    uint32_t input_size = max_windows_ * max_depth_per_window_ * MAX_SEQUENCE_SIZE;
    cudaHostAlloc((void**) &inputs_h_, input_size * sizeof(uint8_t),
                  cudaHostAllocDefault);
    cudaHostAlloc((void**) &num_sequences_per_window_h_, max_windows * sizeof(uint8_t),
            cudaHostAllocDefault);
    cudaHostAlloc((void**) &sequence_lengths_h_, max_windows * max_depth_per_window_* sizeof(uint8_t),
            cudaHostAllocDefault);

    cudaMallocPitch((void**) &inputs_d_,
                    &input_pitch_,
                    sizeof(uint8_t) * MAX_SEQUENCE_SIZE,
                    max_windows_ * max_depth_per_window_);

    cudaMalloc((void**)&num_sequences_per_window_d_, max_windows * sizeof(uint8_t));
    cudaMalloc((void**)&sequence_lengths_d_, max_windows * max_depth_per_window_ * sizeof(uint8_t));

    cudaCheckError();
    std::cout << "Allocated input buffers of size " << input_size << std::endl;

    // Output buffers.
    input_size = max_windows_ * MAX_SEQUENCE_SIZE;
    consensus_h_.reset(new uint8_t[input_size]);

    cudaMallocPitch((void**) &consensus_d_,
                    &consensus_pitch_,
                    sizeof(uint8_t) * MAX_SEQUENCE_SIZE,
                    max_windows_);
    cudaCheckError();
    std::cout << "Allocated output buffers of size " << input_size << std::endl;
}

CUDABatchProcessor::~CUDABatchProcessor()
{
    // Destroy CUDA stream.
    cudaStreamDestroy(stream_);
    cudaCheckError();

    // Free all the host and CUDA memory.
    cudaFree(consensus_d_);
    cudaCheckError();

    std::cout << "Destroyed buffers." << std::endl;
}

bool CUDABatchProcessor::doesWindowFit(std::shared_ptr<Window> window) const
{
    // Checks if adding new window will go over either the MAX_SEQUENCES
    // or max_windows_ count of the batch.
    return (window->sequences_.size() + sequence_count_ < max_windows_ * max_depth_per_window_) &&
        (windows_.size() + 1 < max_windows_);
}

bool CUDABatchProcessor::addWindow(std::shared_ptr<Window> window)
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

bool CUDABatchProcessor::hasWindows() const
{
    return (windows_.size() != 0);
}

void CUDABatchProcessor::generateMemoryMap()
{
    // Fill host/cuda memory with sequence information.
    for(uint32_t i = 0; i < windows_.size(); i++)
    {
        auto window = windows_.at(i);
        uint32_t input_window_offset = i * max_depth_per_window_;
        int num_seqs = 0;
        for(uint32_t j = 0; j < window->sequences_.size(); j++)
        {
            uint32_t input_sequence_offset = input_window_offset + j;
            auto seq = window->sequences_.at(j);

            memcpy(&(inputs_h_[input_sequence_offset * MAX_SEQUENCE_SIZE]),
                   seq.first,
                   seq.second);

            num_seqs++;
            sequence_lengths_h_[i * max_depth_per_window_ + j] = seq.second;
        }
        num_sequences_per_window_h_[i] = num_seqs;
    }

    cudaMemcpy2DAsync(inputs_d_, input_pitch_,
                      inputs_h_.get(), MAX_SEQUENCE_SIZE,
                      MAX_SEQUENCE_SIZE, max_windows_ * max_depth_per_window_,
                      cudaMemcpyHostToDevice, stream_);

    cudaMemcpyAsync(num_sequences_per_window_d_, num_sequences_per_window_h_,
                    max_windows_ * sizeof(uint8_t), cudaMemcpyHostToDevice, stream_);

    cudaMemcpyAsync(sequence_lengths_d_, sequence_lengths_h_,
                    max_depth_per_window_ * max_windows_ * sizeof(uint8_t), cudaMemcpyHostToDevice, stream_);


    cudaCheckError();
}

void CUDABatchProcessor::generatePOA()
{
    // Launch kernel to run 1 POA per thread in thread block.
    uint32_t NUM_THREADS = 32;
    uint32_t num_blocks = (windows_.size() / NUM_THREADS) + 1;
    nvidia::cudapoa::generatePOA(consensus_d_,
                                 inputs_d_,
                                 MAX_SEQUENCE_SIZE,
                                 num_sequences_per_window_d_,
                                 sequence_lengths_d_,
                                 max_depth_per_window_,
                                 windows_.size(),
                                 NUM_THREADS,
                                 num_blocks,
                                 stream_);
    cudaCheckError();
}

void CUDABatchProcessor::getConsensus()
{
    cudaMemcpy2DAsync(consensus_h_.get(),
                      MAX_SEQUENCE_SIZE,
                      consensus_d_,
                      consensus_pitch_,
                      MAX_SEQUENCE_SIZE,
                      max_windows_,
                      cudaMemcpyDeviceToHost,
                      stream_);
    cudaCheckError();

    cudaStreamSynchronize(stream_);
    cudaCheckError();

    for(uint32_t i = 0; i < windows_.size(); i++)
    {
        char* c = reinterpret_cast<char *>(&consensus_h_[i * MAX_SEQUENCE_SIZE]);
        windows_.at(i)->consensus_ = std::string(c);
    }
}

bool CUDABatchProcessor::generateConsensus()
{
    // Generate consensus for all windows in the batch
    generateMemoryMap();
    generatePOA();
    getConsensus();

    return true;
}

void CUDABatchProcessor::reset()
{
    windows_.clear();
    sequence_count_ = 0;

    // Clear host and device memory.
    memset(&inputs_h_[0], 0, max_windows_ * max_depth_per_window_ * MAX_SEQUENCE_SIZE);
    cudaMemsetAsync(inputs_d_, 0, max_windows_ * max_depth_per_window_ * input_pitch_, stream_);
    cudaCheckError();
}

} // namespace racon
