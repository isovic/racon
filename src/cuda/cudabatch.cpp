/*!
 * @file cudabatch.cpp
 *
 * @brief CUDABatch class source file
 */

#include <string>
#include <iostream>
#include <cstring>
#include <algorithm>

#include "cudabatch.hpp"
#include "cudautils.hpp"

namespace racon {

uint32_t CUDABatchProcessor::batches = 0;

std::unique_ptr<CUDABatchProcessor> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth, uint32_t device)
{
    return std::unique_ptr<CUDABatchProcessor>(new CUDABatchProcessor(max_windows, max_window_depth, device));
}

CUDABatchProcessor::CUDABatchProcessor(uint32_t max_windows, uint32_t max_window_depth, uint32_t device)
    : windows_()
    , max_windows_(max_windows)
    , cudapoa_batch_(max_windows, max_window_depth)
{
    bid_ = CUDABatchProcessor::batches++;

    //Set the device
    CU_CHECK_ERR(cudaSetDevice(device));
    cudapoa_batch_.set_device_id(device);

    // Create new CUDA stream.
    CU_CHECK_ERR(cudaStreamCreate(&stream_));
    cudapoa_batch_.set_cuda_stream(stream_);
}

CUDABatchProcessor::~CUDABatchProcessor()
{
    // Destroy CUDA stream.
    CU_CHECK_ERR(cudaStreamDestroy(stream_));
}

bool CUDABatchProcessor::doesWindowFit(std::shared_ptr<Window> window) const
{
    // Checks if adding new window will go over either the MAX_SEQUENCES
    // or max_windows_ count of the batch.
    return (windows_.size() + 1 <= max_windows_);
}

bool CUDABatchProcessor::addWindow(std::shared_ptr<Window> window)
{
    if (doesWindowFit(window))
    {
        windows_.push_back(window);
        return true;
    }
    else
    {
        return false;
    }
}

bool CUDABatchProcessor::hasWindows() const
{
    return (windows_.size() != 0);
}

void CUDABatchProcessor::generateMemoryMap()
{
    auto num_windows = windows_.size();
    for(uint32_t i = 0; i < num_windows; i++)
    {
        // Add new poa
        nvidia::cudapoa::status s = cudapoa_batch_.add_poa();
        if (s != nvidia::cudapoa::CUDAPOA_SUCCESS)
        {
            fprintf(stderr, "Failed to add new to batch %d.\n",
                    cudapoa_batch_.batch_id());
            exit(1);
        }

        std::shared_ptr<Window> window = windows_.at(i);
        uint32_t num_seqs = window->sequences_.size();
        for(uint32_t j = 0; j < num_seqs; j++)
        {
            // Add sequences to latest poa in batch.
            std::pair<const char*, uint32_t> seq = window->sequences_.at(j);
            nvidia::cudapoa::status s = cudapoa_batch_.add_seq_to_poa(seq.first, seq.second);
            if (s != nvidia::cudapoa::CUDAPOA_SUCCESS)
            {
                fprintf(stderr, "Could not add sequence to POA in batch %d.\n",
                        cudapoa_batch_.batch_id());
                exit(1);
            }
        }
    }
}

void CUDABatchProcessor::generatePOA()
{
    // call generate poa function
    cudapoa_batch_.generate_poa();
}

void CUDABatchProcessor::getConsensus()
{
    const std::vector<std::string>& consensuses = cudapoa_batch_.get_consensus();

    for(uint32_t i = 0; i < windows_.size(); i++)
    {
        auto window = windows_.at(i);

        // This is a special case borrowed from the CPU version.
        // TODO: We still run this case through the GPU, but could take it out.
        if (window->sequences_.size() < 3)
        {
            window->consensus_ = std::string(window->sequences_.front().first,
                                             window->sequences_.front().second);

            // This status is borrowed from the CPU version which considers this
            // a failed consensus. All other cases are true.
            window_consensus_status_.emplace_back(false);
        }
        else
        {
            window->consensus_ = consensuses.at(i);
            window_consensus_status_.emplace_back(true);
#ifdef DEBUG
            printf("%s\n", window->consensus_.c_str());
#endif
        }
    }
}

const std::vector<bool>& CUDABatchProcessor::generateConsensus()
{
    // Generate consensus for all windows in the batch
    generateMemoryMap();
    generatePOA();
    getConsensus();

    return window_consensus_status_;
}

void CUDABatchProcessor::reset()
{
    windows_.clear();
    window_consensus_status_.clear();
    cudapoa_batch_.reset();
}

} // namespace racon
