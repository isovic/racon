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
    : max_windows_(max_windows)
    , cudapoa_batch_(max_windows, max_window_depth)
    , windows_()
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

bool CUDABatchProcessor::addWindow(std::shared_ptr<Window> window)
{
    if (windows_.size() < max_windows_)
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
    for(uint32_t w = 0; w < num_windows; w++)
    {
        // Add new poa
        nvidia::cudapoa::status s = cudapoa_batch_.add_poa();
        if (s != nvidia::cudapoa::CUDAPOA_SUCCESS)
        {
            fprintf(stderr, "Failed to add new to batch %d.\n",
                    cudapoa_batch_.batch_id());
            exit(1);
        }

        std::shared_ptr<Window> window = windows_.at(w);
        uint32_t num_seqs = window->sequences_.size();

        // Add first sequence as backbone to graph.
        std::pair<const char*, uint32_t> seq = window->sequences_.front();
        cudapoa_batch_.add_seq_to_poa(seq.first, seq.second);

        // Add the rest of the sequences in sorted order of starting positions.
        std::vector<uint32_t> rank;
        rank.reserve(window->sequences_.size());

        for (uint32_t i = 0; i < num_seqs; ++i) {
            rank.emplace_back(i);
        }

        std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
                return window->positions_[lhs].first < window->positions_[rhs].first; });

        // Start from index 1 since first sequence has already been added as backbone.
        for(uint32_t j = 1; j < num_seqs; j++)
        {
            uint32_t i = rank.at(j);
            // Add sequences to latest poa in batch.
            seq = window->sequences_.at(i);
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
