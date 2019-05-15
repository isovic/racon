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

#include <cudautils/cudautils.hpp>

namespace racon {

std::atomic<uint32_t> CUDABatchProcessor::batches;

std::unique_ptr<CUDABatchProcessor> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth, uint32_t device, int8_t gap, int8_t mismatch, int8_t match, bool cuda_banded_alignment)
{
    return std::unique_ptr<CUDABatchProcessor>(new CUDABatchProcessor(max_windows, max_window_depth, device, gap, mismatch, match, cuda_banded_alignment));
}

CUDABatchProcessor::CUDABatchProcessor(uint32_t max_windows, uint32_t max_window_depth, uint32_t device, int8_t gap, int8_t mismatch, int8_t match, bool cuda_banded_alignment)
    : max_windows_(max_windows)
    , cudapoa_batch_(genomeworks::cudapoa::create_batch(max_windows, max_window_depth, device, gap, mismatch, match, cuda_banded_alignment))
    , windows_()
{
    bid_ = CUDABatchProcessor::batches++;
    
    // Create new CUDA stream.
    GW_CU_CHECK_ERR(cudaStreamCreate(&stream_));
    cudapoa_batch_->set_cuda_stream(stream_);
}

CUDABatchProcessor::~CUDABatchProcessor()
{
    // Destroy CUDA stream.
    GW_CU_CHECK_ERR(cudaStreamDestroy(stream_));
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
        genomeworks::cudapoa::StatusType s = cudapoa_batch_->add_poa();
        if (s != genomeworks::cudapoa::StatusType::success)
        {
            fprintf(stderr, "Failed to add new to batch %d.\n",
                    cudapoa_batch_->batch_id());
            exit(1);
        }

        std::shared_ptr<Window> window = windows_.at(w);
        uint32_t num_seqs = window->sequences_.size();

        // Add first sequence as backbone to graph.
        std::pair<const char*, uint32_t> seq = window->sequences_.front();
        cudapoa_batch_->add_seq_to_poa(seq.first, seq.second);

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
            genomeworks::cudapoa::StatusType s = cudapoa_batch_->add_seq_to_poa(seq.first, seq.second);
            if (s == genomeworks::cudapoa::StatusType::exceeded_maximum_sequence_size)
            {
                fprintf(stderr, "Sequence length exceeds allowed size, discarding this sequence.\n");
                continue;
            } 
            else if (s == genomeworks::cudapoa::StatusType::exceeded_maximum_sequences_per_poa) 
            {
                fprintf(stderr, "More sequences than allowed size, discarding the extra sequences.\n");
                break;
            } 
            else if (s != genomeworks::cudapoa::StatusType::success)
            {
                fprintf(stderr, "Could not add sequence to POA in batch %d.\n",
                        cudapoa_batch_->batch_id());
                exit(1);
            }
        }
    }
}

void CUDABatchProcessor::generatePOA()
{
    // call generate poa function
    cudapoa_batch_->generate_poa();
}

void CUDABatchProcessor::getConsensus()
{
    std::vector<std::string> consensuses;
    std::vector<std::vector<uint16_t>> coverages;
    cudapoa_batch_->get_consensus(consensuses, coverages);

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
            if (window->type_ ==  WindowType::kTGS)
            {
                uint32_t average_coverage = (window->sequences_.size() - 1) / 2;

                int32_t begin = 0, end =  window->consensus_.size() - 1;
                for (; begin < static_cast<int32_t>( window->consensus_.size()); ++begin) {
                    if (coverages.at(i).at(begin) >= average_coverage) {
                        break;
                    }
                }
                for (; end >= 0; --end) {
                    if (coverages.at(i).at(end) >= average_coverage) {
                        break;
                    }
                }

                if (begin >= end) {
                    fprintf(stderr, "[CUDABatchProcessor] warning: "
                            "contig might be chimeric in window %lu!\n", window->id_);
                } else {
                    window->consensus_ =  window->consensus_.substr(begin, end - begin + 1);
                }
            }
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
    cudapoa_batch_->reset();
}

} // namespace racon
