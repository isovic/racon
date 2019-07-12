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

#include "spoa/spoa.hpp"
#include <cudautils/cudautils.hpp>

namespace racon {

std::atomic<uint32_t> CUDABatchProcessor::batches;

std::unique_ptr<CUDABatchProcessor> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth, uint32_t device, int8_t gap, int8_t mismatch, int8_t match, bool cuda_banded_alignment)
{
    return std::unique_ptr<CUDABatchProcessor>(new CUDABatchProcessor(max_windows, max_window_depth, device, gap, mismatch, match, cuda_banded_alignment));
}

CUDABatchProcessor::CUDABatchProcessor(uint32_t max_windows, uint32_t max_window_depth, uint32_t device, int8_t gap, int8_t mismatch, int8_t match, bool cuda_banded_alignment)
    : max_windows_(max_windows)
    , cudapoa_batch_(claragenomics::cudapoa::create_batch(max_windows, max_window_depth, device, claragenomics::cudapoa::OutputType::consensus, gap, mismatch, match, cuda_banded_alignment))
    , windows_()
    , seqs_added_per_window_()
{
    bid_ = CUDABatchProcessor::batches++;
    
    // Create new CUDA stream.
    CGA_CU_CHECK_ERR(cudaStreamCreate(&stream_));
    cudapoa_batch_->set_cuda_stream(stream_);
}

CUDABatchProcessor::~CUDABatchProcessor()
{
    // Destroy CUDA stream.
    CGA_CU_CHECK_ERR(cudaStreamDestroy(stream_));
}

bool CUDABatchProcessor::addWindow(std::shared_ptr<Window> window)
{
    if (windows_.size() < max_windows_)
    {
        windows_.push_back(window);
        seqs_added_per_window_.push_back(0);
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

void CUDABatchProcessor::convertPhredQualityToWeights(const char* qual,
                                                      uint32_t qual_length,
                                                      std::vector<int8_t>& weights)
{
    weights.clear();
    for(uint32_t i = 0; i < qual_length; i++)
    {
        weights.push_back(static_cast<uint8_t>(qual[i]) - 33); // PHRED quality
    }
}

claragenomics::cudapoa::StatusType CUDABatchProcessor::addSequenceToPoa(std::pair<const char*, uint32_t>& seq,
                                                                      std::pair<const char*, uint32_t>& qualities)
{
    // Add sequences to latest poa in batch.
    std::vector<int8_t> weights;
    claragenomics::cudapoa::StatusType status = claragenomics::cudapoa::StatusType::success;
    if (qualities.first == nullptr)
    {
        status = cudapoa_batch_->add_seq_to_poa(seq.first, nullptr, seq.second);
    }
    else
    {
        convertPhredQualityToWeights(qualities.first, qualities.second, weights);
        status = cudapoa_batch_->add_seq_to_poa(seq.first, weights.data(), seq.second);
    }
    return status;
}

void CUDABatchProcessor::generateMemoryMap()
{
    auto num_windows = windows_.size();
    for(uint32_t w = 0; w < num_windows; w++)
    {
        // Add new poa
        claragenomics::cudapoa::StatusType status = cudapoa_batch_->add_poa();
        if (status != claragenomics::cudapoa::StatusType::success)
        {
            fprintf(stderr, "Failed to add new POA to batch %d.\n",
                    cudapoa_batch_->batch_id());
            exit(1);
        }

        std::shared_ptr<Window> window = windows_.at(w);
        uint32_t num_seqs = window->sequences_.size();
        std::vector<uint8_t> weights;

        // Add first sequence as backbone to graph.
        std::pair<const char*, uint32_t> seq = window->sequences_.front();
        std::pair<const char*, uint32_t> qualities = window->qualities_.front();
        status = addSequenceToPoa(seq, qualities);
        if (status != claragenomics::cudapoa::StatusType::success)
        {
            fprintf(stderr, "Could not add backbone to window. Fatal error.\n");
            exit(1);
        }

        // Add the rest of the sequences in sorted order of starting positions.
        std::vector<uint32_t> rank;
        rank.reserve(window->sequences_.size());

        for (uint32_t i = 0; i < num_seqs; ++i) {
            rank.emplace_back(i);
        }

        std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
                return window->positions_[lhs].first < window->positions_[rhs].first; });

        // Start from index 1 since first sequence has already been added as backbone.
        uint32_t long_seq = 0;
        uint32_t skipped_seq = 0;
        for(uint32_t j = 1; j < num_seqs; j++)
        {
            uint32_t i = rank.at(j);
            seq = window->sequences_.at(i);
            qualities = window->qualities_.at(i);
            // Add sequences to latest poa in batch.
            status = addSequenceToPoa(seq, qualities);
            if (status == claragenomics::cudapoa::StatusType::exceeded_maximum_sequence_size)
            {
                long_seq++;
                continue;
            } 
            else if (status == claragenomics::cudapoa::StatusType::exceeded_maximum_sequences_per_poa)
            {
                skipped_seq++;
                continue;
            } 
            else if (status != claragenomics::cudapoa::StatusType::success)
            {
                fprintf(stderr, "Could not add sequence to POA in batch %d.\n",
                        cudapoa_batch_->batch_id());
                exit(1);
            }

            seqs_added_per_window_[w] = seqs_added_per_window_[w] + 1;
        }
#ifndef NDEBUG
        if (long_seq > 0)
        {
            fprintf(stderr, "Too long (%d / %d)\n", long_seq, num_seqs);
        }
        if (skipped_seq > 0)
        {
            fprintf(stderr, "Skipped (%d / %d)\n", skipped_seq, num_seqs);
        }
#endif
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
    std::vector<claragenomics::cudapoa::StatusType> output_status;
    cudapoa_batch_->get_consensus(consensuses, coverages, output_status);

    for(uint32_t i = 0; i < windows_.size(); i++)
    {
        auto window = windows_.at(i);
        if (output_status.at(i) != claragenomics::cudapoa::StatusType::success)
        {
            // leave the failure cases to CPU polisher
            window_consensus_status_.emplace_back(false);
        }
        else
        {
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
                    uint32_t num_seqs_in_window = seqs_added_per_window_[i];
                    uint32_t average_coverage = num_seqs_in_window / 2;

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
            }
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
    seqs_added_per_window_.clear();
    cudapoa_batch_->reset();
}

} // namespace racon
