/*!
 * @file cudaaligner.cpp
 *
 * @brief CUDABatchAligner class source file
 */

#include <claraparabricks/genomeworks/utils/cudautils.hpp>

#include "cudaaligner.hpp"

namespace racon {

using namespace claraparabricks::genomeworks::cudaaligner;

std::atomic<uint32_t> CUDABatchAligner::batches;

std::unique_ptr<CUDABatchAligner> createCUDABatchAligner(uint32_t max_bandwidth,
                                                         uint32_t device_id,
                                                         int64_t max_gpu_memory)
{
    return std::unique_ptr<CUDABatchAligner>(new CUDABatchAligner(max_bandwidth,
                                                                  device_id,
                                                                  max_gpu_memory));
}

CUDABatchAligner::CUDABatchAligner(uint32_t max_bandwidth,
                                   uint32_t device_id,
                                   int64_t max_gpu_memory)
    : overlaps_()
    , stream_(0)
{
    bid_ = CUDABatchAligner::batches++;

    GW_CU_CHECK_ERR(cudaSetDevice(device_id));

    GW_CU_CHECK_ERR(cudaStreamCreate(&stream_));

    aligner_ = create_aligner(AlignmentType::global_alignment,
                              max_bandwidth,
                              stream_,
                              device_id,
                              max_gpu_memory);
}

CUDABatchAligner::~CUDABatchAligner()
{
    aligner_.reset();
    GW_CU_CHECK_ERR(cudaStreamDestroy(stream_));
}

bool CUDABatchAligner::addOverlap(Overlap* overlap, std::vector<std::unique_ptr<Sequence>>& sequences)
{
    const char* q = !overlap->strand_ ? &(sequences[overlap->q_id_]->data()[overlap->q_begin_]) :
        &(sequences[overlap->q_id_]->reverse_complement()[overlap->q_length_ - overlap->q_end_]);
    int32_t q_len = overlap->q_end_ - overlap->q_begin_;
    const char* t = &(sequences[overlap->t_id_]->data()[overlap->t_begin_]);
    int32_t t_len = overlap->t_end_ - overlap->t_begin_;

    // NOTE: The cudaaligner API for adding alignments is the opposite of edlib. Hence, what is
    // treated as target in edlib is query in cudaaligner and vice versa.
    StatusType s = aligner_->add_alignment(t, t_len,
                                                                       q, q_len);
    if (s == StatusType::exceeded_max_alignments)
    {
        return false;
    }
    else if (s == StatusType::exceeded_max_alignment_difference
             || s == StatusType::exceeded_max_length)
    {
        // Do nothing as this case will be handled by CPU aligner.
    }
    else if (s != StatusType::success)
    {
        fprintf(stderr, "Unknown error in cuda aligner!\n");
    }
    else
    {
        overlaps_.push_back(overlap);
    }
    return true;
}

void CUDABatchAligner::alignAll()
{
    aligner_->align_all();
}

void CUDABatchAligner::generate_cigar_strings()
{
    aligner_->sync_alignments();

    const std::vector<std::shared_ptr<Alignment>>& alignments = aligner_->get_alignments();
    // Number of alignments should be the same as number of overlaps.
    if (overlaps_.size() != alignments.size())
    {
        throw std::runtime_error("Number of alignments doesn't match number of overlaps in cudaaligner.");
    }
    for(std::size_t a = 0; a < alignments.size(); a++)
    {
        overlaps_[a]->cigar_ = alignments[a]->convert_to_cigar();
    }
}

void CUDABatchAligner::reset()
{
    overlaps_.clear();
    cpu_overlap_data_.clear();
    aligner_->reset();
}

}
