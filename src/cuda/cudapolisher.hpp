/*!
 * @file cudapolisher.hpp
 *
 * @brief CUDA Polisher class header file
 */

#pragma once

#include <mutex>

#include "polisher.hpp"
#include "cudabatch.hpp"

namespace racon {

class CUDAPolisher : public Polisher {
public:
    ~CUDAPolisher();

    virtual void polish(std::vector<std::unique_ptr<Sequence>>& dst,
        bool drop_unpolished_sequences) override;

    friend std::unique_ptr<Polisher> createPolisher(const std::string& sequences_path,
        const std::string& overlaps_path, const std::string& target_path,
        PolisherType type, uint32_t window_length, double quality_threshold,
        double error_threshold, int8_t match, int8_t mismatch, int8_t gap,
        uint32_t num_threads, uint32_t cuda_batches);

protected:
    CUDAPolisher(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
        std::unique_ptr<bioparser::Parser<Overlap>> oparser,
        std::unique_ptr<bioparser::Parser<Sequence>> tparser,
        PolisherType type, uint32_t window_length, double quality_threshold,
        double error_threshold, int8_t match, int8_t mismatch, int8_t gap,
        uint32_t num_threads, uint32_t cuda_batches);
    CUDAPolisher(const CUDAPolisher&) = delete;
    const CUDAPolisher& operator=(const CUDAPolisher&) = delete;

    // Insert new windows into the batch referred to by batch_id. Return
    // the range of windows added to the batch. Interval closed in front, open
    // at the end.
    std::pair<uint32_t, uint32_t> fillNextBatchOfWindows(uint32_t batch_id);

    // Generate POA for all windows in the batch.
    void processBatch(uint32_t batch_id);

    // Vector of batches. Generated during construction time.
    std::vector<std::unique_ptr<CUDABatchProcessor>> batch_processors_;

    // Vector of bool indicating consensus generation status for each window.
    std::vector<bool> window_consensus_status_;

    // Mutex for accessing the vector of windows.
    std::mutex mutex_windows_;

    // Index of next window to be added to a batch.
#ifdef DEBUG
    uint32_t next_window_index_ = 5000;
#else
    uint32_t next_window_index_ = 0;
#endif

    // Number of batches for POA processing.
    uint32_t cuda_batches_;

    // Number of GPU devices to run with.
    int32_t num_devices_;
};

}
