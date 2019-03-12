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
        uint32_t num_threads, bool use_cuda);

protected:
    CUDAPolisher(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
        std::unique_ptr<bioparser::Parser<Overlap>> oparser,
        std::unique_ptr<bioparser::Parser<Sequence>> tparser,
        PolisherType type, uint32_t window_length, double quality_threshold,
        double error_threshold, int8_t match, int8_t mismatch, int8_t gap,
        uint32_t num_threads);
    CUDAPolisher(const CUDAPolisher&) = delete;
    const CUDAPolisher& operator=(const CUDAPolisher&) = delete;

    // Insert new windows into the batch referred to by batch_id.
    void fillNextBatchOfWindows(uint32_t batch_id);

    // Generate POA for all windows in the batch.
    bool processBatch(uint32_t batch_id);

    // Vector of batches. Generated during construction time.
    std::vector<std::unique_ptr<CUDABatchProcessor>> batch_processors_;

    // Mutex for accessing the vector of windows.
    std::mutex mutex_windows_;

    // Index of next window to be added to a batch.
    uint32_t next_window_index_ = 0;
};

}
