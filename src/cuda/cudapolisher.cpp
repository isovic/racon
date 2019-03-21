/*!
 * @file cudapolisher.cpp
 *
 * @brief CUDA Polisher class source file
 */

#include <future>
#include <iostream>
#include <chrono>
#include <cuda_profiler_api.h>

#include "sequence.hpp"
#include "cudapolisher.hpp"

#include "bioparser/bioparser.hpp"

#pragma message("TODO: Include logger/logger.hpp in cudapolisher.cpp")

namespace racon {

CUDAPolisher::CUDAPolisher(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
    std::unique_ptr<bioparser::Parser<Overlap>> oparser,
    std::unique_ptr<bioparser::Parser<Sequence>> tparser,
    PolisherType type, uint32_t window_length, double quality_threshold,
    double error_threshold, int8_t match, int8_t mismatch, int8_t gap,
    uint32_t num_threads)
        : Polisher(std::move(sparser), std::move(oparser), std::move(tparser),
                type, window_length, quality_threshold,
                error_threshold, match, mismatch, gap, num_threads)
{
    std::cerr << "[CUDAPolisher] Constructed." << std::endl;

    const uint32_t MAX_WINDOWS = 256;
    const uint32_t MAX_DEPTH_PER_WINDOW = 500;

    int32_t num_devices;
    CU_CHECK_ERR(cudaGetDeviceCount(&num_devices));


#ifdef DEBUG
    std::cerr << "In DEBUG mode. Using window size of 200." << std::endl;
    window_length_ = 200;
    for(uint32_t i = 0; i < 1; i++)
#else
    for(uint32_t i = 0; i < 6; i++) //TODO: Make the number of batch processors a CLI arg
#endif
    {
        uint32_t device = i % num_devices;
        batch_processors_.emplace_back(createCUDABatch(MAX_WINDOWS, MAX_DEPTH_PER_WINDOW, device));
    }
}

CUDAPolisher::~CUDAPolisher()
{
    cudaDeviceSynchronize();
    cudaProfilerStop();
}

std::pair<uint32_t, uint32_t> CUDAPolisher::fillNextBatchOfWindows(uint32_t batch_id)
{
    batch_processors_.at(batch_id)->reset();

    // Use mutex to read the vector containing windows in a threadsafe manner.
    std::lock_guard<std::mutex> guard(mutex_windows_);

    // TODO: Reducing window wize by 10 for debugging.
    uint32_t initial_count = next_window_index_;
#ifdef DEBUG
    uint32_t count = 5001;//windows_.size();
#else
    uint32_t count = windows_.size();
#endif
    while(next_window_index_ < count)
    {
        if (batch_processors_.at(batch_id)->addWindow(windows_.at(next_window_index_)))
        {
            next_window_index_++;
        }
        else
        {
            break;
        }
    }
    if (next_window_index_ - initial_count > 0)
    {
        fprintf(stderr, "Processing windows %d - %d (of %d) in batch %d\n",
                initial_count,
                next_window_index_ ,
                count,
                batch_processors_.at(batch_id)->getBatchID());
    }

    return std::pair<uint32_t, uint32_t>(initial_count, next_window_index_);
}

void CUDAPolisher::processBatch(uint32_t batch_id)
{
    while(true)
    {
        std::pair<uint32_t, uint32_t> range = fillNextBatchOfWindows(batch_id);
        if (batch_processors_.at(batch_id)->hasWindows())
        {
            // Launch workload.
            const std::vector<bool>& results = batch_processors_.at(batch_id)->generateConsensus();

            // Check if the number of batches processed is same as the range of
            // of windows that were added.
            if (results.size() != (range.second - range.first))
            {
                throw std::runtime_error("Windows processed doesn't match \
                        range of windows passed to batch\n");
            }

            // Copy over the results from the batch into the per window
            // result vector of the CUDAPolisher.
            for(uint32_t i = 0; i < results.size(); i++)
            {
                window_consensus_status_.at(range.first + i) = results.at(i);
            }
        }
        else
        {
            break;
        }
    }
}

void CUDAPolisher::polish(std::vector<std::unique_ptr<Sequence>>& dst,
    bool drop_unpolished_sequences)
{
    // Dummy polish function.
    std::cerr << "Starting CUDA polish" << std::endl;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // Initialize window consensus statuses.
    window_consensus_status_.resize(windows_.size(), false);

    // Process each of the batches in a separate thread.
    std::vector<std::future<void>> thread_futures;
    for(uint32_t i = 0; i < batch_processors_.size(); i++)
    {
        thread_futures.emplace_back(std::async(std::launch::async,
                                               &CUDAPolisher::processBatch,
                                               this,
                                               i)
                                   );
    }

    // Wait for threads to finish, and collect their results.
    for(uint32_t i = 0; i < thread_futures.size(); i++)
    {
        thread_futures.at(i).wait();
    }

    // Collect results from all windows into final output.
    std::string polished_data = "";
    uint32_t num_polished_windows = 0;

    for (uint64_t i = 0; i < windows_.size(); ++i) {

        num_polished_windows += window_consensus_status_.at(i) == true ? 1 : 0;
        polished_data += windows_[i]->consensus();

        if (i == windows_.size() - 1 || windows_[i + 1]->rank() == 0) {
            double polished_ratio = num_polished_windows /
                static_cast<double>(windows_[i]->rank() + 1);
            //double polished_ratio = 1.0f;

            if (!drop_unpolished_sequences || polished_ratio > 0) {
                std::string tags = type_ == PolisherType::kF ? "r" : "";
                tags += " LN:i:" + std::to_string(polished_data.size());
                tags += " RC:i:" + std::to_string(targets_coverages_[windows_[i]->id()]);
                tags += " XC:f:" + std::to_string(polished_ratio);
                dst.emplace_back(createSequence(sequences_[windows_[i]->id()]->name() +
                    tags, polished_data));
            }

            num_polished_windows = 0;
            polished_data.clear();
        }
        windows_[i].reset();
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cerr << "[CUDAPolisher] Polished in " << ((double)duration / 1000) << " s." << std::endl;
}

}
