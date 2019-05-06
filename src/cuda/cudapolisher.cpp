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
#include <cudautils/cudautils.hpp>

#include "bioparser/bioparser.hpp"

#pragma message("TODO: Include logger/logger.hpp in cudapolisher.cpp")

namespace racon {

CUDAPolisher::CUDAPolisher(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
    std::unique_ptr<bioparser::Parser<Overlap>> oparser,
    std::unique_ptr<bioparser::Parser<Sequence>> tparser,
    PolisherType type, uint32_t window_length, double quality_threshold,
    double error_threshold, int8_t match, int8_t mismatch, int8_t gap,
    uint32_t num_threads, uint32_t cuda_batches)
        : Polisher(std::move(sparser), std::move(oparser), std::move(tparser),
                type, window_length, quality_threshold,
                error_threshold, match, mismatch, gap, num_threads)
        , cuda_batches_(cuda_batches)
{
#ifdef DEBUG
    window_length_ = 200;
    std::cerr << "In DEBUG mode. Using window size of " << window_length_ << std::endl;
#endif

    genomeworks::cudapoa::Init();

    GW_CU_CHECK_ERR(cudaGetDeviceCount(&num_devices_));

    if (num_devices_ < 1)
    {
        throw std::runtime_error("No GPU devices found.");
    }

    std::cerr << "Using " << num_devices_ << " GPU(s) to perform polishing" << std::endl;

    // Run dummy call on each device to initialize CUDA context.
    for(int32_t dev_id = 0; dev_id < num_devices_; dev_id++)
    {
        std::cerr << "Initialize device " << dev_id << std::endl;
        GW_CU_CHECK_ERR(cudaSetDevice(dev_id));
        GW_CU_CHECK_ERR(cudaFree(0));
    }

    std::cerr << "[CUDAPolisher] Constructed." << std::endl;
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
    std::cerr << "[CUDAPolisher] Allocating memory on GPUs." << std::endl;
    std::chrono::high_resolution_clock::time_point alloc_start = std::chrono::high_resolution_clock::now();
    const uint32_t MAX_WINDOWS = 256;
    const uint32_t MAX_DEPTH_PER_WINDOW = 500;

    // Bin batches into each GPU.
    std::vector<uint32_t> batches_per_gpu(num_devices_, 0);

#ifdef DEBUG
    for(uint32_t i = 0; i < 1; i++)
#else
    for(uint32_t i = 0; i < cuda_batches_; i++)
#endif
    {
        uint32_t device = i % num_devices_;
        batches_per_gpu.at(device) = batches_per_gpu.at(device) + 1;
    }

    for(int32_t device = 0; device < num_devices_; device++)
    {
        for(uint32_t batch = 0; batch < batches_per_gpu.at(device); batch++)
        {
            batch_processors_.emplace_back(createCUDABatch(MAX_WINDOWS, MAX_DEPTH_PER_WINDOW, device));
        }
    }

    std::chrono::high_resolution_clock::time_point polish_start = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(polish_start - alloc_start).count();
    std::cerr << "[CUDAPolisher] Allocated memory in " << ((double)duration / 1000) << " s." << std::endl;

    std::cerr << "[CUDAPolisher] Starting CUDA polish." << std::endl;

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

    std::chrono::high_resolution_clock::time_point polish_end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(polish_end - polish_start).count();
    std::cerr << "[CUDAPolisher] Polished in " << ((double)duration / 1000) << " s." << std::endl;

    // Clear POA processors.
    batch_processors_.clear();
}

}
