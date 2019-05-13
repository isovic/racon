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
        , gap_(gap)
        , mismatch_(mismatch)
        , match_(match)
{
#ifdef DEBUG
    window_length_ = 200;
    std::cerr << "In DEBUG mode. Using window size of " << window_length_ << std::endl;
#endif

    genomeworks::cudapoa::Init();
    genomeworks::cudaaligner::Init();

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

void CUDAPolisher::find_overlap_breaking_points(std::vector<std::unique_ptr<Overlap>>& overlaps)
{
    std::cerr << "[CUDAPolisher] Finding overlaps on GPUs." << std::endl;
    std::chrono::high_resolution_clock::time_point overlap_start = std::chrono::high_resolution_clock::now();

    std::mutex mutex_overlaps;
    uint32_t next_overlap_index = 0;

    // Lambda expression for filling up next batch of alignments.
    auto fill_next_batch = [&mutex_overlaps, &next_overlap_index, &overlaps, this](CUDABatchAligner* batch) -> void {
        batch->reset();

        // Use mutex to read the vector containing windows in a threadsafe manner.
        std::lock_guard<std::mutex> guard(mutex_overlaps);

        uint32_t initial_count = next_overlap_index;
        uint32_t count = overlaps.size();
        while(next_overlap_index < count)
        {
            if (batch->addOverlap(overlaps.at(next_overlap_index).get(), sequences_))
            {
                next_overlap_index++;
            }
            else
            {
                break;
            }
        }
        if (next_overlap_index - initial_count > 0)
        {
            fprintf(stderr, "Processing overlaps %d - %d (of %d) in batch %d\n",
                    initial_count,
                    next_overlap_index ,
                    count,
                    batch->getBatchID());
        }
    };

    // Lambda expression for processing a batch of alignments.
    auto process_batch = [&fill_next_batch, this](CUDABatchAligner* batch) -> void {
        while(true)
        {
            fill_next_batch(batch);
            if (batch->hasOverlaps())
            {
                // Launch workload.
                batch->alignAll();
                batch->find_breaking_points(window_length_);
            }
            else
            {
                break;
            }
        }
    };

    // Create batches based on arguments provided to program.
    for(uint32_t batch = 0; batch < cuda_batches_; batch++)
    {
        batch_aligners_.emplace_back(createCUDABatchAligner(20000, 20000, 1000, 0));
    }

    // Run batched alignment.
    std::vector<std::future<void>> thread_futures;
    for(uint32_t i = 0; i < batch_aligners_.size(); i++)
    {
        thread_futures.emplace_back(std::async(std::launch::async,
                                               process_batch,
                                               batch_aligners_.at(i).get())
                                   );
    }

    // Wait for threads to finish, and collect their results.
    for (const auto& future : thread_futures) {
        future.wait();
    }

    batch_aligners_.clear();

    std::chrono::high_resolution_clock::time_point overlap_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(overlap_end - overlap_start).count();
    std::cerr << "[CUDAPolisher] Polished in " << ((double)duration / 1000) << " s." << std::endl;

    // TODO: Kept CPU overlap alignment right now while GPU is a dummy implmentation.
    Polisher::find_overlap_breaking_points(overlaps);
}

void CUDAPolisher::polish(std::vector<std::unique_ptr<Sequence>>& dst,
    bool drop_unpolished_sequences)
{
    // Mutex for accessing the vector of windows.
    std::mutex mutex_windows;

    // Index of next window to be added to a batch.
#ifdef DEBUG
    uint32_t next_window_index = 5000;
#else
    uint32_t next_window_index = 0;
#endif

    // Lambda function for adding windows to batches.
    auto fill_next_batch = [&mutex_windows, &next_window_index, this](CUDABatchProcessor* batch) -> std::pair<uint32_t, uint32_t> {
        batch->reset();

        // Use mutex to read the vector containing windows in a threadsafe manner.
        std::lock_guard<std::mutex> guard(mutex_windows);

        // TODO: Reducing window wize by 10 for debugging.
        uint32_t initial_count = next_window_index;
#ifdef DEBUG
        uint32_t count = 5001;//windows_.size();
#else
        uint32_t count = windows_.size();
#endif
        while(next_window_index < count)
        {
            if (batch->addWindow(windows_.at(next_window_index)))
            {
                next_window_index++;
            }
            else
            {
                break;
            }
        }
        if (next_window_index - initial_count > 0)
        {
            fprintf(stderr, "Processing windows %d - %d (of %d) in batch %d\n",
                    initial_count,
                    next_window_index ,
                    count,
                    batch->getBatchID());
        }

        return std::pair<uint32_t, uint32_t>(initial_count, next_window_index);
    };

    // Lambda function for processing each batch.
    auto process_batch = [&fill_next_batch, this](CUDABatchProcessor* batch) -> void {
        while(true)
        {
            std::pair<uint32_t, uint32_t> range = fill_next_batch(batch);
            if (batch->hasWindows())
            {
                // Launch workload.
                const std::vector<bool>& results = batch->generateConsensus();

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
    };

    // Creation and use of batches.
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
            batch_processors_.emplace_back(createCUDABatch(MAX_WINDOWS, MAX_DEPTH_PER_WINDOW, device, gap_, mismatch_, match_));
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
                                               process_batch,
                                               batch_processors_.at(i).get())
                                   );
    }

    // Wait for threads to finish, and collect their results.
    for (const auto& future : thread_futures) {
        future.wait();
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
