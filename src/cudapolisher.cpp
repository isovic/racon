/*!
 * @file cudapolisher.cpp
 *
 * @brief CUDA Polisher class source file
 */

#include <future>
#include <iostream>

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
    std::cout << "[CUDAPolisher] Constructed." << std::endl;

    for(uint32_t i = 0; i < num_threads; i++)
    {
        batches_.emplace_back(createCUDABatch());
    }
}

CUDAPolisher::~CUDAPolisher()
{
    // Empty destructor for now.
}

void CUDAPolisher::fillNextBatchOfWindows(uint32_t batch_id)
{
    batches_.at(batch_id)->reset();

    // Use mutex to read the vector containing windows in a threadsafe manner.
    std::lock_guard<std::mutex> guard(mutex_windows_);

    std::cout << "Processing batches for " << batch_id << std::endl;
    while(next_window_index_ < windows_.size())
    {
        if (batches_.at(batch_id)->addWindow(windows_.at(next_window_index_)))
        {
            std::cout << "    Added window " << next_window_index_ << std::endl;
            next_window_index_++;
        }
        else
        {
            break;
        }
    }
}

bool CUDAPolisher::processBatch(uint32_t batch_id)
{
    bool result = true;
    while(result)
    {
        fillNextBatchOfWindows(batch_id);
        if (batches_.at(batch_id)->hasWindows())
        {
            // Create new CUDA stream.
            cudaStream_t new_stream;
            cudaStreamCreate(&new_stream);
            batches_.at(batch_id)->setCUDAStream(new_stream);

            // Launch workload.
            result = batches_.at(batch_id)->generateConsensus();
            // Stub sleep command to simulate workload.
            usleep(1000);
        }
        else
        {
            break;
        }
    }
    return result;
}

void CUDAPolisher::polish(std::vector<std::unique_ptr<Sequence>>& dst,
    bool drop_unpolished_sequences)
{
    // Dummy polish function.
    std::cout << "Starting CUDA polish" << std::endl;

    // Process each of the batches in a separate thread.
    std::vector<std::future<bool>> thread_futures;
    for(uint32_t i = 0; i < batches_.size(); i++)
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
        if (!thread_futures.at(i).get())
        {
            std::cerr << "Batch " << i << " had issues \
                processing its windows." << std::endl;
        }
    }

    (void) dst;
    (void) drop_unpolished_sequences;
    std::cout << "[CUDAPolisher] Polished." << std::endl;
}

}
