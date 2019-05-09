/*!
* @file cudabatch.hpp
 *
 * @brief CUDA batch class header file
 */

#pragma once

#include <memory>
#include <cuda_runtime_api.h>
#include <atomic>

#include "window.hpp"
#include "cudapoa/batch.hpp"

namespace racon {

class Window;

class CUDABatchProcessor;
std::unique_ptr<CUDABatchProcessor> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth, uint32_t device);

class CUDABatchProcessor
{
public:
    ~CUDABatchProcessor();

    /**
     * @brief Add a new window to the batch.
     *
     * @param[in] window : The window to add to the batch.
     *
     * @return True of window could be added to the batch.
     */
    bool addWindow(std::shared_ptr<Window> window);

    /**
     * @brief Checks if batch has any windows to process.
     */
    bool hasWindows() const;

    /**
     * @brief Runs the core computation to generate consensus for
     *        all windows in the batch.
     *
     * @return Vector of bool indicating succesful generation of consensus
     *         for each window in the batch.
     */
    const std::vector<bool>& generateConsensus();

    /**
     * @brief Resets the state of the object, which includes
     *        resetting buffer states and counters.
     */
    void reset();

    /**
     * @brief Get batch ID.
     */
    uint32_t getBatchID() const { return bid_; }

    // Builder function to create a new CUDABatchProcessor object.
    friend std::unique_ptr<CUDABatchProcessor>
    createCUDABatch(uint32_t max_windows, uint32_t max_window_depth, uint32_t device);

protected:
    /**
     * @brief Constructor for CUDABatch class.
     *
     * @param[in] max_windows      : Maximum number windows in batch
     * @param[in] max_window_depth : Maximum number of sequences per window
     */
    CUDABatchProcessor(uint32_t max_windows, uint32_t max_window_depth, uint32_t device);
    CUDABatchProcessor(const CUDABatchProcessor&) = delete;
    const CUDABatchProcessor& operator=(const CUDABatchProcessor&) = delete;

    /*
     * @brief Process all the windows and re-map them into
     *        memory for more efficient processing in the CUDA
     *        kernels.
     */
    void generateMemoryMap();

    /*
     * @brief Run the CUDA kernel for generating POA on the batch.
     *        This call is asynchronous.
     */
    void generatePOA();

    /*
     * @brief Wait for execution to complete and grab the output
     *        consensus from the device.
     */
    void getConsensus();

protected:
    // Static batch count used to generate batch IDs.
    static std::atomic<uint32_t> batches;

    // Batch ID.
    uint32_t bid_ = 0;

    // Maximum windows allowed in batch.
    uint32_t max_windows_;

    // CUDA-POA library object that manages POA batch.
    std::unique_ptr<genomeworks::cudapoa::Batch> cudapoa_batch_;

    // Stream for running POA batch.
    cudaStream_t stream_;
    // Windows belonging to the batch.
    std::vector<std::shared_ptr<Window>> windows_;

    // Consensus generation status for each window.
    std::vector<bool> window_consensus_status_;
};

} // namespace racon
