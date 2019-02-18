/*!
* @file cudabatch.hpp
 *
 * @brief CUDA batch class header file
 */

#pragma once

#include <memory>
#include <cuda_runtime_api.h>

#include "../window.hpp"

namespace racon {

class Window;

class CUDABatch;
std::unique_ptr<CUDABatch> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth);

class CUDABatch
{
    const uint32_t MAX_SEQUENCE_SIZE = 2048;

public:
    ~CUDABatch();

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
     * @return Boolean indicating succesful generation of consensus
     *         for all windows.
     */
    bool generateConsensus();

    /**
     * @brief Resets the state of the object, which includes
     *        resetting buffer states and counters.
     */
    void reset();

    // Builder function to create a new CUDABatch object.
    friend std::unique_ptr<CUDABatch> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth);

protected:
    /**
     * @brief Constructor for CUDABatch class.
     *
     * @param[in] max_windows      : Maximum number windows in batch
     * @param[in] max_window_depth : Maximum number of sequences per window
     */
    CUDABatch(uint32_t max_windows, uint32_t max_window_depth);
    CUDABatch(const CUDABatch&) = delete;
    const CUDABatch& operator=(const CUDABatch&) = delete;

    /*
     * @brief Checks if a new window can fit into the batch.
     *
     * @param[in] window : Window to check for fit.
     *
     * @return Whether it can be fit or not.
     */
    bool doesWindowFit(std::shared_ptr<Window>) const;

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

    // Data limits.
    uint32_t max_windows_;
    uint32_t max_depth_per_window_;

    // Windows belonging to the batch.
    std::vector<std::shared_ptr<Window>> windows_;

    // Current count of sequences across all windows in batch.
    uint32_t sequence_count_;

    // CUDA stream for launching kernels.
    cudaStream_t stream_;

    // Buffers for storing results.
    std::unique_ptr<uint8_t[]> consensus_h_;
    uint8_t *consensus_d_;
    size_t consensus_pitch_;

    // Buffers for input data.
    std::unique_ptr<uint8_t[]> inputs_h_;
    uint8_t *inputs_d_;
    size_t input_pitch_;
};

} // namespace racon
