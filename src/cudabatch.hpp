/*!
* @file cudabatch.hpp
 *
 * @brief CUDA batch class header file
 */

#pragma once

#include <memory>
#include <cuda_runtime_api.h>

#include "window.hpp"

namespace racon {

class Window;

class CUDABatch;
std::unique_ptr<CUDABatch> createCUDABatch();

class CUDABatch
{
    // Maximum number of seqeunces and targets.
    // Used to pre-allocate host an cuda buffers.
    const uint32_t MAX_SEQUENCES = 10000;
    const uint32_t MAX_WINDOWS = 100;

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
     * @brief Set CUDA stream for the batch.
     *
     * @param[in] stream : CUDA stream to use for cuda processing.
     */
    void setCUDAStream(cudaStream_t stream);

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
    friend std::unique_ptr<CUDABatch> createCUDABatch();

protected:
    CUDABatch();
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
     */
    void generatePOA();

    // Windows belonging to the batch.
    std::vector<std::shared_ptr<Window>> windows_;

    // Current count of sequences across all windows in batch.
    uint32_t sequence_count_;

    // CUDA stream for launching kernels.
    cudaStream_t stream_;

};

} // namespace racon
