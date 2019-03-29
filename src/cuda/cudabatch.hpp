/*!
* @file cudabatch.hpp
 *
 * @brief CUDA batch class header file
 */

#pragma once

#include <memory>
#include <cuda_runtime_api.h>

#include "window.hpp"
#include "cudapoa_kernels.cuh"

namespace racon {

class Window;

class CUDABatchProcessor;
std::unique_ptr<CUDABatchProcessor> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth, uint32_t device);

class CUDABatchProcessor
{
    const uint32_t NUM_THREADS = 32;

public:
    ~CUDABatchProcessor();

    uint32_t NUM_BLOCKS = 1;
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

    // GPU Device ID
    uint32_t device_id_;

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
    std::vector<bool> window_consensus_status_;

    // Buffers for input data.
    uint8_t *inputs_h_;
    uint8_t *inputs_d_;
    uint8_t * num_sequences_per_window_h_;
    uint16_t * sequence_lengths_h_;
    uint16_t * sequence_lengths_d_;
    nvidia::cudapoa::WindowDetails * window_details_d_;
    nvidia::cudapoa::WindowDetails * window_details_h_;

    // Buffers for temp data.
    int16_t* scores_d_;
    int16_t* traceback_i_d_;
    int16_t* traceback_j_d_;

    // Buffer for temp graph data.
    uint8_t* nodes_d_;
    uint16_t* node_alignments_d_;
    uint16_t* node_alignment_count_d_;
    uint16_t* incoming_edges_d_;
    uint16_t* incoming_edge_count_d_;
    uint16_t* outgoing_edges_d_;
    uint16_t* outgoing_edge_count_d_;
    uint16_t* incoming_edges_weights_d_;
    uint16_t* outoing_edges_weights_d_;
    uint16_t* sorted_poa_d_;
    uint16_t* sorted_poa_node_map_d_;
    uint16_t* sorted_poa_local_edge_count_d_;
    int32_t* consensus_scores_d_;
    int16_t* consensus_predecessors_d_;


    static uint32_t batches;
    uint32_t bid_ = 0;

};

} // namespace racon