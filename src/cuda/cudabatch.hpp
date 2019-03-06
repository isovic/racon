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

class CUDABatchProcessor;
std::unique_ptr<CUDABatchProcessor> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth);

class CUDABatchProcessor
{
    const uint32_t MAX_SEQUENCE_SIZE = 1024;
    const uint32_t NUM_THREADS = 128;
    const uint32_t NUM_BLOCKS = 1;

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
     * @return Boolean indicating succesful generation of consensus
     *         for all windows.
     */
    bool generateConsensus();

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
    friend std::unique_ptr<CUDABatchProcessor> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth);

protected:
    /**
     * @brief Constructor for CUDABatch class.
     *
     * @param[in] max_windows      : Maximum number windows in batch
     * @param[in] max_window_depth : Maximum number of sequences per window
     */
    CUDABatchProcessor(uint32_t max_windows, uint32_t max_window_depth);
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
    uint16_t * num_sequences_per_window_h_;
    uint16_t * sequence_lengths_h_;
    uint16_t * num_sequences_per_window_d_;
    uint16_t * sequence_lengths_d_;

    // Buffers for temp data.
    int32_t* scores_d_;
    int32_t* scores_h_;
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


    static uint32_t batches;
    uint32_t bid_ = 0;

};

} // namespace racon
