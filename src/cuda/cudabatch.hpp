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
    const uint32_t NUM_THREADS = 64;

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

    // Maximum windows to process in batch.
    uint32_t max_windows_;

    // Maximum sequences per window.
    uint32_t max_depth_per_window_;

    // GPU Device ID
    uint32_t device_id_;

    // Windows belonging to the batch.
    std::vector<std::shared_ptr<Window>> windows_;

    // CUDA stream for launching kernels.
    cudaStream_t stream_;

    // Host buffer for storing consensus.
    std::unique_ptr<uint8_t[]> consensus_h_;

    // Device buffer pointer for storing consensus.
    uint8_t *consensus_d_;

    // Pitch of consensus device buffer since it's a 2D array.
    size_t consensus_pitch_;

    // Consensus generation status for each window.
    std::vector<bool> window_consensus_status_;

    // Host and device buffer pointer for input data.
    uint8_t *inputs_h_;
    uint8_t *inputs_d_;

    // Host buffer pointfer number of sequences per window.
    uint8_t * num_sequences_per_window_h_;

    // Host and device buffer for sequence lengths.
    uint16_t * sequence_lengths_h_;
    uint16_t * sequence_lengths_d_;

    // Host and device buffer pointers that hold Window Details struct.
    nvidia::cudapoa::WindowDetails * window_details_d_;
    nvidia::cudapoa::WindowDetails * window_details_h_;

    // Device buffer for the scoring matrix for all windows.
    int16_t* scores_d_;

    // Device buffers for alignment backtrace.
    // i for graph
    // j for sequence
    int16_t* alignment_graph_d_;
    int16_t* alignment_read_d_;

    // Device buffer to store nodes of the graph. The node itself is the base
    // (A, T, C, G) and the id of the node is it's position in the buffer.
    uint8_t* nodes_d_;

    // Device buffer to store the list of nodes aligned to a 
    // specific node in the graph.
    uint16_t* node_alignments_d_;
    uint16_t* node_alignment_count_d_;

    // Device buffer to store incoming edges to a node.
    uint16_t* incoming_edges_d_;
    uint16_t* incoming_edge_count_d_;

    // Device buffer to store outgoing edges from a node.
    uint16_t* outgoing_edges_d_;
    uint16_t* outgoing_edge_count_d_;

    // Devices buffers to store incoming and outgoing edge weights.
    uint16_t* incoming_edges_weights_d_;
    uint16_t* outoing_edges_weights_d_;

    // Device buffer to store the topologically sorted graph. Each element
    // of this buffer is an ID of the node.
    uint16_t* sorted_poa_d_;

    // Device buffer that maintains a mapping between the node ID and its
    // position in the topologically sorted graph.
    uint16_t* sorted_poa_node_map_d_;

    // Device buffer used during topological sort to store incoming
    // edge counts for nodes.
    uint16_t* sorted_poa_local_edge_count_d_;

    // Device buffer to store scores calculated during traversal
    // of graph for consensus generation.
    int32_t* consensus_scores_d_;

    // Device buffer to store the predecessors of nodes during
    // graph traversal.
    int16_t* consensus_predecessors_d_;

    // Static batch count used to generate batch IDs.
    static uint32_t batches;

    // Batch ID.
    uint32_t bid_ = 0;

};

} // namespace racon
