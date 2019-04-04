#pragma once

#include <memory>
#include <vector>
#include <string>
#include <iostream>

#include <cuda_runtime_api.h>

#include "cudapoa_kernels.cuh"

namespace nvidia {

namespace cudapoa {

class Batch
{
    const uint32_t NUM_THREADS = 64;

public:
    enum output_type
    {
        consensus,
        multi_sequence_alignment
    };

    Batch(uint32_t max_poas, uint32_t max_sequences_per_poa);
    ~Batch();

    void add_poa();

    void add_seq_to_poa(const char* seq, uint32_t seq_len);

    uint32_t get_total_poas() const;

    void generate_poa();

    const std::vector<std::string>& get_consensus();

    void reset();

    void set_device_id(uint32_t);

    void set_cuda_stream(cudaStream_t stream);

    uint32_t batch_id() const;

protected:
    // Maximum POAs to process in batch.
    uint32_t max_poas_ = 0;

    // Maximum sequences per POA.
    uint32_t max_sequences_per_poa_ = 0;

    // GPU Device ID
    uint32_t device_id_ = 0;

    // CUDA stream for launching kernels.
    cudaStream_t stream_;

    // Host buffer for storing consensus.
    std::unique_ptr<uint8_t[]> consensus_h_;

    // Device buffer pointer for storing consensus.
    uint8_t *consensus_d_;

    // Pitch of consensus device buffer since it's a 2D array.
    size_t consensus_pitch_;

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

    // Device buffer to store node marks when performing spoa accurate topsort.
    uint8_t* node_marks_d_;

    // Device buffer to store check for aligned nodes.
    bool* check_aligned_nodes_d_;

    // Device buffer to store stack for nodes to be visited.
    uint16_t* nodes_to_visit_d_;

    // Static batch count used to generate batch IDs.
    static uint32_t batches;

    // Batch ID.
    uint32_t bid_ = 0;

    // Total POAs added.
    uint32_t poa_count_ = 0;

    // Number of nucleotides already already inserted.
    uint32_t num_nucleotides_copied_ = 0;

    // Global sequence index.
    uint32_t global_sequence_idx_ = 0;
    
    // Vector of consensus results.
    std::vector<std::string> consensus_strings_;

    uint32_t NUM_BLOCKS = 1;

};

}

}
