// Header for for CUDA POA host kernel wrappers.

#pragma once

#include <stdint.h>
#include <cuda_runtime_api.h>
#include <stdio.h>

#define CUDAPOA_MAX_NODE_EDGES 50
#define CUDAPOA_MAX_NODE_ALIGNMENTS 50
#define CUDAPOA_MAX_NODES_PER_WINDOW 2048
#define CUDAPOA_MAX_SEQUENCE_SIZE 1024
#define CUDAPOA_MAX_MATRIX_DIMENSION (CUDAPOA_MAX_NODES_PER_WINDOW + 1)

#define CU_CHECK_ERR(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPU Error:: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

namespace nvidia {

namespace cudapoa {

void generatePOA(uint8_t* consensus_d,
                 size_t conensus_pitch,
                 uint8_t* sequences_d,
                 size_t sequences_pitch,
                 uint32_t max_sequence_size,
                 uint16_t * num_sequences_per_window_d,
                 uint16_t * sequence_lengths_d,
                 uint32_t max_depth_per_window,
                 uint32_t total_windows,
                 uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream,
                 int16_t* scores, int16_t* ti, int16_t* tj,
                 uint8_t* nodes,
                 uint16_t* incoming_edges, uint16_t* incoming_edge_count,
                 uint16_t* outgoing_edges, uint16_t* outgoing_edge_count,
                 uint16_t* incoming_edge_w, uint16_t* outgoing_edge_w,
                 uint16_t* sorted_poa, uint16_t* node_id_to_pos,
                 uint16_t* node_alignments, uint16_t* node_alignment_count,
                 uint16_t* sorted_poa_local_edge_count);

void topologicalSort(uint16_t* sorted_poa_d,
                     uint16_t* sorted_poa_node_map_d,
                     uint16_t node_count,
                     uint16_t* incoming_edge_count_d,
                     uint16_t* outoing_edges_d,
                     uint16_t* outgoing_edge_count_d,
                     uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream);

void needlemanWunsch(uint8_t* nodes,
                    uint16_t* graph,
                    uint16_t* graph_node_map,
                    uint16_t graph_count,
                    uint16_t* incoming_edge_count,
                    uint16_t* incoming_edges,
                    uint16_t* outgoing_edge_count,
                    uint16_t* outgoing_edges,
                    uint8_t* read,
                    uint16_t read_count,
                    int16_t* scores,
                    int16_t* traceback_i,
                    int16_t* traceback_j,
                    uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream);

void addAlignmentToGraphHost(uint8_t* nodes,
                         uint16_t node_count,
                         uint16_t* node_alignments, uint16_t* node_alignment_count,
                         uint16_t* incoming_edges, uint16_t* incoming_edge_count,
                         uint16_t* outgoing_edges, uint16_t* outgoing_edge_count,
                         uint16_t* incoming_edge_w, uint16_t* outgoing_edge_w,
                         uint16_t alignment_length,
                         uint16_t* graph,
                         int16_t* alignment_graph,
                         uint8_t* read,
                         int16_t* alignment_read, cudaStream_t stream);

}

}
