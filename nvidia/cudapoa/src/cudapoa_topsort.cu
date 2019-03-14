#include "cudapoa_kernels.cuh"
#include <stdio.h>

namespace nvidia {

namespace cudapoa {

// Device function for running topoligical sort on graph.
__device__
void topologicalSortDeviceUtil(uint16_t* sorted_poa,
                               uint16_t* sorted_poa_node_map,
                               uint16_t node_count,
                               uint16_t* incoming_edge_count,
                               uint16_t* outgoing_edges,
                               uint16_t* outgoing_edge_count,
                               uint16_t* local_incoming_edge_count)
{
    //printf("Running top sort\n");
    // Clear the incoming edge count for each node.
    //__shared__ int16_t local_incoming_edge_count[CUDAPOA_MAX_NODES_PER_WINDOW];
    //memset(local_incoming_edge_count, -1, CUDAPOA_MAX_NODES_PER_WINDOW);
    uint16_t sorted_poa_position = 0;

    // Iterate through node IDs (since nodes are from 0
    // through node_count -1, a simple loop works) and fill 
    // out the incoming edge count.
    for(uint16_t n = 0; n < node_count; n++)
    {
        local_incoming_edge_count[n] = incoming_edge_count[n];
        // If we find a node ID has 0 incoming edges, add it to sorted nodes list.
        if (local_incoming_edge_count[n] == 0)
        {
            sorted_poa_node_map[n] = sorted_poa_position;
            sorted_poa[sorted_poa_position++] = n;
        }
    }

    // Loop through set of node IDs with no incoming edges,
    // then iterate through their children. For each child decrement their 
    // incoming edge count. If incoming edge count of child == 0, 
    // add its node ID to the sorted order list.
    for(uint16_t n = 0; n < sorted_poa_position; n++)
    {
        uint16_t node = sorted_poa[n];
        for(uint16_t edge = 0; edge < outgoing_edge_count[node]; edge++)
        {
            uint16_t out_node = outgoing_edges[node * CUDAPOA_MAX_NODE_EDGES + edge];
            //printf("%d\n", out_node);
            local_incoming_edge_count[out_node]--;
            if (local_incoming_edge_count[out_node] == 0)
            {
                sorted_poa_node_map[out_node] = sorted_poa_position;
                sorted_poa[sorted_poa_position++] = out_node;
            }
        }
    }

    // sorted_poa will have final ordering of node IDs.
}

// Kernel for running topological independently.
__global__
void topologicalSortKernel(uint16_t* sorted_poa_d,
                           uint16_t* sorted_poa_node_map_d,
                           uint16_t node_count,
                           uint16_t* incoming_edge_count_d,
                           uint16_t* outgoing_edges_d,
                           uint16_t* outgoing_edge_count_d)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
    {
        topologicalSortDeviceUtil(sorted_poa_d,
                                  sorted_poa_node_map_d,
                                  node_count,
                                  incoming_edge_count_d,
                                  outgoing_edges_d,
                                  outgoing_edge_count_d,
                                  NULL);
    }
}

// Host function for running topological sort kernel.
void topologicalSort(uint16_t* sorted_poa_d,
                     uint16_t* sorted_poa_node_map_d,
                     uint16_t node_count,
                     uint16_t* incoming_edge_count_d,
                     uint16_t* outgoing_edges_d,
                     uint16_t* outgoing_edge_count_d,
                     uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream)
{
    topologicalSortKernel<<<num_blocks, num_threads, 0, stream>>>(sorted_poa_d,
                                                                  sorted_poa_node_map_d,
                                                                  node_count,
                                                                  incoming_edge_count_d,
                                                                  outgoing_edges_d,
                                                                  outgoing_edge_count_d);
}

}

}
