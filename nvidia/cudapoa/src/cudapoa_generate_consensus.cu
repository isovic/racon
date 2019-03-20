
#include "cudapoa_kernels.cuh"
#include <stdio.h>

namespace nvidia {

namespace cudapoa {

__device__
uint16_t branchCompletion(uint16_t max_score_id_pos,
                         uint8_t* nodes,
                         uint16_t node_count,
                         uint16_t* graph,
                         uint16_t* incoming_edges,
                         uint16_t* incoming_edge_count,
                         uint16_t* outgoing_edges,
                         uint16_t* outgoing_edge_count,
                         uint16_t* incoming_edge_w,
                         int32_t* scores,
                         int16_t* predecessors)
{
    uint16_t node_id = graph[max_score_id_pos];
    uint16_t out_edges = outgoing_edge_count[node_id];
    for(uint16_t oe = 0; oe < out_edges; oe++)
    {
        uint16_t out_node_id = outgoing_edges[node_id * CUDAPOA_MAX_NODE_EDGES + oe];
        uint16_t out_node_in_edges = incoming_edge_count[out_node_id];
        for(uint16_t ie = 0; ie < out_node_in_edges; ie++)
        {
            uint16_t id = incoming_edges[out_node_in_edges * CUDAPOA_MAX_NODE_EDGES + ie];
            if (id != node_id)
            {
                scores[node_id] = -1;
            }
        }
    }

    int32_t max_score = 0;
    uint16_t max_score_id = 0;
    for(uint16_t graph_pos = max_score_id_pos + 1; graph_pos < node_count; graph_pos++)
    {
        node_id = graph[graph_pos];
        predecessors[node_id] = -1;

        int32_t score_node_id = -1;

        uint16_t in_edges = incoming_edge_count[node_id];
        for(uint16_t e = 0; e < in_edges; e++)
        {
            uint16_t begin_node_id = incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + e];
            if (scores[begin_node_id] == -1)
            {
                continue;
            }

            uint16_t edge_w = incoming_edge_w[node_id * CUDAPOA_MAX_NODE_EDGES + e];
            if (score_node_id < edge_w ||
                    (score_node_id == edge_w &&
                     scores[predecessors[node_id]] <= scores[begin_node_id]))
            {
                score_node_id = edge_w;
                predecessors[node_id] = begin_node_id;
            }
        }

        if (predecessors[node_id] != -1)
        {
            score_node_id += scores[predecessors[node_id]];
        }

        if (max_score < score_node_id)
        {
            max_score = score_node_id;
            max_score_id = node_id;
        }

        scores[node_id] = score_node_id;
    }

    return max_score_id;
}

__device__
void generateConsensus(uint8_t* nodes,
                         uint16_t node_count,
                         uint16_t* graph,
                         uint16_t* node_id_to_pos,
                         uint16_t* incoming_edges,
                         uint16_t* incoming_edge_count,
                         uint16_t* outgoing_edges,
                         uint16_t* outgoing_edge_count,
                         uint16_t* incoming_edge_w,
                         int16_t* predecessors,
                         int32_t* scores,
                         uint8_t* consensus)
{
    for(uint16_t i = 0; i < node_count; i++)
    {
        predecessors[i] = -1;
        scores[i] = -1;
    }

    uint16_t max_score_id = 0;
    int32_t max_score = -1;

    for(uint16_t graph_pos = 0; graph_pos < node_count; graph_pos++)
    {
        //printf("%d\n", max_score_id);
        uint16_t node_id = graph[graph_pos];
        //printf("node id %d\n", node_id);
        uint16_t in_edges = incoming_edge_count[node_id];

        int32_t score_node_id = scores[node_id];
        for(uint16_t e = 0; e < in_edges; e++)
        {
            uint16_t edge_w = incoming_edge_w[node_id * CUDAPOA_MAX_NODE_EDGES + e];
            uint16_t begin_node_id = incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + e];
            if (score_node_id < edge_w ||
                    (score_node_id == edge_w &&
                     scores[predecessors[node_id]] <= scores[begin_node_id]))
            {
                score_node_id = edge_w;
                predecessors[node_id] = begin_node_id;
            }
        }

        if (predecessors[node_id] != -1)
        {
            score_node_id += scores[predecessors[node_id]];
        }

        if (max_score < score_node_id)
        {
            max_score_id = node_id;
            max_score = score_node_id;
        }

        scores[node_id] = score_node_id;
    }

    if (outgoing_edge_count[max_score_id] != 0)
    {
        while(outgoing_edge_count[max_score_id] != 0)
        {
            max_score_id = branchCompletion(node_id_to_pos[max_score_id],
                                             nodes,
                                             node_count,
                                             graph,
                                             incoming_edges,
                                             incoming_edge_count,
                                             outgoing_edges,
                                             outgoing_edge_count,
                                             incoming_edge_w,
                                             scores,
                                             predecessors);
        }
    }

    uint16_t consensus_pos = 0;
    while(predecessors[max_score_id] != -1)
    {
        consensus[consensus_pos] = nodes[max_score_id];
        max_score_id = predecessors[max_score_id];
        consensus_pos++;
    }
    consensus[consensus_pos] = nodes[max_score_id];
    consensus++;
    consensus[consensus_pos] = '\0';
}

}

}
