
#include "cudapoa_kernels.cuh"
#include <stdio.h>

namespace nvidia {

namespace cudapoa {

// Device function for running Needleman-Wunsch dynamic programming loop.
__device__
uint16_t runNeedlemanWunsch(uint8_t* nodes,
                        uint16_t* graph,
                        uint16_t* node_id_to_pos,
                        uint16_t graph_count,
                        uint16_t* incoming_edge_count,
                        uint16_t* incoming_edges,
                        uint16_t* outgoing_edge_count,
                        uint16_t* outgoing_edges,
                        uint8_t* read,
                        uint16_t read_count,
                        int32_t* scores,
                        int16_t* traceback_i,
                        int16_t* traceback_j)
{
    //printf("Running NW\n");
    // Set gap/mismatch penalty. Currently acquired from default racon settings.
    // TODO: Pass scores from arguments.
    const int32_t GAP = -8;
    const int32_t MISMATCH = -6;
    const int32_t MATCH = 8;

    //printf("graph len %d, read len %d\n", graph_count, read_count);

    // Init vertical boundary (graph).
    for(uint16_t graph_pos = 0; graph_pos < graph_count; graph_pos++)
    {
        uint16_t node_id = graph[graph_pos];
        uint16_t i = graph_pos + 1;
        uint16_t pred_count = incoming_edge_count[node_id];
        if (pred_count == 0)
        {
            scores[i * MAX_DIMENSION] = GAP;
        }
        else
        {
            int32_t penalty = INT_MIN;
            for(uint16_t p = 0; p < pred_count; p++)
            {
                uint16_t pred_node_id = incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + p];
                uint16_t pred_node_graph_pos = node_id_to_pos[pred_node_id] + 1;
                //printf("pred score %d at pos %d\n", 
                //        scores[pred_node_graph_pos * MAX_DIMENSION],
                //        pred_node_graph_pos);
                //printf("node id %d parent id %d\n", node_id, pred_node_id);
                penalty = max(penalty, scores[pred_node_graph_pos * MAX_DIMENSION]);
            }
            scores[i * MAX_DIMENSION] = penalty + GAP;
        }
        //printf("node %c, score %d\n", nodes[node_id], scores[(graph_pos+1) * MAX_DIMENSION]);
    }

    // Init horizonal boundary conditions (read).
    for(uint16_t j = 1; j < read_count + 1; j++)
    {
        scores[j] = j * GAP;
    }

    // Run DP loop for calculating scores.
    int32_t max_score = INT_MIN;
    int16_t max_i = -1;
    int16_t max_j = -1;

    //printf("iterate through grf\n");

    // Iterate through nodes in graph.
    for(uint16_t graph_pos = 0; graph_pos < graph_count; graph_pos++)
    {
        uint16_t node_id = graph[graph_pos];
        uint16_t i = graph_pos + 1;
        int32_t* scores_i = &scores[i * MAX_DIMENSION];
        //printf("%c%03d ", nodes[node_id], node_id);

        uint16_t pred_count = incoming_edge_count[node_id];

        uint16_t pred_i = (pred_count == 0 ? 0 :
                node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES]] + 1);
        int32_t* scores_pred_i = &scores[pred_i * MAX_DIMENSION];

        // Iterate through bases in sequence and fill out score for diagonal move
        // and vertical move.
        for(uint16_t read_pos = 0; read_pos < read_count; read_pos++)
        {
            int32_t char_profile = (nodes[node_id] == read[read_pos] ? MATCH : MISMATCH);
            // Index into score matrix.
            uint16_t j = read_pos + 1;
            scores_i[j] = max(scores_pred_i[j-1] + char_profile,
                    scores_pred_i[j] + GAP);
            //printf("%d ", scores_i[j]);
            //printf("%d %d\n", i, j);
        }
        //printf("\n");

        // Perform same score updates as above, but for rest of predecessors.
        for (uint16_t p = 1; p < pred_count; p++)
        {
            uint16_t pred_i = node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + p]] + 1;
            scores_pred_i = &scores[pred_i * MAX_DIMENSION];

            for(uint16_t read_pos = 0; read_pos < read_count; read_pos++)
            {
                int32_t char_profile = (nodes[node_id] == read[read_pos] ? MATCH : MISMATCH);
                uint16_t j = read_pos + 1;
                scores_i[j] = max(scores_pred_i[j - 1] + char_profile,
                        max(scores_i[j], scores_pred_i[j] + GAP));
            }
        }

        // Perform score updates for horizontal moves.
        for(uint16_t read_pos = 0; read_pos < read_count; read_pos++)
        {
            // Index into score matrix.
            uint16_t j = read_pos + 1;

            scores_i[j] = max(scores_i[j-1] + GAP, scores_i[j]);
            //printf("%d ", scores_i[j]);

            //if (j == read_count)
            //{
            //    printf("oe %d ie %d for pos %d id %d\n", 
            //            outgoing_edge_count[node_id], 
            //            incoming_edge_count[node_id], 
            //            i, node_id);
            //}

            // Once last column of read is reached, update max score
            // and positions.
            if (j == read_count && outgoing_edge_count[node_id] == 0)
            {
                if (max_score < scores_i[j])
                {
                    max_score = scores_i[j];
                    max_i = i;
                    max_j = j;
                }
            }
        }
        //printf("\n");
    }
        //printf("\n");

    //for(uint32_t i = 0; i < graph_count + 1; i++)
    //{
    //    printf("%d\n", scores[i * MAX_DIMENSION]);
    //}
    //for(uint32_t i = 0; i < read_count + 1; i++)
    //{
    //    printf("%d ", scores[i]);
    //}

    // Fill in backtrace
    int16_t i = max_i;
    int16_t j = max_j;
    int16_t prev_i = 0;
    int16_t prev_j = 0;

    //printf("maxi %d maxj %d\n", i, j);

    // Trace back from maximum score position to generate alignment.
    // Trace back is done by re-calculating the score at each cell
    // along the path to see which preceding cell the move could have
    // come from. This seems computaitonally more expensive, but doesn't
    // require storing any traceback buffer during alignment.
    uint16_t aligned_nodes = 0;
    while(!(i == 0 && j == 0))
    {
        //printf("%d %d\n", i, j);
        int32_t scores_ij = scores[i * MAX_DIMENSION + j];
        bool pred_found = false;

        // Check if move is diagonal.
        if (i != 0 && j != 0)
        {
            uint16_t node_id = graph[i - 1];
            int32_t match_cost = (nodes[node_id] == read[j-1] ? MATCH : MISMATCH);

            uint16_t pred_count = incoming_edge_count[node_id];
            uint16_t pred_i = (pred_count == 0 ? 0 :
                    node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES]] + 1);

            //printf("j %d\n", j-1);
            if (scores_ij == scores[pred_i * MAX_DIMENSION + (j - 1)] + match_cost)
            {
                prev_i = pred_i;
                prev_j = j - 1;
                pred_found = true;
            }

            if (!pred_found)
            {
                for(uint16_t p = 1; p < pred_count; p++)
                {
                    pred_i = node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + p]] + 1;

                    if (scores_ij == scores[pred_i * MAX_DIMENSION + (j - 1)] + match_cost)
                    {
                        prev_i = pred_i;
                        prev_j = j - 1;
                        pred_found = true;
                        break;
                    }
                }
            }
        }

        // Check if move is vertical.
        if (!pred_found && i != 0)
        {
            uint16_t node_id = graph[i - 1];
            uint16_t pred_count = incoming_edge_count[node_id];
            uint16_t pred_i = (pred_count == 0 ? 0 :
                    node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES]] + 1);

            if (scores_ij == scores[pred_i * MAX_DIMENSION + j] + GAP)
            {
                prev_i = pred_i;
                prev_j = j;
                pred_found = true;
            }

            if (!pred_found)
            {
                for(uint16_t p = 1; p < pred_count; p++)
                {
                    pred_i = node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + p]] + 1;

                    if (scores_ij == scores[pred_i * MAX_DIMENSION + j] + GAP)
                    {
                        prev_i = pred_i;
                        prev_j = j;
                        pred_found = true;
                        break;
                    }
                }
            }
        }

        // Check if move is horizontal.
        if (!pred_found && scores_ij == scores[i * MAX_DIMENSION + (j - 1)] + GAP)
        {
            prev_i = i;
            prev_j = j - 1;
            pred_found = true;
        }

        //printf("(%c, %c)\n", (i == prev_i ? '-' : nodes[graph[i-1]]),
        //        (j == prev_j ? '-' : read[j-1]));
        traceback_i[aligned_nodes] = (i == prev_i ? -1 : graph[i-1]);
        traceback_j[aligned_nodes] = (j == prev_j ? -1 : j-1);
        aligned_nodes++;

        //printf("%d %d\n", traceback_i[aligned_nodes - 1], traceback_j[aligned_nodes-1]);

        i = prev_i;
        j = prev_j;
    }

    //printf("aligned nodes %d\n", aligned_nodes);

    //for(int16_t pos = aligned_nodes - 1; pos >= aligned_nodes - 20; pos--)
    //for(int16_t pos = 0; pos < 20; pos++)
    //for(int16_t pos = aligned_nodes - 1; pos >= 0; pos--)
    //{
    //    printf("(%c, %c)\n", 
    //            (traceback_i[pos] == -1 ? '-' : nodes[traceback_i[pos]]), 
    //            (traceback_j[pos] == -1 ? '-' : read[traceback_j[pos]]));
    //    printf("(%d, %d) ", traceback_i[pos], traceback_j[pos]);
    //}
    //printf("\n");

    //if (aligned_nodes == 0)
    //{
    //    printf("Found an alignment with 0 length\n");
    //}

    return aligned_nodes;
}


// Kernel for running Needleman-Wunsch independently.
__global__
void needlemanWunschKernel(uint8_t* nodes,
                           uint16_t* graph,
                           uint16_t* node_id_to_pos,
                           uint16_t graph_count,
                           uint16_t* incoming_edge_count,
                           uint16_t* incoming_edges,
                           uint16_t* outgoing_edge_count,
                           uint16_t* outgoing_edges,
                           uint8_t* read,
                           uint16_t read_count,
                           int32_t* scores,
                           int16_t* traceback_i,
                           int16_t* traceback_j)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
    {
    runNeedlemanWunsch(nodes,
                       graph,
                       node_id_to_pos,
                       graph_count,
                       incoming_edge_count,
                       incoming_edges,
                       outgoing_edge_count,
                       outgoing_edges,
                       read,
                       read_count,
                       scores,
                       traceback_i,
                       traceback_j);
    }
}

// Host function for needleman wunsch kernel.
void needlemanWunsch(uint8_t* nodes,
                    uint16_t* graph,
                    uint16_t* node_id_to_pos,
                    uint16_t graph_count,
                    uint16_t* incoming_edge_count,
                    uint16_t* incoming_edges,
                    uint16_t* outgoing_edge_count,
                    uint16_t* outgoing_edges,
                    uint8_t* read,
                    uint16_t read_count,
                    int32_t* scores,
                    int16_t* traceback_i,
                    int16_t* traceback_j,
                    uint32_t num_threads, uint32_t num_blocks, cudaStream_t stream)
{
    needlemanWunschKernel<<<num_blocks, num_threads, 0, stream>>>(nodes,
                                                                  graph,
                                                                  node_id_to_pos,
                                                                  graph_count,
                                                                  incoming_edge_count,
                                                                  incoming_edges,
                                                                  outgoing_edge_count,
                                                                  outgoing_edges,
                                                                  read,
                                                                  read_count,
                                                                  scores,
                                                                  traceback_i,
                                                                  traceback_j);
}

}

}

