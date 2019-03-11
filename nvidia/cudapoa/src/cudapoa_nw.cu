
#include "cudapoa_kernels.cuh"
#include <stdio.h>

namespace nvidia {

namespace cudapoa {

// Device function for running Needleman-Wunsch dynamic programming loop.
__device__
uint16_t runNeedlemanWunsch(uint8_t* nodes,
        //uint8_t* nodes_global,
                        uint16_t* graph,
                        uint16_t* node_id_to_pos_global,
                        //uint16_t* node_id_to_pos,
                        uint16_t graph_count,
                        uint16_t* incoming_edge_count,
                        uint16_t* incoming_edges,
                        //uint16_t* incoming_edges_global,
                        uint16_t* outgoing_edge_count,
                        uint16_t* outgoing_edges,
                        //uint8_t* read_global,
                        uint8_t* read,
                        uint16_t read_count,
                        int16_t* scores,
                        int16_t* traceback_i,
                        int16_t* traceback_j)
{
    //printf("Running NW\n");
    // Set gap/mismatch penalty. Currently acquired from default racon settings.
    // TODO: Pass scores from arguments.
    const int16_t GAP = -8;
    const int16_t MISMATCH = -6;
    const int16_t MATCH = 8;

    __shared__ int16_t max_scores[256];
    __shared__ int16_t max_is[256];
    __shared__ int16_t max_js[256];

    __shared__ uint16_t node_id_to_pos[1024];
    __shared__ uint16_t local_incoming_edges[256 * CUDAPOA_MAX_NODE_EDGES];

    uint32_t thread_idx = threadIdx.x;

    max_scores[thread_idx] = SHRT_MIN;
    max_is[thread_idx] = -1;
    max_js[thread_idx] = -1;

    long long int start = clock64();
    long long int init = 0;

    if (thread_idx == 0)
    {
        //memcpy(incoming_edges, incoming_edges_global, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * 50);
        //printf("graph %d, read %d\n", graph_count, read_count);

        //memset(scores, INT_MIN, sizeof(int32_t) * 
        //        CUDAPOA_MAX_MATRIX_DIMENSION * CUDAPOA_MAX_MATRIX_DIMENSION);

        // Init vertical boundary (graph).
        for(uint16_t graph_pos = 0; graph_pos < graph_count; graph_pos++)
        {
            node_id_to_pos[graph_pos] = node_id_to_pos_global[graph_pos];
            //nodes[graph_pos] = nodes_global[graph_pos];

            scores[0] = 0;
            uint16_t node_id = graph[graph_pos];
            uint16_t i = graph_pos + 1;
            uint16_t pred_count = incoming_edge_count[node_id];
            if (pred_count == 0)
            {
                scores[i * CUDAPOA_MAX_MATRIX_DIMENSION] = GAP;
            }
            else
            {
                int16_t penalty = SHRT_MIN;
                for(uint16_t p = 0; p < pred_count; p++)
                {
                    uint16_t pred_node_id = incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + p];
                    uint16_t pred_node_graph_pos = node_id_to_pos_global[pred_node_id] + 1;
                    //printf("pred score %d at pos %d\n", 
                    //        scores[pred_node_graph_pos * CUDAPOA_MAX_MATRIX_DIMENSION],
                    //        pred_node_graph_pos);
                    //printf("node id %d parent id %d\n", node_id, pred_node_id);
                    penalty = max(penalty, scores[pred_node_graph_pos * CUDAPOA_MAX_MATRIX_DIMENSION]);
                }
                scores[i * CUDAPOA_MAX_MATRIX_DIMENSION] = penalty + GAP;
            }
            //printf("node %c, score %d\n", nodes[node_id], scores[(graph_pos+1) * CUDAPOA_MAX_MATRIX_DIMENSION]);
        }

        // Init horizonal boundary conditions (read).
        for(uint16_t j = 1; j < read_count + 1; j++)
        {
            scores[j] = j * GAP;
            //read[j-1] = read_global[j-1];
        }
        init = clock64() - start;

    }

    __syncthreads();

    //if (graph_count == 479 && read_count == 13)
    //{
        //for(uint32_t i = 0; i < graph_count; i++)
        //{
        //    printf("node-%d pos %d %d %d, ", i, /*node_id_to_pos[i],*/
        //            scores[(node_id_to_pos[i] + 1) * CUDAPOA_MAX_MATRIX_DIMENSION],
        //            incoming_edge_count[i],
        //            outgoing_edge_count[i]);
        //    for(uint16_t j  = 0; j < incoming_edge_count[i]; j++)
        //    {
        //        printf("%d ", incoming_edges[i * CUDAPOA_MAX_NODE_EDGES + j]);
        //    }
        //    printf(", ");
        //    for(uint16_t j  = 0; j < outgoing_edge_count[i]; j++)
        //    {
        //        printf("%d ", outgoing_edges[i * CUDAPOA_MAX_NODE_EDGES + j]);
        //    }
        //    printf("\n");
        //}
        //for(uint32_t i = 0; i < read_count + 1; i++)
        //{
        //    printf("%d ", scores[i]);
        //}
    //}

    // Run DP loop for calculating scores.
    //printf("iterate through grf\n");
    uint16_t max_diag = min(graph_count, read_count);

    start = clock64();

    for(uint16_t c = 0; c < (graph_count + read_count -1); c++)
    {
        uint16_t i_0 = min(c, graph_count - 1);
        uint16_t j_0 = c - i_0;
        uint16_t diag = min(max_diag, min(c + 1, graph_count + read_count - (c + 1)));

        //printf("i0 %d j0 %d diag %d\n", i_0, j_0, diag);

        for(uint32_t t = thread_idx; t < diag; t += blockDim.x)
        {
            uint16_t i = i_0 - t + 1;
            uint16_t j = j_0 + t + 1;

            //printf("%d %d\n", i, j);

            uint16_t graph_pos = i - 1;
            uint16_t read_pos = j - 1;
            uint16_t node_id = graph[graph_pos];

            //int16_t* scores_i = &scores[i * CUDAPOA_MAX_MATRIX_DIMENSION];
            int16_t score_ij = SHRT_MIN;
            //printf("%c%03d ", nodes[node_id], node_id);

            uint16_t pred_count = incoming_edge_count[node_id];
            memcpy(&local_incoming_edges[thread_idx * CUDAPOA_MAX_NODE_EDGES],
                   &incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES],
                   sizeof(uint16_t) * pred_count);

            uint16_t pred_i = (pred_count == 0 ? 0 :
                    node_id_to_pos[local_incoming_edges[thread_idx * CUDAPOA_MAX_NODE_EDGES]] + 1);
            //int16_t* scores_pred_i = &scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION];

            // Iterate through bases in sequence and fill out score for diagonal move
            // and vertical move.
            int16_t char_profile = (nodes[node_id] == read[read_pos] ? MATCH : MISMATCH);
            // Index into score matrix.
            //scores_i[j] = max(scores_pred_i[j-1] + char_profile,
            //        scores_pred_i[j] + GAP);

            //score_ij = max(scores_pred_i[j-1] + char_profile,
            //        scores_pred_i[j] + GAP);
            score_ij = max(scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION + j-1] + char_profile,
                    scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION + j] + GAP);

            //printf("%d ", scores_i[j]);
            //printf("%d %d\n", i, j);
            //printf("\n");

            // Perform same score updates as above, but for rest of predecessors.
            for (uint16_t p = 1; p < pred_count; p++)
            {
                uint16_t pred_i = node_id_to_pos[local_incoming_edges[thread_idx * CUDAPOA_MAX_NODE_EDGES + p]] + 1;
                //scores_pred_i = &scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION];

                //int32_t char_profile = (nodes[node_id] == read[read_pos] ? MATCH : MISMATCH);
                score_ij = max(scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION + j - 1] + char_profile,
                        max(score_ij, scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION + j] + GAP));
            }

            // Perform score updates for horizontal moves.
            // Index into score matrix.

            //score_ij = max(scores_i[j-1] + GAP, score_ij);
            score_ij = max(scores[i * CUDAPOA_MAX_MATRIX_DIMENSION + (j-1)] + GAP, score_ij);
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
                if (max_scores[thread_idx] < score_ij)
                {
                    max_scores[thread_idx] = score_ij;
                    max_is[thread_idx] = i;
                    max_js[thread_idx] = j;
                }
            }

            //scores_i[j] = score_ij;
            scores[i * CUDAPOA_MAX_MATRIX_DIMENSION + j] = score_ij;
        }

        __syncthreads();
    }

    long long int nw = clock() - start;

    //int32_t max_score = INT_MIN;
    //int16_t max_i = -1;
    //int16_t max_j = -1;

    // Iterate through nodes in graph.
    //for(uint16_t graph_pos = 0; graph_pos < graph_count; graph_pos++)
    //{
    //    uint16_t node_id = graph[graph_pos];
    //    uint16_t i = graph_pos + 1;
    //    int32_t* scores_i = &scores[i * CUDAPOA_MAX_MATRIX_DIMENSION];
    //    //printf("%c%03d ", nodes[node_id], node_id);

    //    uint16_t pred_count = incoming_edge_count[node_id];

    //    uint16_t pred_i = (pred_count == 0 ? 0 :
    //            node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES]] + 1);
    //    int32_t* scores_pred_i = &scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION];

    //    // Iterate through bases in sequence and fill out score for diagonal move
    //    // and vertical move.
    //    for(uint16_t read_pos = 0; read_pos < read_count; read_pos++)
    //    {
    //        int32_t char_profile = (nodes[node_id] == read[read_pos] ? MATCH : MISMATCH);
    //        // Index into score matrix.
    //        uint16_t j = read_pos + 1;
    //        scores_i[j] = max(scores_pred_i[j-1] + char_profile,
    //                scores_pred_i[j] + GAP);
    //        //printf("%d ", scores_i[j]);
    //        //printf("%d %d\n", i, j);
    //    }
    //    //printf("\n");

    //    // Perform same score updates as above, but for rest of predecessors.
    //    for (uint16_t p = 1; p < pred_count; p++)
    //    {
    //        uint16_t pred_i = node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + p]] + 1;
    //        scores_pred_i = &scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION];

    //        for(uint16_t read_pos = 0; read_pos < read_count; read_pos++)
    //        {
    //            int32_t char_profile = (nodes[node_id] == read[read_pos] ? MATCH : MISMATCH);
    //            uint16_t j = read_pos + 1;
    //            scores_i[j] = max(scores_pred_i[j - 1] + char_profile,
    //                    max(scores_i[j], scores_pred_i[j] + GAP));
    //        }
    //    }

    //    // Perform score updates for horizontal moves.
    //    for(uint16_t read_pos = 0; read_pos < read_count; read_pos++)
    //    {
    //        // Index into score matrix.
    //        uint16_t j = read_pos + 1;

    //        scores_i[j] = max(scores_i[j-1] + GAP, scores_i[j]);
    //        //printf("%d ", scores_i[j]);

    //        //if (j == read_count)
    //        //{
    //        //    printf("oe %d ie %d for pos %d id %d\n", 
    //        //            outgoing_edge_count[node_id], 
    //        //            incoming_edge_count[node_id], 
    //        //            i, node_id);
    //        //}

    //        // Once last column of read is reached, update max score
    //        // and positions.
    //        if (j == read_count && outgoing_edge_count[node_id] == 0)
    //        {
    //            if (max_score < scores_i[j])
    //            {
    //                max_score = scores_i[j];
    //                max_i = i;
    //                max_j = j;
    //            }
    //        }
    //    }
    //    //printf("\n");
    //}
        //printf("\n");

    //if (graph_count == 479 && read_count == 13)
    //{
    //    for(uint16_t i = 0; i < graph_count; i++)
    //    {
    //    uint16_t     a = node_id_to_pos[i] + 1;
    //        for(uint16_t j = 1; j < read_count+1; j++)
    //        {
    //            printf("%d ", scores[a*CUDAPOA_MAX_MATRIX_DIMENSION + j]);
    //        }
    //        printf("\n");
    //    }
    //    return;
    //}
    start = clock64();

    long long int tb = 0;

    uint16_t aligned_nodes = 0;
    if (thread_idx == 0)
    {
        int16_t i = -1;
        int16_t j = -1;
        int16_t max_score = SHRT_MIN;
        for(uint32_t t = 0; t < blockDim.x; t++)
        {
            //printf("%d %d\n", t, max_scores[t]);
            if (max_scores[t] >= max_score)
            {
                i = max_is[t];
                j = max_js[t];
                max_score = max_scores[t];
            }
        }
        // Fill in backtrace

        int16_t prev_i = 0;
        int16_t prev_j = 0;

        //printf("maxi %d maxj %d\n", i, j);

        // Trace back from maximum score position to generate alignment.
        // Trace back is done by re-calculating the score at each cell
        // along the path to see which preceding cell the move could have
        // come from. This seems computaitonally more expensive, but doesn't
        // require storing any traceback buffer during alignment.
        while(!(i == 0 && j == 0))
        {
            //printf("%d %d\n", i, j);
            int16_t scores_ij = scores[i * CUDAPOA_MAX_MATRIX_DIMENSION + j];
            bool pred_found = false;
            //if (graph_count == 479 && read_count == 13)
            //{
            //    printf("%d %d node %d\n", i, j, graph[i-1]);
            //}

            // Check if move is diagonal.
            if (i != 0 && j != 0)
            {
                uint16_t node_id = graph[i - 1];
                int16_t match_cost = (nodes[node_id] == read[j-1] ? MATCH : MISMATCH);

                uint16_t pred_count = incoming_edge_count[node_id];
                uint16_t pred_i = (pred_count == 0 ? 0 :
                        (node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES]] + 1));

                //printf("j %d\n", j-1);
                if (scores_ij == (scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION + (j - 1)] + match_cost))
                {
                    prev_i = pred_i;
                    prev_j = j - 1;
                    pred_found = true;
                }

                if (!pred_found)
                {
                    for(uint16_t p = 1; p < pred_count; p++)
                    {
                        pred_i = (node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + p]] + 1);

                        if (scores_ij == (scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION + (j - 1)] + match_cost))
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

                if (scores_ij == scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION + j] + GAP)
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

                        if (scores_ij == scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION + j] + GAP)
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
            if (!pred_found && scores_ij == scores[i * CUDAPOA_MAX_MATRIX_DIMENSION + (j - 1)] + GAP)
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

        tb = clock64() - start;
    }

    //if (thread_idx == 0)
    //{
    //    long long int total = init + nw + tb;
    //    printf("Total time of init is %lf %\n", ((double)init / total) * 100.f);
    //    printf("Total time of nw is %lf %\n", ((double)nw / total) * 100.f);
    //    printf("Total time of tb is %lf %\n", ((double)tb / total) * 100.f);
    //}


    //for(int16_t pos = aligned_nodes - 1; pos >= aligned_nodes - 20; pos--)
    //for(int16_t pos = 0; pos < 20; pos++)
    //if (graph_count == 479 && read_count == 13)
    //{
    //    printf("maxi %d maxj %d\n", max_i, max_j);
    //    for(int16_t pos = 0; pos < aligned_nodes; pos++)
    //    {
    //        //printf("(%c, %c)\n", 
    //        //        (traceback_i[pos] == -1 ? '-' : nodes[traceback_i[pos]]), 
    //        //        (traceback_j[pos] == -1 ? '-' : read[traceback_j[pos]]));
    //        printf("(%d, %d) ", traceback_i[pos], traceback_j[pos]);
    //    }
    //    printf("\n");
    //}

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
                           int16_t* scores,
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
                    int16_t* scores,
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

