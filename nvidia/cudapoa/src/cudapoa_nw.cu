
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
                        //uint16_t* incoming_edge_count_global,
                        uint16_t* incoming_edge_count,
                        uint16_t* incoming_edges,
                        //uint16_t* incoming_edges_global,
                        uint16_t* outgoing_edge_count_global,
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

    //__shared__ uint16_t local_incoming_edges[CUDAPOA_MAX_NODE_EDGES];
    __shared__ uint16_t node_id_to_pos[1024];
    //__shared__ uint16_t incoming_edge_count[1024];
    //__shared__ uint16_t outgoing_edge_count[1024];

    __shared__ int16_t score_i[1024];
    //__shared__ int16_t score_prev_i[1024];

    uint32_t thread_idx = threadIdx.x;

    long long int start = clock64();
    long long int init = 0;

    for(uint16_t graph_pos = thread_idx; graph_pos < graph_count; graph_pos += blockDim.x)
    {
        node_id_to_pos[graph_pos] = node_id_to_pos_global[graph_pos];
        //incoming_edge_count[graph_pos] = incoming_edge_count_global[graph_pos];
    }

    __syncthreads();

    if (thread_idx == 0)
    {
        //printf("graph %d, read %d\n", graph_count, read_count);

        // Init vertical boundary (graph).
        for(uint16_t graph_pos = 0; graph_pos < graph_count; graph_pos++)
        {
            //node_id_to_pos[graph_pos] = node_id_to_pos_global[graph_pos];
            //incoming_edge_count[graph_pos] = incoming_edge_count_global[graph_pos];
            //outgoing_edge_count[graph_pos]= outgoing_edge_count_global[graph_pos];
            //nodes[graph_pos] = nodes_global[graph_pos];

            scores[0] = 0;
            uint16_t node_id = graph[graph_pos];
            uint16_t i = graph_pos + 1;
            //uint16_t pred_count = incoming_edge_count_global[node_id];
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
                    //uint16_t pred_node_graph_pos = node_id_to_pos_global[pred_node_id] + 1;
                    uint16_t pred_node_graph_pos = node_id_to_pos[pred_node_id] + 1;
                    //printf("pred score %d at pos %d\n", 
                    //        scores[pred_node_graph_pos * CUDAPOA_MAX_MATRIX_DIMENSION],
                    //        pred_node_graph_pos);
                    //printf("node id %d parent id %d\n", node_id, pred_node_id);
                    penalty = max(penalty, scores[pred_node_graph_pos * CUDAPOA_MAX_MATRIX_DIMENSION]);
                }
                scores[i * CUDAPOA_MAX_MATRIX_DIMENSION] = penalty + GAP;
            }

            //printf("%d \n", scores[i * CUDAPOA_MAX_MATRIX_DIMENSION]);
            //printf("node %c, score %d\n", nodes[node_id], scores[(graph_pos+1) * CUDAPOA_MAX_MATRIX_DIMENSION]);
        }

        //score_prev_i[0] = 0;
        // Init horizonal boundary conditions (read).
        for(uint16_t j = 1; j < read_count + 1; j++)
        {
            //score_prev_i[j] = j * GAP;
            scores[j] = j * GAP;
            //printf("%d ", scores[j]);
        }
        //printf("\n");

        init = clock64() - start;

    }

    __syncthreads();

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

    start = clock64();

    int16_t max_score = SHRT_MIN;
    int16_t max_i = -1;
    int16_t max_j = -1;

    long long int serial = 0;

    // Run DP loop for calculating scores. Process each row at a time, and
    // compute vertical and diagonal values in parallel.
    for(uint16_t graph_pos = 0; graph_pos < graph_count; graph_pos++)
    {
        uint16_t node_id = graph[graph_pos];
        uint16_t i = graph_pos + 1;
        //int32_t* scores_i = &scores[i * CUDAPOA_MAX_MATRIX_DIMENSION];
                int16_t init_score = scores[i * CUDAPOA_MAX_MATRIX_DIMENSION];
                uint16_t out_edge_count = outgoing_edge_count_global[node_id];

        uint16_t pred_count = incoming_edge_count[node_id];

        uint16_t pred_i_1 = (pred_count == 0 ? 0 :
                node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES]] + 1);
                //node_id_to_pos[local_incoming_edges[0]] + 1);
        int16_t* scores_pred_i_1 = &scores[pred_i_1 * CUDAPOA_MAX_MATRIX_DIMENSION];

        //if (pred_i_1 == i - 1)
        //{
        //    //printf("using optim 1\n");
        //    scores_pred_i_1 = &score_prev_i[0];
        //}
        //else
        //{
        //    scores_pred_i_1 = &scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION];
        //}

        uint8_t n = nodes[node_id];

        int16_t prev_score = init_score;
        uint16_t max_cols = max(read_count, (((read_count - 1) / blockDim.x) + 1) * blockDim.x);
        for(uint16_t read_pos = thread_idx; read_pos < max_cols; read_pos += blockDim.x)
        {
            if (read_pos < read_count)
            {
                //printf("updating vertical for pos %d thread %d\n", read_pos, thread_idx);
                int32_t char_profile = (n == read[read_pos] ? MATCH : MISMATCH);
                // Index into score matrix.
                uint16_t j = read_pos + 1;
                int16_t score = max(scores_pred_i_1[j-1] + char_profile,
                        scores_pred_i_1[j] + GAP);

                // Perform same score updates as above, but for rest of predecessors.
                for (uint16_t p = 1; p < pred_count; p++)
                {
                    int16_t pred_i_2 = node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + p]] + 1;
                    int16_t* scores_pred_i_2 = &scores[pred_i_2 * CUDAPOA_MAX_MATRIX_DIMENSION];

                    //int16_t pred_i_2 = node_id_to_pos[local_incoming_edges[p]] + 1;
                    //if (pred_i_2 == i - 1)
                    //{
                    //    scores_pred_i_2 = score_prev_i;
                    //}
                    //else
                    //{
                    //    scores_pred_i_2 = &scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION];
                    //}
                    //scores_pred_i = &scores[pred_i * CUDAPOA_MAX_MATRIX_DIMENSION];

                    score = max(scores_pred_i_2[j - 1] + char_profile,
                            max(score, scores_pred_i_2[j] + GAP));
                }

                score_i[j] = score;
            }
            //printf("%d \n", j);

            __syncthreads();

            // Calculate horizontal max scores in single thread, using data in shared memory only.
            if (thread_idx == 0)
            {
                long long int temp = clock64();
                //int16_t init_score = scores[i * CUDAPOA_MAX_MATRIX_DIMENSION];
                //score_i[0] = init_score;
                //int16_t prev_score = init_score;
                //score_prev_i[0] = init_score;

                // Perform score updates for horizontal moves.
                //printf("horizontal update for %d to %d\n", read_pos, min(read_pos + blockDim.x, read_count));
                uint16_t pos;
                for(pos = read_pos; pos < min(read_pos + blockDim.x, read_count); pos++)
                {
                    // Index into score matrix.
                    uint16_t jpos = pos + 1;
                    int16_t score = score_i[jpos];

                    //score = max(score_i[j-1] + GAP, score);
                    score = max(prev_score + GAP, score);

                    //if (j == read_count)
                    //{
                    //    printf("oe %d ie %d for pos %d id %d\n", 
                    //            outgoing_edge_count[node_id], 
                    //            incoming_edge_count[node_id], 
                    //            i, node_id);
                    //}

                    prev_score = score;
                    score_i[jpos] = score;
                    //score_prev_i[j] = score;
                }

                // Once last column of read is reached, update max score
                // and positions. So process last column data separately.
                //uint16_t out_edge_count = outgoing_edge_count_global[node_id];
                //uint16_t j = read_count;
                //int16_t score = score_i[j];
                ////score = max(score_i[j-1] + GAP, score);
                //score = max(prev_score + GAP, score);
                if (pos == read_count && out_edge_count == 0)
                {
                    if (max_score < prev_score)
                    {
                        max_score = prev_score;
                        max_i = i;
                        max_j = pos;
                        //printf("updating max\n");
                    }
                }
                //score_i[j] = score;
                //printf("\n");

                //uint16_t next_id = graph[min(graph_pos + 1, graph_count - 1)];
                //for(uint16_t e = 0; e < incoming_edge_count[next_id]; e++)
                //{
                //    local_incoming_edges[e] = incoming_edges[next_id * CUDAPOA_MAX_NODE_EDGES + e];
                //}
                serial += (clock64() - temp);
                //printf("horizontal done\n");
            }
            //__syncthreads();
        }

        __syncthreads();

        //printf("writing scores\n");

        // Write scores back to global memory.
        for(uint32_t read_pos = thread_idx; read_pos < read_count; read_pos += blockDim.x)
        {
            uint16_t j = read_pos + 1;
            scores[i * CUDAPOA_MAX_MATRIX_DIMENSION + j] = score_i[j];
        }
        __syncthreads();

    }

//    if (thread_idx == 0)
//    {
//        for(uint32_t i = 0; i < graph_count + 1; i++)
//        {
//            for(uint32_t j = 0; j < read_count + 1; j++)
//            {
//                printf("%05d\n", scores[i * CUDAPOA_MAX_MATRIX_DIMENSION + j]);
//            }
//        }
//    }

    long long int nw = clock64() - start;
    long long int tb = 0;

    int16_t i = max_i;
    int16_t j = max_j;

    start = clock64();

    uint16_t aligned_nodes = 0;
    if (thread_idx == 0)
    {
        // Fill in backtrace

        int16_t prev_i = 0;
        int16_t prev_j = 0;

        //printf("maxi %d maxj %d score %d\n", i, j, scores[i * CUDAPOA_MAX_MATRIX_DIMENSION + j]);

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
            // printf("%d %d node %d\n", i, j, graph[i-1]);

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

            //printf("loop %d %d\n",i, j);
        }
        //printf("aligned nodes %d\n", aligned_nodes);

        tb = clock64() - start;
    }

    //if (thread_idx == 0)
    //{
    //    long long int total = init + nw + tb;
    //    printf("Total time of init is %lf %\n", ((double)init / total) * 100.f);
    //    printf("Total time of serial is %lf %\n", ((double)serial / total) * 100.f);
    //    printf("Total time of nw is %lf %\n", ((double)(nw - serial) / total) * 100.f);
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

