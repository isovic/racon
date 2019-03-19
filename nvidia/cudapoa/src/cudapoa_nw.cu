
#include "cudapoa_kernels.cuh"
#include <stdio.h>

#define WARP_SIZE 32

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
#pragma message("TODO: Pass match/gap/mismatch scores into NW kernel as parameters.")
    const int16_t GAP = -8;
    const int16_t MISMATCH = -6;
    const int16_t MATCH = 8;

    //__shared__ uint16_t local_incoming_edges[CUDAPOA_MAX_NODE_EDGES];
    //__shared__ uint16_t node_id_to_pos[2048];
    //__shared__ uint16_t incoming_edge_count[1024];
    //__shared__ uint16_t outgoing_edge_count[1024];

    //__shared__ int16_t score_i[1024];
    //__shared__ int16_t score_prev_i[1024];

    __shared__ int16_t prev_score[1];

    uint32_t thread_idx = threadIdx.x;

    long long int start = clock64();
    long long int init = 0;

    //for(uint16_t graph_pos = thread_idx; graph_pos < graph_count; graph_pos += blockDim.x)
    //{
    //    //node_id_to_pos[graph_pos] = node_id_to_pos_global[graph_pos];
    //    //incoming_edge_count[graph_pos] = incoming_edge_count_global[graph_pos];
    //}

    // Init horizonal boundary conditions (read).
    for(uint16_t j = thread_idx + 1; j < read_count + 1; j += blockDim.x)
    {
        //score_prev_i[j] = j * GAP;
        scores[j] = j * GAP;
    }


    if (thread_idx == 0)
    {
#ifdef DEBUG
        printf("graph %d, read %d\n", graph_count, read_count);
#endif

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
                scores[i * CUDAPOA_MAX_SEQUENCE_SIZE] = GAP;
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
                    penalty = max(penalty, scores[pred_node_graph_pos * CUDAPOA_MAX_SEQUENCE_SIZE]);
                }
                scores[i * CUDAPOA_MAX_SEQUENCE_SIZE] = penalty + GAP;
            }

            //printf("%d \n", scores[i * CUDAPOA_MAX_MATRIX_DIMENSION]);
            //printf("node %c, score %d\n", nodes[node_id], scores[(graph_pos+1) * CUDAPOA_MAX_MATRIX_DIMENSION]);
        }

        //score_prev_i[0] = 0;
        //for(uint16_t j = 1; j < read_count + 1; j++)
        //{
        //    //printf("%d ", scores[j]);
        //}
        //printf("\n");

        init = clock64() - start;

    }


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

    long long int serial = 0;

    // Run DP loop for calculating scores. Process each row at a time, and
    // compute vertical and diagonal values in parallel.
    for(uint16_t graph_pos = 0; graph_pos < graph_count; graph_pos++)
    {
        uint16_t node_id = graph[graph_pos];
        uint16_t i = graph_pos + 1;
        //int32_t* scores_i = &scores[i * CUDAPOA_MAX_MATRIX_DIMENSION];

        int16_t init_score = scores[i * CUDAPOA_MAX_SEQUENCE_SIZE];
        uint16_t out_edge_count = outgoing_edge_count_global[node_id];

        uint16_t pred_count = incoming_edge_count[node_id];

        uint16_t pred_i_1 = (pred_count == 0 ? 0 :
                node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES]] + 1);
                //node_id_to_pos[local_incoming_edges[0]] + 1);
        int16_t* scores_pred_i_1 = &scores[pred_i_1 * CUDAPOA_MAX_SEQUENCE_SIZE];

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

        //int16_t prev_score = init_score;
        if (thread_idx == 0)
        {
            prev_score[0] = init_score;
        }

        //__syncthreads();

        // max_cols is the first warp boundary multiple beyond read_count. This is done
        // so all threads in the warp enter the loop.
        uint16_t max_cols = max(read_count, (((read_count - 1) / (blockDim.x * 4)) + 1) * (blockDim.x * 4));
        //if (thread_idx == 0)
        //{
        //    printf("max cols %d\n", max_cols);
        //}
        for(uint16_t read_pos = thread_idx * 4; read_pos < max_cols; read_pos += blockDim.x * 4)
        {
            //printf("updating vertical for pos %d thread %d\n", read_pos, thread_idx);
            int16_t char_profile0 = (n == read[read_pos + 0] ? MATCH : MISMATCH);
            int16_t char_profile1 = (n == read[read_pos + 1] ? MATCH : MISMATCH);
            int16_t char_profile2 = (n == read[read_pos + 2] ? MATCH : MISMATCH);
            int16_t char_profile3 = (n == read[read_pos + 3] ? MATCH : MISMATCH);
            // Index into score matrix.
            uint16_t j0 = read_pos + 1;
            uint16_t j1 = read_pos + 2;
            uint16_t j2 = read_pos + 3;
            uint16_t j3 = read_pos + 4;
            //printf("thread idx %d locations %d %d %d %d\n", thread_idx, j0, j1, j2, j3);
            int16_t score0 = max(scores_pred_i_1[j0-1] + char_profile0,
                    scores_pred_i_1[j0] + GAP);
            int16_t score1 = max(scores_pred_i_1[j1-1] + char_profile1,
                    scores_pred_i_1[j1] + GAP);
            int16_t score2 = max(scores_pred_i_1[j2-1] + char_profile2,
                    scores_pred_i_1[j2] + GAP);
            int16_t score3 = max(scores_pred_i_1[j3-1] + char_profile3,
                    scores_pred_i_1[j3] + GAP);

            // Perform same score updates as above, but for rest of predecessors.
            for (uint16_t p = 1; p < pred_count; p++)
            {
                int16_t pred_i_2 = node_id_to_pos[incoming_edges[node_id * CUDAPOA_MAX_NODE_EDGES + p]] + 1;
                int16_t* scores_pred_i_2 = &scores[pred_i_2 * CUDAPOA_MAX_SEQUENCE_SIZE];

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

                score0 = max(scores_pred_i_2[j0 - 1] + char_profile0,
                        max(score0, scores_pred_i_2[j0] + GAP));
                score1 = max(scores_pred_i_2[j1 - 1] + char_profile1,
                        max(score1, scores_pred_i_2[j1] + GAP));
                score2 = max(scores_pred_i_2[j2 - 1] + char_profile2,
                        max(score2, scores_pred_i_2[j2] + GAP));
                score3 = max(scores_pred_i_2[j3 - 1] + char_profile3,
                        max(score3, scores_pred_i_2[j3] + GAP));
            }

            //__syncthreads();

            long long int temp = clock64();

            //const uint32_t WARP_SIZE = 32;
            //uint8_t num_warps = ((blockDim.x - 1) / WARP_SIZE) + 1; // Assuming warp size of 32.
            //if (thread_idx < 32)
            //for(uint8_t warp = 0; warp < num_warps; warp++)
            //{
            // For the first thread, calculate the score from score of the last thread in the
            // last warp. If j == 1, then prev_score is initialized from score matrix.
            //uint32_t first_warp_thread_idx = warp * WARP_SIZE;
            //uint32_t last_warp_thread_idx = (warp + 1) * WARP_SIZE - 1;

            if (thread_idx == 0)
            {
                //printf("last score for thread 0 is %d\n", init_score);
                //score0 = max(prev_score[0] + GAP, score);
                score0 = max(init_score + GAP, score0);
                score1 = max(score0 + GAP, score1);
                score2 = max(score1 + GAP, score2);
                score3 = max(score2 + GAP, score3);
                //printf("score for thread %d location %d is %d with prev_score %d\n", thread_idx, j, score, prev_score);
            }

            // While there are changes to the horizontal score values, keep updating the matrix.
            // So loop will only run the number of time there are corrections in the matrix.
            // The any_sync warp primitive lets us easily check if any of the threads had an update.
            //bool loop = true;
            uint32_t loop = 0xffffffff >> 1;
            while(loop)
            {
                //loop = false;
                // The shfl_up lets us grab a value from the lane below.
                //printf("thread %d line %d\n", thread_idx, __LINE__);
                int16_t last_score = __shfl_up_sync(0xffffffff << 1, score3, 1);
                //printf("thread %d line %d\n", thread_idx, __LINE__);
                score0 = max(last_score + GAP, score0);
                score1 = max(score0 + GAP, score1);
                score2 = max(score1 + GAP, score2);
                score3 = max(score2 + GAP, score3);
                //if (thread_idx == 1)
                //{
                //    printf("thread idx %d scores %d, %d, %d, %d\n", thread_idx, score0, score1, score2, score3);
                //}
                //if (thread_idx == 2)
                //{
                //    printf("thread idx %d scores %d, %d, %d, %d\n", thread_idx, score0, score1, score2, score3);
                //}
                //int16_t new_score = max(__shfl_up_sync(0xffffffff << 1, score, 1) + GAP, score);
                //if (new_score > score)
                //{
                //    loop = true;
                //    score = new_score;
                //}
                loop = loop >> 1;
                //printf("thread %d line %d\n", thread_idx, __LINE__);
            }

            //printf("score after shuffle operations for thread %d location %d is %d\n", thread_idx, j, score);

            // Move the score of the last thread in the warp to a variable
            // into the first thread, so in the next warp the the first thread has a valid
            // prev value.
            //prev_score = __shfl_down_sync(0x1, score, 31);//score;

            //if (thread_idx == (WARP_SIZE - 1))
            //{
            //    prev_score[0] = score3;
            //}

                //printf("thread %d line %d\n", thread_idx, __LINE__);
            init_score = __shfl_sync(0xffffffff, score3, 31);
                //printf("thread %d line %d\n", thread_idx, __LINE__);

            if (j0 == read_count)
            {
                printf("last score is %d\n", score0);
            }
            if (j1 == read_count)
            {
                printf("last score is %d\n", score1);
            }
            if (j2 == read_count)
            {
                printf("last score is %d\n", score2);
            }
            if (j3 == read_count)
            {
                printf("last score is %d\n", score3);
            }

            //if (thread_idx == 0)
            //{
            //    printf("previous score for thread %d is %d\n", thread_idx, prev_score);
            //}

            // Index into score matrix.
            scores[i * CUDAPOA_MAX_SEQUENCE_SIZE + j0] = score0;
            scores[i * CUDAPOA_MAX_SEQUENCE_SIZE + j1] = score1;
            scores[i * CUDAPOA_MAX_SEQUENCE_SIZE + j2] = score2;
            scores[i * CUDAPOA_MAX_SEQUENCE_SIZE + j3] = score3;
            serial += (clock64() - temp);

            __syncthreads();
        }

        //__syncthreads();

    }

    //__syncthreads();

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

    start = clock64();

    uint16_t aligned_nodes = 0;
    if (thread_idx == 0)
    {
        // Find location of the maximum score in the matrix.
        int16_t i = 0;
        int16_t j = read_count;
        int16_t mscore = SHRT_MIN;

        for (int16_t idx = 1; idx <= graph_count; idx++)
        {
            if (outgoing_edge_count_global[graph[idx - 1]] == 0)
            {
                int16_t s = scores[idx * CUDAPOA_MAX_SEQUENCE_SIZE + j];
                if (mscore < s)
                {
                    mscore = s;
                    i = idx;
                }
            }
        }
        // Fill in backtrace

        int16_t prev_i = 0;
        int16_t prev_j = 0;

        //printf("maxi %d maxj %d score %d\n", i, j, scores[i * CUDAPOA_MAX_SEQUENCE_SIZE + j]);

        // Trace back from maximum score position to generate alignment.
        // Trace back is done by re-calculating the score at each cell
        // along the path to see which preceding cell the move could have
        // come from. This seems computaitonally more expensive, but doesn't
        // require storing any traceback buffer during alignment.
        while(!(i == 0 && j == 0))
        {
            //printf("%d %d\n", i, j);
            int16_t scores_ij = scores[i * CUDAPOA_MAX_SEQUENCE_SIZE + j];
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
                if (scores_ij == (scores[pred_i * CUDAPOA_MAX_SEQUENCE_SIZE + (j - 1)] + match_cost))
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

                        if (scores_ij == (scores[pred_i * CUDAPOA_MAX_SEQUENCE_SIZE + (j - 1)] + match_cost))
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

                if (scores_ij == scores[pred_i * CUDAPOA_MAX_SEQUENCE_SIZE + j] + GAP)
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

                        if (scores_ij == scores[pred_i * CUDAPOA_MAX_SEQUENCE_SIZE + j] + GAP)
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
            if (!pred_found && scores_ij == scores[i * CUDAPOA_MAX_SEQUENCE_SIZE + (j - 1)] + GAP)
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
#ifdef DEBUG
        printf("aligned nodes %d\n", aligned_nodes);
#endif

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

