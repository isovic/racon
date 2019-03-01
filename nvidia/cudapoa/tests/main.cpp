// Writing some basic tests for CUDA kernels.

#include <stdio.h>
#include <iostream>
#include <assert.h>

#include "cudapoa_kernels.cuh"

void cudaCheckError(cudaError_t error, const char* msg, uint32_t line, const char* file)
{
    if (error != cudaSuccess)
    {
        std::cerr << msg << " (" << cudaGetErrorString(error) << ") " << \
            "on " << file << ":" << line << std::endl;
        exit(-1);
    }
}

void testTopologicalSort_1()
{
    std::cout << "Running " << __func__ << std::endl;
    /*
             |----->F
             |
             |
       E---->A----->B----->D
                    ^
                    |
             C------|
    */
    uint32_t num_nodes = 6;
    uint8_t nodes[CUDAPOA_MAX_NODES_PER_WINDOW];
    uint16_t sorted[CUDAPOA_MAX_NODES_PER_WINDOW];
    uint16_t sorted_node_map[CUDAPOA_MAX_NODES_PER_WINDOW];
    nodes[0] = 'A';
    nodes[1] = 'B';
    nodes[2] = 'C';
    nodes[3] = 'D';
    nodes[4] = 'E';
    nodes[5] = 'F';

    cudaError_t error = cudaSuccess;

    // Allocate device buffer for results.
    uint16_t* sorted_d;
    error = cudaMalloc((void**) &sorted_d,
                       sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
    cudaCheckError(error, "", __LINE__, __FILE__);
    uint16_t* sorted_node_map_d;
    error = cudaMalloc((void**) &sorted_node_map_d,
                       sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Allocate host and device buufer for incoming edge counts.
    uint16_t incoming_edge_count[CUDAPOA_MAX_NODES_PER_WINDOW];
    incoming_edge_count[0] = 1; 
    incoming_edge_count[1] = 2;
    incoming_edge_count[2] = 0;
    incoming_edge_count[3] = 1;
    incoming_edge_count[4] = 0;
    incoming_edge_count[5] = 1;
    uint16_t* incoming_edge_count_d;
    error = cudaMalloc((void**) &incoming_edge_count_d,
                        sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(incoming_edge_count_d, incoming_edge_count,
               sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW,
               cudaMemcpyHostToDevice);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Allocate device and host buffer for outgoing edges.
    uint16_t outgoing_edges[CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
    outgoing_edges[0 * CUDAPOA_MAX_NODE_EDGES + 0] = 5;
    outgoing_edges[0 * CUDAPOA_MAX_NODE_EDGES + 1] = 1;
    outgoing_edges[1 * CUDAPOA_MAX_NODE_EDGES + 0] = 3;
    outgoing_edges[2 * CUDAPOA_MAX_NODE_EDGES + 0] = 1;
    outgoing_edges[4 * CUDAPOA_MAX_NODE_EDGES + 0] = 0;
    uint16_t* outgoing_edges_d;
    error = cudaMalloc((void**) &outgoing_edges_d, 
                       sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(outgoing_edges_d, outgoing_edges,
               sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES,
               cudaMemcpyHostToDevice);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Allocate device and host buffer for outgoing edge count.
    uint16_t outgoing_edge_count[CUDAPOA_MAX_NODES_PER_WINDOW];
    outgoing_edge_count[0] = 2;
    outgoing_edge_count[1] = 1;
    outgoing_edge_count[2] = 1;
    outgoing_edge_count[4] = 1;
    uint16_t* outgoing_edge_count_d;
    error = cudaMalloc((void**) &outgoing_edge_count_d,
                        sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(outgoing_edge_count_d, outgoing_edge_count,
               sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW,
               cudaMemcpyHostToDevice);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Run topological sort.
    nvidia::cudapoa::topologicalSort(sorted_d, sorted_node_map_d, num_nodes, incoming_edge_count_d,
                    outgoing_edges_d, outgoing_edge_count_d,
                    1, 1, 0);
    error = cudaDeviceSynchronize();
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Copy results back to memory and print them.
    error = cudaMemcpy(sorted, sorted_d, sizeof(uint16_t) * num_nodes,
                       cudaMemcpyDeviceToHost);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(sorted_node_map, sorted_node_map_d, sizeof(uint16_t) * num_nodes,
                       cudaMemcpyDeviceToHost);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Verify results.
    assert (sorted[0] == 2); // C
    assert (sorted[1] == 4); // E
    assert (sorted[2] == 0); // A
    assert (sorted[3] == 5); // F
    assert (sorted[4] == 1); // B
    assert (sorted[5] == 3); // D

    assert (sorted_node_map[0] == 2); // C
    assert (sorted_node_map[1] == 4); // E
    assert (sorted_node_map[2] == 0); // A
    assert (sorted_node_map[3] == 5); // F
    assert (sorted_node_map[4] == 1); // B
    assert (sorted_node_map[5] == 3); // D

    // Uncomment to print final sorted list.
    for(uint32_t i = 0; i < num_nodes; i++)
    {
        std::cout << nodes[sorted[i]] << ", " << std::endl;
    }
}

void testNW_1()
{
    std::cout << "Running " << __func__ << std::endl;
    /*
             |----->F------|
             |             |
             |             \/
       E---->A----->B----->D
                    ^
                    |
             C------|
    */
    uint32_t num_nodes = 6;
    uint8_t nodes[CUDAPOA_MAX_NODES_PER_WINDOW];
    uint16_t sorted[CUDAPOA_MAX_NODES_PER_WINDOW];
    uint16_t sorted_node_map[CUDAPOA_MAX_NODES_PER_WINDOW];
    nodes[0] = 'A';
    nodes[1] = 'B';
    nodes[2] = 'C';
    nodes[3] = 'D';
    nodes[4] = 'E';
    nodes[5] = 'F';

    cudaError_t error = cudaSuccess;
    uint8_t* nodes_d;
    error = cudaMalloc((void**) &nodes_d,
                       sizeof(uint8_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
    cudaCheckError(error, "", __LINE__, __FILE__);

    error = cudaMemcpy(nodes_d, nodes, sizeof(uint8_t) * CUDAPOA_MAX_NODES_PER_WINDOW,
            cudaMemcpyHostToDevice);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Allocate device buffer for results.
    uint16_t* sorted_d;
    error = cudaMalloc((void**) &sorted_d,
                       sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
    cudaCheckError(error, "", __LINE__, __FILE__);
    uint16_t* sorted_node_map_d;
    error = cudaMalloc((void**) &sorted_node_map_d,
                       sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Allocate device and host buffer for incoming edges.
    uint16_t incoming_edges[CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
    incoming_edges[0 * CUDAPOA_MAX_NODE_EDGES + 0] = 4;
    incoming_edges[1 * CUDAPOA_MAX_NODE_EDGES + 0] = 0;
    incoming_edges[1 * CUDAPOA_MAX_NODE_EDGES + 1] = 2;
    incoming_edges[3 * CUDAPOA_MAX_NODE_EDGES + 0] = 1;
    incoming_edges[3 * CUDAPOA_MAX_NODE_EDGES + 1] = 5;
    incoming_edges[5 * CUDAPOA_MAX_NODE_EDGES + 0] = 0;
    uint16_t* incoming_edges_d;
    error = cudaMalloc((void**) &incoming_edges_d, 
                       sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(incoming_edges_d, incoming_edges,
               sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES,
               cudaMemcpyHostToDevice);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Allocate host and device buufer for incoming edge counts.
    uint16_t incoming_edge_count[CUDAPOA_MAX_NODES_PER_WINDOW];
    incoming_edge_count[0] = 1; 
    incoming_edge_count[1] = 2;
    incoming_edge_count[2] = 0;
    incoming_edge_count[3] = 2;
    incoming_edge_count[4] = 0;
    incoming_edge_count[5] = 1;
    uint16_t* incoming_edge_count_d;
    error = cudaMalloc((void**) &incoming_edge_count_d,
                        sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(incoming_edge_count_d, incoming_edge_count,
               sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW,
               cudaMemcpyHostToDevice);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Allocate device and host buffer for outgoing edges.
    uint16_t outgoing_edges[CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
    outgoing_edges[0 * CUDAPOA_MAX_NODE_EDGES + 0] = 5;
    outgoing_edges[0 * CUDAPOA_MAX_NODE_EDGES + 1] = 1;
    outgoing_edges[1 * CUDAPOA_MAX_NODE_EDGES + 0] = 3;
    outgoing_edges[2 * CUDAPOA_MAX_NODE_EDGES + 0] = 1;
    outgoing_edges[4 * CUDAPOA_MAX_NODE_EDGES + 0] = 0;
    outgoing_edges[5 * CUDAPOA_MAX_NODE_EDGES + 0] = 3;
    uint16_t* outgoing_edges_d;
    error = cudaMalloc((void**) &outgoing_edges_d, 
                       sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(outgoing_edges_d, outgoing_edges,
               sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES,
               cudaMemcpyHostToDevice);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Allocate device and host buffer for outgoing edge count.
    uint16_t outgoing_edge_count[CUDAPOA_MAX_NODES_PER_WINDOW];
    outgoing_edge_count[0] = 2;
    outgoing_edge_count[1] = 1;
    outgoing_edge_count[2] = 1;
    outgoing_edge_count[4] = 1;
    outgoing_edge_count[5] = 1;
    uint16_t* outgoing_edge_count_d;
    error = cudaMalloc((void**) &outgoing_edge_count_d,
                        sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(outgoing_edge_count_d, outgoing_edge_count,
               sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW,
               cudaMemcpyHostToDevice);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Run topological sort.
    nvidia::cudapoa::topologicalSort(sorted_d, sorted_node_map_d, num_nodes, incoming_edge_count_d,
                    outgoing_edges_d, outgoing_edge_count_d,
                    1, 1, 0);
    error = cudaDeviceSynchronize();
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Copy results back to memory and print them.
    error = cudaMemcpy(sorted, sorted_d, sizeof(uint16_t) * num_nodes,
                       cudaMemcpyDeviceToHost);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(sorted_node_map, sorted_node_map_d, sizeof(uint16_t) * num_nodes,
                       cudaMemcpyDeviceToHost);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Verify results.
    assert (sorted[0] == 2); // C
    assert (sorted[1] == 4); // E
    assert (sorted[2] == 0); // A
    assert (sorted[3] == 5); // F
    assert (sorted[4] == 1); // B
    assert (sorted[5] == 3); // D

    assert (sorted_node_map[0] == 2); // C
    assert (sorted_node_map[1] == 4); // E
    assert (sorted_node_map[2] == 0); // A
    assert (sorted_node_map[3] == 5); // F
    assert (sorted_node_map[4] == 1); // B
    assert (sorted_node_map[5] == 3); // D

    // Uncomment to print final sorted list.
    for(uint32_t i = 0; i < num_nodes; i++)
    {
        std::cout << nodes[sorted[i]] << ", " << std::endl;
    }

    uint8_t read[4];
    read[0] = 'E';
    read[1] = 'C';
    read[2] = 'B';
    read[3] = 'D';
    uint8_t *read_d;
    cudaMalloc((void**) &read_d, sizeof(uint8_t) * 4);
    cudaMemcpy(read_d, read, 4, cudaMemcpyHostToDevice);

    int32_t* scores_d;
    int16_t *traceback_i_d, *traceback_j_d;
    cudaMalloc((void**) &scores_d, sizeof(int32_t) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1));
    cudaMalloc((void**) &traceback_i_d, sizeof(int16_t) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1));
    cudaMalloc((void**) &traceback_j_d, sizeof(int16_t) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1));

    // Run needleman wunsch.
    nvidia::cudapoa::needlemanWunsch(nodes_d, sorted_d, sorted_node_map_d, num_nodes, 
                                     incoming_edge_count_d,incoming_edges_d, 
                                     outgoing_edge_count_d,outgoing_edges_d, 
                                     read_d, 4,
                                     scores_d, traceback_i_d, traceback_j_d,
                                     1, 1, 0);

    error = cudaDeviceSynchronize();
    int16_t traceback_i[(CUDAPOA_MAX_NODES_PER_WINDOW + 1) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1)];
    int16_t traceback_j[(CUDAPOA_MAX_NODES_PER_WINDOW + 1) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1)];
    int32_t scores[(CUDAPOA_MAX_NODES_PER_WINDOW + 1) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1)];

    cudaMemcpy(traceback_i, traceback_i_d, sizeof(int16_t) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1), cudaMemcpyDeviceToHost);
    cudaMemcpy(traceback_j, traceback_j_d, sizeof(int16_t) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1), cudaMemcpyDeviceToHost);
    cudaMemcpy(scores, scores_d, sizeof(int32_t) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) * (CUDAPOA_MAX_NODES_PER_WINDOW + 1), cudaMemcpyDeviceToHost);

    for(int i = 0; i < num_nodes + 1; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            //printf("(%d, %d, %d) ", traceback_i[i * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) + j],
            //                        traceback_j[i * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) + j],
            //                        scores[i * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) + j]);
            printf("(%d, %d) ", traceback_i[i * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) + j],
                                    traceback_j[i * (CUDAPOA_MAX_NODES_PER_WINDOW + 1) + j]);
        }
        printf("\n");
    }

}

int main()
{
    //testTopologicalSort_1();
    testNW_1();
    std::cout << "All tests passed successfully." << std::endl;
    return 0;
}
