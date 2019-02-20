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
    uint8_t sorted[CUDAPOA_MAX_NODES_PER_WINDOW];
    nodes[0] = 'A';
    nodes[1] = 'B';
    nodes[2] = 'C';
    nodes[3] = 'D';
    nodes[4] = 'E';
    nodes[5] = 'F';

    cudaError_t error = cudaSuccess;

    // Allocate device buffer for results.
    uint8_t* sorted_d;
    error = cudaMalloc((void**) &sorted_d,
                       sizeof(uint8_t) * CUDAPOA_MAX_NODES_PER_WINDOW);
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
    uint16_t outoing_edges[CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES];
    outoing_edges[0 * CUDAPOA_MAX_NODE_EDGES + 0] = 5;
    outoing_edges[0 * CUDAPOA_MAX_NODE_EDGES + 1] = 1;
    outoing_edges[1 * CUDAPOA_MAX_NODE_EDGES + 0] = 3;
    outoing_edges[2 * CUDAPOA_MAX_NODE_EDGES + 0] = 1;
    outoing_edges[4 * CUDAPOA_MAX_NODE_EDGES + 0] = 0;
    uint16_t* outgoing_edges_d;
    error = cudaMalloc((void**) &outgoing_edges_d, 
                       sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES);
    cudaCheckError(error, "", __LINE__, __FILE__);
    error = cudaMemcpy(outgoing_edges_d, outoing_edges,
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
    nvidia::cudapoa::topologicalSort(sorted_d, num_nodes, incoming_edge_count_d,
                    outgoing_edges_d, outgoing_edge_count_d,
                    1, 1, 0);
    error = cudaDeviceSynchronize();
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Copy results back to memory and print them.
    error = cudaMemcpy(sorted, sorted_d, sizeof(uint8_t) * num_nodes,
                       cudaMemcpyDeviceToHost);
    cudaCheckError(error, "", __LINE__, __FILE__);

    // Verify results.
    assert (sorted[0] == 2); // C
    assert (sorted[1] == 4); // E
    assert (sorted[2] == 0); // A
    assert (sorted[3] == 5); // F
    assert (sorted[4] == 1); // B
    assert (sorted[5] == 3); // D

    // Uncomment to print final sorted list.
    //for(uint32_t i = 0; i < num_nodes; i++)
    //{
    //    std::cout << nodes[sorted[i]] << ", " << std::endl;
    //}
}

int main()
{
    testTopologicalSort_1();
    std::cout << "All tests passed successfully." << std::endl;
    return 0;
}
