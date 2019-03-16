/*!
 * @file cudabatch.cpp
 *
 * @brief CUDABatch class source file
 */

#include <string>
#include <iostream>
#include <cstring>

#include "cudapoa_kernels.cuh"

#include "cudabatch.hpp"
#include "cudautils.hpp"

#ifndef TABS
#define TABS printTabs(bid_)
#endif

namespace racon {

uint32_t CUDABatchProcessor::batches = 0;

inline std::string printTabs(uint32_t tab_count)
{
    std::string s;
    for(uint32_t i = 0; i < tab_count; i++)
    {
        s += "\t";
    }
    return s;
}

std::unique_ptr<CUDABatchProcessor> createCUDABatch(uint32_t max_windows, uint32_t max_window_depth)
{
    return std::unique_ptr<CUDABatchProcessor>(new CUDABatchProcessor(max_windows, max_window_depth));
}

CUDABatchProcessor::CUDABatchProcessor(uint32_t max_windows, uint32_t max_window_depth)
    : max_windows_(max_windows)
    , max_depth_per_window_(max_window_depth)
    , windows_()
    , sequence_count_()
    , stream_()
    , consensus_pitch_()
{
    bid_ = CUDABatchProcessor::batches++;

    NUM_BLOCKS = max_windows_;

    // Verify that thread block size is in multiples of warp size.
    if (NUM_THREADS % 32 != 0)
    {
        std::cerr << "Thread block size needs to be in multiples of 32." << std::endl;
        exit(-1);
    }

    // Create new CUDA stream.
    CU_CHECK_ERR(cudaStreamCreate(&stream_));

    // Allocate host memory and CUDA memory based on max sequence and target counts.

    uint32_t input_size = max_windows_ * max_depth_per_window_ * CUDAPOA_MAX_SEQUENCE_SIZE; //TODO how big does this need to be

    // Input buffers.
    CU_CHECK_ERR(cudaHostAlloc((void**) &inputs_h_, input_size * sizeof(uint8_t),
                  cudaHostAllocDefault));
    CU_CHECK_ERR(cudaHostAlloc((void**) &num_sequences_per_window_h_, max_windows * sizeof(uint16_t),
            cudaHostAllocDefault));
    CU_CHECK_ERR(cudaHostAlloc((void**) &sequence_lengths_h_, max_windows * max_depth_per_window_* sizeof(uint16_t),
            cudaHostAllocDefault));
    CU_CHECK_ERR(cudaHostAlloc((void**) &window_details_h_, max_windows * sizeof(nvidia::cudapoa::WindowDetails),
            cudaHostAllocDefault));

    //device allocations
    CU_CHECK_ERR(cudaMalloc((void**)&inputs_d_, input_size * sizeof(uint8_t)));
    CU_CHECK_ERR(cudaMalloc((void**)&sequence_lengths_d_, max_windows * max_depth_per_window_ * sizeof(uint16_t)));
    CU_CHECK_ERR(cudaMalloc((void**)&window_details_d_, max_windows * sizeof(nvidia::cudapoa::WindowDetails)));

    std::cout << TABS << bid_ << " Allocated input buffers of size " << (static_cast<float>(input_size)  / (1024 * 1024)) << "MB" << std::endl;

    // Output buffers.
    input_size = max_windows_ * CUDAPOA_MAX_SEQUENCE_SIZE;
    CU_CHECK_ERR(cudaHostAlloc((void**) &consensus_h_, input_size * sizeof(uint8_t),
                  cudaHostAllocDefault));

    CU_CHECK_ERR(cudaMallocPitch((void**) &consensus_d_,
                    &consensus_pitch_,
                    sizeof(uint8_t) * CUDAPOA_MAX_SEQUENCE_SIZE,
                    max_windows_));
    std::cout << TABS << bid_ << " Allocated output buffers of size " << (static_cast<float>(input_size)  / (1024 * 1024)) << "MB" << std::endl;

    // Buffers for storing NW scores and backtrace.
    CU_CHECK_ERR(cudaMalloc((void**) &scores_d_, sizeof(int16_t) * CUDAPOA_MAX_MATRIX_DIMENSION * CUDAPOA_MAX_SEQUENCE_SIZE * NUM_BLOCKS));
    CU_CHECK_ERR(cudaMalloc((void**) &traceback_i_d_, sizeof(int16_t) * CUDAPOA_MAX_MATRIX_DIMENSION * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &traceback_j_d_, sizeof(int16_t) * CUDAPOA_MAX_MATRIX_DIMENSION * NUM_BLOCKS ));

    // Debug print for size allocated.
    uint32_t temp_size = (sizeof(int16_t) * CUDAPOA_MAX_MATRIX_DIMENSION * CUDAPOA_MAX_SEQUENCE_SIZE * NUM_BLOCKS );
    temp_size += 2 * (sizeof(int16_t) * CUDAPOA_MAX_MATRIX_DIMENSION * NUM_BLOCKS );
    std::cout << TABS << bid_ << " Allocated temp buffers of size " << (static_cast<float>(temp_size)  / (1024 * 1024)) << "MB" << std::endl;

    // Allocate graph buffers.
    CU_CHECK_ERR(cudaMalloc((void**) &nodes_d_, sizeof(uint8_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &node_alignments_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_ALIGNMENTS * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &node_alignment_count_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &incoming_edges_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &incoming_edge_count_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &outgoing_edges_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &outgoing_edge_count_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &incoming_edges_weights_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &outoing_edges_weights_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &sorted_poa_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &sorted_poa_node_map_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMalloc((void**) &sorted_poa_local_edge_count_d_, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));

    CU_CHECK_ERR(cudaMemset(nodes_d_, 0, sizeof(uint8_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(node_alignments_d_, 0, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_ALIGNMENTS * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(node_alignment_count_d_, 0, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(incoming_edges_d_,0,  sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(incoming_edge_count_d_,0,  sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(outgoing_edges_d_, 0, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(outgoing_edge_count_d_, 0, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(incoming_edges_weights_d_, 0, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(outoing_edges_weights_d_, 0, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(sorted_poa_d_, 0, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));
    CU_CHECK_ERR(cudaMemset(sorted_poa_local_edge_count_d_, 0, sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ));

    // Debug print for size allocated.
    temp_size = sizeof(uint8_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_ALIGNMENTS * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * CUDAPOA_MAX_NODE_EDGES * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ;
    temp_size += sizeof(uint16_t) * CUDAPOA_MAX_NODES_PER_WINDOW * NUM_BLOCKS ;
    std::cout << TABS << bid_ << " Allocated temp buffers of size " << (static_cast<float>(temp_size)  / (1024 * 1024)) << "MB" << std::endl;
}

CUDABatchProcessor::~CUDABatchProcessor()
{
    // Destroy CUDA stream.
    CU_CHECK_ERR(cudaStreamDestroy(stream_));

    // Free all the host and CUDA memory.
    CU_CHECK_ERR(cudaFree(consensus_d_));

    CU_CHECK_ERR(cudaFree(scores_d_));
    CU_CHECK_ERR(cudaFree(traceback_i_d_));
    CU_CHECK_ERR(cudaFree(traceback_j_d_));

    std::cout << TABS << "Destroyed buffers." << std::endl;

    CU_CHECK_ERR(cudaFree(nodes_d_));
    CU_CHECK_ERR(cudaFree(node_alignments_d_));
    CU_CHECK_ERR(cudaFree(node_alignment_count_d_));
    CU_CHECK_ERR(cudaFree(incoming_edges_d_));
    CU_CHECK_ERR(cudaFree(incoming_edge_count_d_));
    CU_CHECK_ERR(cudaFree(outgoing_edges_d_));
    CU_CHECK_ERR(cudaFree(outgoing_edge_count_d_));
    CU_CHECK_ERR(cudaFree(incoming_edges_weights_d_));
    CU_CHECK_ERR(cudaFree(outoing_edges_weights_d_));
    CU_CHECK_ERR(cudaFree(sorted_poa_d_));
}

bool CUDABatchProcessor::doesWindowFit(std::shared_ptr<Window> window) const
{
    if (window->sequences_.size() > max_depth_per_window_)
    {
        std::cerr << TABS << bid_ << " Number of sequences in window is greater than max depth!\n";
    }
    // Checks if adding new window will go over either the MAX_SEQUENCES
    // or max_windows_ count of the batch.
    return (windows_.size() + 1 <= max_windows_);
}

bool CUDABatchProcessor::addWindow(std::shared_ptr<Window> window)
{
    if (doesWindowFit(window))
    {
        windows_.push_back(window);
        sequence_count_ += std::min(max_depth_per_window_, (uint32_t) window->sequences_.size());
        return true;
    }
    else
    {
        return false;
    }
}

bool CUDABatchProcessor::hasWindows() const
{
    return (windows_.size() != 0);
}



void CUDABatchProcessor::generateMemoryMap()
{
    // Fill host/cuda memory with sequence information.
    uint32_t global_sequence_idx = 0; // count total number of sequences across all windows
    uint32_t num_nucleotides_copied = 0; // count total number of bases across all windows

    auto num_windows = windows_.size();
    for(uint32_t i = 0; i < num_windows; i++)
    {
        nvidia::cudapoa::WindowDetails window_details;
        window_details.seq_starts = num_nucleotides_copied;
        window_details.seq_len_buffer_offset = global_sequence_idx;

        auto window = windows_.at(i);
        auto num_seqs = (uint16_t) window->sequences_.size();
        for(uint32_t j = 0; j < num_seqs; j++){
            auto seq = window->sequences_.at(j);

            memcpy(&(inputs_h_[num_nucleotides_copied]),
                   seq.first,
                   seq.second);

            auto seq_len = (uint16_t) seq.second;
            num_nucleotides_copied += seq_len;

            sequence_lengths_h_[global_sequence_idx] = seq_len;
            global_sequence_idx++;
        }

        window_details.num_seqs = num_seqs;
        window_details_h_ [i] = window_details;
    }

    //Copy sequencecs, sequence lengths and window details to device
    CU_CHECK_ERR(cudaMemcpyAsync(inputs_d_, inputs_h_,
                                 num_nucleotides_copied * sizeof(uint8_t), cudaMemcpyHostToDevice, stream_));
    CU_CHECK_ERR(cudaMemcpyAsync(window_details_d_, window_details_h_,
                                 num_windows * sizeof(nvidia::cudapoa::WindowDetails), cudaMemcpyHostToDevice, stream_));
    CU_CHECK_ERR(cudaMemcpyAsync(sequence_lengths_d_, sequence_lengths_h_,
                                 global_sequence_idx * sizeof(uint16_t), cudaMemcpyHostToDevice, stream_));

}


void CUDABatchProcessor::generatePOA()
{
    // Launch kernel to run 1 POA per thread in thread block.
    std::cout << TABS << bid_ << " Launching kernel for " << windows_.size() << std::endl;
    nvidia::cudapoa::generatePOA(consensus_d_,
                                 inputs_d_,
                                 sequence_lengths_d_,
                                 window_details_d_,
                                 windows_.size(),
                                 NUM_THREADS,
                                 NUM_BLOCKS,
                                 stream_,
                                 scores_d_,
                                 traceback_i_d_,
                                 traceback_j_d_,
                                 nodes_d_,
                                 incoming_edges_d_,
                                 incoming_edge_count_d_,
                                 outgoing_edges_d_,
                                 outgoing_edge_count_d_,
                                 incoming_edges_weights_d_,
                                 outoing_edges_weights_d_,
                                 sorted_poa_d_,
                                 sorted_poa_node_map_d_,
                                 node_alignments_d_,
                                 node_alignment_count_d_,
                                 sorted_poa_local_edge_count_d_);
    CU_CHECK_ERR(cudaPeekAtLastError());
    std::cout << TABS << bid_ << " Launched kernel" << std::endl;
}

void CUDABatchProcessor::getConsensus()
{
    std::cout << TABS << bid_ << " Launching memcpy D2H" << std::endl;
    CU_CHECK_ERR(cudaMemcpy2DAsync(consensus_h_.get(),
				   CUDAPOA_MAX_SEQUENCE_SIZE,
				   consensus_d_,
				   consensus_pitch_,
				   CUDAPOA_MAX_SEQUENCE_SIZE,
				   max_windows_,
				   cudaMemcpyDeviceToHost,
				   stream_));
    CU_CHECK_ERR(cudaStreamSynchronize(stream_));

    std::cout << TABS << bid_ << " Finished memcpy D2H" << std::endl;

    for(uint32_t i = 0; i < windows_.size(); i++)
    {
        char* c = reinterpret_cast<char *>(&consensus_h_[i * CUDAPOA_MAX_SEQUENCE_SIZE]);
        windows_.at(i)->consensus_ = std::string(c);
    }
}

bool CUDABatchProcessor::generateConsensus()
{
    // Generate consensus for all windows in the batch
    generateMemoryMap();
    generatePOA();
    getConsensus();

    return true;
}

void CUDABatchProcessor::reset()
{
    windows_.clear();
    sequence_count_ = 0;

    // Clear host and device memory.
    memset(&inputs_h_[0], 0, max_windows_ * max_depth_per_window_ * CUDAPOA_MAX_SEQUENCE_SIZE);
    CU_CHECK_ERR(cudaMemsetAsync(inputs_d_, 0, max_windows_ * max_depth_per_window_ * CUDAPOA_MAX_SEQUENCE_SIZE, stream_));
}

} // namespace racon
