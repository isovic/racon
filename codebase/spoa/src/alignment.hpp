/*!
 * @file alignment.hpp
 *
 * @brief Alignment class header file
 */

#pragma once

#include <assert.h>
#include <string>
#include <memory>
#include <vector>

enum class AlignmentType {
    kSW, // Smith Waterman
    kNW, // Needleman Wunsch
    kOV // Overlap
};

class AlignmentParams {
public:

    AlignmentParams(int16_t match, int16_t mismatch, int16_t gap_open,
        int16_t gap_extend, AlignmentType type);
    AlignmentParams(int16_t match, int16_t mismatch, int16_t insertion_open,
        int16_t insertion_extend, int16_t deletion_open, int16_t deletion_extend,
        AlignmentType type);
    ~AlignmentParams();

    int16_t match;
    int16_t mismatch;
    int16_t insertion_open;
    int16_t insertion_extend;
    int16_t deletion_open;
    int16_t deletion_extend;
    AlignmentType type;
};

class Graph;
using GraphSharedPtr = std::shared_ptr<Graph>;

class Alignment;
std::unique_ptr<Alignment> createAlignment(const std::string& sequence,
    GraphSharedPtr graph, AlignmentParams params);

class Alignment {
public:

    ~Alignment();

    void align_sequence_to_graph();

    void backtrack();

    int32_t alignment_score() const;

    const std::vector<int32_t>& alignment_node_ids() const {
        assert(is_backtracked_ == true && "No backtrack done!");
        return alignment_node_ids_;
    }

    const std::vector<int32_t>& alignment_seq_ids() const {
        assert(is_backtracked_ == true && "No backtrack done!");
        return alignment_seq_ids_;
    }

    friend std::unique_ptr<Alignment> createAlignment(const std::string& sequence,
        GraphSharedPtr graph, AlignmentParams params);

private:

    Alignment(const std::string& sequence, GraphSharedPtr graph,
        AlignmentParams params);
    Alignment(const Alignment&) = delete;
    const Alignment& operator=(const Alignment&) = delete;

    class MatrixElement {
    public:

        MatrixElement(int32_t score, int32_t prev_i, int32_t prev_j,
            int16_t insertion_cost, int16_t deletion_cost);
        ~MatrixElement();

        int32_t score;
        int32_t prev_i;
        int32_t prev_j;
        int16_t deletion_cost;
        int16_t insertion_cost;
    };

    class MatrixMove {
    public:

        MatrixMove(int32_t score, int32_t i, int32_t j, int32_t type);
        ~MatrixMove();

        int32_t score;
        int32_t i;
        int32_t j;
        int32_t type; // 0 - diagonal, 1 - insertion, 2 - deletion
    };

    inline MatrixElement& matrix(uint32_t i, uint32_t j) {
        assert(i < matrix_height_);
        assert(j < matrix_width_);
        return matrix_[i * matrix_width_ + j];
    }

    inline const MatrixElement& matrix(uint32_t i, uint32_t j) const {
        assert(i < matrix_height_);
        assert(j < matrix_width_);
        return matrix_[i * matrix_width_ + j];
    }

    void print_matrix();

    std::string sequence_;
    GraphSharedPtr graph_;
    AlignmentParams params_;

    uint32_t matrix_width_;
    uint32_t matrix_height_;
    std::vector<MatrixElement> matrix_;

    bool is_aligned_;
    int32_t max_i_;
    int32_t max_j_;
    int32_t max_score_;
    std::vector<uint32_t> node_id_to_graph_id_;

    bool is_backtracked_;
    std::vector<int32_t> alignment_node_ids_;
    std::vector<int32_t> alignment_seq_ids_;
};
