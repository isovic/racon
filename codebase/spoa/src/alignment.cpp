/*!
 * @file alignment.cpp
 *
 * @brief Alignment class source file
 */

#include <limits>
#include <algorithm>

#include "node.hpp"
#include "edge.hpp"
#include "graph.hpp"
#include "alignment.hpp"


AlignmentParams::AlignmentParams(int16_t m, int16_t mm, int16_t gap_opn,
    int16_t gap_ext, AlignmentType t) :
        match(m), mismatch(mm), insertion_open(gap_opn), insertion_extend(gap_ext),
        deletion_open(gap_opn), deletion_extend(gap_ext), type(t) {
}

AlignmentParams::AlignmentParams(int16_t m, int16_t mm, int16_t ins_opn,
    int16_t ins_ext, int16_t del_opn, int16_t del_ext, AlignmentType t) :
        match(m), mismatch(mm), insertion_open(ins_opn), insertion_extend(ins_ext),
        deletion_open(del_opn), deletion_extend(del_ext), type(t) {
}

AlignmentParams::~AlignmentParams() {
}

Alignment::MatrixElement::MatrixElement(int32_t s, int32_t p_i, int32_t p_j,
    int16_t del, int16_t ins) :
        score(s), prev_i(p_i), prev_j(p_j), deletion_cost(del), insertion_cost(ins) {
}

Alignment::MatrixElement::~MatrixElement() {
}

Alignment::MatrixMove::MatrixMove(int32_t s, int32_t k, int32_t l, int32_t t) :
        score(s), i(k), j(l), type(t) {
}

Alignment::MatrixMove::~MatrixMove() {
}

std::unique_ptr<Alignment> createAlignment(const std::string& sequence,
    GraphSharedPtr graph, AlignmentParams params) {

    return std::unique_ptr<Alignment>(new Alignment(sequence, graph,
        std::move(params)));
}

Alignment::Alignment(const std::string& sequence, GraphSharedPtr graph,
    AlignmentParams params) :
        sequence_(sequence), graph_(graph), params_(std::move(params)) {

    assert(sequence_.size() != 0);

    matrix_width_ = sequence_.size() + 1;
    matrix_height_ = graph_->nodes().size() + 1;

    matrix_.resize(matrix_width_ * matrix_height_, MatrixElement(0, -1, -1,
        params_.insertion_open, params_.deletion_open));

    is_aligned_ = false;
    max_i_ = -1;
    max_j_ = -1;
    max_score_ = params.type == AlignmentType::kNW ? std::numeric_limits<int32_t>::min() : 0;

    graph_->topological_sort();
    const auto& sorted_nodes_ids = graph_->sorted_nodes_ids();
    node_id_to_graph_id_.resize(sorted_nodes_ids.size());

    for (uint32_t i = 0; i < sorted_nodes_ids.size(); ++i) {
        node_id_to_graph_id_[sorted_nodes_ids[i]] = i;
    }

    is_backtracked_ = false;

    if (params_.type == AlignmentType::kNW) {

        for (uint32_t j = 1; j < matrix_width_; ++j) {
            matrix(0, j).score = params_.insertion_open + (j - 1) * params_.insertion_extend;
            matrix(0, j).prev_i = 0;
            matrix(0, j).prev_j = j - 1;
        }

        for (uint32_t node_id: sorted_nodes_ids) {
            const auto& node = graph_->node(node_id);
            uint32_t i = node_id_to_graph_id_[node_id] + 1;

            if (node->in_edges().size() == 0) {
                matrix(i, 0).score = params_.deletion_open;
                matrix(i, 0).prev_i = 0;
                matrix(i, 0).prev_j = 0;
            } else {
                matrix(i, 0).score = std::numeric_limits<int32_t>::min();
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_i = node_id_to_graph_id_[edge->begin_node_id()] + 1;
                    if (matrix(i, 0).score < matrix(pred_i, 0).score) {
                        matrix(i, 0).score = matrix(pred_i, 0).score;
                        matrix(i, 0).prev_i = pred_i;
                    }
                }
                matrix(i, 0).prev_j = 0;
                matrix(i, 0).score += params_.deletion_extend;
            }
        }
    }

    //print_matrix();
}

Alignment::~Alignment() {
}

void Alignment::align_sequence_to_graph() {

    if (is_aligned_ == true) {
        return;
    }

    graph_->topological_sort();
    const auto& sorted_nodes_ids = graph_->sorted_nodes_ids();

    std::vector<MatrixMove> possible_moves;

    for (uint32_t node_id: sorted_nodes_ids) {
        const auto& node = graph_->node(node_id);
        char graph_letter = node->letter();
        uint32_t i = node_id_to_graph_id_[node_id] + 1;

        for (uint32_t j = 1; j < matrix_width_; ++j) {
            int32_t match_cost = graph_letter == sequence_[j - 1] ? params_.match : params_.mismatch;

            possible_moves.clear();

            if (node->in_edges().size() == 0) {
                // match/mismatch
                possible_moves.emplace_back(matrix(0, j - 1).score + match_cost,
                    0, j - 1, 0);
                // insertion to sequence
                possible_moves.emplace_back(matrix(0, j).score + matrix(0, j).insertion_cost,
                    0, j, 1);

            } else {
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_i = node_id_to_graph_id_[edge->begin_node_id()] + 1;

                    // match/mismatch
                    possible_moves.emplace_back(matrix(pred_i, j - 1).score + match_cost,
                        pred_i, j - 1, 0);
                    // insertion to sequence
                    possible_moves.emplace_back(matrix(pred_i, j).score + matrix(pred_i, j).insertion_cost,
                        pred_i, j, 1);
                }
            }

            // deletion from graph
            possible_moves.emplace_back(matrix(i, j - 1).score + matrix(i, j - 1).deletion_cost,
                i, j - 1, 2);

            // find best move
            int32_t max_idx = 0;
            for (uint32_t idx = 1; idx < possible_moves.size(); ++idx) {
                if (possible_moves[max_idx].score < possible_moves[idx].score) {
                    max_idx = idx;
                }
            }

            // update matrix field (i, j)
            matrix(i, j).score = possible_moves[max_idx].score;
            matrix(i, j).prev_i = possible_moves[max_idx].i;
            matrix(i, j).prev_j = possible_moves[max_idx].j;

            if (possible_moves[max_idx].type == 1) {
                matrix(i, j).insertion_cost = params_.insertion_extend;
            } else if (possible_moves[max_idx].type == 2) {
                matrix(i, j).deletion_cost = params_.deletion_extend;
            }

            if (params_.type == AlignmentType::kSW) {

                if (matrix(i, j).score < 0) {
                    matrix(i, j).score = 0;
                    matrix(i, j).prev_i = -1;
                    matrix(i, j).prev_j = -1;
                }

                if (max_score_ < matrix(i, j).score) {
                    max_score_ = matrix(i, j).score;
                    max_i_ = i;
                    max_j_ = j;
                }

            } else if (params_.type == AlignmentType::kNW) {

                if (j == matrix_width_ - 1 && node->out_edges().size() == 0) {
                    if (max_score_ < matrix(i, j).score) {
                        max_score_ = matrix(i, j).score;
                        max_i_ = i;
                        max_j_ = j;
                    }
                }

            } else if (params_.type == AlignmentType::kOV) {

                if (j == matrix_width_ - 1 || node->out_edges().size() == 0) {
                    if (max_score_ < matrix(i, j).score) {
                        max_score_ = matrix(i, j).score;
                        max_i_ = i;
                        max_j_ = j;
                    }
                }
            }
        }
    }

    is_aligned_ = true;

    //print_matrix();
}

int32_t Alignment::alignment_score() const {
    assert(is_aligned_ == true && "No alignment done!");
    return max_score_;
}

void Alignment::backtrack() {

    if (is_backtracked_ == true) {
        return;
    }
    assert(is_aligned_  == true && "No alignment done!");

    if (max_i_ == -1 && max_j_ == -1) { // no alignment found
        is_backtracked_ = true;
        return;
    }

    uint32_t i = max_i_;
    uint32_t j = max_j_;
    // fprintf(stderr, "Score, i, j = %d, %d, %d\n", max_score_, i, j);

    auto sw_condition = [&]() { return (matrix(i, j).score == 0) ? false : true; };
    auto nw_condition = [&]() { return (i == 0 && j == 0) ? false : true; };
    auto ov_condition = [&]() { return (i == 0 || j == 0) ? false : true; };

    const auto& graph_id_to_node_id = graph_->sorted_nodes_ids();

    //while ((params_.type != AlignmentType::kSW || matrix(i, j).score > 0)
    //    && !(i == 0 && j == 0)) {
    while ((params_.type == AlignmentType::kSW && sw_condition()) ||
        (params_.type == AlignmentType::kNW && nw_condition()) ||
        (params_.type == AlignmentType::kOV && ov_condition())) {

        uint32_t prev_i = matrix(i, j).prev_i;
        uint32_t prev_j = matrix(i, j).prev_j;

        alignment_node_ids_.emplace_back(i == prev_i ? -1 : graph_id_to_node_id[i - 1]);
        alignment_seq_ids_.emplace_back(j == prev_j ? -1 : j - 1);

        i = prev_i;
        j = prev_j;
    }

    std::reverse(alignment_node_ids_.begin(), alignment_node_ids_.end());
    std::reverse(alignment_seq_ids_.begin(), alignment_seq_ids_.end());

    is_backtracked_ = true;
}

void Alignment::print_matrix() {
    for (uint32_t i = 0; i < matrix_height_; ++i) {
        for (uint32_t j = 0; j < matrix_width_; ++j) {
            printf("(%3d, %3d, %3d) ", matrix(i, j).score, matrix(i, j).prev_i,
                matrix(i, j).prev_j);
            }
        printf("\n");
    }
    printf("\n");
}
