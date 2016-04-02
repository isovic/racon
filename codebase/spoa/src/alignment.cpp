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

    assert(type == AlignmentType::kNW || type == AlignmentType::kSW || type == AlignmentType::kOV);
}

AlignmentParams::AlignmentParams(int16_t m, int16_t mm, int16_t ins_opn,
    int16_t ins_ext, int16_t del_opn, int16_t del_ext, AlignmentType t) :
        match(m), mismatch(mm), insertion_open(ins_opn), insertion_extend(ins_ext),
        deletion_open(del_opn), deletion_extend(del_ext), type(t) {

    assert(type == AlignmentType::kNW || type == AlignmentType::kSW || type == AlignmentType::kOV);
}

AlignmentParams::~AlignmentParams() {
}

std::unique_ptr<Alignment> createAlignment(const std::string& sequence,
    std::shared_ptr<Graph> graph, AlignmentParams params) {

    assert(sequence.size() != 0);

    return std::unique_ptr<Alignment>(new Alignment(sequence, graph,
        std::move(params)));
}

Alignment::Alignment(const std::string& sequence, std::shared_ptr<Graph> graph,
    AlignmentParams params) :
        sequence_profile_(256), graph_(graph), params_(std::move(params)),
        matrix_width_(sequence.size() + 1),
        matrix_height_(graph->nodes().size() + 1),
        H_(matrix_width_ * matrix_height_, 0),
        F_(matrix_width_ * matrix_height_, 0),
        E_(matrix_width_ * matrix_height_, 0),
        is_aligned_(false),
        max_i_(-1), max_j_(-1), max_score_(0),
        node_id_to_graph_id_(),
        is_backtracked_(false),
        alignment_node_ids_(),
        alignment_seq_ids_() {

    for (const auto& c: graph->alphabet()) {
        sequence_profile_[c].reserve(sequence.size());
        for (const auto& s: sequence) {
            sequence_profile_[c].push_back(c == s ? params_.match : params_.mismatch);
        }
    }

    int32_t big_negative_value = std::numeric_limits<int32_t>::min() + 1000;
    for (uint32_t j = 1; j < matrix_width_; ++j) {
        F_[j] = big_negative_value;
    }
    for (uint32_t i = 1; i < matrix_height_; ++i) {
        E_[i * matrix_width_] = big_negative_value;
    }

    graph_->topological_sort();
    const auto& sorted_nodes_ids = graph_->sorted_nodes_ids();

    node_id_to_graph_id_.resize(sorted_nodes_ids.size());
    for (uint32_t i = 0; i < sorted_nodes_ids.size(); ++i) {
        node_id_to_graph_id_[sorted_nodes_ids[i]] = i;
    }

    if (params_.type == AlignmentType::kNW) {
        max_score_ = big_negative_value;

        for (uint32_t j = 1; j < matrix_width_; ++j) {
            H_[j] = params_.insertion_open + (j - 1) * params_.insertion_extend;
            E_[j] = H_[j];
        }

        for (uint32_t node_id: sorted_nodes_ids) {
            const auto& node = graph_->node(node_id);
            uint32_t i = node_id_to_graph_id_[node_id] + 1;

            if (node->in_edges().size() == 0) {
                H_[i * matrix_width_] = params_.deletion_open;
            } else {
                int32_t max_score = big_negative_value;
                for (const auto& edge: node->in_edges()) {
                    uint32_t pred_i = node_id_to_graph_id_[edge->begin_node_id()] + 1;
                    max_score = std::max(max_score, H_[pred_i * matrix_width_]);
                }
                H_[i * matrix_width_] = max_score + params_.deletion_extend;
            }
            F_[i * matrix_width_] = H_[i * matrix_width_];
        }
    }

    // print_matrix();
}

Alignment::~Alignment() {
}

void Alignment::align_sequence_to_graph() {

    if (is_aligned_ == true) {
        return;
    }

    graph_->topological_sort();
    const auto& sorted_nodes_ids = graph_->sorted_nodes_ids();

    for (uint32_t node_id: sorted_nodes_ids) {
        const auto& node = graph_->node(node_id);
        const auto& char_profile = sequence_profile_[node->letter()];
        uint32_t i = node_id_to_graph_id_[node_id] + 1;

        int32_t* H_row = &H_[i * matrix_width_];
        int32_t* F_row = &F_[i * matrix_width_];

        uint32_t pred_i = node->in_edges().empty() ? 0 :
            node_id_to_graph_id_[node->in_edges().front()->begin_node_id()] + 1;

        int32_t* H_pred_row = &H_[pred_i * matrix_width_];
        int32_t* F_pred_row = &F_[pred_i * matrix_width_];

        for (uint32_t j = 1; j < matrix_width_; ++j) {
            // update F
            F_row[j] = std::max(H_pred_row[j] + params_.insertion_open, F_pred_row[j] + params_.insertion_extend);
            // update H
            H_row[j] = std::max(H_pred_row[j - 1] + char_profile[j - 1], F_row[j]);
        }

        // check other predeccessors
        for (uint32_t p = 1; p < node->in_edges().size(); ++p) {
            pred_i = node_id_to_graph_id_[node->in_edges()[p]->begin_node_id()] + 1;

            H_pred_row = &H_[pred_i * matrix_width_];
            F_pred_row = &F_[pred_i * matrix_width_];

            for (uint32_t j = 1; j < matrix_width_; ++j) {
                // update F
                F_row[j] = std::max(F_row[j], std::max(H_pred_row[j] + params_.insertion_open, F_pred_row[j] + params_.insertion_extend));
                // update H
                H_row[j] = std::max(H_row[j], std::max(H_pred_row[j - 1] + char_profile[j - 1], F_row[j]));
            }
        }

        int32_t* E_row = &E_[i * matrix_width_];

        for (uint32_t j = 1; j < matrix_width_; ++j) {
            // update E
            E_row[j] = std::max(H_row[j - 1] + params_.deletion_open, E_row[j - 1] + params_.deletion_extend);
            // update H
            H_row[j] = std::max(H_row[j], E_row[j]);

            if (params_.type == AlignmentType::kSW) {
                H_row[j] = std::max(H_row[j], 0);
                update_max_score(H_row, i, j);
            } else if (params_.type == AlignmentType::kNW) {
                if (j == matrix_width_ - 1 && node->out_edges().size() == 0) {
                    update_max_score(H_row, i, j);
                }
            } else if (params_.type == AlignmentType::kOV) {
                if (j == matrix_width_ - 1 || node->out_edges().size() == 0) {
                    update_max_score(H_row, i, j);
                }
            }
        }
    }

    is_aligned_ = true;

    // print_matrix();
}

int32_t Alignment::score() const {
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

    auto sw_condition = [&]() { return (H_[i * matrix_width_ + j] == 0) ? false : true; };
    auto nw_condition = [&]() { return (i == 0 && j == 0) ? false : true; };
    auto ov_condition = [&]() { return (i == 0 || j == 0) ? false : true; };

    const auto& graph_id_to_node_id = graph_->sorted_nodes_ids();

    uint32_t prev_i = 0, prev_j = 0;

    while ((params_.type == AlignmentType::kSW && sw_condition()) ||
        (params_.type == AlignmentType::kNW && nw_condition()) ||
        (params_.type == AlignmentType::kOV && ov_condition())) {

        // bloody backtrack
        auto H_ij = H_[i * matrix_width_ + j];
        bool predecessor_found = false;

        if (i != 0) {
            const auto& node = graph_->node(graph_id_to_node_id[i - 1]);
            int32_t match_cost = j != 0 ? sequence_profile_[node->letter()][j - 1] : 0;

            uint32_t pred_i = node->in_edges().empty() ? 0 :
                node_id_to_graph_id_[node->in_edges().front()->begin_node_id()] + 1;

            if (j != 0 && H_ij == H_[pred_i * matrix_width_ + (j - 1)] + match_cost) {
                prev_i = pred_i;
                prev_j = j - 1;
                predecessor_found = true;
            } else if ((H_ij == F_[pred_i * matrix_width_ + j] + params_.insertion_extend) ||
                (H_ij == H_[pred_i * matrix_width_ + j] + params_.insertion_open)) {
                prev_i = pred_i;
                prev_j = j;
                predecessor_found = true;
            }

            if (!predecessor_found) {
                const auto& edges = node->in_edges();
                for (uint32_t p = 1; p < edges.size(); ++p) {
                    uint32_t pred_i = node_id_to_graph_id_[edges[p]->begin_node_id()] + 1;
                    if (j != 0 && H_ij == H_[pred_i * matrix_width_ + (j - 1)] + match_cost) {
                        prev_i = pred_i;
                        prev_j = j - 1;
                        predecessor_found = true;
                        break;
                    }
                    if ((H_ij == F_[pred_i * matrix_width_ + j] + params_.insertion_extend) ||
                        (H_ij == H_[pred_i * matrix_width_ + j] + params_.insertion_open)){
                        prev_i = pred_i;
                        prev_j = j;
                        predecessor_found = true;
                        break;
                    }
                }
            }
        }

        if (!predecessor_found && H_ij == E_[i * matrix_width_ + j]) {
            prev_i = i;
            prev_j = j - 1;
        }

        alignment_node_ids_.emplace_back(i == prev_i ? -1 : graph_id_to_node_id[i - 1]);
        alignment_seq_ids_.emplace_back(j == prev_j ? -1 : j - 1);

        i = prev_i;
        j = prev_j;
    }

    std::reverse(alignment_node_ids_.begin(), alignment_node_ids_.end());
    std::reverse(alignment_seq_ids_.begin(), alignment_seq_ids_.end());

    is_backtracked_ = true;
}

inline void Alignment::update_max_score(int32_t* H_row, uint32_t i, uint32_t j) {
    if (max_score_ < H_row[j]) {
        max_score_ = H_row[j];
        max_i_ = i;
        max_j_ = j;
    }
}

void Alignment::print_matrix() {
    for (uint32_t i = 0; i < matrix_height_; ++i) {
        for (uint32_t j = 0; j < matrix_width_; ++j) {
            printf("(%3d %3d %3d) ", H_[i * matrix_width_ + j],
                E_[i * matrix_width_ + j], F_[i * matrix_width_ + j]);
        }
        printf("\n");
    }
    printf("\n");
}
