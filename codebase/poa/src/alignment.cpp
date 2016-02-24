/**
 * @file alignment.cpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Alignment class implementation file 
 * @details Implementation file for Alignment class used for
 * calculating local alignment between graph and sequence. Class is
 * based on https://github.com/ljdursi/poapy/blob/master/seqgraphalignment.py
 * python implementation
 */
#include <algorithm>
#include <vector>
#include <string>
#include <tuple>
#include <limits>
#include <iostream>
#include <cassert>

#include "./alignment.hpp"
#include "./graph.hpp"

using std::get;
using std::make_pair;
using std::make_tuple;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::string;
using std::tuple;
using std::vector;

namespace POA {

    int Alignment::match_score_ = 4;
    int Alignment::mismatch_score_ = -2;
    int Alignment::open_gap_score_ = -4;
    int Alignment::extend_gap_score_ = -2;


    Alignment::Alignment(const string& sequence,
                         const Graph& graph): sequence_(sequence),
                                              graph_(graph) {
        assert(sequence.size() > 0);
    }

    /**
     * @brief returns move score
     * @details returns move score
     *
     * @param move move
     * @return move score
     */
    int move_score(const move& move) {
        return get<0>(move);
    }


    /**
     * @brief returns move graph_idx
     * @details returns move graph_idx
     *
     * @param move move
     * @return move graph_index
     */
    int move_graph_idx(const move& move) {
        return get<1>(move);
    }


    /**
     * @brief returns move seq_idx
     * @details returns move seq_idx
     *
     * @param move move
     * @return move seq_idx
     */
    int move_seq_idx(const move& move) {
        return get<2>(move);
    }


    /**
     * @brief returns move type
     * @details returns move type
     *
     * @param move move
     * @return move type
     */
    char move_type(const move& move) {
        return get<3>(move);
    }


    void Alignment::backtrack() {
        while (dp(max_i_, max_j_).score > 0 && !(max_i_ == 0 && max_j_ == 0)) {
            int next_i = dp(max_i_, max_j_).prev_graph_idx;
            int next_j = dp(max_i_, max_j_).prev_seq_idx;
            uint32_t seq_idx = max_j_ - 1;
            uint32_t node_id = index_to_nodeID_[max_i_ - 1];

            seq_idxs_.emplace_front(next_j != max_j_ ? seq_idx : -1);
            node_ids_.emplace_front(next_i != max_i_ ? node_id : -1);

            max_i_ = next_i;
            max_j_ = next_j;
        }
    }


    void Alignment::align() {
        align_banded_starting_at(0, -1);
    }

    void Alignment::align_starting_at(const uint32_t pos) {
        align_banded_starting_at(pos, -1);
    }

    void Alignment::align_banded_starting_at(const uint32_t min_pos, const int band_width) {
        max_score_ = numeric_limits<int>::min();
        max_i_ = -1;
        max_j_ = -1;

        init_dp(min_pos, band_width);

        auto match_score = [](char seq_base, char node_base) -> int {
            if (seq_base == node_base) {
                return match_score_;
            } else {
                return mismatch_score_;
            }
        };

        // modified smith-waterman(SW)
        // operation candidates
        // char values:
        // -> 'I' - insertion
        // -> 'D' - deletion
        // -> 'M' - match or mismatch
        for (uint32_t i = 0; i < valid_nodes_num_; ++i) {
            const auto node_id = index_to_nodeID_[i];
            auto& node = graph_.getNode(node_id);
            const char base = node->base();

            // calculate sequence substring that should be part of DP
            int lo = seq_limits_[i + 1].first, hi = seq_limits_[i + 1].second;
            for (int j = lo; j < hi; ++j) {
                // insertion from sequence, unchanged as in SW
                move best_candidate(dp(i + 1, j).score + dp(i + 1, j).insert_cost,
                                        i + 1, j, 'I');

                if (node->getPredecessorsIds().size()) {
                  // for every other operation I have to check for all predeccesors of
                  // current node
                  for (auto node_id : node->getPredecessorsIds()) {
                    auto prev_index = index_from_node_id(node_id);

                    // match/mismatch
                    move mm(dp(prev_index + 1, j).score +
                        match_score(sequence_[j], base),
                        prev_index + 1, j, 'M');
                    if (move_score(mm) >= move_score(best_candidate)) {
                      best_candidate = mm;
                    }

                    // insertion from graph to sequence / deletion
                    move ins(dp(prev_index + 1, j + 1).score +
                        dp(prev_index + 1, j + 1).delete_cost,
                        prev_index + 1, j + 1, 'D');
                    if (move_score(ins) >= move_score(best_candidate)) {
                      best_candidate = ins;
                    }
                  }
                } else {
                  int prev_index = -1;
                  // match/mismatch
                  move mm(dp(prev_index + 1, j).score +
                      match_score(sequence_[j], base),
                      prev_index + 1, j, 'M');
                  if (move_score(mm) >= move_score(best_candidate)) {
                    best_candidate = mm;
                  }

                  // insertion from graph to sequence / deletion
                  move ins(dp(prev_index + 1, j + 1).score +
                      dp(prev_index + 1, j + 1).delete_cost,
                      prev_index + 1, j + 1, 'D');
                  if (move_score(ins) >= move_score(best_candidate)) {
                    best_candidate = ins;
                  }
                }

                // update dp and backtrack tables
                dp(i + 1, j + 1).score = move_score(best_candidate);
                dp(i + 1, j + 1).prev_graph_idx = move_graph_idx(best_candidate);
                dp(i + 1, j + 1).prev_seq_idx = move_seq_idx(best_candidate);

                char move = move_type(best_candidate);
                if (move == 'I') {
                    dp(i + 1, j + 1).insert_cost = extend_gap_score_;
                } else if (move == 'D') {
                    dp(i + 1, j + 1).delete_cost = extend_gap_score_;
                }


                // because of local alignment minimum score should be zero
                if (dp(i + 1, j + 1).score < 0) {
                    dp(i + 1, j + 1).score = 0;
                    dp(i + 1, j + 1).prev_graph_idx = -1;
                    dp(i + 1, j + 1).prev_seq_idx = -1;
                }


                // update max score and its position
                if (dp(i + 1, j + 1).score >= max_score_) {
                    max_i_ = i + 1;
                    max_j_ = j + 1;
                    max_score_ = dp(i + 1, j + 1).score;
                }
            }
        }

        backtrack();
    }


    void Alignment::init_dp(const int min_pos, const int band_width) {
        valid_nodes_num_ = 0;

        max_valid_node_id_ = numeric_limits<int>::min();
        min_valid_node_id_ = numeric_limits<int>::max();

        const vector<uint32_t>& nodes_ids = const_cast<Graph&>(graph_).getNodesIds();
        for (int node_id: nodes_ids) {
          if (graph_.max_node_distance(node_id) >= min_pos) {
            valid_nodes_num_++;
            max_valid_node_id_ = max(max_valid_node_id_, node_id);
            min_valid_node_id_ = min(min_valid_node_id_, node_id);
          }
        }

        index_to_nodeID_.resize(valid_nodes_num_);
        nodeID_to_index_.resize(max_valid_node_id_ - min_valid_node_id_ + 1, -1);
        seq_limits_.resize(valid_nodes_num_ + 1);
        dp_width_ = 0;

        uint32_t idx = 0;
        for (int node_id: nodes_ids) {
            if (graph_.max_node_distance(node_id) < min_pos) {
                continue;
            }

            nodeID_to_index_[node_id - min_valid_node_id_] = idx;
            index_to_nodeID_[idx] = node_id;

            // calculate sequence_ substring that has to be aligned to current node
            int lo = 0, hi = sequence_.length();
            if (band_width >= 0) {
                lo = max(0, (int) (graph_.min_node_distance(node_id) - band_width - min_pos));
                hi = min((int) sequence_.length(), (int) (graph_.max_node_distance(node_id) + band_width - min_pos + 1));
            }
            seq_limits_[idx + 1] = make_pair(lo, hi);

            dp_width_ = max((int) dp_width_, hi - lo);

            idx++;
        }

        // init first row limits (belongs to no node)
        seq_limits_[0] = seq_limits_[1];

        // init dynamic programming smith waterman matrix
        dp_.resize((valid_nodes_num_ + 1)*(dp_width_ + 1), dp_el(0, -1, -1, open_gap_score_, open_gap_score_));
    }

    int Alignment::index_from_node_id(const uint32_t node_id) const {
        if ((int) node_id < min_valid_node_id_ || (int) node_id > max_valid_node_id_) {
            return -1;
        }
        return nodeID_to_index_[(int) node_id - min_valid_node_id_];
    }

    inline dp_el& Alignment::dp(const uint32_t i, const uint32_t j) {
        /**
         * since we do not store the whole matrix in memory (because we do not need it),
         * we use precalculated values of sequence edges for each node to store/fetch
         * just these parts (marked with x below).
         * _______________________
         * |xxxxxxxxxxxx.........|
         * |xxxxxxxxxxxxx........|
         * |.xxxxxxxxxxxxx.......|
         * |..xxxxxxxxxxxxx......|
         * |...xxxxxxxxxxxxx.....|
         * |.....xxxxxxxxxxxxx...|
         * |......xxxxxxxxxxxxx..|
         * |.......xxxxxxxxxxxxx.|
         * |........xxxxxxxxxxxxx|
         * |.........xxxxxxxxxxxx|
         * |..........xxxxxxxxxxx|
         * -----------------------
         */
        const int lo = seq_limits_[i].first;
        return dp_[i*(dp_width_ + 1) + j - lo];
    }
}
