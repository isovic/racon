/**
 * @file graph.cpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Graph class implementation file
 * @details Implementation file for Graph class used for
 * partial order alignment algorithm. Class is
 * based on https://github.com/ljdursi/poapy/blob/master/poagraph.py
 * Graph class python implementation
 */
#include <cassert>
#include <string>
#include <vector>
#include <memory>
#include <stack>
#include <algorithm>
#include <deque>
#include <utility>
#include <unordered_map>
#include <limits>
#include <iostream>

#include "./graph.hpp"
#include "./node.hpp"

using std::vector;
using std::string;
using std::unordered_map;
using std::shared_ptr;
using std::stack;
using std::reverse;
using std::deque;
using std::pair;
using std::make_pair;
using std::max;
using std::min;
using std::numeric_limits;

namespace POA {

    Graph::Graph(const string& seq, const string& label): next_id_(0),
                                                          is_sorted(true) {
        addUnmatchedSequence(seq, label, true);
        topological_sort();
    }

    void Graph::addNode(char base) {
        shared_ptr<Node> node(new Node(next_id_, base, *this));
        nodes_.emplace_back(node);
        ++next_id_;
        is_sorted = false;
    }


    void Graph::addEdge(uint32_t id_A, uint32_t id_B, const string& label) {
        if (id_A >= nodes_.size() || id_B >= nodes_.size()) {
            fprintf(stderr, "Error: Adding an edge with non existing nodes\n");
            exit(1);
        }
        nodes_[id_A]->addOutEdge(id_B, label);
        nodes_[id_B]->addInEdge(id_A, label);
        is_sorted = false;
    }


    const shared_ptr<Node>& Graph::getNode(uint32_t id) const {
        return nodes_[id];
    }


    const vector<uint32_t>& Graph::getNodesIds() {
        if (!is_sorted) {
            topological_sort();
        }
        return nodes_order_;
    }

    const vector<shared_ptr< Node >>& Graph::getNodes() const {
        return nodes_;
    }


    uint32_t Graph::getNodesNum() const {
        return nodes_.size();
    }


    pair<uint32_t, uint32_t> Graph::addUnmatchedSequence(const string& sequence,
                                                         const string& label,
                                                         bool update) {
        if (sequence.empty()) {
            return make_pair(-1, -1);
        }

        int start_node_id = -1;
        int last_node_id = -1;

        for (auto& base : sequence) {
            addNode(base);
            int added_node_id = next_id_ - 1;
            if (start_node_id == -1) {
                start_node_id = added_node_id;
            }
            if (last_node_id != -1 && added_node_id != -1) {
                addEdge(last_node_id, added_node_id, label);
            }
            last_node_id = added_node_id;
        }

        // is_sorted = true;
        if (update) {
            sequences_.emplace_back(sequence);
            labels_.emplace_back(label);
            start_ids_.emplace_back(start_node_id);
        }
        return make_pair(start_node_id, last_node_id);
    }


    void Graph::insertSequenceAlignment(const Alignment& alignment,
                                        const string& seq,
                                        const string& label) {
        const string& aln_sequence = alignment.sequence();
        const deque<int>& seq_idxs = alignment.seq_idxs();
        const deque<int>& node_ids = alignment.node_ids();

        int first_id = -1;
        int head_id = -1;
        int tail_id = -1;
        pair<uint32_t, uint32_t> prefix_end_ids;
        pair<uint32_t, uint32_t> suffix_end_ids;

        // because of local alignment prefix or sufix of sequence can be unaligned
        // so we add it directly to graph
        deque<uint32_t> valid_seq_idxs;
        for (auto idx : seq_idxs) {
            if (idx != -1) {
                valid_seq_idxs.emplace_back(idx);
            }
        }
        uint32_t aln_seq_start_idx = valid_seq_idxs.front();
        uint32_t aln_seq_end_idx = valid_seq_idxs.back();

        if (aln_seq_start_idx > 0) {
            prefix_end_ids = addUnmatchedSequence(
                                    aln_sequence.substr(0, aln_seq_start_idx),
                                    label,
                                    false);
            first_id = prefix_end_ids.first;
            head_id = prefix_end_ids.second;
        }

        if (aln_seq_end_idx < aln_sequence.length()) {
            suffix_end_ids = addUnmatchedSequence(
                                    aln_sequence.substr(aln_seq_end_idx + 1),
                                    label,
                                    false);
            tail_id = suffix_end_ids.first;
        }

        // aligned part of sequence
        uint32_t size = max(seq_idxs.size(), node_ids.size());
        for (uint32_t i = 0; i < size; ++i) {
            auto& seq_idx = seq_idxs[i];
            auto& match_id = node_ids[i];

            if (seq_idx == -1) {
                continue;
            }

            int node_id = -1;
            char base = aln_sequence[seq_idx];
            if (match_id == -1) {
                // if sequence base unmatched with graph node add new node
                addNode(base);
                node_id = next_id_ - 1;
            } else if (nodes_[match_id]->base() == base) {
                // if sequence base matched to a node with same base
                node_id = match_id;
            } else {
                // if sequence base matched to a node with different base
                // which is aligned to a node with same base
                int found_node_id = -1;
                for (auto id : nodes_[match_id]->getAlignedIds()) {
                    if (nodes_[id]->base() == base) {
                        found_node_id = id;
                        break;
                    }
                }

                if (found_node_id == -1) {
                    // we didn't find aligned node with same base
                    addNode(base);
                    node_id = next_id_ - 1;
                    // add all aligned to nodes to newly created node
                    for (auto id : nodes_[match_id]->getAlignedIds()) {
                        nodes_[node_id]->addAlignedNode(id);
                    }
                    nodes_[node_id]->addAlignedNode(match_id);

                    // to nodes aligned to newly created node add this node
                    // as aligned to
                    for (auto id : nodes_[node_id]->getAlignedIds()) {
                        nodes_[id]->addAlignedNode(node_id);
                    }
                } else {
                    // node id is found node id
                    node_id = found_node_id;
                }
            }

            if (head_id != -1 && node_id != -1) {
                addEdge(head_id, node_id, label);
            }
            head_id = node_id;
            if (first_id == -1) {
                first_id = head_id;
            }
        }

        // connect aligned part with unaligned suffix
        if (head_id != -1 && tail_id != -1) {
            addEdge(head_id, tail_id, label);
        }

        // resort nodes order
        topological_sort();

        sequences_.emplace_back(seq);
        labels_.emplace_back(label);
        start_ids_.emplace_back(first_id);
    }


    void Graph::DFS(uint32_t start_id) {
        unordered_set<uint32_t> started;
        stack<uint32_t> s;
        s.emplace(start_id);
        while (!s.empty()) {
            auto id = s.top();
            s.pop();

            if (visited_.find(id) != visited_.end()) {
                continue;
            }

            if (started.find(id) != started.end()) {
                visited_.insert(id);
                nodes_order_.emplace_back(id);
                started.erase(id);
                continue;
            }

            // necessary because i want to put node
            // on top of order stack only when all successors
            // of this node are processed
            started.insert(id);
            s.emplace(id);
            // push all successors on top of stack
            for (auto next_id : getNode(id)->getSuccessorsIds()) {
                if (visited_.find(next_id) == visited_.end()) {
                    s.emplace(next_id);
                }
            }
        }
    }


    void Graph::topological_sort() {
        visited_.clear();
        nodes_order_.clear();
        for (auto& node : nodes_) {
            if (visited_.find(node->id()) == visited_.end()) {
                DFS(node->id());
            }
        }
        // expensive - added linear time complexity
        // figure out a better way
        reverse(nodes_order_.begin(), nodes_order_.end());
        is_sorted = true;

        // calculate max distances
        calc_nodes_distances();
    }


    // method finds maximum-weight path in graph
    void Graph::generate_consensus(string *pconsensus) {
        auto& consensus = *pconsensus;

        vector<int> next_node(nodes_.size(), -1);
        vector<int> dp(nodes_.size(), -1);
        int max_weight = numeric_limits<int>::min();
        int start_node_id = -1;

        // topological sort on reverse graph
        if (!is_sorted) {
            topological_sort();
        }
        reverse(nodes_order_.begin(), nodes_order_.end());

        // now no node is visited before all its children are visited
        // when iterating over nodes
        for (auto id : nodes_order_) {
            auto& node = nodes_[id];

            const Edges& out_edges = node->getOutEdges();
            if (out_edges.empty()) {
                // if the node has no outgoing edges, weight of heaviest
                // path starting at that node is zero: d(u) = 0
                dp[id] = 0;
                next_node[id] = id;
            } else {
                // otherwise, for each outgoing edge(u, v) compute w(e) + d(v)
                // and set d(u) to be the largest value attained this way
                dp[id] = 0;
                int max_edge = -1;
                for (const auto& e: out_edges) {
                    int weight = e.second->getLabels().size();
                    /*if (weight + dp[e.first] > dp[id]) {
                        dp[id] = weight + dp[e.first];
                        next_node[id] = e.first;
                    }*/
                    if (weight > dp[id] || (weight == dp[id] && dp[e.first] > dp[id])) {
                        dp[id] = weight;
                        max_edge = e.first;
                    }
                }
                dp[id] += dp[max_edge];
                next_node[id] = max_edge;

                // update max
                if (dp[id] > max_weight) {
                    max_weight = dp[id];
                    start_node_id = id;
                }
            }
        }

        // generate consensus sequence
        int curr_id = start_node_id;
        while (curr_id != next_node[curr_id]) {
            consensus += nodes_[curr_id]->base();
            curr_id = next_node[curr_id];
        }
        consensus += nodes_[curr_id]->base();
    }


    void Graph::alignment_strings() {
        unordered_map<int, int> node_to_colID;
        int curr_col = 0;

        const vector<uint32_t>& node_ids = getNodesIds();
        for (auto& id : node_ids) {
            int found_idx = -1;
            for (auto aligned_id : nodes_[id]->getAlignedIds()) {
                if (node_to_colID.find(aligned_id) != node_to_colID.end()) {
                    found_idx = node_to_colID[aligned_id];
                    break;
                }
            }
            if (found_idx == -1) {
                found_idx = curr_col;
                ++curr_col;
            }
            node_to_colID[id] = found_idx;
        }

        int num_columns = curr_col;

        vector<string> seq_names;
        vector<string> aln_strings;
        for (size_t i = 0; i < labels_.size(); ++i) {
            auto& label = labels_[i];
            int curr_node_id = start_ids_[i];
            string aln_string(num_columns, '-');
            while (curr_node_id != -1) {
                int col_idx = node_to_colID[curr_node_id];
                aln_string[col_idx] = nodes_[curr_node_id]->base();
                curr_node_id = nodes_[curr_node_id]->nextNode(label);
            }
            std::cout << aln_string << std::endl;
        }
    }

    void Graph::calc_nodes_distances() {
      assert(is_sorted == true);

      max_node_distance_.clear();
      max_node_distance_.resize(nodes_order_.size(), 0);

      min_node_distance_.clear();
      min_node_distance_.resize(nodes_order_.size(), 0);
      for (auto curr_id : nodes_order_) {
          min_node_distance_[curr_id] = getNode(curr_id)->getPredecessorsIds().size() > 0 ? numeric_limits<int>::max() : 0;
          for (auto prev_id : getNode(curr_id)->getPredecessorsIds()) {
              max_node_distance_[curr_id] = max(max_node_distance_[prev_id] + 1, max_node_distance_[curr_id]);
              min_node_distance_[curr_id] = min(min_node_distance_[prev_id] + 1, min_node_distance_[curr_id]);
          }
      }
    }

    int Graph::max_node_distance(const int node_id) const {
      return max_node_distance_[node_id];
    }

    int Graph::min_node_distance(const int node_id) const {
      return min_node_distance_[node_id];
    }
}
