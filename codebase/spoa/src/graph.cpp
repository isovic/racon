/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <set>
#include <algorithm>
#include <math.h>

#include "node.hpp"
#include "edge.hpp"
#include "alignment.hpp"
#include "graph.hpp"

std::unique_ptr<Graph> createGraph(const std::string& sequence, float weight) {
    std::vector<float> weights(sequence.size(), weight);
    return createGraph(sequence, weights);
}

std::unique_ptr<Graph> createGraph(const std::string& sequence, const std::string& quality) {
    std::vector<float> weights;
    for (const auto& q: quality) {
        weights.emplace_back(1. - pow(10., (float) (q - 33) / (-10.)));
    }
    return createGraph(sequence, weights);
}

std::unique_ptr<Graph> createGraph(const std::string& sequence, const std::vector<float>& weights) {
    assert(sequence.size() != 0);
    assert(sequence.size() == weights.size());
    return std::unique_ptr<Graph>(new Graph(sequence, weights));
}


Graph::Graph(const std::string& sequence, const std::vector<float>& weights) :
        num_sequences_(), num_nodes_(), nodes_(), alphabet_(), is_sorted_(false),
        sorted_nodes_ids_(), sequences_start_nodes_ids_(), consensus_() {

    for (const auto& c: sequence) {
        alphabet_.insert(c);
    }

    int32_t start_node_id = this->add_sequence(sequence, weights, 0, sequence.size());

    sequences_start_nodes_ids_.emplace_back(start_node_id);
    ++num_sequences_;
}

Graph::~Graph() {
}

uint32_t Graph::add_node(char letter) {
    nodes_.emplace_back(createNode(num_nodes_, letter));
    return num_nodes_++;
}

void Graph::add_edge(uint32_t begin_node_id, uint32_t end_node_id, float weight) {
    assert(begin_node_id < num_nodes_ && end_node_id < num_nodes_);

    for (const auto& edge: nodes_[begin_node_id]->out_edges()) {
        if (edge->end_node_id() == end_node_id) {
            edge->add_sequence(num_sequences_, weight);
            return;
        }
    }

    std::shared_ptr<Edge> edge = createEdge(begin_node_id, end_node_id, num_sequences_, weight);
    nodes_[begin_node_id]->add_out_edge(edge);
    nodes_[end_node_id]->add_in_edge(edge);
}

void Graph::topological_sort() {

    if (is_sorted_) {
        return;
    }
    sorted_nodes_ids_.clear();

    // 0 - unmarked, 1 - temporarily marked, 2 - permanently marked
    std::vector<uint8_t> marks(num_nodes_, 0);

    uint32_t i = 0;
    while (true) {
        for (; i < num_nodes_; ++i) {
            if (marks[nodes_[i]->id()] == 0) {
                break;
            }
        }
        if (i == nodes_.size()) {
            break;
        }

        this->visit_node(i, marks);
    }

    assert(this->is_topologically_sorted() == true);
    is_sorted_ = true;
}

void Graph::visit_node(uint32_t node_id, std::vector<uint8_t>& marks) {
    assert(marks[node_id] != 1 && "Graph is not a DAG!");

    if (marks[node_id] == 0) {
        marks[node_id] = 1;

        const auto& in_edges = nodes_[node_id]->in_edges();
        for (const auto& edge: in_edges) {
            this->visit_node(edge->begin_node_id(), marks);
        }

        marks[node_id] = 2;
        sorted_nodes_ids_.emplace_back(node_id);
    }
}

bool Graph::is_topologically_sorted() const {
    assert(nodes_.size() == sorted_nodes_ids_.size());

    std::set<uint32_t> visited_nodes;
    for (uint32_t node_id: sorted_nodes_ids_) {
        for (const auto& edge: nodes_[node_id]->in_edges()) {
            if (visited_nodes.count(edge->begin_node_id()) == 0) {
                return false;
            }
        }
        visited_nodes.insert(node_id);
    }

    return true;
}

void Graph::add_alignment(std::shared_ptr<Alignment> alignment, const std::string& sequence,
    float weight) {

    std::vector<float> weights(sequence.size(), weight);
    this->add_alignment(alignment, sequence, weights);
}

void Graph::add_alignment(std::shared_ptr<Alignment> alignment, const std::string& sequence,
    const std::string& quality) {

    std::vector<float> weights;
    for (const auto& q: quality) {
        weights.emplace_back(1. - pow(10., (float) (q - 33) / (-10.)));
    }
    this->add_alignment(alignment, sequence, weights);
}

void Graph::add_alignment(std::shared_ptr<Alignment> alignment, const std::string& sequence,
    const std::vector<float>& weights) {

    assert(sequence.size() != 0);
    assert(sequence.size() == weights.size());

    const auto& node_ids = alignment->node_ids();
    const auto& seq_ids = alignment->seq_ids();

    assert(node_ids.size() == seq_ids.size());

    for (const auto& c: sequence) {
        alphabet_.insert(c);
    }

    if (seq_ids.size() == 0) { // no local alignment!
        int32_t start_node_id = this->add_sequence(sequence, weights, 0, sequence.size());
        ++num_sequences_;
        sequences_start_nodes_ids_.emplace_back(start_node_id);

        is_sorted_ = false;
        this->topological_sort();
        return;
    }

    std::vector<uint32_t> valid_seq_ids;
    for (const auto& id: seq_ids) {
        if (id != -1) {
            valid_seq_ids.emplace_back(id);
        }
    }

    uint32_t tmp = num_nodes_;
    int32_t start_node_id = this->add_sequence(sequence, weights, 0, valid_seq_ids.front());
    int32_t head_node_id = tmp == num_nodes_ ? -1 : num_nodes_ - 1;

    int32_t tail_node_id = this->add_sequence(sequence, weights, valid_seq_ids.back() + 1, sequence.size());

    int32_t new_node_id = -1;

    for (uint32_t i = 0; i < seq_ids.size(); ++i) {
        if (seq_ids[i] == -1) {
            continue;
        }

        char letter = sequence[seq_ids[i]];
        if (node_ids[i] == -1) {
            new_node_id = this->add_node(letter);

        } else {
            auto node = nodes_[node_ids[i]];
            if (node->letter() == letter) {
                new_node_id = node_ids[i];

            } else {
                int32_t aligned_to_node_id = -1;
                for (const auto& aid: node->aligned_nodes_ids()) {
                    if (nodes_[aid]->letter() == letter) {
                        aligned_to_node_id = aid;
                        break;
                    }
                }

                if (aligned_to_node_id == -1) {
                    new_node_id = this->add_node(letter);

                    for (const auto& aid: node->aligned_nodes_ids()) {
                        nodes_[new_node_id]->add_aligned_node_id(aid);
                        nodes_[aid]->add_aligned_node_id(new_node_id);
                    }

                    nodes_[new_node_id]->add_aligned_node_id(node_ids[i]);
                    node->add_aligned_node_id(new_node_id);

                } else {
                    new_node_id = aligned_to_node_id;
                }
            }
        }

        if (start_node_id == -1) {
            start_node_id = new_node_id;
        }

        if (head_node_id != -1) {
            this->add_edge(head_node_id, new_node_id, weights[seq_ids[i]]);
        }

        head_node_id = new_node_id;
    }

    if (tail_node_id != -1) {
        this->add_edge(head_node_id, tail_node_id, weights[valid_seq_ids.back()+1]);
    }

    ++num_sequences_;
    sequences_start_nodes_ids_.emplace_back(start_node_id);

    is_sorted_ = false;
    this->topological_sort();
}

int32_t Graph::add_sequence(const std::string& sequence, const std::vector<float>& weights,
    uint32_t begin, uint32_t end) {

    if (begin == end) {
        return -1;
    }

    assert(begin < sequence.size() && end <= sequence.size());

    int32_t first_node_id = this->add_node(sequence[begin]);

    uint32_t node_id;
    for (uint32_t i = begin + 1; i < end; ++i) {
        node_id = this->add_node(sequence[i]);
        this->add_edge(node_id - 1, node_id, weights[i]);
    }

    return first_node_id;
}

void Graph::generate_msa(std::vector<std::string>& dst) {

    std::vector<int32_t> msa_node_ids(num_nodes_, -1);

    int32_t base_counter = 0;
    for (const auto& node_id: sorted_nodes_ids_) {

        if (nodes_[node_id]->aligned_nodes_ids().size() == 0) {
            msa_node_ids[node_id] = base_counter;
        } else {
            int32_t min_id = base_counter;
            for (const auto& aid: nodes_[node_id]->aligned_nodes_ids()) {
                if (msa_node_ids[aid] != -1 && min_id > msa_node_ids[aid]) {
                    min_id = msa_node_ids[aid];
                }
            }
            msa_node_ids[node_id] = min_id;
        }

        if (msa_node_ids[node_id] == base_counter) {
            ++base_counter;
        }
    }

    for (uint32_t i = 0; i < num_sequences_; ++i) {
        std::string alignment_str(base_counter, '-');
        uint32_t curr_node_id = sequences_start_nodes_ids_[i];

        while (true) {
            alignment_str[msa_node_ids[curr_node_id]] = nodes_[curr_node_id]->letter();

            uint32_t prev_node_id = curr_node_id;
            for (const auto& edge: nodes_[prev_node_id]->out_edges()) {
                for (const auto& label: edge->sequence_labels()) {
                    if (label == i) {
                        curr_node_id = edge->end_node_id();
                        break;
                    }
                }
                if (prev_node_id != curr_node_id) {
                    break;
                }
            }

            if (prev_node_id == curr_node_id) {
                break;
            }
        }

        dst.emplace_back(alignment_str);
    }

    // do the same for consensus sequence
    this->traverse_heaviest_bundle();

    std::string alignment_str(base_counter, '-');
    for (const auto& id: consensus_) {
        alignment_str[msa_node_ids[id]] = nodes_[id]->letter();
    }
    dst.emplace_back(alignment_str);
}

std::string Graph::generate_consensus() {

    this->traverse_heaviest_bundle();
    std::string consensus_str = "";
    for (const auto& node_id: consensus_) {
        consensus_str += nodes_[node_id]->letter();
    }

    return consensus_str;
}

void Graph::traverse_heaviest_bundle() {

    this->topological_sort();

    std::vector<int32_t> predecessors(num_nodes_, -1);
    std::vector<float> scores(num_nodes_, 0);

    uint32_t max_score_id = 0;
    for (const auto& id: sorted_nodes_ids_) {
        for (const auto& edge: nodes_[id]->in_edges()) {
            if (scores[id] < edge->total_weight() ||
                (scores[id] == edge->total_weight() &&
                scores[predecessors[id]] <= scores[edge->begin_node_id()])) {

                scores[id] = edge->total_weight();
                predecessors[id] = edge->begin_node_id();
            }
        }

        if (predecessors[id] != -1) {
            scores[id] += scores[predecessors[id]];
        }

        if (scores[max_score_id] < scores[id]) {
            max_score_id = id;
        }
    }

    if (nodes_[max_score_id]->out_edges().size() != 0) {
        fprintf(stderr, "MAX SCORE HAS OUT EDGES (%d)!\n", max_score_id);

        std::vector<uint32_t> node_id_to_rank(num_nodes_, 0);
        for (uint32_t i = 0; i < num_nodes_; ++i) {
            node_id_to_rank[sorted_nodes_ids_[i]] = i;
        }

        while (nodes_[max_score_id]->out_edges().size() != 0) {
            max_score_id = this->branch_completion(scores, predecessors,
                node_id_to_rank[max_score_id]);
            fprintf(stderr, "Max score id = %d\n", max_score_id);
        }
    }

    // traceback
    consensus_.clear();
    while (predecessors[max_score_id] != -1) {
        consensus_.emplace_back(max_score_id);
        max_score_id = predecessors[max_score_id];
    }
    consensus_.emplace_back(max_score_id);

    std::reverse(consensus_.begin(), consensus_.end());
}

uint32_t Graph::branch_completion(std::vector<float>& scores,
    std::vector<int32_t>& predecessors, uint32_t rank) {

    uint32_t node_id = sorted_nodes_ids_[rank];
    for (const auto& edge: nodes_[node_id]->out_edges()) {
        for (const auto& o_edge: nodes_[edge->end_node_id()]->in_edges()) {
            if (o_edge->begin_node_id() != node_id) {
                scores[o_edge->begin_node_id()] = -1;
            }
        }
    }

    float max_score = 0;
    uint32_t max_score_id = 0;
    for (uint32_t i = rank + 1; i < sorted_nodes_ids_.size(); ++i) {

        uint32_t id = sorted_nodes_ids_[i];
        scores[id] = -1;
        predecessors[id] = -1;

        for (const auto& edge: nodes_[id]->in_edges()) {
            if (scores[edge->begin_node_id()] == -1) {
                continue;
            }

            if (scores[id] < edge->total_weight() ||
                (scores[id] == edge->total_weight() &&
                scores[predecessors[id]] <= scores[edge->begin_node_id()])) {

                scores[id] = edge->total_weight();
                predecessors[id] = edge->begin_node_id();
            }
        }

        if (predecessors[id] != -1) {
            scores[id] += scores[predecessors[id]];
        }

        if (max_score < scores[id]) {
            max_score = scores[id];
            max_score_id = id;
        }
    }

    return max_score_id;
}

void Graph::print() const {
    printf("digraph %d {\n", num_sequences_);
    printf("    graph [rankdir=LR]\n");
    for (uint32_t i = 0; i < num_nodes_; ++i) {
        printf("    %d [label = \"%d|%c\"]\n", i, i, nodes_[i]->letter());
        for (const auto& edge: nodes_[i]->out_edges()) {
            printf("    %d -> %d [label = \"%.3f\"]\n", i, edge->end_node_id(),
                edge->total_weight());
        }
        for (const auto& aid: nodes_[i]->aligned_nodes_ids()) {
            if (aid > i) {
                printf("    %d -> %d [style = dotted, arrowhead = none]\n", i, aid);
            }
        }
    }
    printf("}\n");
}
