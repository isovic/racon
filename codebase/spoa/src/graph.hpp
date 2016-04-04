/*!
 * @file graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <assert.h>
#include <memory>
#include <string>
#include <vector>
#include <unordered_set>

class Node;
class Alignment;

class Graph;
std::unique_ptr<Graph> createGraph(const std::string& sequence, float weight = 1.0);
std::unique_ptr<Graph> createGraph(const std::string& sequence, const std::string& quality);
std::unique_ptr<Graph> createGraph(const std::string& sequence, const std::vector<float>& weights);

class Graph {
public:

    ~Graph();

    uint32_t num_sequences() const {
        return num_sequences_;
    }

    const std::shared_ptr<Node> node(uint32_t id) const {
        assert(id < num_nodes_);
        return nodes_[id];
    }

    const std::vector<std::shared_ptr<Node>>& nodes() const {
        return nodes_;
    }

    const std::vector<uint32_t>& sorted_nodes_ids() const {
        assert(is_sorted_ == true);
        return sorted_nodes_ids_;
    }

    const std::unordered_set<uint8_t>& alphabet() const {
        return alphabet_;
    }

    void topological_sort();

    void add_alignment(std::shared_ptr<Alignment> alignment, const std::string& sequence,
        float weight = 1.0);
    void add_alignment(std::shared_ptr<Alignment> alignment, const std::string& sequence,
        const std::string& quality);
    void add_alignment(std::shared_ptr<Alignment> alignment, const std::string& sequence,
        const std::vector<float>& weights);

    void generate_msa(std::vector<std::string>& dst);

    std::string generate_consensus();

    void print() const;

    friend std::unique_ptr<Graph> createGraph(const std::string& sequence,
        const std::vector<float>& weights);

private:

    Graph(const std::string& sequence, const std::vector<float>& weights);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    uint32_t add_node(char letter);

    void add_edge(uint32_t begin_node_id, uint32_t end_node_id, float weight);

    void visit_node(uint32_t node_id, std::vector<uint8_t>& marks);

    bool is_topologically_sorted() const;

    int32_t add_sequence(const std::string& sequence, const std::vector<float>& weights,
        uint32_t begin, uint32_t end);

    void traverse_heaviest_bundle();

    uint32_t branch_completion(std::vector<float>& scores, std::vector<int32_t>& predecessors,
        uint32_t rank);

    uint32_t num_sequences_;
    uint32_t num_nodes_;
    std::vector<std::shared_ptr<Node>> nodes_;

    std::unordered_set<uint8_t> alphabet_;

    bool is_sorted_;
    std::vector<uint32_t> sorted_nodes_ids_;

    std::vector<uint32_t> sequences_start_nodes_ids_;

    std::vector<uint32_t> consensus_;
};
