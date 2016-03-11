/*!
 * @file node.hpp
 *
 * @brief Node class header file
 */

#pragma once

#include <vector>
#include <memory>

class Node;
std::unique_ptr<Node> createNode(uint32_t id, char letter);

class Edge;
using EdgeSharedPtr = std::shared_ptr<Edge>;

class Node {
public:

    ~Node();

    uint32_t id() const {
        return id_;
    }

    char letter() const {
        return letter_;
    }

    const std::vector<EdgeSharedPtr>& in_edges() const {
        return in_edges_;
    }

    void add_in_edge(EdgeSharedPtr edge) {
        in_edges_.push_back(edge);
    }

    const std::vector<EdgeSharedPtr>& out_edges() const {
        return out_edges_;
    }

    void add_out_edge(EdgeSharedPtr edge) {
        out_edges_.emplace_back(edge);
    }

    const std::vector<uint32_t>& aligned_nodes_ids() const {
        return aligned_nodes_ids_;
    }

    void add_aligned_node_id(uint32_t id) {
        aligned_nodes_ids_.emplace_back(id);
    }

    friend std::unique_ptr<Node> createNode(uint32_t id, char letter);

private:

    Node(uint32_t id, char letter);
    Node(const Node&) = delete;
    const Node& operator=(const Node&) = delete;

    uint32_t id_;
    char letter_;

    std::vector<EdgeSharedPtr> in_edges_;
    std::vector<EdgeSharedPtr> out_edges_;

    std::vector<uint32_t> aligned_nodes_ids_;
};
