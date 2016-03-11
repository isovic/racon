/*!
 * @file edge.hpp
 *
 * @brief Edge class header file
 */

#pragma once

#include <assert.h>
#include <stdint.h>
#include <vector>
#include <memory>


class Edge;
std::unique_ptr<Edge> createEdge(uint32_t begin_node_id, uint32_t end_node_id,
    uint32_t label, float weight = 1.0);

class Edge {
public:

    ~Edge();

    uint32_t begin_node_id() const {
        return begin_node_id_;
    }

    uint32_t end_node_id() const {
        return end_node_id_;
    }

    const std::vector<uint32_t>& sequence_labels() const {
        return sequence_labels_;
    }

    const std::vector<float>& sequence_weights() const {
        return sequence_weights_;
    }

    float total_weight() const {
        return total_weight_;
    }

    void add_sequence(uint32_t label, float weight = 1.0);

    friend std::unique_ptr<Edge> createEdge(uint32_t begin_node_id,
        uint32_t end_node_id, uint32_t label, float weight);

private:

    Edge(uint32_t begin_node_id, uint32_t end_node_id, uint32_t label, float weight);
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    uint32_t begin_node_id_;
    uint32_t end_node_id_;
    std::vector<uint32_t> sequence_labels_;
    std::vector<float> sequence_weights_;
    float total_weight_;
};
