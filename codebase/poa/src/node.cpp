/**
 * @file node.cpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Node class implementation file 
 * @details Implemenation file for Node class used for
 * storing nodes in graph. Class is based on
 * https://github.com/ljdursi/poapy/blob/master/poagraph.py
 * Node class python implementation
 */
#include "./node.hpp"
#include <string>
#include <memory>
#include <vector>
#include <algorithm>

using std::string;
using std::shared_ptr;
using std::vector;
using std::sort;

namespace POA {

    Node::Node(uint32_t id, char base, const Graph& graph): id_(id),
                                                            base_(base),
                                                            graph_(graph) {}


    void Node::addInEdge(uint32_t A_id, const string& label) {
        if (in_edges.find(A_id) == in_edges.end()) {
            shared_ptr<Edge> edge(new Edge(A_id, id_, label, graph_));
            in_edges[A_id] = edge;
            predecessors.emplace_back(A_id);
        } else {
            in_edges[A_id]->addLabel(label);
        }
    }


    void Node::addOutEdge(uint32_t B_id, const string& label) {
        if (out_edges.find(B_id) == out_edges.end()) {
            shared_ptr<Edge> edge(new Edge(id_, B_id, label, graph_));
            out_edges[B_id] = edge;
            successors.emplace_back(B_id);
        } else {
            out_edges[B_id]->addLabel(label);
        }
    }

    const vector<uint32_t>& Node::getSuccessorsIds() const {
        return successors;
    }

    const vector<uint32_t>& Node::getPredecessorsIds() const {
        return predecessors;
    }

    const vector<uint32_t>& Node::getAlignedIds() const {
        return aligned_to_ids;
    }


    void Node::addAlignedNode(uint32_t id) {
        aligned_to_ids.emplace_back(id);
    }

    int Node::nextNode(const string& label) {
        int node_id = -1;
        for (auto& e : out_edges) {
            for (auto& l : (e.second)->getLabels()) {
                if (label == l) {
                    node_id = e.first;
                }
            }
        }
        return node_id;
    }
}
