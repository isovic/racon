/**
 * @file node.hpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Node class header file 
 * @details Header file with declaration of Node class used for
 * storing nodes in graph. Class is based on
 * https://github.com/ljdursi/poapy/blob/master/poagraph.py
 * Node class python implementation
 */
#ifndef NODE_H
#define NODE_H

#include <string>
#include <unordered_map>
#include <map>
#include <memory>
#include <vector>

#include "./graph.hpp"

using std::string;
using std::unordered_map;
using std::map;
using std::shared_ptr;
using std::vector;

namespace POA {

    class Graph;
    class Edge;

    typedef map<uint32_t, shared_ptr< Edge >> Edges;


    class Node {
     public:
        /**
         * @brief Node constructor
         * 
         * @param id node id
         * @param base node base(nucleotid)
         * @param graph graph node belongs to
         */
        Node(uint32_t id, char base, const Graph& graph);


        /**
         * @brief Getter for id
         * @return node id
         */
        uint32_t id() const { return id_; }


        /**
         * @brief Getter for base
         * @return node base
         */
        char base() const { return base_; }


        /**
         * @brief Getter for out edges
         * @return All node outgoing edges and
         * edges end ids (@see Edges typedef)
         */
        const Edges& getOutEdges() const { return out_edges; }


        /**
         * @brief Getter for node in edges
         * @return All node ingoing edges and
         * edges start ids (@see Edges typedef)
         */
        const Edges& getInEdges() const { return in_edges; }


        /**
         * @brief Getter for node predecessors
         * @details Node predecessors are all nodes
         * which outgoing edges end in this node. Equivalently
         * said, these are start nodes of all ingoing edges
         * of this node
         * @return Predecessors ids
         */
        const vector<uint32_t>& getPredecessorsIds() const;


        /**
         * @brief Getter for node successors
         * @details Node successors are all nodes
         * which ingoing edges start in this node. Equivalently
         * said, these are end nodes of all outgoing edges
         * of this node
         * @return [description]
         */
        const vector<uint32_t>& getSuccessorsIds() const;


        /**
         * @brief Getter for aligned nodes ids
         * @details Aligned nodes are nodes which represent
         * same position in alignment
         * @return Ids of aligned nodes
         */
        const vector<uint32_t>& getAlignedIds() const;


        /**
         * @brief Add ingoing edge to node
         * 
         * @param A_id start node id
         * @param label sequence label associated with edge
         */
        void addInEdge(uint32_t A_id, const string& label);


        /**
         * @brief Add outgoind edge to node
         * 
         * @param B_id end node id
         * @param label sequence label associated with edge
         */
        void addOutEdge(uint32_t B_id, const string& label);


        /**
         * @brief Add aligned node to this node
         * 
         * @param id aligned node id
         */
        void addAlignedNode(uint32_t id);


        /**
         * @brief Getter for next node id
         * @details Next node is end node of outgoging
         * edge associated with sequence label given
         * as parameter
         * 
         * @param label sequence label
         * @return next node id if exists, -1 otherwise
         */
        int nextNode(const string& label);

     private:
        // node id
        uint32_t id_;
        // node base
        char base_;
        // graph node belongs to
        const Graph& graph_;
        // ingoing edges
        Edges in_edges;
        // outgoing edges
        Edges out_edges;

        // ids of predecessors
        vector<uint32_t> predecessors;
        // ids of successsors
        vector<uint32_t> successors;
        // ids of aligned nodes
        vector<uint32_t> aligned_to_ids;
    };
}

#endif  // NODE_H
