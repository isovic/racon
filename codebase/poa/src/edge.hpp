/**
 * @file edge.hpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Edge class header file 
 * @details Header file with declaration of Edge class which
 * represents edge in directed acyclic graph. Class is
 * based on https://github.com/ljdursi/poapy/blob/master/poagraph.py
 * python implemented edge class.
 */
#ifndef EDGE_H
#define EDGE_H

#include <string>
#include <vector>
#include <memory>

#include "./graph.hpp"
#include "./node.hpp"

using std::shared_ptr;
using std::vector;
using std::string;

namespace POA {

    class Graph;
    class Node;

    /** 
     * Class represents edge in DAG
     */
    class Edge {
     public:
        /**
         * @brief Edge constructor
         * 
         * @param id_A start node id
         * @param id_B end node id
         * @param label sequence label
         * @param graph graph containing this edge
         */
        Edge(uint32_t id_A, uint32_t id_B,
             const string& label, const Graph& graph);


        /**
         * @brief Adds label to edge.
         * 
         * @param label sequence label
         */
        void addLabel(const string& label);


        /**
         * @brief Getter for edge labels.
         * @return vector of labels
         */
        const vector<string>& getLabels() { return labels; }

     private:
        // Graph containing this edge
        const Graph& graph_;
        // Edge end nodes
        shared_ptr<Node> A_, B_;
        // Sequence labels - represent sequence which
        // are contained in this edge. Number of labels
        // is weight of this edge.
        vector<string> labels;
    };
}

#endif  // EDGE_H
