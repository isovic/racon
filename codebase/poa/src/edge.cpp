/**
 * @file edge.cpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Edge class implementation file
 * @details Implementation file for Edge class which
 * represents edge in directed acyclic graph. Class is
 * based on https://github.com/ljdursi/poapy/blob/master/poagraph.py
 * python implemented edge class.
 */
#include "./edge.hpp"
#include <string>
#include <vector>
#include <memory>

using std::shared_ptr;
using std::vector;
using std::string;

namespace POA {

    Edge::Edge(uint32_t id_A, uint32_t id_B,
               const string& label, const Graph& graph): graph_(graph),
                                                         A_(graph.getNode(id_A)),
                                                         B_(graph.getNode(id_B)) {
        labels.emplace_back(label);
    }

    void Edge::addLabel(const string& label) {
        // check for duplicates, if so use set instead of vector
        labels.emplace_back(label);
    }
}
