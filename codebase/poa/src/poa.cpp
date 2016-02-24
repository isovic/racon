// Copyright 2015 @mculinovic
#include <vector>
#include <string>

#include "./graph.hpp"
#include "./alignment.hpp"
#include "./poa.hpp"

using std::string;
using std::to_string;
using std::vector;

namespace POA {

    string poa_consensus(const vector<string>& sequences) {
        if (sequences.empty()) {
            return "";
        }

        Graph graph(sequences[0], "seq0");

        for (size_t i = 1; i < sequences.size(); ++i) {
            Alignment aln(const_cast<string&>(sequences[i]), graph);
            aln.align();
            graph.insertSequenceAlignment(aln, sequences[i], "seq" + to_string(i));
        }

        string consensus;
        graph.generate_consensus(&consensus);
        return consensus;
    }
}
