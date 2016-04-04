/*!
 * @file poa.cpp
 *
 * @brief Poa source file which encapsulates the implementation
 */

#include <stdlib.h>
#include <stdint.h>
#include <algorithm>

#include "graph.hpp"
#include "spoa.hpp"

void prepare_indices(std::vector<uint32_t>& dst, const std::vector<std::string>& sequences, bool sort) {
    dst.resize(sequences.size());
    std::iota(dst.begin(), dst.end(), static_cast<uint32_t>(0));

    if (sort) {
        std::sort(dst.begin(), dst.end(),
            [&](uint32_t lhs, uint32_t rhs) {
                return sequences[lhs].size() > sequences[rhs].size();
            }
        );
    }
}

namespace SPOA {

std::string generate_consensus(const std::vector<std::string>& sequences,
    AlignmentParams params, bool sorted) {

    std::vector<uint32_t> indices;
    prepare_indices(indices, sequences, sorted);

    std::shared_ptr<Graph> graph = createGraph(sequences[indices.front()]);
    graph->topological_sort();

    for (uint32_t i = 1; i < sequences.size(); ++i) {
        auto alignment = createAlignment(sequences[indices[i]], graph, params);
        alignment->align_sequence_to_graph();
        alignment->backtrack();
        graph->add_alignment(std::move(alignment), sequences[indices[i]]);
    }

    return graph->generate_consensus();
}

std::string generate_consensus(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, AlignmentParams params, bool sorted) {

    std::vector<uint32_t> indices;
    prepare_indices(indices, sequences, sorted);

    std::shared_ptr<Graph> graph = createGraph(sequences[indices.front()], qualities[indices.front()]);
    graph->topological_sort();

    for (uint32_t i = 1; i < sequences.size(); ++i) {
        auto alignment = createAlignment(sequences[indices[i]], graph, params);
        alignment->align_sequence_to_graph();
        alignment->backtrack();
        graph->add_alignment(std::move(alignment), sequences[indices[i]], qualities[indices[i]]);
    }

    return graph->generate_consensus();
}

}
