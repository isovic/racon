/*!
 * @file window.cpp
 *
 * @brief Window class source file
 */

#include <algorithm>

#include "window.hpp"

#include "spoa/spoa.hpp"

namespace racon {

std::unique_ptr<Window> createWindow(uint32_t id, uint32_t rank, WindowType type,
    const char* backbone, uint32_t backbone_length, const char* quality,
    uint32_t quality_length) {

    if (backbone_length == 0 || backbone_length != quality_length) {
        fprintf(stderr, "racon::createWindow error: "
            "empty backbone sequence/unequal quality length!\n");
        exit(1);
    }

    return std::unique_ptr<Window>(new Window(id, rank, type, backbone,
        backbone_length, quality, quality_length));
}

Window::Window(uint32_t id, uint32_t rank, WindowType type, const char* backbone,
    uint32_t backbone_length, const char* quality, uint32_t quality_length)
        : id_(id), rank_(rank), type_(type), consensus_(), sequences_(),
        qualities_(), positions_() {

    sequences_.emplace_back(backbone, backbone_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(0, 0);
}

Window::~Window() {
}

void Window::add_layer(const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length, uint32_t begin, uint32_t end) {

    if (sequence_length != quality_length) {
        fprintf(stderr, "racon::Window::add_layer error: "
            "unequal quality size!\n");
        exit(1);
    }
    if (begin >= end || begin > sequences_.front().size() || end > sequences_.front().size()) {
        fprintf(stderr, "racon::Window::add_layer error: "
            "layer begin and end positions are invalid!\n");
        exit(1);
    }

    sequences_.emplace_back(sequence, sequence_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(begin, end);
}

bool Window::generate_consensus(std::shared_ptr<spoa::AlignmentEngine> alignment_engine) {

    if (sequences_.size() < 3) {
        consensus_ = sequences_.front();
        return false;
    }

    auto graph = spoa::createGraph();
    graph->add_alignment(spoa::Alignment(), sequences_.front(), qualities_.front());

    std::vector<uint32_t> rank;
    rank.reserve(sequences_.size());
    for (uint32_t i = 0; i < sequences_.size(); ++i) {
        rank.emplace_back(i);
    }

    std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
        return positions_[lhs].first < positions_[rhs].first; });

    uint32_t offset = 0.01 * sequences_.front().size();
    for (uint32_t j = 1; j < sequences_.size(); ++j) {
        uint32_t i = rank[j];
        if (positions_[i].first < offset && positions_[i].second >
            sequences_.front().size() - offset) {

            auto alignment = alignment_engine->align_sequence_with_graph(
                sequences_[i], graph);
            graph->add_alignment(alignment, sequences_[i], qualities_[i]);
        } else {
            std::vector<int32_t> mapping;
            auto subgraph = graph->subgraph(positions_[i].first,
                positions_[i].second, mapping);
            auto alignment = alignment_engine->align_sequence_with_graph(
                sequences_[i], subgraph);
            subgraph->update_alignment(alignment, mapping);
            graph->add_alignment(alignment, sequences_[i], qualities_[i]);
        }
    }

    std::vector<uint32_t> coverages;
    consensus_ = graph->generate_consensus(coverages);

    if (type_ == WindowType::kTGS) {
        uint32_t average_coverage = (sequences_.size() - 1) / 2;

        int32_t begin = 0, end = consensus_.size() - 1;
        for (; begin < static_cast<int32_t>(consensus_.size()); ++begin) {
            if (coverages[begin] >= average_coverage) {
                break;
            }
        }
        for (; end >= 0; --end) {
            if (coverages[end] >= average_coverage) {
                break;
            }
        }

        if (begin >= end) {
            fprintf(stderr, "racon::Window::generate_consensus error: "
                "window %d %d is broken!\n", id_, rank_);
            exit(1);
        }
        consensus_ = consensus_.substr(begin, end - begin + 1);
    }

    return true;
}

}
