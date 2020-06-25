/*!
 * @file util.hpp
 *
 * @brief Utility functions used throughout the code.
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <tuple>
#include "IntervalTree.h"

using IntervalTreeInt64 = IntervalTree<int64_t, size_t>;
using IntervalVectorInt64 = IntervalTreeInt64::interval_vector;
using IntervalInt64 = IntervalTreeInt64::interval;

namespace racon {

template<typename T>
bool transmuteId(const std::unordered_map<T, uint64_t>& t_to_id, const T& t,
    uint64_t& id) {

    auto it = t_to_id.find(t);
    if (it == t_to_id.end()) {
        return false;
    }
    id = it->second;
    return true;
}

std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> find_breaking_points_from_cigar(
            const std::string& cigar, std::vector<IntervalInt64> windows);

}
