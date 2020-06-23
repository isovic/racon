/*!
 * @file util.hpp
 *
 * @brief Utility functions used throughout the code.
 */

#pragma once

#include <cstdint>
#include <unordered_map>

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

}
