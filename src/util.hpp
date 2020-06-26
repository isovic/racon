/*!
 * @file util.hpp
 *
 * @brief Utility functions used throughout the code.
 */

#pragma once

#include <algorithm>
#include <cstdint>
#include <ostream>
#include <vector>
#include <unordered_map>
#include <tuple>
#include "IntervalTree.h"

using IntervalTreeInt64 = IntervalTree<int64_t, size_t>;
using IntervalVectorInt64 = IntervalTreeInt64::interval_vector;
using IntervalInt64 = IntervalTreeInt64::interval;

namespace racon {

class WindowInterval {
public:
    WindowInterval() = default;

    WindowInterval(int64_t _q_start, int64_t _q_end, int64_t _t_start, int64_t _t_end, int64_t _window_id)
        : q_start(_q_start)
        , q_end(_q_end)
        , t_start(_t_start)
        , t_end(_t_end)
        , window_id(_window_id)
    {}

    bool operator==(const WindowInterval& rhs) const
    {
        return q_start == rhs.q_start && q_end == rhs.q_end && t_start == rhs.t_start
                    && t_end == rhs.t_end && window_id == rhs.window_id;
    }

    friend std::ostream& operator<<(::std::ostream& os, const WindowInterval& a);

    int64_t q_start = 0;
    int64_t q_end = 0;
    int64_t t_start = 0;
    int64_t t_end = 0;
    int64_t window_id = -1;
};

inline std::ostream& operator<<(::std::ostream& os, const WindowInterval& a) {
    os << "q_start = " << a.q_start << ", q_end = " << a.q_end
        << ", t_start = " << a.t_start << ", t_end = " << a.t_end
        << ", window_id = " << a.window_id;
    return os;
}

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

std::vector<WindowInterval> generate_window_breakpoints(
            const std::string& cigar, int64_t q_start, int64_t t_start,
            std::vector<std::tuple<int64_t, int64_t, int64_t>> windows);

}
