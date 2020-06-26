/*!
 * @file util.cpp
 *
 * @brief Utility source file.
 */

#include "util.hpp"

namespace racon {

// Tuple of: <start, end, window_id>
// using WindowInterval = std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>;

// std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> find_breaking_points_from_cigar(
//             const std::string& cigar, int32_t q_start, int32_t t_start,
//             std::vector<IntervalInt64> windows)
// {
// }

/*
 * \param windows Vector of tuples, where each tuple element is: <start, end, window_id>.
*/
std::vector<WindowInterval> generate_window_breakpoints(
            const std::string& cigar, int64_t q_start, int64_t t_start,
            std::vector<std::tuple<int64_t, int64_t, int64_t>> windows)
{
    std::vector<WindowInterval> ret;

    if (windows.empty()) {
        return ret;
    }

    std::sort(windows.begin(), windows.end(),
                [](const std::tuple<int64_t, int64_t, int64_t>& a,
                    const std::tuple<int64_t, int64_t, int64_t>& b) {
        return std::get<0>(a) < std::get<0>(b);
    });

    // Avoid using the operators.
    const char* cigar_str = cigar.c_str();
    int64_t num_windows = static_cast<int64_t>(windows.size());

    bool found_first_match = false;
    std::tuple<int64_t, int64_t, int64_t> first_match(0, 0, 0);
    std::tuple<int64_t, int64_t, int64_t> last_match(0, 0, 0);

    // The "-1" is because of the while loops below (increment is done at the top).
    int64_t q_pos = q_start - 1;
    int64_t t_pos = t_start - 1;
    int64_t curr_w = 0;

    while (curr_w < num_windows && std::get<1>(windows[curr_w]) <= t_start) {
        ++curr_w;
    }

    for (uint32_t i = 0, j = 0; i < cigar.size(); ++i) {
        if (cigar_str[i] == 'M' || cigar_str[i] == '=' || cigar_str[i] == 'X') {
            uint32_t num_bases = atoi(&cigar_str[j]);
            j = i + 1;
            for (uint32_t k = 0; k < num_bases; ++k) {
                ++q_pos;
                ++t_pos;

                while (curr_w < num_windows && std::get<1>(windows[curr_w]) <= t_pos) {
                    ++curr_w;
                }
                if (curr_w >= num_windows) {
                    break;
                }

                const auto& win_start = std::get<0>(windows[curr_w]);
                const auto& win_end = std::get<1>(windows[curr_w]);
                const auto& win_id = std::get<2>(windows[curr_w]);

                if (t_pos < win_start) {
                    continue;
                }

                if (!found_first_match) {
                    found_first_match = true;
                    first_match = std::make_tuple(t_pos, q_pos, win_id);
                }
                last_match = std::make_tuple(t_pos + 1, q_pos + 1, win_id);
                if (t_pos == (win_end - 1)) {
                    if (found_first_match) {
                        WindowInterval wi(
                                std::get<1>(first_match),
                                std::get<1>(last_match),
                                std::get<0>(first_match),
                                std::get<0>(last_match),
                                std::get<2>(last_match));
                        ret.emplace_back(std::move(wi));
                    }
                    found_first_match = false;
                    ++curr_w;
                    if (curr_w >= num_windows) {
                        break;
                    }
                }
            }
        } else if (cigar_str[i] == 'I') {
            q_pos += atoi(&cigar_str[j]);
            j = i + 1;
        } else if (cigar_str[i] == 'D' || cigar_str[i] == 'N') {
            uint32_t num_bases = atoi(&cigar_str[j]);
            j = i + 1;
            for (uint32_t k = 0; k < num_bases; ++k) {
                ++t_pos;

                while (curr_w < num_windows && std::get<1>(windows[curr_w]) <= t_pos) {
                    ++curr_w;
                }
                if (curr_w >= num_windows) {
                    break;
                }

                const auto& win_start = std::get<0>(windows[curr_w]);
                const auto& win_end = std::get<1>(windows[curr_w]);

                if (t_pos < win_start) {
                    continue;
                }

                if (t_pos == (win_end - 1)) {
                    if (found_first_match) {
                        WindowInterval wi(
                                std::get<1>(first_match),
                                std::get<1>(last_match),
                                std::get<0>(first_match),
                                std::get<0>(last_match),
                                std::get<2>(last_match));
                        ret.emplace_back(std::move(wi));
                    }
                    found_first_match = false;
                    ++curr_w;
                    if (curr_w >= num_windows) {
                        break;
                    }
                }
            }
        } else if (cigar_str[i] == 'S' || cigar_str[i] == 'H' || cigar_str[i] == 'P') {
            j = i + 1;
        }
    }

    // Add the final piece.
    if (curr_w < static_cast<int64_t>(windows.size())) {
        const auto& win_start = std::get<0>(windows[curr_w]);
        const auto& win_end = std::get<1>(windows[curr_w]);
        const auto& final_t_start = std::get<0>(first_match);
        const auto& final_t_end = std::get<0>(last_match);
        if (found_first_match && final_t_start >= win_start && final_t_end <= win_end) {
            WindowInterval wi(
                    std::get<1>(first_match),
                    std::get<1>(last_match),
                    std::get<0>(first_match),
                    std::get<0>(last_match),
                    std::get<2>(last_match));
            ret.emplace_back(std::move(wi));
        }
    }

    return ret;
}

}
