/*!
 * @file util.cpp
 *
 * @brief Utility source file.
 */

#include "util.hpp"

namespace racon {

std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> find_breaking_points_from_cigar(
            const std::string& cigar, std::vector<IntervalInt64> windows)
{
    std::sort(windows.begin(), windows.end(), [](const IntervalInt64& a, const IntervalInt64& b) {
        return a.start < b.start;
    });

    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> ret;

    // // find breaking points from cigar
    // std::vector<int32_t> window_ends;
    // uint32_t w_offset = 0;
    // for (uint32_t i = 0; i < t_end_; i += window_length) {
    //     if (i > t_begin_) {
    //         window_ends.emplace_back(i - 1);
    //     } else {
    //         ++w_offset;
    //     }
    // }
    // window_ends.emplace_back(t_end_ - 1);

    // uint32_t w = 0;
    // bool found_first_match = false;
    // std::tuple<uint32_t, uint32_t, uint32_t> first_match = {0, 0, 0}, last_match = {0, 0, 0};

    // int32_t q_ptr = (strand_ ? (q_length_ - q_end_) : q_begin_) - 1;
    // int32_t t_ptr = t_begin_ - 1;

    // for (uint32_t i = 0, j = 0; i < cigar_.size(); ++i) {
    //     if (cigar_[i] == 'M' || cigar_[i] == '=' || cigar_[i] == 'X') {
    //         uint32_t k = 0, num_bases = atoi(&cigar_[j]);
    //         j = i + 1;
    //         while (k < num_bases) {
    //             ++q_ptr;
    //             ++t_ptr;

    //             if (!found_first_match) {
    //                 found_first_match = true;
    //                 first_match = std::make_tuple(t_ptr, q_ptr, w + w_offset);
    //             }
    //             last_match = std::make_tuple(t_ptr + 1, q_ptr + 1, w + w_offset);
    //             if (t_ptr == window_ends[w]) {
    //                 if (found_first_match) {
    //                     breaking_points_.emplace_back(first_match);
    //                     breaking_points_.emplace_back(last_match);
    //                 }
    //                 found_first_match = false;
    //                 ++w;
    //             }

    //             ++k;
    //         }
    //     } else if (cigar_[i] == 'I') {
    //         q_ptr += atoi(&cigar_[j]);
    //         j = i + 1;
    //     } else if (cigar_[i] == 'D' || cigar_[i] == 'N') {
    //         uint32_t k = 0, num_bases = atoi(&cigar_[j]);
    //         j = i + 1;
    //         while (k < num_bases) {
    //             ++t_ptr;
    //             if (t_ptr == window_ends[w]) {
    //                 if (found_first_match) {
    //                     breaking_points_.emplace_back(first_match);
    //                     breaking_points_.emplace_back(last_match);
    //                 }
    //                 found_first_match = false;
    //                 ++w;
    //             }
    //             ++k;
    //         }
    //     } else if (cigar_[i] == 'S' || cigar_[i] == 'H' || cigar_[i] == 'P') {
    //         j = i + 1;
    //     }
    // }
}

}
