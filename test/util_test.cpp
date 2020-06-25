/*!
 * @file util.cpp
 *
 * @brief Unit tests for the utils.
 */

#include "gtest/gtest.h"
#include "racon_test_config.h"
#include "util.hpp"

#include <iostream>
#include <string>
#include <tuple>
#include <vector>

TEST(Utils, find_breaking_points_from_cigar) {
    // Tuple: test_name, CIGAR string, windows as tuples, expected breakpoints split into windows.
    std::vector<std::tuple<std::string, std::string, int64_t,
        int64_t, std::vector<std::tuple<int64_t, int64_t, int64_t>>, std::vector<racon::WindowInterval>>> test_data = {
        {"Empty input", "", 0, 0, {}, {}},

        {"Simple windowing", "1000M", 0, 0,
            {
                {100, 500, 0},
                {700, 800, 1},
                {900, 950, 2},
            },
            {
                racon::WindowInterval(100, 500, 100, 500, 0),
                racon::WindowInterval(700, 800, 700, 800, 1),
                racon::WindowInterval(900, 950, 900, 950, 2),
            }
        },

        {"Simple windowing with indels", "200=100I300=50D450=", 0, 0,
            {
                {100, 500, 0},
                {500, 700, 1},
                {700, 800, 2},
                {900, 950, 3},
            },
            {
                racon::WindowInterval(100, 600, 100, 500, 0),
                racon::WindowInterval(600, 750, 550, 700, 1),
                racon::WindowInterval(750, 850, 700, 800, 2),
                racon::WindowInterval(950, 1000, 900, 950, 3),
            }
        },

        {"Intra-window alignments", "100=", 0, 250,
            {
                {0, 200, 0},
                {200, 400, 1},
                {400, 600, 2},
                {600, 800, 3},
                {800, 1000, 4},
            },
            {
                racon::WindowInterval(0, 100, 250, 350, 1),
            }
        },

        {"Cross-window alignments", "200=", 0, 250,
            {
                {0, 200, 0},
                {200, 400, 1},
                {400, 600, 2},
                {600, 800, 3},
                {800, 1000, 4},
            },
            {
                racon::WindowInterval(0, 150, 250, 400, 1),
                racon::WindowInterval(150, 200, 400, 450, 2),
            }
        },

        {"Flanking insertions", "10I180=10I", 0, 200,
            {
                {0, 200, 0},
                {200, 400, 1},
                {400, 600, 2},
                {600, 800, 3},
                {800, 1000, 4},
            },
            {
                racon::WindowInterval(10, 190, 200, 380, 1),
            }
        },

        {"Flanking insertions", "10D180=10D", 0, 200,
            {
                {0, 200, 0},
                {200, 400, 1},
                {400, 600, 2},
                {600, 800, 3},
                {800, 1000, 4},
            },
            {
                racon::WindowInterval(0, 180, 210, 390, 1),
            }
        },

        {"Alignment outside of the windows, front", "50=", 0, 50,
            {
                {200, 400, 0},
                {400, 600, 1},
                {600, 800, 2},
                {800, 1000, 3},
            },
            { }
        },

        {"Alignment outside of the windows, back", "50=", 0, 1050,
            {
                {200, 400, 0},
                {400, 600, 1},
                {600, 800, 2},
                {800, 1000, 3},
            },
            { }
        },

        {"Alignment in between windows (not overlapping)", "50=", 0, 450,
            {
                {0, 200, 0},
                {200, 400, 1},
                {600, 800, 3},
                {800, 1000, 4},
            },
            { }
        },

        {"Out of order windows, should still work because internal sort", "200=100I300=50D450=", 0, 0,
            {
                {900, 950, 3},
                {500, 700, 1},
                {100, 500, 0},
                {700, 800, 2},
            },
            {
                racon::WindowInterval(100, 600, 100, 500, 0),
                racon::WindowInterval(600, 750, 550, 700, 1),
                racon::WindowInterval(750, 850, 700, 800, 2),
                racon::WindowInterval(950, 1000, 900, 950, 3),
            }
        },

    };

    for (const auto& single_test: test_data) {
        const std::string& test_name = std::get<0>(single_test);
        const std::string& cigar = std::get<1>(single_test);
        int64_t aln_q_start = std::get<2>(single_test);
        int64_t aln_t_start = std::get<3>(single_test);
        const std::vector<std::tuple<int64_t, int64_t, int64_t>>& windows = std::get<4>(single_test);
        const std::vector<racon::WindowInterval>& expected = std::get<5>(single_test);

        SCOPED_TRACE(test_name);

        auto result = racon::generate_window_breakpoints(cigar, aln_q_start, aln_t_start, windows);

        // std::cerr << test_name << "\n";
        // std::cerr << cigar << "\n";
        // for (const auto& val: result) {
        //     std::cerr << "q_start = " << val.q_start << ", q_end = " << val.q_end
        //     << ", t_start = " << val.t_start << ", t_end = " << val.t_end
        //     << ", window_id = " << val.window_id << "\n";
        //     // std::cerr << val << "\n";
        // }
        // std::cerr << "\n";

        EXPECT_EQ(expected, result);
    }
}
