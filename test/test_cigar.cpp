/*!
 * @file test_cigar.cpp
 *
 * @brief Unit tests for the CIGAR tools.
 */

#include "gtest/gtest.h"
#include "racon_test_config.h"
#include "cigar.hpp"

#include <string>
#include <tuple>
#include <vector>

TEST(Cigar, ParseCigarString_AllTests) {
    // Tuple: test_name, input_cigar_string, expected_cigar, expected_throw
    using TestTupleType = std::tuple<std::string, std::string, racon::Cigar, bool>;
    std::vector<TestTupleType> test_data{
        TestTupleType("Empty input", "", racon::Cigar(), false),
        TestTupleType("Single op", "10=", racon::Cigar{{'=', 10}}, false),
        TestTupleType("Multiple ops", "10=1X1I1D10=", racon::Cigar{{'=', 10}, {'X', 1}, {'I', 1}, {'D', 1}, {'=', 10}}, false),
        TestTupleType("Bad CIGAR, no op char", "10", racon::Cigar{}, true),
        TestTupleType("Bad CIGAR, no length value in middle", "10=I", racon::Cigar{}, true),
        TestTupleType("Bad CIGAR, no length value at front", "X", racon::Cigar{}, true),
    };

    for (const auto& single_test: test_data) {
        const std::string& test_name = std::get<0>(single_test);
        const std::string& in_cigar_str = std::get<1>(single_test);
        const racon::Cigar& expected = std::get<2>(single_test);
        const bool expected_throw = std::get<3>(single_test);

        SCOPED_TRACE(test_name);

        if (expected_throw) {
            EXPECT_THROW(
                {
                    racon::Cigar results = racon::ParseCigarString(in_cigar_str);
                },
                std::runtime_error);
        } else {
            racon::Cigar results = racon::ParseCigarString(in_cigar_str);
            EXPECT_EQ(expected, results);
        }
    }
}
