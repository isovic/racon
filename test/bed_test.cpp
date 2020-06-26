/*!
 * @file bed_test.cpp
 *
 * @brief Unit tests for the BED parser.
 */

#include "gtest/gtest.h"
#include "racon_test_config.h"
#include "bed.hpp"

#include <string>
#include <tuple>
#include <vector>

TEST(BedFile, DeserializeTests) {
    // Tuple: test_name, input line, expected record, expected return value, expected to throw.
    using TestTupleType = std::tuple<std::string, std::string, racon::BedRecord, bool, bool>;
    std::vector<std::tuple<std::string, std::string, racon::BedRecord, bool, bool>> test_data{
        TestTupleType("Empty input", "", racon::BedRecord(), false,false),
        TestTupleType("Three columns", "chr01 0 1000", racon::BedRecord("chr01", 0, 1000), true, false),
        TestTupleType("Multiple columns", "chr01 1000 2000 some other columns that are ignored", racon::BedRecord("chr01", 1000, 2000), true, false),
        TestTupleType("Invalid BED line", "bad_bed", racon::BedRecord(), true, true),
    };

    for (const auto& single_test: test_data) {
        const std::string& test_name = std::get<0>(single_test);
        const std::string& in_line = std::get<1>(single_test);
        const racon::BedRecord& expected_record = std::get<2>(single_test);
        const bool expected_rv = std::get<3>(single_test);
        const bool expected_throw = std::get<4>(single_test);

        SCOPED_TRACE(test_name);

        if (expected_throw) {
            EXPECT_THROW(
                {
                    racon::BedRecord record;
                    racon::BedFile::Deserialize(in_line, record);
                },
                std::runtime_error);
        } else {
            racon::BedRecord record;
            const bool rv = racon::BedFile::Deserialize(in_line, record);
            EXPECT_EQ(expected_record, record);
            EXPECT_EQ(expected_rv, rv);
        }

    }
}

TEST(BedReader, AllTests) {
    // Tuple: test_name, input line, expected record, expected return value.
    using TestTupleType = std::tuple<std::string, std::string, std::vector<racon::BedRecord>, bool>;
    std::vector<std::tuple<std::string, std::string, std::vector<racon::BedRecord>, bool>> test_data {
        TestTupleType("Empty input", "", {}, false),

        TestTupleType(
            "Normal BED file",
            R"(chr01 0 1000
chr02 1000 2000
chr03 2000 3000
)",
            std::vector<racon::BedRecord>{
                racon::BedRecord("chr01", 0, 1000),
                racon::BedRecord("chr02", 1000, 2000),
                racon::BedRecord("chr03", 2000, 3000),
            },
            false
        ),
        TestTupleType(
            "Normal BED file",
            R"(chr01 0 1000
bad_line
)",
            std::vector<racon::BedRecord>(), true
        ),

    };

    for (const auto& single_test: test_data) {
        const std::string& test_name = std::get<0>(single_test);
        const std::string in_line = std::get<1>(single_test);
        const std::vector<racon::BedRecord>& expected_records = std::get<2>(single_test);
        const bool expected_throw = std::get<3>(single_test);

        SCOPED_TRACE(test_name);

        std::istringstream iss(in_line);

        if (expected_throw) {
            EXPECT_THROW(
                {
                    std::vector<racon::BedRecord> records = racon::BedReader::ReadAll(iss);
                },
                std::runtime_error);
        } else {
            std::vector<racon::BedRecord> records = racon::BedReader::ReadAll(iss);
            EXPECT_EQ(expected_records, records);
        }
    }
}
