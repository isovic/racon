/*!
 * @file test_cigar.cpp
 *
 * @brief Unit tests for the CIGAR tools.
 */

#include "gtest/gtest.h"
#include "racon_test_config.h"
#include "cigar.hpp"
#include "vcf.hpp"

#include <string>
#include <tuple>
#include <vector>

TEST(VCF, ExtractVCFEventsFromCigarString_AllTests) {
    // Tuple: test_name, in_cigar_string, in_query_seq, in_target_seq, expected_diffs, expected_throw
    using TestTupleType = std::tuple<std::string, std::string, std::string, std::string, bool, std::vector<racon::VcfDiff>>;
    std::vector<TestTupleType> test_data{
        TestTupleType("Empty input", "", "", "", false, {}),
        TestTupleType("Exact match, no diffs", "4=", "ACTG", "ACTG", false, {}),

        /*
         * Note: the single mismatch test includes one match base before the mismatch too.
         * Racon reports that because it's better to be consistent for all diffs.
         * Typically, the preceding base would only be reported for indels, and not for
         * mismatches.
        */
        TestTupleType("Single mismatch", "4=1X4=", "ACTGAACTG", "ACTGCACTG", false,
                        {
                            // pos, ref, alt
                            {3, "GC", "GA"},
                        }),
        TestTupleType("Single insertion.", "4=1I4=", "ACTGAACTG", "ACTGACTG", false,
                        {
                            // pos, ref, alt
                            {3, "G", "GA"},
                        }),
        TestTupleType("Single deletion.", "4=1D4=", "ACTGACTG", "ACTGCACTG", false,
                        {
                            // pos, ref, alt
                            {3, "GC", "G"},
                        }),
        TestTupleType("Multiple non-neighboring events.", "4=1I1=1D1=1X1=", "ACTGAATCG", "ACTGACTGG", false,
                        // R: ACTG-ACTGG
                        // Q: ACTGAA-TCG
                        {
                            // pos, ref, alt
                            {3, "G", "GA"},
                            {4, "AC", "A"},
                            {6, "TG", "TC"},
                        }),
        TestTupleType("Multiple neighboring events.", "4=1I1D1X1=", "ACTGACG", "ACTGCGG", false,
                        // R: ACTG-CGG
                        // Q: ACTGA-CG
                        {
                            // pos, ref, alt
                            {3, "GCG", "GAC"},
                        }),

        TestTupleType("Leading mismatch", "1X4=", "AACTG", "CACTG", false,
                        {
                            // pos, ref, alt
                            {0, "CA", "AA"},
                        }),
        TestTupleType("Leading insertions", "3I4=", "AAACTGG", "CTGG", false,
                        // R: ---CTGG
                        // Q: AAACTGG
                        {
                            // pos, ref, alt
                            {0, "C", "AAAC"},
                        }),
        TestTupleType("Leading deletions", "3D4=", "CTGG", "AAACTGG", false,
                        // R: AAACTGG
                        // Q: ---CTGG
                        {
                            // pos, ref, alt
                            {0, "AAAC", "C"},
                        }),
        TestTupleType("Leading multiple events", "2I3D1X3=", "GGTTGG", "AAACTGG", false,
                        // R: --AAACTGG
                        // Q: GG---TTGG
                        {
                            // pos, ref, alt
                            {0, "AAACT", "GGTT"},
                        }),

        TestTupleType("Trailing mismatch", "4=1X", "ACTGA", "ACTGC", false,
                        {
                            // pos, ref, alt
                            {3, "GC", "GA"},
                        }),
        TestTupleType("Trailing insertion", "4=1I", "ACTGA", "ACTG", false,
                        {
                            // pos, ref, alt
                            {3, "G", "GA"},
                        }),
        TestTupleType("Trailing deletion", "4=1D", "ACTG", "ACTGC", false,
                        {
                            // pos, ref, alt
                            {3, "GC", "G"},
                        }),
        TestTupleType("Trailing multiple events", "4=1X1D1I", "ACTGTC", "ACTGAC", false,
                        // R: ACTGAC-
                        // Q: ACTGT-C
                        {
                            // pos, ref, alt
                            {3, "GAC", "GTC"},
                        }),

        TestTupleType("Invalid input, entire CIGAR is a set of diffs", "1X1D1I", "TC", "AC", true,
                        {}),
        TestTupleType("Invalid input, the sequence lengths do not match the CIGAR", "4=1I", "ACTG", "ACTG", true,
                        {}),
    };

    for (const auto& single_test: test_data) {
        const std::string& test_name = std::get<0>(single_test);
        const std::string& in_cigar_str = std::get<1>(single_test);
        const std::string& in_qseq = std::get<2>(single_test);
        const std::string& in_tseq = std::get<3>(single_test);
        const bool expected_throw = std::get<4>(single_test);
        const std::vector<racon::VcfDiff>& expected = std::get<5>(single_test);

        SCOPED_TRACE(test_name);

        // std::cerr << "TestName: " << test_name << "\n";

        const racon::Cigar cigar = racon::ParseCigarString(in_cigar_str);

        if (expected_throw) {
            EXPECT_THROW(
                {
                    std::vector<racon::VcfDiff> results = racon::ExtractVCFEventsFromCigarString(cigar, in_qseq, in_tseq);
                },
                std::runtime_error);
        } else {
            std::vector<racon::VcfDiff> results = racon::ExtractVCFEventsFromCigarString(cigar, in_qseq, in_tseq);
            // std::cerr << "Results:\n";
            // for (size_t i = 0; i < results.size(); ++i) {
            //     std::cerr << "[i = " << i << "] pos = " << results[i].pos << ", ref: '" << results[i].ref << "', alt: '" << results[i].alt << "'\n";
            // }
            EXPECT_EQ(expected, results);
        }
    }
}
