/*!
 * @file racon_test.cpp
 *
 * @brief Racon unit test source file
 */

#include "racon_test_config.h"

#include "sequence.hpp"
#include "polisher.hpp"

#include "edlib.h"
#include "bioparser/bioparser.hpp"
#include "gtest/gtest.h"

uint32_t calculateEditDistance(const std::string& query, const std::string& target) {

    EdlibAlignResult result = edlibAlign(query.c_str(), query.size(), target.c_str(),
        target.size(), edlibDefaultAlignConfig());

    uint32_t edit_distance = result.editDistance;
    edlibFreeAlignResult(result);

    return edit_distance;
}

class RaconPolishingTest: public ::testing::Test {
public:
    void SetUp(const std::string& sequences_path, const std::string& overlaps_path,
        const std::string& target_path, racon::PolisherType type,
        uint32_t window_length, double quality_threshold, double error_threshold,
        int8_t match,int8_t mismatch, int8_t gap) {

        polisher = createPolisher(sequences_path, overlaps_path, target_path,
            type, window_length, quality_threshold, error_threshold, match,
            mismatch, gap, 1);
    }

    void TearDown() {}

    void initialize() {
        polisher->initialize();
    }

    void polish(std::vector<std::unique_ptr<racon::Sequence>>& dst,
        bool drop_unpolished_sequences) {

        return polisher->polish(dst, drop_unpolished_sequences);
    }

    std::unique_ptr<racon::Polisher> polisher;
};

TEST(RaconTest, PolisherTypeError) {
    EXPECT_DEATH((createPolisher("", "", "", static_cast<racon::PolisherType>(3),
        0, 0, 0, 0, 0, 0, 0)), ".racon::createPolisher. error: invalid polisher"
        " type!");
}

TEST(RaconTest, WindowLengthError) {
    EXPECT_DEATH((createPolisher("", "", "", racon::PolisherType::kC, 0, 0, 0, 0,
        0, 0, 0)), ".racon::createPolisher. error: invalid window length!");
}

TEST(RaconTest, SequencesPathExtensionError) {
    EXPECT_DEATH((createPolisher("", "", "", racon::PolisherType::kC, 500, 0,
        0, 0, 0, 0, 0)), ".racon::createPolisher. error: file  has unsupported "
        "format extension .valid extensions: .fasta, .fa, .fastq, .fq.!");
}

TEST(RaconTest, OverlapsPathExtensionError) {
    EXPECT_DEATH((createPolisher(racon_test_data_path + "sample_reads.fastq", "",
        "", racon::PolisherType::kC, 500, 0, 0, 0, 0, 0, 0)),
        ".racon::createPolisher. error: file  has unsupported format extension "
        ".valid extensions: .mhap, .paf, .sam.!");
}

TEST(RaconTest, TargetPathExtensionError) {
    EXPECT_DEATH((createPolisher(racon_test_data_path + "sample_reads.fastq",
        racon_test_data_path + "sample_overlaps.paf", "", racon::PolisherType::kC,
        500, 0, 0, 0, 0, 0, 0)), ".racon::createPolisher. error: file  has "
        "unsupported format extension .valid extensions: .fasta, .fa, .fastq, .fq.!");
}

TEST_F(RaconPolishingTest, ConsensusWithQualities) {
    SetUp(racon_test_data_path + "sample_reads.fastq", racon_test_data_path +
        "sample_overlaps.paf", racon_test_data_path + "sample_layout.fasta",
        racon::PolisherType::kC, 500, 10, 0.3, 5, -4, -8);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, true);
    EXPECT_EQ(polished_sequences.size(), 1);

    polished_sequences[0]->create_reverse_complement();

    auto parser = bioparser::createParser<bioparser::FastaParser, racon::Sequence>(
        racon_test_data_path + "sample_reference.fasta");
    parser->parse_objects(polished_sequences, -1);
    EXPECT_EQ(polished_sequences.size(), 2);

    EXPECT_EQ(calculateEditDistance(polished_sequences[0]->reverse_complement(),
        polished_sequences[1]->data()), 1312);
}

TEST_F(RaconPolishingTest, ConsensusWithoutQualities) {
    SetUp(racon_test_data_path + "sample_reads.fasta", racon_test_data_path +
        "sample_overlaps.paf", racon_test_data_path + "sample_layout.fasta",
        racon::PolisherType::kC, 500, 10, 0.3, 5, -4, -8);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, true);
    EXPECT_EQ(polished_sequences.size(), 1);

    polished_sequences[0]->create_reverse_complement();

    auto parser = bioparser::createParser<bioparser::FastaParser, racon::Sequence>(
        racon_test_data_path + "sample_reference.fasta");
    parser->parse_objects(polished_sequences, -1);
    EXPECT_EQ(polished_sequences.size(), 2);

    EXPECT_EQ(calculateEditDistance(polished_sequences[0]->reverse_complement(),
        polished_sequences[1]->data()), 1566);
}

TEST_F(RaconPolishingTest, ConsensusWithQualitiesAndAlignments) {
    SetUp(racon_test_data_path + "sample_reads.fastq", racon_test_data_path +
        "sample_overlaps.sam", racon_test_data_path + "sample_layout.fasta",
        racon::PolisherType::kC, 500, 10, 0.3, 5, -4, -8);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, true);
    EXPECT_EQ(polished_sequences.size(), 1);

    polished_sequences[0]->create_reverse_complement();

    auto parser = bioparser::createParser<bioparser::FastaParser, racon::Sequence>(
        racon_test_data_path + "sample_reference.fasta");
    parser->parse_objects(polished_sequences, -1);
    EXPECT_EQ(polished_sequences.size(), 2);

    EXPECT_EQ(calculateEditDistance(polished_sequences[0]->reverse_complement(),
        polished_sequences[1]->data()), 1317);
}

TEST_F(RaconPolishingTest, ConsensusWithoutQualitiesAndWithAlignments) {
    SetUp(racon_test_data_path + "sample_reads.fasta", racon_test_data_path +
        "sample_overlaps.sam", racon_test_data_path + "sample_layout.fasta",
        racon::PolisherType::kC, 500, 10, 0.3, 5, -4, -8);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, true);
    EXPECT_EQ(polished_sequences.size(), 1);

    polished_sequences[0]->create_reverse_complement();

    auto parser = bioparser::createParser<bioparser::FastaParser, racon::Sequence>(
        racon_test_data_path + "sample_reference.fasta");
    parser->parse_objects(polished_sequences, -1);
    EXPECT_EQ(polished_sequences.size(), 2);

    EXPECT_EQ(calculateEditDistance(polished_sequences[0]->reverse_complement(),
        polished_sequences[1]->data()), 1770);
}

TEST_F(RaconPolishingTest, ConsensusWithQualitiesLargerWindow) {
    SetUp(racon_test_data_path + "sample_reads.fastq", racon_test_data_path +
        "sample_overlaps.paf", racon_test_data_path + "sample_layout.fasta",
        racon::PolisherType::kC, 1000, 10, 0.3, 5, -4, -8);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, true);
    EXPECT_EQ(polished_sequences.size(), 1);

    polished_sequences[0]->create_reverse_complement();

    auto parser = bioparser::createParser<bioparser::FastaParser, racon::Sequence>(
        racon_test_data_path + "sample_reference.fasta");
    parser->parse_objects(polished_sequences, -1);
    EXPECT_EQ(polished_sequences.size(), 2);

    EXPECT_EQ(calculateEditDistance(polished_sequences[0]->reverse_complement(),
        polished_sequences[1]->data()), 1289);
}

TEST_F(RaconPolishingTest, ConsensusWithQualitiesEditDistance) {
    SetUp(racon_test_data_path + "sample_reads.fastq", racon_test_data_path +
        "sample_overlaps.paf", racon_test_data_path + "sample_layout.fasta",
        racon::PolisherType::kC, 500, 10, 0.3, 1, -1, -1);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, true);
    EXPECT_EQ(polished_sequences.size(), 1);

    polished_sequences[0]->create_reverse_complement();

    auto parser = bioparser::createParser<bioparser::FastaParser, racon::Sequence>(
        racon_test_data_path + "sample_reference.fasta");
    parser->parse_objects(polished_sequences, -1);
    EXPECT_EQ(polished_sequences.size(), 2);

    EXPECT_EQ(calculateEditDistance(polished_sequences[0]->reverse_complement(),
        polished_sequences[1]->data()), 1321);
}

TEST_F(RaconPolishingTest, FragmentCorrectionWithQualities) {
    SetUp(racon_test_data_path + "sample_reads.fastq", racon_test_data_path +
        "sample_ava_overlaps.paf", racon_test_data_path + "sample_reads.fastq",
        racon::PolisherType::kC, 500, 10, 0.3, 1, -1, -1);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, true);
    EXPECT_EQ(polished_sequences.size(), 34);

    uint32_t total_length = 0;
    for (const auto& it: polished_sequences) {
        total_length += it->data().size();
    }
    EXPECT_EQ(total_length, 334619);
}

TEST_F(RaconPolishingTest, FragmentCorrectionWithQualitiesFull) {
    SetUp(racon_test_data_path + "sample_reads.fastq", racon_test_data_path +
        "sample_ava_overlaps.paf", racon_test_data_path + "sample_reads.fastq",
        racon::PolisherType::kF, 500, 10, 0.3, 1, -1, -1);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, false);
    EXPECT_EQ(polished_sequences.size(), 236);

    uint32_t total_length = 0;
    for (const auto& it: polished_sequences) {
        total_length += it->data().size();
    }
    EXPECT_EQ(total_length, 1658339);
}

TEST_F(RaconPolishingTest, FragmentCorrectionWithoutQualitiesFull) {
    SetUp(racon_test_data_path + "sample_reads.fasta", racon_test_data_path +
        "sample_ava_overlaps.paf", racon_test_data_path + "sample_reads.fasta",
        racon::PolisherType::kF, 500, 10, 0.3, 1, -1, -1);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, false);
    EXPECT_EQ(polished_sequences.size(), 236);

    uint32_t total_length = 0;
    for (const auto& it: polished_sequences) {
        total_length += it->data().size();
    }
    EXPECT_EQ(total_length, 1664087);
}

TEST_F(RaconPolishingTest, FragmentCorrectionWithQualitiesFullMhap) {
    SetUp(racon_test_data_path + "sample_reads.fastq", racon_test_data_path +
        "sample_ava_overlaps.mhap", racon_test_data_path + "sample_reads.fastq",
        racon::PolisherType::kF, 500, 10, 0.3, 1, -1, -1);

    initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polish(polished_sequences, false);
    EXPECT_EQ(polished_sequences.size(), 236);

    uint32_t total_length = 0;
    for (const auto& it: polished_sequences) {
        total_length += it->data().size();
    }
    EXPECT_EQ(total_length, 1658339);
}
