/*!
 * @file polisher.cpp
 *
 * @brief Polisher class source file
 */

#include <unordered_map>

#include "overlap.hpp"
#include "sequence.hpp"
#include "window.hpp"
#include "polisher.hpp"

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "spoa/spoa.hpp"

namespace racon {

constexpr uint32_t kChunkSize = 1024 * 1024 * 1024; // ~ 1GB

std::unique_ptr<Polisher> createPolisher(const std::string& sequences_path,
    const std::string& overlaps_path, const std::string& target_path,
    PolisherType type, uint32_t window_length, double quality_threshold,
    double error_threshold, int8_t match, int8_t mismatch, int8_t gap,
    uint32_t num_threads) {

    std::unique_ptr<bioparser::Parser<Sequence>> sparser = nullptr,
        tparser = nullptr;
    std::unique_ptr<bioparser::Parser<Overlap>> oparser = nullptr;

    auto extension = sequences_path.substr(std::min(sequences_path.rfind('.'),
        sequences_path.size()));
    if (extension == ".fasta" || extension == ".fa") {
        sparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            sequences_path);
    } else if (extension == ".fastq" || extension == ".fq") {
        sparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            sequences_path);
    } else {
        fprintf(stderr, "racon::createPolisher error:"
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fa, .fastq, .fq)!", sequences_path.c_str());
        exit(1);
    }

    extension = overlaps_path.substr(std::min(overlaps_path.rfind('.'),
        overlaps_path.size()));
    if (overlaps_path == ".mhap") {
        oparser = bioparser::createParser<bioparser::MhapParser, Overlap>(
            overlaps_path);
    } else if (overlaps_path == ".paf") {
        oparser = bioparser::createParser<bioparser::PafParser, Overlap>(
            overlaps_path);
    } else if (overlaps_path == ".sam") {
        oparser = bioparser::createParser<bioparser::SamParser, Overlap>(
            overlaps_path);
    } else {
        fprintf(stderr, "racon::createPolisher error:"
            "file %s has unsupported format extension (valid extensions: "
            ".mhap, .paf, .sam)!", overlaps_path.c_str());
        exit(1);
    }

    extension = target_path.substr(std::min(target_path.rfind('.'),
        target_path.size()));
    if (extension == ".fasta" || extension == ".fa") {
        tparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            target_path);
    } else if (extension == ".fastq" || extension == ".fq") {
        tparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            target_path);
    } else {
        fprintf(stderr, "racon::createPolisher error:"
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fa, .fastq, .fq)!", target_path.c_str());
        exit(1);
    }

    return std::unique_ptr<Polisher>(new Polisher(std::move(sparser),
        std::move(oparser), std::move(tparser), type, window_length,
        quality_threshold, error_threshold, match, mismatch, gap,
        num_threads));
}

Polisher::Polisher(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
    std::unique_ptr<bioparser::Parser<Overlap>> oparser,
    std::unique_ptr<bioparser::Parser<Sequence>> tparser,
    PolisherType type, uint32_t window_length, double quality_threshold,
    double error_threshold, int8_t match, int8_t mismatch, int8_t gap,
    uint32_t num_threads)
        : sparser_(std::move(sparser)), oparser_(std::move(oparser)),
        tparser_(std::move(tparser)), type_(type), quality_threshold_(
        quality_threshold), error_threshold_(error_threshold),
        alignment_engines_(), window_length_(window_length), windows_(),
        thread_pool(thread_pool::createThreadPool(num_threads)) {

    for (uint32_t i = 0; i < num_threads; ++i) {
        alignment_engines_.emplace_back(spoa::createAlignmentEngine(
            spoa::AlignmentType::kNW, match, mismatch, gap));
    }
}

Polisher::~Polisher() {
}

void Polisher::initialize() {

    std::vector<std::unique_ptr<Sequence>> sequences;
    std::vector<std::unique_ptr<Overlap>> overlaps;

}

}
