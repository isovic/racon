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

template<class T>
void shrinkToFit(std::vector<std::unique_ptr<T>>& src, uint64_t begin) {

    uint64_t i = begin;
    for (uint64_t j = begin; i < src.size(); ++i) {
        if (src[i] != nullptr) {
            continue;
        }

        j = std::max(j, i);
        while (j < src.size() && src[j] == nullptr) {
            ++j;
        }

        if (j >= src.size()) {
            break;
        } else if (i != j) {
            src[i].swap(src[j]);
        }
    }
    if (i < src.size()) {
        src.resize(i);
    }
}

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
        fprintf(stderr, "racon::createPolisher error: "
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fa, .fastq, .fq)!\n", sequences_path.c_str());
        exit(1);
    }

    extension = overlaps_path.substr(std::min(overlaps_path.rfind('.'),
        overlaps_path.size()));
    if (extension == ".mhap") {
        oparser = bioparser::createParser<bioparser::MhapParser, Overlap>(
            overlaps_path);
    } else if (extension == ".paf") {
        oparser = bioparser::createParser<bioparser::PafParser, Overlap>(
            overlaps_path);
    } else if (extension == ".sam") {
        oparser = bioparser::createParser<bioparser::SamParser, Overlap>(
            overlaps_path);
    } else {
        fprintf(stderr, "racon::createPolisher error: "
            "file %s has unsupported format extension (valid extensions: "
            ".mhap, .paf, .sam)!\n", overlaps_path.c_str());
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
        fprintf(stderr, "racon::createPolisher error: "
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fa, .fastq, .fq)!\n", target_path.c_str());
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
    std::unordered_map<std::string, uint32_t> name_to_nid;
    std::unordered_map<uint32_t, uint32_t> id_to_nid;

    sparser_->reset();
    sparser_->parse_objects(sequences, -1);

    uint32_t sequences_size = sequences.size();
    for (uint32_t i = 0; i < sequences.size(); ++i) {
        name_to_nid[sequences[i]->name()] = i;
        id_to_nid[i << 1 | 0] = i;
    }

    tparser_->reset();
    uint32_t t = sequences_size, j = sequences_size;
    while (true) {
        uint32_t l = sequences.size();
        auto status = tparser_->parse_objects(sequences, kChunkSize);

        for (uint32_t i = l; i < sequences.size(); ++i, ++j) {
            auto it = name_to_nid.find(sequences[i]->name());
            if (it != name_to_nid.end()) {
                sequences[i].reset();
                id_to_nid[(j - sequences_size) << 1 | 1] = it->second;
            } else {
                name_to_nid[sequences[i]->name()] = t;
                id_to_nid[(j - sequences_size) << 1 | 1] = t;
                ++t;
            }
        }

        shrinkToFit(sequences, l);

        if (!status) {
            break;
        }
    }

    std::vector<std::unique_ptr<Overlap>> overlaps;

    auto remove_invalid_overlaps = [&](uint64_t begin, uint64_t end) -> void {
        for (uint64_t i = begin; i < end; ++i) {
            if (overlaps[i] == nullptr) {
                continue;
            }
            if (overlaps[i]->error() > error_threshold_) {
                overlaps[i].reset();
                continue;
            }
            if (type_ == PolisherType::kC) {
                for (uint64_t j = i + 1; j < end; ++j) {
                    if (overlaps[j] == nullptr) {
                        continue;
                    }
                    if (overlaps[i]->error() < overlaps[j]->error()) {
                        overlaps[j].reset();
                    } else {
                        overlaps[i].reset();
                        break;
                    }
                }
            }
        }
    };

    oparser_->reset();
    uint64_t l = 0;
    while (true) {
        auto status = oparser_->parse_objects(overlaps, kChunkSize);

        uint64_t c = l;
        for (uint64_t i = l; i < overlaps.size(); ++i) {
            if (!overlaps[i]->is_valid()) {
                overlaps[i].reset();
                continue;
            }
            overlaps[i]->transmute(name_to_nid, id_to_nid);

            while (overlaps[c] == nullptr) {
                ++c;
            }
            if (overlaps[c]->q_id() != overlaps[i]->q_id()) {
                remove_invalid_overlaps(c, i);
                c = i;
            } else if (!status && i + 1 == overlaps.size()) {
                remove_invalid_overlaps(c, i + 1);
            }
        }

        uint64_t n = 0;
        for (uint64_t i = l; i < c; ++i) {
            if (overlaps[i] == nullptr) {
                ++n;
            }
        }
        shrinkToFit(overlaps, l);
        l = c - n;

        if (!status) {
            break;
        }
    }
}

}
