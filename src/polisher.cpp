/*!
 * @file polisher.cpp
 *
 * @brief Polisher class source file
 */

#include <algorithm>
#include <unordered_set>

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
uint64_t shrinkToFit(std::vector<std::unique_ptr<T>>& src, uint64_t begin) {

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
    uint64_t num_deletions = src.size() - i;
    if (i < src.size()) {
        src.resize(i);
    }
    return num_deletions;
}

std::unique_ptr<Polisher> createPolisher(const std::string& sequences_path,
    const std::string& overlaps_path, const std::string& target_path,
    PolisherType type, uint32_t window_length, double quality_threshold,
    double error_threshold, int8_t match, int8_t mismatch, int8_t gap,
    uint32_t num_threads) {

    if (type != PolisherType::kC && type != PolisherType::kF) {
        fprintf(stderr, "[racon::createPolisher] error: invalid polisher type!\n");
        exit(1);
    }

    if (window_length == 0) {
        fprintf(stderr, "[racon::createPolisher] error: invalid window length!\n");
        exit(1);
    }

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
        fprintf(stderr, "[racon::createPolisher] error: "
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
        fprintf(stderr, "[racon::createPolisher] error: "
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
        fprintf(stderr, "[racon::createPolisher] error: "
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
        alignment_engines_(), sequences_(),
        dummy_backbone_quality_(window_length * 2, '!'),
        dummy_layer_quality_(window_length * 4, '"'),
        window_length_(window_length), windows_(),
        thread_pool_(thread_pool::createThreadPool(num_threads)),
        thread_to_id_() {

    uint32_t id = 0;
    for (const auto& it: thread_pool_->thread_identifiers()) {
        thread_to_id_[it] = id++;
    }

    for (uint32_t i = 0; i < num_threads; ++i) {
        alignment_engines_.emplace_back(spoa::createAlignmentEngine(
            spoa::AlignmentType::kNW, match, mismatch, gap));
        alignment_engines_.back()->prealloc(window_length_, 5);
    }
}

Polisher::~Polisher() {
}

void Polisher::initialize() {

    if (!windows_.empty()) {
        fprintf(stderr, "[racon::Polisher::initialize] warning: "
            "object already initialized!\n");
        return;
    }

    tparser_->reset();
    tparser_->parse_objects(sequences_, -1);

    uint64_t targets_size = sequences_.size();
    if (targets_size == 0) {
        fprintf(stderr, "[racon::Polisher::initialize] error: "
            "empty target sequences set!\n");
        exit(1);
    }

    std::unordered_map<std::string, uint64_t> name_to_id;
    std::unordered_map<uint64_t, uint64_t> id_to_id;
    for (uint64_t i = 0; i < targets_size; ++i) {
        name_to_id[sequences_[i]->name()] = i;
        id_to_id[i << 1 | 1] = i;
    }

    std::vector<bool> has_name(targets_size, true);
    std::vector<bool> has_data(targets_size, true);
    std::vector<bool> has_reverse_data(targets_size, false);

    fprintf(stderr, "[racon::Polisher::initialize] loaded target sequences\n");

    std::unordered_set<uint64_t> duplicate_sequences;
    uint64_t sequences_size = 0, total_sequences_length = 0;

    sparser_->reset();
    while (true) {
        uint64_t l = sequences_.size();
        auto status = sparser_->parse_objects(sequences_, kChunkSize);

        uint64_t n = 0;
        for (uint64_t i = l; i < sequences_.size(); ++i, ++sequences_size) {
            total_sequences_length += sequences_[i]->data().size();

            auto it = name_to_id.find(sequences_[i]->name());
            if (it != name_to_id.end()) {
                uint64_t j = it->second;
                if (j >= targets_size) {
                    fprintf(stderr, "[racon::Polisher::initialize] error: "
                        "duplicate sequence %s in file\n",
                        sequences_[i]->name().c_str());
                    exit(1);
                }
                if (sequences_[i]->data().size() != sequences_[j]->data().size() ||
                    sequences_[i]->quality().size() != sequences_[j]->quality().size()) {

                    fprintf(stderr, "[racon::Polisher::initialize] error: "
                        "duplicate sequence %s with unequal data\n",
                        sequences_[i]->name().c_str());
                    exit(1);
                }

                duplicate_sequences.insert(j);
                sequences_[i].reset();
                ++n;
                id_to_id[sequences_size << 1 | 0] = j;
            } else {
                name_to_id[sequences_[i]->name()] = i - n;
                id_to_id[sequences_size << 1 | 0] = i - n;
            }
        }

        shrinkToFit(sequences_, l);

        if (!status) {
            break;
        }
    }

    if (sequences_size == 0) {
        fprintf(stderr, "[racon::Polisher::initialize] error: "
            "empty sequences set!\n");
        exit(1);
    }

    has_name.resize(sequences_.size(), false);
    has_data.resize(sequences_.size(), false);
    has_reverse_data.resize(sequences_.size(), false);

    WindowType window_type = static_cast<double>(total_sequences_length) /
        sequences_size <= 1000 ? WindowType::kNGS : WindowType::kTGS;

    fprintf(stderr, "[racon::Polisher::initialize] loaded sequences\n");

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
                    if (overlaps[i]->length() > overlaps[j]->length()) {
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
            overlaps[i]->transmute(name_to_id, id_to_id);

            while (overlaps[c] == nullptr) {
                ++c;
            }
            if (overlaps[c]->q_id() != overlaps[i]->q_id()) {
                remove_invalid_overlaps(c, i);
                c = i;
            }
        }
        if (!status) {
            remove_invalid_overlaps(c, overlaps.size());
            c = overlaps.size();
        }

        for (uint64_t i = l; i < c; ++i) {
            if (overlaps[i] == nullptr) {
                continue;
            }

            if (overlaps[i]->strand()) {
                has_reverse_data[overlaps[i]->q_id()] = true;
                if (type_ == PolisherType::kF &&
                    duplicate_sequences.find(overlaps[i]->q_id()) != duplicate_sequences.end() &&
                    duplicate_sequences.find(overlaps[i]->t_id()) != duplicate_sequences.end()) {

                    has_reverse_data[overlaps[i]->t_id()] = true;
                }
            } else {
                has_data[overlaps[i]->q_id()] = true;
            }
        }

        fprintf(stderr, "[racon::Polisher::initialize] loaded batch of overlaps\n");

        uint64_t n = shrinkToFit(overlaps, l);
        l = c - n;

        if (!status) {
            break;
        }
    }

    std::unordered_map<std::string, uint64_t>().swap(name_to_id);
    std::unordered_map<uint64_t, uint64_t>().swap(id_to_id);

    if (overlaps.empty()) {
        fprintf(stderr, "[racon::Polisher::initialize] error: "
            "empty overlap set!\n");
        exit(1);
    }

    std::vector<std::future<void>> thread_futures;
    for (uint64_t i = 0; i < sequences_.size(); ++i) {
        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint64_t j) -> void {
                sequences_[j]->transmute(has_name[j], has_data[j], has_reverse_data[j]);
            }, i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }

    thread_futures.clear();
    for (uint64_t i = 0; i < overlaps.size(); ++i) {
        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint64_t j) -> void {
                overlaps[j]->find_breaking_points(sequences_, window_length_);
            }, i));
    }
    for (uint64_t i = 0; i < thread_futures.size(); ++i) {
        thread_futures[i].wait();
        fprintf(stderr, "[racon::Polisher::initialize] aligned overlap %lu/%lu\r",
            i + 1, overlaps.size());
    }
    fprintf(stderr, "\n");

    if (type_ == PolisherType::kF && !duplicate_sequences.empty()) {
        uint64_t num_overlaps = overlaps.size();
        for (uint64_t i = 0; i < num_overlaps; ++i) {
            if (duplicate_sequences.find(overlaps[i]->q_id()) != duplicate_sequences.end() &&
                duplicate_sequences.find(overlaps[i]->t_id()) != duplicate_sequences.end()) {

                overlaps.emplace_back(overlaps[i]->dual_overlap());
            }
        }

        std::sort(overlaps.begin(), overlaps.end(),
            [](const std::unique_ptr<Overlap>& lhs, const std::unique_ptr<Overlap>& rhs) {
                return lhs->q_id() < rhs->q_id();
            });
    }

    std::vector<uint64_t> id_to_first_window_id(targets_size + 1, 0);
    for (uint64_t i = 0; i < targets_size; ++i) {
        uint32_t k = 0;
        for (uint32_t j = 0; j < sequences_[i]->data().size(); j += window_length_, ++k) {

            uint32_t length = std::min(j + window_length_,
                static_cast<uint32_t>(sequences_[i]->data().size())) - j;

            windows_.emplace_back(createWindow(i, k, window_type,
                &(sequences_[i]->data()[j]), length,
                sequences_[i]->quality().empty() ? &(dummy_backbone_quality_[0]) :
                &(sequences_[i]->quality()[j]), length));
        }

        id_to_first_window_id[i + 1] = id_to_first_window_id[i] + k;
    }

    for (uint64_t i = 0; i < overlaps.size(); ++i) {

        const auto& sequence = sequences_[overlaps[i]->q_id()];
        const auto& breaking_points = overlaps[i]->breaking_points();

        for (uint32_t j = 0; j < breaking_points.size(); j += 2) {
            if (breaking_points[j + 1].second - breaking_points[j].second < 0.02 * window_length_) {
                continue;
            }

            if (!sequence->quality().empty() ||
                !sequence->reverse_quality().empty()) {

                const auto& quality = overlaps[i]->strand() ?
                    sequence->reverse_quality() : sequence->quality();
                double average_quality = 0;
                for (uint32_t k = breaking_points[j].second; k < breaking_points[j + 1].second; ++k) {
                    average_quality += static_cast<uint32_t>(quality[k]) - 33;
                }
                average_quality /= breaking_points[j + 1].second - breaking_points[j].second;

                if (average_quality < quality_threshold_) {
                    continue;
                }
            }

            uint64_t window_id = id_to_first_window_id[overlaps[i]->t_id()] +
                breaking_points[j].first / window_length_;
            uint32_t window_start = (breaking_points[j].first / window_length_) *
                window_length_;

            const char* data = overlaps[i]->strand() ?
                &(sequence->reverse_complement()[breaking_points[j].second]) :
                &(sequence->data()[breaking_points[j].second]);

            const char* quality = overlaps[i]->strand() ?
                (sequence->reverse_quality().empty() ?
                    &(dummy_layer_quality_[0]) :
                    &(sequence->reverse_quality()[breaking_points[j].second]))
                :
                (sequence->quality().empty() ?
                    &(dummy_layer_quality_[0]) :
                    &(sequence->quality()[breaking_points[j].second]));

            uint32_t length = breaking_points[j + 1].second -
                breaking_points[j].second;

            windows_[window_id]->add_layer(data, length, quality, length,
                breaking_points[j].first - window_start,
                breaking_points[j + 1].first - window_start - 1);
        }

        overlaps[i].reset();
    }

    fprintf(stderr, "[racon::Polisher::initialize] transformed data into windows\n");
}

void Polisher::polish(std::vector<std::unique_ptr<Sequence>>& dst,
    bool drop_unpolished_sequences) {

    std::vector<std::future<bool>> thread_futures;
    for (uint64_t i = 0; i < windows_.size(); ++i) {
        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint64_t j) -> bool {
                auto it = thread_to_id_.find(std::this_thread::get_id());
                if (it == thread_to_id_.end()) {
                    fprintf(stderr, "[racon::Polisher::polish] error: "
                        "thread identifier not present!\n");
                    exit(1);
                }
                return windows_[j]->generate_consensus(
                    alignment_engines_[it->second]);
            }, i));
    }

    std::string polished_data = "";
    uint32_t num_polished_windows = 0;

    for (uint64_t i = 0; i < thread_futures.size(); ++i) {
        thread_futures[i].wait();

        num_polished_windows += thread_futures[i].get() == true ? 1 : 0;
        polished_data += windows_[i]->consensus();

        if (i == windows_.size() - 1 || windows_[i + 1]->rank() == 0) {
            double polished_ratio = num_polished_windows /
                static_cast<double>(windows_[i]->rank() + 1);

            if (!drop_unpolished_sequences || polished_ratio > 0) {
                std::string ratio_str = "_C:" + std::to_string(polished_ratio);
                dst.emplace_back(createSequence(sequences_[windows_[i]->id()]->name() +
                    ratio_str, polished_data));
            }

            num_polished_windows = 0;
            polished_data.clear();
        }
        windows_[i].reset();

        fprintf(stderr, "[racon::Polisher::polish] generated consensus for window %lu/%lu\r",
            i + 1, thread_futures.size());
    }
    fprintf(stderr, "\n");

    std::vector<std::unique_ptr<Window>>().swap(windows_);
    std::vector<std::unique_ptr<Sequence>>().swap(sequences_);
}

}
