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
        alignment_engines_(), window_length_(window_length), windows_(),
        target_names_(), thread_pool_(thread_pool::createThreadPool(num_threads)),
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

    std::vector<std::unique_ptr<Sequence>> targets;
    std::unordered_map<std::string, uint64_t> t_name_to_id;

    tparser_->reset();
    tparser_->parse_objects(targets, -1);

    if (targets.empty()) {
        fprintf(stderr, "[racon::Polisher::initialize] error: "
            "empty target sequences set!\n");
        exit(1);
    }

    for (uint64_t i = 0; i < targets.size(); ++i) {
        target_names_.emplace_back(targets[i]->name());
        t_name_to_id[targets[i]->name()] = i;
    }

    fprintf(stderr, "[racon::Polisher::initialize] loaded target sequences\n");

    std::vector<std::unique_ptr<Sequence>> sequences;
    std::unordered_map<std::string, uint64_t> q_name_to_id;
    std::unordered_set<uint64_t> shared_sequences;

    uint64_t num_sequences = 0, total_sequences_length = 0;

    sparser_->reset();
    while (true) {
        auto status = sparser_->parse_objects(sequences, kChunkSize);

        for (uint64_t i = 0; i < sequences.size(); ++i, ++num_sequences) {
            q_name_to_id[sequences[i]->name()] = num_sequences;
            total_sequences_length += sequences[i]->data().size();
            if (t_name_to_id.find(sequences[i]->name()) != t_name_to_id.end()) {
                shared_sequences.emplace(num_sequences);
            }
        }
        sequences.clear();

        if (!status) {
            break;
        }
    }

    if (num_sequences == 0) {
        fprintf(stderr, "[racon::Polisher::initialize] error: "
            "empty sequences set!\n");
        exit(1);
    }

    WindowType window_type = static_cast<double>(total_sequences_length) /
        num_sequences <= 1000 ? WindowType::kNGS : WindowType::kTGS;

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
            overlaps[i]->transmute(q_name_to_id, t_name_to_id);

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

        fprintf(stderr, "[racon::Polisher::initialize] loaded batch of overlaps\n");

        uint64_t n = shrinkToFit(overlaps, l);
        l = c - n;

        if (!status) {
            break;
        }
    }

    if (overlaps.empty()) {
        fprintf(stderr, "[racon::Polisher::initialize] error: "
            "empty overlap set!\n");
        exit(1);
    }

    std::sort(overlaps.begin(), overlaps.end(),
        [](const std::unique_ptr<Overlap>& lhs, const std::unique_ptr<Overlap>& rhs) {
            return lhs->q_id() < rhs->q_id();
        });

    auto prepare_batch_of_sequences = [&]() -> bool {
        uint64_t l = sequences.size();
        auto status = sparser_->parse_objects(sequences, kChunkSize / 2);

        std::vector<std::future<void>> thread_futures;
        for (uint64_t i = l; i < sequences.size(); ++i) {
            thread_futures.emplace_back(thread_pool_->submit_task(
                [&](uint64_t j) -> void {
                    sequences[j]->create_reverse_complement();
                }, i));
        }
        for (const auto& it: thread_futures) {
            it.wait();
        }

        return status;
    };

    l = 0;
    sparser_->reset();
    while (true) {
        uint64_t ls = sequences.size();
        auto status = prepare_batch_of_sequences();

        std::vector<std::future<void>> thread_futures;
        for (uint64_t i = l; i < overlaps.size(); ++i) {
            if (overlaps[i]->q_id() > sequences.size() - 1) {
                break;
            }

            thread_futures.emplace_back(thread_pool_->submit_task(
                [&](uint64_t j) -> void {
                    overlaps[j]->find_breaking_points(window_length_, sequences,
                        targets);
                }, i));
        }

        for (uint64_t i = 0; i < thread_futures.size(); ++i, ++l) {
            thread_futures[i].wait();
            fprintf(stderr, "[racon::Polisher::initialize] aligned overlap %lu/%lu\r",
                l + 1, overlaps.size());
        }

        for (uint64_t i = ls; i < sequences.size(); ++i) {
            sequences[i].reset();
        }

        if (!status) {
            sequences.clear();
            break;
        }
    }
    fprintf(stderr, "\n");

    if (type_ == PolisherType::kF && !shared_sequences.empty()) {
        uint64_t num_overlaps = overlaps.size();
        for (uint64_t i = 0; i < num_overlaps; ++i) {
            if (shared_sequences.find(overlaps[i]->q_id()) != shared_sequences.end() &&
                q_name_to_id.find(targets[overlaps[i]->t_id()]->name()) != q_name_to_id.end()) {

                overlaps.emplace_back(overlaps[i]->dual_overlap());
            }
        }

        std::sort(overlaps.begin(), overlaps.end(),
            [](const std::unique_ptr<Overlap>& lhs, const std::unique_ptr<Overlap>& rhs) {
                return lhs->q_id() < rhs->q_id();
            });
    }

    std::string dummy_backbone_quality(window_length_ * 2, '!');
    std::vector<uint64_t> id_to_first_window_id(targets.size() + 1, 0);
    for (uint64_t i = 0; i < targets.size(); ++i) {
        uint32_t k = 0;
        for (uint32_t j = 0; j < targets[i]->data().size(); j += window_length_, ++k) {

            uint32_t length = std::min(j + window_length_,
                static_cast<uint32_t>(targets[i]->data().size())) - j;

            windows_.emplace_back(createWindow(i, k, window_type,
                &(targets[i]->data()[j]), length,
                targets[i]->quality().empty() ? &(dummy_backbone_quality[0]) :
                &(targets[i]->quality()[j]), length));
        }

        id_to_first_window_id[i + 1] = id_to_first_window_id[i] + k;
    }

    std::string dummy_layer_quality(window_length_ * 4, '"');
    l = 0;
    sparser_->reset();
    while (true) {
        uint64_t ls = sequences.size();
        auto status = prepare_batch_of_sequences();

        std::vector<std::future<void>> thread_futures;
        for (uint64_t i = l; i < overlaps.size(); ++i, ++l) {
            if (overlaps[i]->q_id() > sequences.size() - 1) {
                break;
            }

            const auto& breaking_points = overlaps[i]->breaking_points();
            for (uint32_t j = 0; j < breaking_points.size(); j += 2) {
                if (breaking_points[j + 1].second - breaking_points[j].second < 0.02 * window_length_) {
                    continue;
                }

                if (!sequences[overlaps[i]->q_id()]->quality().empty()) {
                    const auto& quality = overlaps[i]->strand() ?
                        sequences[overlaps[i]->q_id()]->reverse_quality() :
                        sequences[overlaps[i]->q_id()]->quality();
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

                const char* sequence = overlaps[i]->strand() ?
                    &(sequences[overlaps[i]->q_id()]->reverse_complement()[breaking_points[j].second]) :
                    &(sequences[overlaps[i]->q_id()]->data()[breaking_points[j].second]);

                const char* quality = sequences[overlaps[i]->q_id()]->quality().empty() ?
                    &(dummy_layer_quality[0]) : (overlaps[i]->strand() ?
                    &(sequences[overlaps[i]->q_id()]->reverse_quality()[breaking_points[j].second]) :
                    &(sequences[overlaps[i]->q_id()]->quality()[breaking_points[j].second]));

                uint32_t length = breaking_points[j + 1].second -
                    breaking_points[j].second;

                windows_[window_id]->add_layer(sequence, length, quality, length,
                    breaking_points[j].first - window_start,
                    breaking_points[j + 1].first - window_start - 1);
            }

            overlaps[i].reset();
        }

        if (!status) {
            sequences.clear();
            break;
        }

        for (uint64_t i = ls; i < sequences.size(); ++i) {
            sequences[i].reset();
        }
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
                std::string ratio_str = " C:" + std::to_string(polished_ratio);
                dst.emplace_back(createSequence(target_names_[windows_[i]->id()] +
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
    std::vector<std::string>().swap(target_names_);
}

}
