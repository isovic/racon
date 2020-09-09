/*!
 * @file polisher.cpp
 *
 * @brief Polisher class source file
 */

#include <algorithm>
#include <unordered_set>
#include <iostream>

#include "overlap.hpp"
#include "sequence.hpp"
#include "logger.hpp"
#include "polisher.hpp"
#include "util.hpp"
#ifdef CUDA_ENABLED
#include "cuda/cudapolisher.hpp"
#endif

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "spoa/spoa.hpp"

// #define BED_FEATURE_TEST

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
    const std::string& bed_path,
    PolisherType type, uint32_t window_length, double quality_threshold,
    double error_threshold, bool trim, bool produce_liftover, int8_t match, int8_t mismatch, int8_t gap,
    uint32_t num_threads, uint32_t cudapoa_batches, bool cuda_banded_alignment,
    uint32_t cudaaligner_batches, uint32_t cudaaligner_band_width) {

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

    auto is_suffix = [](const std::string& src, const std::string& suffix) -> bool {
        if (src.size() < suffix.size()) {
            return false;
        }
        return src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    if (is_suffix(sequences_path, ".fasta") || is_suffix(sequences_path, ".fasta.gz") ||
        is_suffix(sequences_path, ".fna") || is_suffix(sequences_path, ".fna.gz") ||
        is_suffix(sequences_path, ".fa") || is_suffix(sequences_path, ".fa.gz")) {
        sparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            sequences_path);
    } else if (is_suffix(sequences_path, ".fastq") || is_suffix(sequences_path, ".fastq.gz") ||
        is_suffix(sequences_path, ".fq") || is_suffix(sequences_path, ".fq.gz")) {
        sparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            sequences_path);
    } else {
        fprintf(stderr, "[racon::createPolisher] error: "
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
            ".fq, .fq.gz)!\n",
            sequences_path.c_str());
        exit(1);
    }

    if (is_suffix(overlaps_path, ".mhap") || is_suffix(overlaps_path, ".mhap.gz")) {
        oparser = bioparser::createParser<bioparser::MhapParser, Overlap>(
            overlaps_path);
    } else if (is_suffix(overlaps_path, ".paf") || is_suffix(overlaps_path, ".paf.gz")) {
        oparser = bioparser::createParser<bioparser::PafParser, Overlap>(
            overlaps_path);
    } else if (is_suffix(overlaps_path, ".sam") || is_suffix(overlaps_path, ".sam.gz")) {
        oparser = bioparser::createParser<bioparser::SamParser, Overlap>(
            overlaps_path);
    } else {
        fprintf(stderr, "[racon::createPolisher] error: "
            "file %s has unsupported format extension (valid extensions: "
            ".mhap, .mhap.gz, .paf, .paf.gz, .sam, .sam.gz)!\n", overlaps_path.c_str());
        exit(1);
    }

    if (is_suffix(target_path, ".fasta") || is_suffix(target_path, ".fasta.gz") ||
        is_suffix(target_path, ".fna") || is_suffix(target_path, ".fna.gz") ||
        is_suffix(target_path, ".fa") || is_suffix(target_path, ".fa.gz")) {
        tparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            target_path);
    } else if (is_suffix(target_path, ".fastq") || is_suffix(target_path, ".fastq.gz") ||
        is_suffix(target_path, ".fq") || is_suffix(target_path, ".fq.gz")) {
        tparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            target_path);
    } else {
        fprintf(stderr, "[racon::createPolisher] error: "
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
            ".fq, .fq.gz)!\n",
            target_path.c_str());
        exit(1);
    }

    bool use_bed = false;
    std::vector<BedRecord> bed_records;
    if (bed_path.size() > 0) {
        use_bed = true;
        bed_records = BedReader::ReadAll(bed_path);
    }
    std::cerr << "[racon::createPolisher] Use bed: " << (use_bed ? "true" : "false") << ", bed records: " << bed_records.size() << "\n";

    if (cudapoa_batches > 0 || cudaaligner_batches > 0)
    {
#ifdef CUDA_ENABLED
        // If CUDA is enabled, return an instance of the CUDAPolisher object.
        return std::unique_ptr<Polisher>(new CUDAPolisher(std::move(sparser),
                    std::move(oparser), std::move(tparser),
                    std::move(bed_records), use_bed,
                    type, window_length,
                    quality_threshold, error_threshold, trim, match, mismatch, gap,
                    num_threads, cudapoa_batches, cuda_banded_alignment, cudaaligner_batches,
                    cudaaligner_band_width));
#else
        fprintf(stderr, "[racon::createPolisher] error: "
                "Attemping to use CUDA when CUDA support is not available.\n"
                "Please check logic in %s:%s\n",
                __FILE__, __func__);
        exit(1);
#endif
    }
    else
    {
        (void) cuda_banded_alignment;
        return std::unique_ptr<Polisher>(new Polisher(std::move(sparser),
                    std::move(oparser), std::move(tparser),
                    std::move(bed_records), use_bed,
                    type, window_length,
                    quality_threshold, error_threshold, trim, produce_liftover, match, mismatch, gap,
                    num_threads));
    }
}

Polisher::Polisher(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
    std::unique_ptr<bioparser::Parser<Overlap>> oparser,
    std::unique_ptr<bioparser::Parser<Sequence>> tparser,
    std::vector<BedRecord> bed_records, bool use_bed,
    PolisherType type, uint32_t window_length, double quality_threshold,
    double error_threshold, bool trim, bool produce_liftover, int8_t match, int8_t mismatch, int8_t gap,
    uint32_t num_threads)
        : sparser_(std::move(sparser)), oparser_(std::move(oparser)),
        tparser_(std::move(tparser)), bed_records_(std::move(bed_records)),
        use_bed_(use_bed), type_(type), quality_threshold_(
        quality_threshold), error_threshold_(error_threshold), trim_(trim),
        produce_liftover_(produce_liftover),
        alignment_engines_(), sequences_(), dummy_quality_(window_length, '!'),
        window_length_(window_length), windows_(),
        thread_pool_(thread_pool::createThreadPool(num_threads)),
        thread_to_id_(), logger_(new Logger()) {

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
    logger_->total("[racon::Polisher::] total =");
}

void Polisher::initialize() {

    if (!windows_.empty()) {
        fprintf(stderr, "[racon::Polisher::initialize] warning: "
            "object already initialized!\n");
        return;
    }

    logger_->log();

    tparser_->reset();
    tparser_->parse(sequences_, -1);

    uint64_t targets_size = sequences_.size();
    if (targets_size == 0) {
        fprintf(stderr, "[racon::Polisher::initialize] error: "
            "empty target sequences set!\n");
        exit(1);
    }

    std::unordered_map<std::string, uint64_t> name_to_id;
    std::unordered_map<uint64_t, uint64_t> id_to_id;
    for (uint64_t i = 0; i < targets_size; ++i) {
        name_to_id[sequences_[i]->name() + "t"] = i;
        id_to_id[i << 1 | 1] = i;
    }

    std::vector<bool> has_name(targets_size, true);
    std::vector<bool> has_data(targets_size, true);
    std::vector<bool> has_reverse_data(targets_size, false);

    logger_->log("[racon::Polisher::initialize] loaded target sequences");
    logger_->log();

    uint64_t sequences_size = 0, total_sequences_length = 0;

    sparser_->reset();
    while (true) {
        uint64_t l = sequences_.size();
        auto status = sparser_->parse(sequences_, kChunkSize);

        uint64_t n = 0;
        for (uint64_t i = l; i < sequences_.size(); ++i, ++sequences_size) {
            total_sequences_length += sequences_[i]->data().size();

            auto it = name_to_id.find(sequences_[i]->name() + "t");
            if (it != name_to_id.end()) {
                if (sequences_[i]->data().size() != sequences_[it->second]->data().size() ||
                    sequences_[i]->quality().size() != sequences_[it->second]->quality().size()) {

                    fprintf(stderr, "[racon::Polisher::initialize] error: "
                        "duplicate sequence %s with unequal data\n",
                        sequences_[i]->name().c_str());
                    exit(1);
                }

                name_to_id[sequences_[i]->name() + "q"] = it->second;
                id_to_id[sequences_size << 1 | 0] = it->second;

                sequences_[i].reset();
                ++n;
            } else {
                name_to_id[sequences_[i]->name() + "q"] = i - n;
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

    logger_->log("[racon::Polisher::initialize] loaded sequences");
    logger_->log();



    //////////////////////////////////
    /// Collect the BED intervals. ///
    //////////////////////////////////
    logger_->log("[racon::Polisher::initialize] building the interval trees");
    logger_->log();
    target_bed_intervals_.clear();
    // Collect target intervals of interest.
    // If the BED file is not specified, construct windows covering full spans of targets.
    if (use_bed_) {
        for (size_t i = 0; i < bed_records_.size(); ++i) {
            const auto& record = bed_records_[i];
            uint64_t t_id = 0;
            if (!transmuteId(name_to_id, record.chrom() + "t", t_id)) {
                throw std::runtime_error("Target sequence '" + record.chrom() +
                            "' specified in the BED file was not found among the target sequences.");
            }
            target_bed_intervals_[t_id].emplace_back(IntervalInt64(record.chrom_start(), record.chrom_end() - 1, i));
        }
    } else {
        for (uint64_t t_id = 0; t_id < targets_size; ++t_id) {
            target_bed_intervals_[t_id].emplace_back(IntervalInt64(0, static_cast<int64_t>(sequences_[t_id]->data().size()) - 1, -1));
        }
    }
    // Sort target intervals.
    for (auto& it: target_bed_intervals_) {
        std::stable_sort(it.second.begin(), it.second.end(),
            [](const IntervalInt64& a, const IntervalInt64& b) { return a.start < b.start; });
    }
    // Construct the trees.
    target_bed_trees_.clear();
    for (const auto& it: target_bed_intervals_) {
        // Make a copy, because the IntervalTree has only the move constructor,
        // and we still need t he intvervals for validation below.
        auto intervals = it.second;
        target_bed_trees_[it.first] = IntervalTreeInt64(std::move(intervals));
    }
    // Validate that there are no overlapping intervals.
    for (const auto& it: target_bed_intervals_) {
        int64_t t_id = it.first;
        for (const auto& interval: it.second) {
            auto foundIntervals = target_bed_trees_[t_id].findOverlapping(interval.start, interval.stop);
            if (foundIntervals.size() != 1 ||
                (foundIntervals.size() == 1 && foundIntervals.front().value != interval.value)) {
                throw std::runtime_error("Invalid BED record: '" +
                    BedFile::Serialize(bed_records_[interval.value]) +
                    "'. It overlaps other BED records, which is not allowed.");
            }
        }
    }
    /////////////////////////////////



    std::vector<std::unique_ptr<Overlap>> overlaps;

    auto remove_invalid_overlaps = [&](uint64_t begin, uint64_t end) -> void {
        for (uint64_t i = begin; i < end; ++i) {
            if (overlaps[i] == nullptr) {
                continue;
            }
            if (overlaps[i]->error() > error_threshold_ ||
                overlaps[i]->q_id() == overlaps[i]->t_id()) {
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
        auto status = oparser_->parse(overlaps, kChunkSize);

        uint64_t c = l;
        for (uint64_t i = l; i < overlaps.size(); ++i) {
            overlaps[i]->transmute(sequences_, name_to_id, id_to_id);

            if (!overlaps[i]->is_valid()) {
                overlaps[i].reset();
                continue;
            }

            // Remove overlaps in regions not specified by BED.
            if (use_bed_) {
                auto foundIntervals = target_bed_trees_[overlaps[i]->t_id()].findOverlapping(
                        static_cast<int64_t>(overlaps[i]->t_begin()),
                        static_cast<int64_t>(overlaps[i]->t_end()) - 1);
                if (foundIntervals.empty()) {
                    overlaps[i].reset();
                    continue;
                }
            }

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
            } else {
                has_data[overlaps[i]->q_id()] = true;
            }
        }

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

    logger_->log("[racon::Polisher::initialize] loaded overlaps");
    logger_->log();

    std::vector<std::future<void>> thread_futures;
    for (uint64_t i = 0; i < sequences_.size(); ++i) {
        thread_futures.emplace_back(thread_pool_->submit(
            [&](uint64_t j) -> void {
                sequences_[j]->transmute(has_name[j], has_data[j], has_reverse_data[j]);
            }, i));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }

    logger_->log("[racon::Polisher::initialize] Constructing windows for BED regions.\n");

    create_and_populate_windows_with_bed(overlaps, targets_size, window_type);

    logger_->log("[racon::Polisher::initialize] transformed data into windows");

// #ifdef BED_FEATURE_TEST
//     int32_t shift = 0;
//     for (int32_t i = 0; i < static_cast<int32_t>(windows_.size()); ++i) {
//         if ((i + 1) % 100 != 0) {
//             windows_[i].reset();
//             ++shift;
//             // std::cerr << "nulling i = " << i << "\n";
//             continue;
//         }
//         std::swap(windows_[i-shift], windows_[i]);
//     }
//     windows_.resize(static_cast<int32_t>(windows_.size()) - shift);
//     std::cerr << "windows_.size() = " << windows_.size() << "\n";
// #endif

}

void Polisher::create_and_populate_windows_with_bed(std::vector<std::unique_ptr<Overlap>>& overlaps,
        uint64_t targets_size, WindowType window_type) {

    // The -1 marks that the target doesn't have any windows.
    id_to_first_window_id_.clear();
    id_to_first_window_id_.resize(targets_size + 1, -1);

    std::unordered_map<int64_t, std::vector<std::tuple<int64_t, int64_t, int64_t>>> windows;

    // Target intervals are sorted ahead of time.
    for (const auto& it: target_bed_intervals_) {
        int64_t t_id = it.first;
        const std::vector<IntervalInt64>& intervals = it.second;

        // Mark the window start.
        id_to_first_window_id_[t_id] = static_cast<int64_t>(windows_.size());

        // Generate windows for each interval separately.
        uint32_t k = 0;
        for (const auto& interval: intervals) {
            for (int64_t win_start = interval.start; win_start < (interval.stop + 1);
                 win_start += window_length_, ++k) {

                int64_t length = std::min(win_start + static_cast<int64_t>(window_length_),
                    (interval.stop + 1)) - win_start;
                int64_t win_id = windows_.size();

                windows_.emplace_back(createWindow(t_id, k, window_type, win_start,
                    &(sequences_[t_id]->data()[win_start]), length,
                    sequences_[t_id]->quality().empty() ? &(dummy_quality_[0]) :
                    &(sequences_[t_id]->quality()[win_start]), length));

                target_window_intervals_[t_id].emplace_back(IntervalInt64(win_start, win_start + length - 1, win_id));
                windows[t_id].emplace_back(win_start, win_start + length, win_id);
            }
        }
    }
    // Construct the trees. The iterator is not const because the IntervalTree only has the
    // move constructor.
    target_window_trees_.clear();
    for (auto& it: target_window_intervals_) {
        target_window_trees_[it.first] = IntervalTreeInt64(std::move(it.second));
    }

// #ifdef BED_FEATURE_TEST
//     for (size_t i = 0; i < windows_.size(); ++i) {
//         std::cerr << "[window " << i << "] " << *windows_[i] << "\n";
//     }
// #endif

    find_overlap_breaking_points(overlaps, windows);

// #ifdef BED_FEATURE_TEST
//     for (uint64_t i = 0; i < overlaps.size(); ++i) {
//         const auto& sequence = sequences_[overlaps[i]->q_id()];
//         const std::vector<WindowInterval>& breaking_points = overlaps[i]->breaking_points();

//         std::cerr << "overlap_id = " << i << "\n";
//         std::cerr << "    " << *overlaps[i] << "\n";
//         std::cerr << "All breaking points:\n";
//         for (uint32_t j = 0; j < breaking_points.size(); ++j) {
//             const auto& bp = breaking_points[j];
//             std::cerr << "[j = " << j << "] bp = " << bp << ", Window: " << *windows_[bp.window_id] << "\n";
//             if (bp.t_start < windows_[bp.window_id]->start() || bp.t_start >= windows_[bp.window_id]->end() ||
//                 bp.t_end < windows_[bp.window_id]->start() || bp.t_end > windows_[bp.window_id]->end()) {
//                 std::cerr << "ERROR! Coordiantes out of bounds!\n";
//                 exit(1);
//             }
//         }
//         std::cerr << "\n";
//     }
// #endif

    assign_sequences_to_windows(overlaps, targets_size);
}

void Polisher::assign_sequences_to_windows(std::vector<std::unique_ptr<Overlap>>& overlaps, uint64_t targets_size) {

    targets_coverages_.resize(targets_size, 0);

    for (uint64_t i = 0; i < overlaps.size(); ++i) {

        ++targets_coverages_[overlaps[i]->t_id()];

        const auto& sequence = sequences_[overlaps[i]->q_id()];
        const std::vector<WindowInterval>& breaking_points = overlaps[i]->breaking_points();

        // std::cerr << "overlap_id = " << i << "\n";
        // std::cerr << "    " << *overlaps[i] << "\n";
        // std::cerr << "All breaking points:\n";
        // for (uint32_t j = 0; j < breaking_points.size(); ++j) {
        //     const auto& bp = breaking_points[j];
        //     std::cerr << "[j = " << j << "] bp = " << bp << "\n";
        // }

        for (uint32_t j = 0; j < breaking_points.size(); ++j) {
            const auto& bp = breaking_points[j];
            // const uint32_t win_t_start = std::get<0>(breaking_points[j]);
            // const uint32_t win_t_end = std::get<0>(breaking_points[j + 1]);
            // const uint32_t win_q_start = std::get<1>(breaking_points[j]);
            // const uint32_t win_q_end = std::get<1>(breaking_points[j + 1]);
            const auto& win_t_start = bp.t_start;
            const auto& win_t_end = bp.t_end;
            const auto& win_q_start = bp.q_start;
            const auto& win_q_end = bp.q_end;

            if ((win_q_end - win_q_start) < 0.02 * window_length_) {
                continue;
            }

            if (!sequence->quality().empty() ||
                !sequence->reverse_quality().empty()) {

                const auto& quality = overlaps[i]->strand() ?
                    sequence->reverse_quality() : sequence->quality();
                double average_quality = 0;
                for (uint32_t k = win_q_start; k < win_q_end; ++k) {
                    average_quality += static_cast<uint32_t>(quality[k]) - 33;
                }
                average_quality /= (win_q_end - win_q_start);

                if (average_quality < quality_threshold_) {
                    continue;
                }
            }

            uint64_t window_id = bp.window_id;
            uint32_t window_start = windows_[bp.window_id]->start();

            const char* data = overlaps[i]->strand() ?
                &(sequence->reverse_complement()[win_q_start]) :
                &(sequence->data()[win_q_start]);
            uint32_t data_length = win_q_end - win_q_start;

            const char* quality = overlaps[i]->strand() ?
                (sequence->reverse_quality().empty() ?
                    nullptr : &(sequence->reverse_quality()[win_q_start]))
                :
                (sequence->quality().empty() ?
                    nullptr : &(sequence->quality()[win_q_start]));
            uint32_t quality_length = quality == nullptr ? 0 : data_length;

            // std::cerr << "[j = " << j << "] Adding layer for bp = " << bp << "\n";

            windows_[window_id]->add_layer(data, data_length,
                quality, quality_length,
                win_t_start - window_start,
                win_t_end - window_start - 1);
        }
        // std::cerr << "\n";

        overlaps[i].reset();
    }

}

void Polisher::find_overlap_breaking_points(std::vector<std::unique_ptr<Overlap>>& overlaps,
    const std::unordered_map<int64_t, std::vector<std::tuple<int64_t, int64_t, int64_t>>>& windows)
{
    std::vector<std::future<void>> thread_futures;
    for (uint64_t i = 0; i < overlaps.size(); ++i) {

        thread_futures.emplace_back(thread_pool_->submit(
            [&](uint64_t j) -> void {
                auto it = windows.find(overlaps[j]->t_id());
                if (it != windows.end()) {
                    overlaps[j]->find_breaking_points(sequences_, it->second);
                }
            }, i));
    }

    uint32_t logger_step = thread_futures.size() / 20;
    for (uint64_t i = 0; i < thread_futures.size(); ++i) {
        thread_futures[i].wait();
        if (logger_step != 0 && (i + 1) % logger_step == 0 && (i + 1) / logger_step < 20) {
            logger_->bar("[racon::Polisher::initialize] aligning overlaps");
        }
    }
    if (logger_step != 0) {
        logger_->bar("[racon::Polisher::initialize] aligning overlaps");
    } else {
        logger_->log("[racon::Polisher::initialize] aligned overlaps");
    }
}

void Polisher::polish(std::vector<std::unique_ptr<Sequence>>& dst,
    bool drop_unpolished_sequences) {

    logger_->log();

    std::vector<std::future<bool>> thread_futures;
    for (uint64_t i = 0; i < windows_.size(); ++i) {
        thread_futures.emplace_back(thread_pool_->submit(
            [&](uint64_t j) -> bool {
                auto it = thread_to_id_.find(std::this_thread::get_id());
                if (it == thread_to_id_.end()) {
                    fprintf(stderr, "[racon::Polisher::polish] error: "
                        "thread identifier not present!\n");
                    exit(1);
                }
                return windows_[j]->generate_consensus(
                    alignment_engines_[it->second], trim_, produce_liftover_);
            }, i));
    }

    std::string polished_data = "";
    uint32_t num_polished_windows = 0;

    uint64_t logger_step = thread_futures.size() / 20;

    uint64_t prev_window_end = 0;

    for (uint64_t i = 0; i < thread_futures.size(); ++i) {
        thread_futures[i].wait();

        num_polished_windows += thread_futures[i].get() == true ? 1 : 0;

        // BED region related: Add the sequence in between windows.
        if (windows_[i]->start() > prev_window_end) {
            uint64_t span = windows_[i]->start() - prev_window_end;
            polished_data += sequences_[windows_[i]->id()]->data().substr(prev_window_end, span);
        }

        // Add the window consensus.
        polished_data += windows_[i]->consensus();

        if (i == windows_.size() - 1 || windows_[i + 1]->rank() == 0) {
            // BED region related: Append the remaining suffix from the last window to the end of the target.
            uint32_t tlen = sequences_[windows_[i]->id()]->data().size();
            if (windows_[i]->end() < tlen) {
                uint64_t suffix_start = windows_[i]->end();
                polished_data += sequences_[windows_[i]->id()]->data().substr(suffix_start);
            }

            double polished_ratio = num_polished_windows /
                static_cast<double>(windows_[i]->rank() + 1);

            if (!drop_unpolished_sequences || polished_ratio > 0) {
                std::string tags = type_ == PolisherType::kF ? "r" : "";
                tags += " LN:i:" + std::to_string(polished_data.size());
                tags += " RC:i:" + std::to_string(targets_coverages_[windows_[i]->id()]);
                tags += " XC:f:" + std::to_string(polished_ratio);
                dst.emplace_back(createSequence(sequences_[windows_[i]->id()]->name() +
                    tags, polished_data));
            }

            num_polished_windows = 0;
            polished_data.clear();
        }
        prev_window_end = windows_[i]->end();
        windows_[i].reset();

        if (logger_step != 0 && (i + 1) % logger_step == 0 && (i + 1) / logger_step < 20) {
            logger_->bar("[racon::Polisher::polish] generating consensus");
        }
    }

    // Write the original sequences if there were no BED windows assigned to them.
    // If the BED is used, then some sequences can have 0 windows, and consensus will not
    // be generated.
    if (use_bed_) {
        // There is an extra element because of legacy code.
        for (int32_t t_id = 0; t_id < (static_cast<int32_t>(id_to_first_window_id_.size()) - 1); ++t_id) {
            if (id_to_first_window_id_[t_id] < 0) {
                std::string tags = type_ == PolisherType::kF ? "r" : "";
                tags += " LN:i:" + std::to_string(sequences_[t_id]->data().size());
                tags += " RC:i:" + std::to_string(targets_coverages_[t_id]);
                tags += " XC:f:0.0";
                dst.emplace_back(createSequence(sequences_[t_id]->name() +
                    tags, sequences_[t_id]->data()));

            }
        }
    }

    if (logger_step != 0) {
        logger_->bar("[racon::Polisher::polish] generating consensus");
    } else {
        logger_->log("[racon::Polisher::polish] generated consensus");
    }

    std::vector<std::shared_ptr<Window>>().swap(windows_);
    std::vector<std::unique_ptr<Sequence>>().swap(sequences_);
}

}
