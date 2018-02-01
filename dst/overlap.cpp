/*!
 * @file overlap.cpp
 *
 * @brief Overlap class source file
 */

#include <algorithm>

#include "sequence.hpp"
#include "overlap.hpp"
#include "edlib.h"

namespace racon {

Overlap::Overlap(uint32_t a_id, uint32_t b_id, double, uint32_t,
    uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
    uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length)
        : q_name_(), q_id_(a_id - 1), q_begin_(a_begin), q_end_(a_end),
        q_length_(a_length), t_name_(), t_id_(b_id - 1), t_begin_(b_begin),
        t_end_(b_end), t_length_(b_length), strand_(a_rc ^ b_rc), length_(),
        error_(), cigar_(), is_valid_(true), is_transmuted_(false),
        breaking_points_(), dual_breaking_points_() {

    length_ = std::max(q_end_ - q_begin_, t_end_ - t_begin_);
    error_ = 1 - std::min(q_end_ - q_begin_, t_end_ - t_begin_) /
        static_cast<double>(length_);
}

Overlap::Overlap(const char* q_name, uint32_t q_name_length, uint32_t q_length,
    uint32_t q_begin, uint32_t q_end, char orientation, const char* t_name,
    uint32_t t_name_length, uint32_t t_length, uint32_t t_begin,
    uint32_t t_end, uint32_t, uint32_t, uint32_t)
        : q_name_(q_name, q_name_length), q_id_(), q_begin_(q_begin),
        q_end_(q_end), q_length_(q_length), t_name_(t_name, t_name_length),
        t_id_(), t_begin_(t_begin), t_end_(t_end), t_length_(t_length),
        strand_(orientation == '-'), length_(), error_(), cigar_(),
        is_valid_(true), is_transmuted_(false), breaking_points_(),
        dual_breaking_points_() {

    length_ = std::max(q_end_ - q_begin_, t_end_ - t_begin_);
    error_ = 1 - std::min(q_end_ - q_begin_, t_end_ - t_begin_) /
        static_cast<double>(length_);
}

Overlap::Overlap(const char* q_name, uint32_t q_name_length, uint32_t flag,
    const char* t_name, uint32_t t_name_length, uint32_t t_begin,
    uint32_t, const char* cigar, uint32_t cigar_length, const char*,
    uint32_t, uint32_t, uint32_t, const char*, uint32_t, const char*,
    uint32_t)
        : q_name_(q_name, q_name_length), q_id_(), q_begin_(0), q_end_(),
        q_length_(0), t_name_(t_name, t_name_length), t_id_(), t_begin_(t_begin - 1),
        t_end_(), t_length_(0), strand_(flag & 0x10), length_(), error_(),
        cigar_(cigar, cigar_length), is_valid_(!(flag & 0x4)),
        is_transmuted_(false), breaking_points_(), dual_breaking_points_() {

    if (cigar_.size() < 2 && is_valid_) {
        fprintf(stderr, "Racon::Overlap::Overlap error: "
            "missing alignment from SAM object!\n");
        exit(1);
    } else {
        for (uint32_t i = 0; i < cigar_.size(); ++i) {
            if (cigar_[i] == 'S' || cigar_[i] == 'H') {
                q_begin_ = atoi(&cigar_[0]);
                break;
            } else if (cigar_[i] == 'M' || cigar_[i] == '=' || cigar_[i] == 'I' ||
                cigar_[i] == 'D' || cigar_[i] == 'N' || cigar_[i] == 'P' ||
                cigar_[i] == 'X') {
                break;
            }
        }

        uint32_t q_alignment_length = 0, q_clip_length = 0, t_alignment_length = 0;
        for (uint32_t i = 0, j = 0; i < cigar_.size(); ++i) {
            if (cigar_[i] == 'M' || cigar_[i] == '=' || cigar_[i] == 'X') {
                auto num_bases = atoi(&cigar_[j]);
                j = i + 1;
                q_alignment_length += num_bases;
                t_alignment_length += num_bases;
            } else if (cigar_[i] == 'I') {
                auto num_bases = atoi(&cigar_[j]);
                j = i + 1;
                q_alignment_length += num_bases;
            } else if (cigar_[i] == 'D' || cigar_[i] == 'N') {
                auto num_bases = atoi(&cigar_[j]);
                j = i + 1;
                t_alignment_length += num_bases;
            } else if (cigar_[i] == 'S' || cigar_[i] == 'H') {
                q_clip_length += atoi(&cigar_[j]);
                j = i + 1;
            } else if (cigar_[i] == 'P') {
                j = i + 1;
            }
        }

        q_end_ = q_begin_ + q_alignment_length;
        q_length_ = q_clip_length + q_alignment_length;
        if (strand_) {
            uint32_t tmp = q_begin_;
            q_begin_ = q_length_ - q_end_;
            q_end_ = q_length_ - tmp;
        }

        t_end_ = t_begin_ + t_alignment_length;

        length_ = std::max(q_alignment_length, t_alignment_length);
        error_ = 1 - std::min(q_alignment_length, t_alignment_length) /
            static_cast<double>(length_);
    }
}

Overlap::Overlap()
        : q_name_(), q_id_(), q_begin_(), q_end_(), q_length_(), t_name_(),
        t_id_(), t_begin_(), t_end_(), t_length_(), strand_(), length_(),
        error_(), cigar_(), is_valid_(true), is_transmuted_(true),
        breaking_points_(), dual_breaking_points_() {
}

template<typename T>
bool transmuteId(const std::unordered_map<T, uint32_t>& hash, const T& key,
    uint32_t& id) {

    auto it = hash.find(key);
    if (it == hash.end()) {
        return false;
    }
    id = it->second;
    return true;
}

void Overlap::transmute(const std::unordered_map<std::string, uint32_t>& name_to_id,
    const std::unordered_map<uint32_t, uint32_t>& id_to_id) {

    if (!is_valid_) {
        fprintf(stderr, "racon::Overlap::transmute error: "
            "overlap is not valid!\n");
        exit(1);
    }

    if (is_transmuted_) {
         return;
    }

    if (!q_name_.empty()) {
        if (!transmuteId(name_to_id, q_name_ + "0", q_id_)) {
            fprintf(stderr, "racon::Overlap::transmute error: "
                "missing sequence with name %s!\n", q_name_.c_str());
            exit(1);
        }
    } else {
        if (!transmuteId(id_to_id, q_id_ << 1 | 0, q_id_)) {
            fprintf(stderr, "racon::Overlap::transmute error: "
                "missing sequence with id %d!\n", q_id_ >> 1);
            exit(1);
        }
    }
    if (!t_name_.empty()) {
        if (!transmuteId(name_to_id, t_name_ + "1", t_id_)) {
            fprintf(stderr, "racon::Overlap::transmute error: "
                "missing target sequence with name %s!\n", t_name_.c_str());
            exit(1);
        }
    } else {
        if (!transmuteId(id_to_id, t_id_ << 1 | 1, t_id_)) {
            fprintf(stderr, "racon::Overlap::transmute error: "
                "missing target sequence with id %d!\n", t_id_ >> 1);
            exit(1);
        }
    }

    is_transmuted_ = true;
}

void Overlap::find_breaking_points(const std::vector<std::unique_ptr<Sequence>>& sequences,
    uint32_t window_length) {

    if (!is_transmuted_) {
        fprintf(stderr, "racon::Overlap::find_breaking_points error: "
            "overlap is not transmuted!\n");
        exit(1);
    }

    if (!breaking_points_.empty()) {
        return;
    }

    if (q_length_ != sequences[q_id_]->data().size()) {
        fprintf(stderr, "racon::overlap::find_breaking_points error: "
            "unmatched query lengths!\n");
        exit(1);
    }
    t_length_ = sequences[t_id_]->data().size();

    if (cigar_.empty()) {
        // align overlaps with edlib
        const char* q = !strand_ ? &(sequences[q_id_]->data()[q_begin_]) :
            &(sequences[q_id_]->reverse_complement()[q_length_ - q_end_]);
        const char* t = &(sequences[t_id_]->data()[t_begin_]);

        EdlibAlignResult result = edlibAlign(q, q_end_ - q_begin_, t, t_end_ -
            t_begin_, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH,
            nullptr, 0));

        if (result.status == EDLIB_STATUS_OK) {
            char* cigar = edlibAlignmentToCigar(result.alignment,
                result.alignmentLength, EDLIB_CIGAR_STANDARD);
            cigar_ = cigar;
            free(cigar);
        } else {
            fprintf(stderr, "racon::Overlap::find_breaking_points error: "
                "edlib unable to align sequences (%u x %u)!\n", q_id_, t_id_);
            exit(1);
        }

        edlibFreeAlignResult(result);
    }

    // find breaking points from cigar
    std::vector<int32_t> window_ends;
    for (uint32_t i = 0; i < t_end_; i += window_length) {
        if (i > t_begin_) {
            window_ends.emplace_back(i - 1);
        }
    }
    window_ends.emplace_back(t_end_ - 1);

    std::vector<int32_t> dual_window_ends;
    if (strand_) {
        dual_window_ends.emplace_back(q_length_ - q_begin_ - 1);
    }
    for (uint32_t i = 0; i < q_end_; i += window_length) {
        if (i > q_begin_) {
            dual_window_ends.emplace_back(strand_ ? q_length_ - i - 1 : i - 1);
        }
    }
    if (strand_) {
        std::reverse(dual_window_ends.begin(), dual_window_ends.end());
    } else {
        dual_window_ends.emplace_back(q_end_ - 1);
    }

    uint32_t w = 0, d_w = 0;
    bool found_first_match = false, found_first_dual_match = false;
    std::pair<uint32_t, uint32_t> first_match = {0, 0}, last_match = {0, 0},
        first_dual_match = {0, 0}, last_dual_match = {0, 0};
    int32_t q_ptr = (strand_ ? (q_length_ - q_end_) : q_begin_) - 1,
        t_ptr = t_begin_ - 1;

    for (uint32_t i = 0, j = 0; i < cigar_.size(); ++i) {
        if (cigar_[i] == 'M' || cigar_[i] == '=' || cigar_[i] == 'X') {
            uint32_t k = 0, num_bases = atoi(&cigar_[j]);
            j = i + 1;
            while (k < num_bases) {
                ++q_ptr;
                ++t_ptr;

                if (!found_first_match) {
                    found_first_match = true;
                    first_match.first = t_ptr;
                    first_match.second = q_ptr;
                }
                last_match.first = t_ptr + 1;
                last_match.second = q_ptr + 1;
                if (t_ptr == window_ends[w]) {
                    if (found_first_match) {
                        breaking_points_.emplace_back(first_match);
                        breaking_points_.emplace_back(last_match);
                    }
                    found_first_match = false;
                    ++w;
                }

                if (!found_first_dual_match) {
                    found_first_dual_match = true;
                    first_dual_match.first = strand_ ? q_length_ - q_ptr : q_ptr;
                    first_dual_match.second = strand_ ? t_length_ - t_ptr : t_ptr;
                }
                last_dual_match.first = strand_ ? q_length_ - q_ptr - 1 : q_ptr + 1;
                last_dual_match.second = strand_ ? t_length_ - t_ptr - 1 : t_ptr + 1;
                if (q_ptr == dual_window_ends[d_w]) {
                    if (found_first_dual_match) {
                        dual_breaking_points_.emplace_back(first_dual_match);
                        dual_breaking_points_.emplace_back(last_dual_match);
                    }
                    found_first_dual_match = false;
                    ++d_w;
                }
                ++k;
            }
        } else if (cigar_[i] == 'I') {
            uint32_t k = 0, num_bases = atoi(&cigar_[j]);
            j = i + 1;
            while (k < num_bases) {
                ++q_ptr;
                if (q_ptr == dual_window_ends[d_w]) {
                    if (found_first_dual_match) {
                        dual_breaking_points_.emplace_back(first_dual_match);
                        dual_breaking_points_.emplace_back(last_dual_match);
                    }
                    found_first_dual_match = false;
                    ++d_w;
                }
                ++k;
            }
        } else if (cigar_[i] == 'D' || cigar_[i] == 'N') {
            uint32_t k = 0, num_bases = atoi(&cigar_[j]);
            j = i + 1;
            while (k < num_bases) {
                ++t_ptr;
                if (t_ptr == window_ends[w]) {
                    if (found_first_match) {
                        breaking_points_.emplace_back(first_match);
                        breaking_points_.emplace_back(last_match);
                    }
                    found_first_match = false;
                    ++w;
                }
                ++k;
            }
        } else if (cigar_[i] == 'S' || cigar_[i] == 'H' || cigar_[i] == 'P') {
            j = i + 1;
        }
    }
    if (strand_) {
        std::reverse(dual_breaking_points_.begin(), dual_breaking_points_.end());
    }

    std::string().swap(cigar_);
}

std::unique_ptr<Overlap> Overlap::dual_overlap() {

    if (!is_transmuted_) {
        fprintf(stderr, "racon::Overlap::dual_overlap error: "
            "overlap is not transmuted!\n");
        exit(1);
    }
    if (dual_breaking_points_.empty()) {
        fprintf(stderr, "racon::Overlap::dual_overlap error: "
            "dual breaking points unavailable!\n");
        exit(1);
    }

    auto other = std::unique_ptr<Overlap>(new Overlap());
    other->q_id_ = t_id_;
    other->t_id_ = q_id_;
    other->strand_ = strand_;
    other->breaking_points_ = dual_breaking_points_;

    std::vector<std::pair<uint32_t, uint32_t>>().swap(dual_breaking_points_);

    return std::move(other);
}

}
