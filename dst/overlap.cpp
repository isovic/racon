/*!
 * @file overlap.cpp
 *
 * @brief Overlap class source file
 */

#include "overlap.hpp"

namespace racon {

Overlap::Overlap(uint32_t a_id, uint32_t b_id, double, uint32_t,
    uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
    uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length)
        : q_name_(), q_id_((a_id - 1) << 1 | 0), q_begin_(a_begin),
        q_end_(a_end), q_length_(a_length), t_name_(),
        t_id_((b_id - 1) << 1 | 1), t_begin_(b_begin), t_end_(b_end),
        t_length_(b_length), strand_(a_rc ^ b_rc), error_(),
        cigar_(), is_valid_(true), is_transmuted_(false), breaking_points_() {

    error_ = 1 - std::min(a_end - a_begin, b_end - b_begin) /
        static_cast<double>(std::max(a_end - a_begin, b_end - b_begin));
}

Overlap::Overlap(const char* q_name, uint32_t q_name_length, uint32_t q_length,
    uint32_t q_begin, uint32_t q_end, char orientation, const char* t_name,
    uint32_t t_name_length, uint32_t t_length, uint32_t t_begin,
    uint32_t t_end, uint32_t, uint32_t, uint32_t)
        : q_name_(q_name, q_name_length), q_id_(), q_begin_(q_begin),
        q_end_(q_end), q_length_(q_length), t_name_(t_name, t_name_length),
        t_id_(), t_begin_(t_begin), t_end_(t_end), t_length_(t_length),
        strand_(orientation == '-'), error_(), cigar_(), is_valid_(true),
        is_transmuted_(false), breaking_points_() {

    error_ = 1 - std::min(q_end - q_begin, t_end - t_begin) /
        static_cast<double>(std::max(q_end - q_begin, t_end - t_begin));
}

Overlap::Overlap(const char* q_name, uint32_t q_name_length, uint32_t flag,
    const char* t_name, uint32_t t_name_length, uint32_t t_begin,
    uint32_t, const char* cigar, uint32_t cigar_length, const char*,
    uint32_t, uint32_t, uint32_t, const char*, uint32_t, const char*,
    uint32_t)
        : q_name_(q_name, q_name_length), q_id_(), q_begin_(), q_end_(),
        q_length_(), t_name_(t_name, t_name_length), t_id_(), t_begin_(t_begin),
        t_end_(), t_length_(), strand_(flag & 0x10), error_(),
        cigar_(cigar, cigar_length), is_valid_(!(flag & 0x4)),
        is_transmuted_(false), breaking_points_() {

    if (cigar_length > 1) {
        uint32_t q_length = 0, t_length = 0;
        for (uint32_t i = 0, j = 0; i < cigar_length; ++i) {
            if (cigar[i] == 'M' || cigar[i] == '=' || cigar[i] == 'X') {
                auto num_bases = atoi(&cigar[j]);
                j = i + 1;
                q_length += num_bases;
                t_length += num_bases;
            } else if (cigar[i] == 'I') {
                auto num_bases = atoi(&cigar[j]);
                j = i + 1;
                q_length += num_bases;
            } else if (cigar[i] == 'D' || cigar[i] == 'N') {
                auto num_bases = atoi(&cigar[j]);
                j = i + 1;
                t_length += num_bases;
            } else if (cigar[i] == 'S' || cigar[i] == 'H' || cigar[i] == 'P') {
                j = i + 1;
            }
        }
        error_ = 1 - std::min(q_length, t_length) / static_cast<double>(std::max(
            q_length, t_length));
    }
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

void Overlap::transmute(const std::unordered_map<std::string, uint32_t>& name_to_nid,
    const std::unordered_map<uint32_t, uint32_t>& id_to_nid) {

    if (is_transmuted_) {
         return;
    }

    if (!q_name_.empty()) {
        if (!transmuteId(name_to_nid, q_name_, q_id_)) {
            fprintf(stderr, "racon::Overlap::transmute error: "
                "missing sequence with name %s!\n", q_name_.c_str());
            exit(1);
        }
    } else {
        if (!transmuteId(id_to_nid, q_id_, q_id_)) {
            fprintf(stderr, "racon::Overlap::transmute error: "
                "missing sequence with id %d!\n", q_id_ >> 1);
            exit(1);
        }
    }
    if (!t_name_.empty()) {
        if (!transmuteId(name_to_nid, t_name_, t_id_)) {
            fprintf(stderr, "racon::Overlap::transmute error: "
                "missing target sequence with name %s!\n", t_name_.c_str());
            exit(1);
        }
    } else {
        if (!transmuteId(id_to_nid, t_id_, t_id_)) {
            fprintf(stderr, "racon::Overlap::transmute error: "
                "missing target sequence with id %d!\n", t_id_ >> 1);
            exit(1);
        }
    }

    is_transmuted_ = true;
}

}
