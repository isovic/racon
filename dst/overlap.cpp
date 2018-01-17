/*!
 * @file overlap.cpp
 *
 * @brief Overlap class source file
 */

#include "overlap.hpp"

namespace racon {

Overlap::Overlap(uint32_t a_id, uint32_t b_id, double accuracy, uint32_t,
    uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
    uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length)
        : q_name_(), q_id_(a_id), q_begin_(a_begin), q_end_(a_end),
        q_length_(a_length), q_sequence_(), q_quality_(), t_name_(), t_id_(b_id),
        t_begin_(b_begin), t_end_(b_end), t_length_(b_length),
        strand_(a_rc ^ b_rc), error_(1.0 - accuracy), breaking_points_() {

}

Overlap::Overlap(const char* q_name, uint32_t q_name_length, uint32_t q_length,
    uint32_t q_begin, uint32_t q_end, char orientation, const char* t_name,
    uint32_t t_name_length, uint32_t t_length, uint32_t t_begin,
    uint32_t t_end, uint32_t matching_bases, uint32_t overlap_length,
    uint32_t)
        : q_name_(q_name, q_name_length), q_id_(), q_begin_(q_begin),
        q_end_(q_end), q_length_(q_length), q_sequence_(), q_quality_(),
        t_name_(t_name, t_name_length), t_id_(), t_begin_(t_begin),
        t_end_(t_end), t_length_(t_length), strand_(orientation == '-'),
        error_(1 - (matching_bases / static_cast<double>(overlap_length))),
        breaking_points_() {

}

Overlap::Overlap(const char* q_name, uint32_t q_name_length, uint32_t flag,
    const char* t_name, uint32_t t_name_length, uint32_t t_begin,
    uint32_t, const char* cigar, uint32_t cigar_length, const char*,
    uint32_t, uint32_t, uint32_t, const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length)
        : q_name_(q_name, q_name_length), q_id_(), q_begin_(), q_end_(),
        q_length_(), q_sequence_(), q_quality_(), t_name_(t_name, t_name_length),
        t_id_(), t_begin_(t_begin), t_end_(), t_length_(), strand_(flag & 0x10),
        error_(), breaking_points_() {

}

}
