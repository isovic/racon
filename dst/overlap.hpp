/*!
 * @file overlap.hpp
 *
 * @brief Overlap class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>
#include <string>
#include <utility>

namespace bioparser {
    template<class T>
    class MhapParser;

    template<class T>
    class PafParser;

    template<class T>
    class SamParser;
}

namespace racon {

class Overlap {
public:

    ~Overlap() = default;

    friend bioparser::MhapParser<Overlap>;
    friend bioparser::PafParser<Overlap>;
    friend bioparser::SamParser<Overlap>;

private:

    Overlap(uint32_t a_id, uint32_t b_id, double accuracy, uint32_t minmers,
        uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
        uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length);
    Overlap(const char* q_name, uint32_t q_name_length, uint32_t q_length,
        uint32_t q_begin, uint32_t q_end, char orientation, const char* t_name,
        uint32_t t_name_length, uint32_t t_length, uint32_t t_begin,
        uint32_t t_end, uint32_t matching_bases, uint32_t overlap_length,
        uint32_t maping_quality);
    Overlap(const char* q_name, uint32_t q_name_length, uint32_t flag,
        const char* t_name, uint32_t t_name_length, uint32_t t_begin,
        uint32_t mapping_quality, const char* cigar, uint32_t cigar_length,
        const char* t_next_name, uint32_t t_next_name_length,
        uint32_t t_next_begin, uint32_t template_length, const char* sequence,
        uint32_t sequence_length, const char* quality, uint32_t quality_length);
    Overlap(const Overlap&) = delete;
    const Overlap& operator=(const Overlap&) = delete;

    std::string q_name_;
    uint32_t q_id_;
    uint32_t q_begin_;
    uint32_t q_end_;
    uint32_t q_length_;
    std::string q_sequence_;
    std::string q_quality_;

    std::string t_name_;
    uint32_t t_id_;
    uint32_t t_begin_;
    uint32_t t_end_;
    uint32_t t_length_;

    uint32_t strand_;
    double error_;

    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> breaking_points_;
};

}
