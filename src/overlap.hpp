/*!
 * @file overlap.hpp
 *
 * @brief Overlap class header file
 */

#pragma once

#include "util.hpp"
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>

namespace bioparser {
    template<class T>
    class MhapParser;

    template<class T>
    class PafParser;

    template<class T>
    class SamParser;
}

namespace racon {

class Sequence;

class Overlap {
public:
    ~Overlap() = default;

    uint32_t q_id() const {
        return q_id_;
    }

    uint32_t t_id() const {
        return t_id_;
    }

    uint32_t t_begin() const {
        return t_begin_;
    }

    uint32_t t_end() const {
        return t_end_;
    }

    uint32_t strand() const {
        return strand_;
    }

    bool is_valid() const {
        return is_valid_;
    }

    void transmute(const std::vector<std::unique_ptr<Sequence>>& sequences,
        const std::unordered_map<std::string, uint64_t>& name_to_id,
        const std::unordered_map<uint64_t, uint64_t>& id_to_id);

    uint32_t length() const {
        return length_;
    }

    double error() const {
        return error_;
    }

    const std::string& cigar() const {
        return cigar_;
    }

    const std::vector<WindowInterval>& breaking_points() const {
        return breaking_points_;
    }

    void find_breaking_points(const std::vector<std::unique_ptr<Sequence>>& sequences,
        std::vector<std::tuple<int64_t, int64_t, int64_t>> windows);

    friend bioparser::MhapParser<Overlap>;
    friend bioparser::PafParser<Overlap>;
    friend bioparser::SamParser<Overlap>;

    friend std::ostream& operator<<(::std::ostream& os, const Overlap& a);

#ifdef CUDA_ENABLED
    friend class CUDABatchAligner;
#endif
private:
    Overlap(uint64_t a_id, uint64_t b_id, double accuracy, uint32_t minmers,
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
    Overlap();
    Overlap(const Overlap&) = delete;
    const Overlap& operator=(const Overlap&) = delete;
    virtual void find_breaking_points_from_cigar(std::vector<std::tuple<int64_t, int64_t, int64_t>> windows);
    virtual void align_overlaps(const char* q, uint32_t q_len, const char* t, uint32_t t_len);

    std::string q_name_;
    uint64_t q_id_;
    uint32_t q_begin_;
    uint32_t q_end_;
    uint32_t q_length_;

    std::string t_name_;
    uint64_t t_id_;
    uint32_t t_begin_;
    uint32_t t_end_;
    uint32_t t_length_;

    uint32_t strand_;
    uint32_t length_;
    double error_;
    std::string cigar_;

    bool is_valid_;
    bool is_transmuted_;
    // std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> breaking_points_;
    std::vector<std::pair<uint32_t, uint32_t>> dual_breaking_points_;

    std::vector<WindowInterval> breaking_points_;
};

inline std::ostream& operator<<(::std::ostream& os, const Overlap& a) {
    std::string delim(" ");
    os << a.q_id_ << delim << a.q_begin_ << delim << a.q_end_ << delim << a.q_length_
        << delim << (a.strand_ ? "-" : "+")
        << delim << a.t_id_ << delim << a.t_begin_ << delim << a.t_end_ << delim << a.t_length_
        << delim << a.cigar_ << delim << a.is_valid_ << delim << a.is_transmuted_;
    return os;
}

}
