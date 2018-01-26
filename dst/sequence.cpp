/*!
 * @file sequence.cpp
 *
 * @brief Sequence class source file
 */

#include "sequence.hpp"

namespace racon {

Sequence::Sequence(const char* name, uint32_t name_length, const char* data,
    uint32_t data_length)
        : name_(name, name_length), data_(data, data_length),
        reverse_complement_(), quality_(), reverse_quality_() {
}

Sequence::Sequence(const char* name, uint32_t name_length, const char* data,
    uint32_t data_length, const char* quality, uint32_t quality_length)
        : name_(name, name_length), data_(data, data_length),
        reverse_complement_(), quality_(quality, quality_length),
        reverse_quality_() {

    uint32_t quality_sum = 0;
    for (const auto& it: quality_) {
        quality_sum += it - '!';
    }

    if (quality_sum == 0) {
        std::string().swap(quality_);
    }
}

void Sequence::create_reverse_complement() {

    reverse_complement_.clear();
    reverse_complement_.reserve(data_.size());

    for (int32_t i = data_.size() - 1; i >= 0; --i) {
        switch (data_[i]) {
            case 'A':
                reverse_complement_ += 'T';
                break;
            case 'T':
                reverse_complement_ += 'A';
                break;
            case 'C':
                reverse_complement_ += 'G';
                break;
            case 'G':
                reverse_complement_ += 'C';
                break;
            default:
                break;
        }
    }

    reverse_quality_.clear();
    reverse_quality_.reserve(quality_.size());

    for (int32_t i = quality_.size() - 1; i >= 0; --i) {
        reverse_quality_ += quality_[i];
    }
}

}
