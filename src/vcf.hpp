/*!
 * @file vcf.hpp
 *
 * @brief Data structures and methods for manipulating the VCF entries.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include "cigar.hpp"

namespace racon {

struct VcfDiff {
    int32_t pos;
    std::string ref;
    std::string alt;
};

std::vector<VcfDiff> ExtractVCFEventsFromCigarString(const racon::Cigar& cigar, const std::string& qseq, const std::string& tseq);

}
