/*!
 * @file cigar.hpp
 *
 * @brief Data structures for storing the CIGAR operations.
 */

#pragma once

#include <cstdint>
#include <vector>

namespace racon {

struct CigarOperation {
    char op;
    int32_t count;
};

using Cigar = std::vector<CigarOperation>;

}
