/*!
 * @file cigar.hpp
 *
 * @brief Data structures for storing the CIGAR operations.
 */

#pragma once

#include <cstdint>
#include <vector>
#include <string>

namespace racon {

struct CigarOperation {
    char op;
    int32_t count;

    bool operator==(const CigarOperation& rhs) const
    {
        return op == rhs.op && count == rhs.count;
    }
};

using Cigar = std::vector<CigarOperation>;

Cigar ParseCigarString(const std::string& cigarStr);

}
