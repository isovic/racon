/*!
 * @file cigar.cpp
 *
 * @brief CIGAR data structures and tools.
 */

#include <cctype>
#include <memory>
#include <string>

#include "cigar.hpp"

namespace racon {

Cigar ParseCigarString(const std::string& cigarStr) {
    Cigar ret;
    for (size_t i = 0, j = 0; i < cigarStr.size(); ++i) {
        if (std::isdigit(cigarStr[i]) == 0) {
            if (j == i) {
                throw std::runtime_error("Invalid CIGAR string passed to ParseCigarString, there is no length value. CIGAR: '" + cigarStr + "'.");
            }
            char op = cigarStr[i];
            int32_t count = atoi(&cigarStr[j]);
            j = i + 1;
            ret.emplace_back(CigarOperation{op, count});
        }
    }
    if (ret.empty() && cigarStr.size() > 0) {
        throw std::runtime_error("Invalid CIGAR string passed to ParseCigarString, there are no ops specified. CIGAR: '" + cigarStr + "'.");
    }
    return ret;
}

}
