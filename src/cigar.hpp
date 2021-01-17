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

class CigarOperation {
public:
    char op;
    int32_t count;

    CigarOperation(char _op, int32_t _count)
        : op(_op)
        , count(_count)
    {
    }

    bool operator==(const CigarOperation& rhs) const
    {
        return op == rhs.op && count == rhs.count;
    }
};

using Cigar = std::vector<CigarOperation>;

Cigar ParseCigarString(const std::string& cigarStr);
std::string CigarToString(const Cigar& cigar);
void AddCigarEvent(Cigar& cigar, char op, int32_t count);
void MergeCigar(Cigar& cigarDest, const Cigar& cigarSrc);

}
