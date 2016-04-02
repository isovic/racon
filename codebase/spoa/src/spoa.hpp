/*!
 * @file poa.hpp
 *
 * @brief Poa header file which encapsulates the implementation
 */

#pragma once

#include <string>
#include <vector>

#include "alignment.hpp"

namespace SPOA {

std::string generate_consensus(const std::vector<std::string>& sequences,
    AlignmentParams params, bool sorted = false);
std::string generate_consensus(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, AlignmentParams params, bool sorted = false);

}
