/**
 * @file poa.hpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief POA programming interface header file 
 */
#include <vector>
#include <string>

using std::string;
using std::vector;

/**
 * @brief Generates consensus using poa
 * @details Method generates consensus sequence
 * of sequences given as parameter
 * 
 * @param sequences sequences for alignment
 * @return consensus sequence
 */
namespace POA{
string poa_consensus(const vector<string>& sequences);
}
