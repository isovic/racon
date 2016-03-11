/*!
 * @file reader.hpp
 *
 * @brief Reader class header file
 */

#pragma once

#include <stdlib.h>
#include <memory>
#include <vector>
#include <string>

class Chain;
class Reader;

using ChainSet = std::vector<std::unique_ptr<Chain>>;

std::unique_ptr<Reader> createReader(const std::string& path);

class Reader {
public:

    ~Reader() = default;

    bool read_chains(ChainSet& dst, size_t max_bytes);

	friend std::unique_ptr<Reader> createReader(const std::string& path);

private:

	Reader(FILE* input_file);
	Reader(const Reader&) = delete;
	const Reader& operator=(const Reader&) = delete;

    std::unique_ptr<FILE, int (*)(FILE*)> input_file_;
    std::vector<char> buffer_;
    uint32_t num_chains_read_;
};
