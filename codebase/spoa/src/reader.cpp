/*!
 * @file reader.cpp
 *
 * @brief Reader class source file
 */

#include <assert.h>

#include "chain.hpp"
#include "reader.hpp"

constexpr uint32_t kBufferSize = 1024 * 1024;
constexpr uint32_t kArraySize = 65000;

std::unique_ptr<Reader> createReader(const std::string& path) {

    auto input_file = fopen(path.c_str(), "r");
    assert(input_file);

    return std::unique_ptr<Reader>(new Reader(input_file));
}

Reader::Reader(FILE* input_file)
        : input_file_(input_file, fclose), buffer_(kBufferSize, '0'),
        num_chains_read_(0) {
}

bool Reader::read_chains(ChainSet& dst, size_t max_bytes) {

    /* Code taken from SW# (author: Matija Korpar) */

    bool status = false;
    size_t bytes_read = 0;
    size_t bytes_over = 0;

    auto input_file = input_file_.get();

    bool is_name = true;
    bool is_end = feof(input_file);

    char name[kArraySize];
    uint32_t name_length = 0;

    char data[kArraySize];
    uint32_t data_length = 0;

    while (!is_end) {

        uint32_t read = fread(buffer_.data(), sizeof(char), kBufferSize, input_file);
        is_end = feof(input_file);

        bytes_read += read;

        if (max_bytes != 0 && bytes_read > max_bytes) {
            fseek(input_file, -(bytes_over + read), SEEK_CUR);
            status = true;
            break;
        }

        for (uint32_t i = 0; i < read; ++i) {

            auto c = buffer_[i];

            if (!is_name && (c == '>' || (is_end && i == read - 1))) {

                bytes_over = 0;
                is_name = true;

                dst.emplace_back(createChain(num_chains_read_++, name, name_length,
                    data, data_length));

                name_length = data_length = 0;
            }

            if (is_name) {
                if (c == '\n') {
                    is_name = false;
                } else if (name_length == kArraySize)  {
                    continue;
                } else if (!(name_length == 0 && (c == '>' || isspace(c)))) {
                    if (c != '\r') {
                        name[name_length++] = c;
                    }
                }
            } else {
                data[data_length++] = c;
            }

            bytes_over++;
        }
    }

    return status;
}
