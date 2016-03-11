/*!
 * @file chain.cpp
 *
 * @brief Chain class source file
 */

#include <assert.h>

#include "reader.hpp"
#include "chain.hpp"

std::vector<char> kCoder = {
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  65,  66,  67,  68,  69,
    70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
    80,  81,  82,  83,  84,  85,  86,  87,  88,  89,
    90,  -1,  -1,  -1,  -1,  -1,  -1,  65,  66,  67,
    68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
    78,  79,  80,  81,  82,  83,  84,  85,  86,  87,
    88,  89,  90,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
    -1,  -1,  -1,  -1,  -1
};

std::unique_ptr<Chain> createChain(uint32_t id, char* name, uint32_t name_length,
    char* data, uint32_t data_length) {

    assert(name_length);
    assert(data_length);

    // remove trailing white spaces
    while (isspace(name[name_length - 1])) {
        --name_length;
    }
    assert(name_length && "name cannot be empty");

    std::string data_;
    data_.reserve(data_length);
    uint32_t data_ptr = 0, valid_data_length = 0;

    for (uint32_t i = 0; i < data_length; ++i) {
        auto c = kCoder[data[i]];
        if (c != -1) {
            data[data_ptr++] = c;
            ++valid_data_length;
        }
    }
    assert(valid_data_length && "no valid chars found");

    return std::unique_ptr<Chain>(new Chain(id, std::string(name, name_length),
        std::string(data, valid_data_length)));
}

void createChainSet(ChainSet& dst, const std::string& path) {

    auto reader = createChainSetPartInitialize(path);
    createChainSetPart(dst, std::move(reader), 0);
}

std::unique_ptr<Reader> createChainSetPartInitialize(const std::string& path) {

    /* maybe in future: check if the chains are cached */
    return createReader(path);
}

bool createChainSetPart(ChainSet& dst, std::shared_ptr<Reader> reader, size_t max_bytes) {

    assert(reader);
    return reader->read_chains(dst, max_bytes);
}

Chain::Chain(uint32_t id, std::string&& name, std::string&& data)
        : id_(id), name_(name), data_(data) {
}
