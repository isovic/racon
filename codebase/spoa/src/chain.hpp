/*!
 * @file chain.hpp
 *
 * @brief Chain class header file
 */

#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <memory>
#include <vector>
#include <string>

class Reader;
class Chain;

using ChainSet = std::vector<std::unique_ptr<Chain>>;

std::unique_ptr<Chain> createChain(uint32_t id, char* name, uint32_t name_length,
    char* data, uint32_t data_length);

void createChainSet(ChainSet& dst, const std::string& path);

std::unique_ptr<Reader> createChainSetPartInitialize(const std::string& path);

bool createChainSetPart(ChainSet& dst, std::shared_ptr<Reader> reader, size_t max_bytes);

class Chain {
public:

    ~Chain() = default;

    uint32_t id() const {
        return id_;
    }

    const std::string& name() const {
        return name_;
    }

    const size_t name_length() const {
        return name_.size();
    }

    const std::string& data() const {
        return data_;
    }

    const uint32_t length() const {
        return data_.size();
    }

    friend std::unique_ptr<Chain> createChain(uint32_t id, char* name,
        uint32_t name_length, char* data, uint32_t data_length);

private:

    Chain(uint32_t id, std::string&& name, std::string&& data);
    Chain(const Chain&) = delete;
    const Chain& operator=(const Chain&) = delete;

    uint32_t id_;
    std::string name_;
    std::string data_;
};
