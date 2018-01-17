/*!
 * @file window.hpp
 *
 * @brief Window class header file
 */

#pragma once

#include <stdlib.h>
#include <vector>
#include <memory>
#include <string>
#include <utility>

namespace spoa {
    class AlignmentEngine;
}

namespace racon {

enum class WindowType {
    kNGS, // Next Generation Sequencing
    kTGS // Third Generation Sequencing
};

class Window;
std::unique_ptr<Window> createWindow(uint32_t id, uint32_t rank, WindowType type,
    const std::string& backbone, const std::string& quality);

class Window {
public:
    ~Window();

    const std::string& consensus() const {
        return consensus_;
    }

    void generate_consensus(std::shared_ptr<spoa::AlignmentEngine> alignment_engine);

    void add_layer(const char* sequence, uint32_t sequence_length,
        const char* quality, uint32_t quality_length, uint32_t begin,
        uint32_t end);

    friend std::unique_ptr<Window> createWindow(uint32_t id, uint32_t rank,
        WindowType type, const std::string& backbone, const std::string& quality);
private:
    Window(uint32_t id, uint32_t rank, WindowType type, const std::string& backbone,
        const std::string& quality);
    Window(const Window&) = delete;
    const Window& operator=(const Window&) = delete;

    uint32_t id_;
    uint32_t rank_;
    WindowType type_;
    std::string consensus_;
    std::vector<std::string> sequences_;
    std::vector<std::string> qualities_;
    std::vector<std::pair<uint32_t, uint32_t>> positions_;
};

}
