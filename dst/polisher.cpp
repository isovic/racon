/*!
 * @file polisher.cpp
 *
 * @brief Polisher class source file
 */

#include "window.hpp"
#include "polisher.hpp"

namespace racon {

std::unique_ptr<Polisher> createPolisher() {
    return std::unique_ptr<Polisher>(new Polisher());
}

Polisher::Polisher()
        : windows_() {
}

Polisher::~Polisher() {
}

}
