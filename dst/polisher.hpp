/*!
 * @file polisher.hpp
 *
 * @brief Polisher class header file
 */

#pragma once

#include <stdlib.h>
#include <vector>
#include <memory>

namespace bioparser {
    template<class T>
    class Reader;
}

namespace racon {

class Window;

class Polisher;
std::unique_ptr<Polisher> createPolisher();

class Polisher {
public:
    ~Polisher();

    void initialize();

    void polish();

    friend std::unique_ptr<Polisher> createPolisher();
private:
    Polisher();
    Polisher(const Polisher&) = delete;
    const Polisher& operator=(const Polisher&) = delete;

    std::vector<std::unique_ptr<Window>> windows_;
};

}
