cmake_minimum_required(VERSION 3.2)
project(racon)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

option(racon_build_tests "Build racon google tests" OFF)

add_executable(racon
    dst/main.cpp
    dst/polisher.cpp)

if (racon_build_tests)
endif(racon_build_tests)

install(TARGETS racon DESTINATION bin)
