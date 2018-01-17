cmake_minimum_required(VERSION 3.2)
project(racon)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(racon_build_tests "Build racon google tests" OFF)

add_executable(racon
    dst/main.cpp
    dst/polisher.cpp
    dst/overlap.cpp
    dst/sequence.cpp
    dst/window.cpp)

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/spoa EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

target_link_libraries(racon bioparser spoa thread_pool pthread)

install(TARGETS racon DESTINATION bin)

if (racon_build_tests)
endif(racon_build_tests)
