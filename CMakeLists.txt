cmake_minimum_required(VERSION 3.2)
project(racon)

# The version number.
set (RACON_VERSION_MAJOR 1)
set (RACON_VERSION_MINOR 3)
set (RACON_VERSION_PATCH 3)
# Get the commit ID.
execute_process(
  COMMAND git log -1 --pretty=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE RACON_VERSION_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Configure the version file.
configure_file("${PROJECT_SOURCE_DIR}/src/version.hpp.in" "${PROJECT_BINARY_DIR}/generated/version.hpp")
# Add the include directory for the generated file.
include_directories(${PROJECT_BINARY_DIR}/generated)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(racon_build_tests "Build racon unit tests" OFF)
option(racon_build_wrapper "Build racon wrapper" OFF)

add_executable(racon
    src/main.cpp
    src/polisher.cpp
    src/overlap.cpp
    src/sequence.cpp
    src/window.cpp)

if (NOT TARGET bioparser)
    add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
endif()
if (NOT TARGET spoa)
    add_subdirectory(vendor/spoa EXCLUDE_FROM_ALL)
endif()
if (NOT TARGET thread_pool)
    add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)
endif()
if (NOT TARGET edlib)
    add_subdirectory(vendor/edlib EXCLUDE_FROM_ALL)
endif()
if (NOT TARGET logger)
    add_subdirectory(vendor/logger EXCLUDE_FROM_ALL)
endif()

target_link_libraries(racon bioparser spoa thread_pool edlib_static logger)

install(TARGETS racon DESTINATION bin)

if (racon_build_tests)
    set(racon_test_data_path ${PROJECT_SOURCE_DIR}/test/data/)
    configure_file("${PROJECT_SOURCE_DIR}/test/racon_test_config.h.in"
        "${PROJECT_BINARY_DIR}/config/racon_test_config.h")
    include_directories(${PROJECT_BINARY_DIR}/config)
    include_directories(${PROJECT_SOURCE_DIR}/src)

    add_executable(racon_test
        test/racon_test.cpp
        src/polisher.cpp
        src/overlap.cpp
        src/sequence.cpp
        src/window.cpp)

    if (NOT TARGET gtest_main)
        add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)
    endif()

    target_link_libraries(racon_test bioparser spoa thread_pool edlib_static logger gtest_main)
endif()

if (racon_build_wrapper)
    set(racon_path ${PROJECT_BINARY_DIR}/bin/racon)
    set(rampler_path ${PROJECT_BINARY_DIR}/vendor/rampler/bin/rampler)
    configure_file(${PROJECT_SOURCE_DIR}/scripts/racon_wrapper.py
        ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/racon_wrapper)
    file(COPY ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/racon_wrapper
        DESTINATION ${PROJECT_BINARY_DIR}/bin
        FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)

    if (NOT TARGET rampler)
        add_subdirectory(vendor/rampler)
    endif()
endif()
