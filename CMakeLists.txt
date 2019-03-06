cmake_minimum_required(VERSION 3.2)
project(racon)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(racon_build_tests "Build racon unit tests" OFF)
option(racon_build_wrapper "Build racon wrapper" OFF)
option(racon_enable_cuda "Build racon with NVIDIA CUDA support" OFF)

# Check CUDA compatibility.
if(racon_enable_cuda)
    find_package(CUDA 9.0 QUIET REQUIRED)
    if(NOT ${CUDA_FOUND})
        message(FATAL_ERROR "CUDA not detected on system. Please install")
    else()
        message(STATUS "Using CUDA ${CUDA_VERSION} from ${CUDA_TOOLKIT_ROOT_DIR}")
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -lineinfo")
    endif()
endif()

include_directories(${PROJECT_SOURCE_DIR}/src)

set(racon_sources
    src/main.cpp
    src/polisher.cpp
    src/overlap.cpp
    src/sequence.cpp
    src/window.cpp)

if(racon_enable_cuda)
    list(APPEND racon_sources src/cuda/cudapolisher.cpp src/cuda/cudabatch.cpp)
    cuda_add_executable(racon ${racon_sources})
    target_compile_definitions(racon PRIVATE CUDA_ENABLED)
else()
    add_executable(racon ${racon_sources})
endif()

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
if (racon_enable_cuda AND NOT TARGET cudapoa)
    add_subdirectory(nvidia/cudapoa EXCLUDE_FROM_ALL)
endif()

target_link_libraries(racon bioparser spoa thread_pool pthread edlib_static logger)
if (racon_enable_cuda)
    target_link_libraries(racon cudapoa)
endif()

install(TARGETS racon DESTINATION bin)

if (racon_build_tests)
    set(racon_test_data_path ${PROJECT_SOURCE_DIR}/test/data/)
    configure_file("${PROJECT_SOURCE_DIR}/test/racon_test_config.h.in"
        "${PROJECT_BINARY_DIR}/config/racon_test_config.h")
    include_directories(${PROJECT_BINARY_DIR}/config)
    include_directories(${PROJECT_SOURCE_DIR}/src)

    set(racon_test_sources
        test/racon_test.cpp
        src/polisher.cpp
        src/overlap.cpp
        src/sequence.cpp
        src/window.cpp)

    if (racon_enable_cuda)
        list(APPEND racon_test_sources src/cuda/cudapolisher.cpp src/cuda/cudabatch.cpp)
        cuda_add_executable(racon_test ${racon_test_sources})
        target_compile_definitions(racon_test PRIVATE CUDA_ENABLED)
    else()
        add_executable(racon_test ${racon_test_sources})
    endif()

    add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)

    target_link_libraries(racon_test bioparser spoa thread_pool pthread
        edlib_static logger gtest_main)
    if (racon_enable_cuda)
        target_link_libraries(racon_test cudapoa)
    endif()
endif(racon_build_tests)

if (racon_build_wrapper)
    set(racon_path ${PROJECT_BINARY_DIR}/bin/racon)
    set(rampler_path ${PROJECT_BINARY_DIR}/vendor/rampler/bin/rampler)
    configure_file(${PROJECT_SOURCE_DIR}/scripts/racon_wrapper.py
        ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/racon_wrapper)
    file(COPY ${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/racon_wrapper
        DESTINATION ${PROJECT_BINARY_DIR}/bin
        FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)

    add_subdirectory(vendor/rampler)
endif(racon_build_wrapper)
