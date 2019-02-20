cmake_minimum_required(VERSION 3.10.2)
project(cudapoa)

find_package(CUDA 9.0 QUIET REQUIRED)
if(NOT ${CUDA_FOUND})
    message(FATAL_ERROR "CUDA not detected on system. Please install")
else()
    message(STATUS "Using CUDA ${CUDA_VERSION} from ${CUDA_TOOLKIT_ROOT_DIR}")
endif()

set(sources
    src/cudapoa_kernels.cu)

cuda_add_library(cudapoa STATIC ${sources})

target_include_directories(cudapoa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add tests folder
add_subdirectory(tests)
