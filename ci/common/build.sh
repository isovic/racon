#!/bin/bash
# Copyright (c) 2018, NVIDIA CORPORATION.
######################################
# racon-gpu CPU/GPU conda build script for CI #
######################################
set -e

# Get commandline arguments
LOCAL_BUILD_DIR=$1

# Logger function for build status output
function logger() {
  echo -e "\n>>>> $@\n"
}

# Set path and build parallel level
export PATH=/conda/bin:/usr/local/cuda/bin:$PATH
export PARALLEL_LEVEL=4

# Set home to the job's workspace
export HOME=$WORKSPACE

# Switch to project root; also root of repo checkout
cd $WORKSPACE

################################################################################
# SETUP - Check environment
################################################################################

logger "Get env..."
env

logger "Check versions..."
gcc --version
g++ --version

# FIX Added to deal with Anancoda SSL verification issues during conda builds
conda config --set ssl_verify False

conda install \
    -c conda-forge \
    -c sarcasm \
    -c bioconda \
    doxygen \
    ninja \
    cmake

CUDA_REL=${CUDA:0:3}
if [ "${CUDA:0:2}" == '10' ]; then
  # CUDA 10 release
  CUDA_REL=${CUDA:0:4}
fi

git clean -xdf

CMAKE_COMMON_VARIABLES="-DCMAKE_BUILD_TYPE=Release -Dracon_build_tests=ON"

if [ "${BUILD_FOR_GPU}" == '1' ]; then
  CMAKE_BUILD_GPU="-Dracon_enable_cuda=ON"
else
  CMAKE_BUILD_GPU="-Dracon_enable_cuda=OFF"
fi

# Use CMake-based build procedure
mkdir --parents ${LOCAL_BUILD_DIR}
cd ${LOCAL_BUILD_DIR}

# configure
cmake $CMAKE_COMMON_VARIABLES ${CMAKE_BUILD_GPU} ..
# build
make -j${PARALLEL_LEVEL} VERBOSE=1 all
