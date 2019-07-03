#!/bin/bash
# Copyright (c) 2018, NVIDIA CORPORATION.
######################################
# cuDF GPU conda build script for CI #
######################################
set -e

# Logger function for build status output
function logger() {
  echo -e "\n>>>> $@\n"
}

cd ${WORKSPACE}

LOCAL_BUILD_DIR=${WORKSPACE}/build

. ci/common/build.sh ${LOCAL_BUILD_DIR}

if [ "${TEST_ON_GPU}" == '1' ]; then
  logger "GPU config..."
  nvidia-smi

  logger "Running GPU-based test..."

  logger "Pulling GPU test data..."
  cd ${WORKSPACE}
  if [ ! -d "ont-racon-data" ]; then
    if [ ! -f "${ont-racon-data.tar.gz}" ]; then
      wget -q -L https://s3.us-east-2.amazonaws.com/racon-data/ont-racon-data.tar.gz
    fi
    tar xvzf ont-racon-data.tar.gz
  fi
  ci/gpu/cuda_test.sh

  logger "Unit test results..."
  cd ${LOCAL_BUILD_DIR}/bin
  ./racon_test --gtest_filter=*CUDA*
fi
