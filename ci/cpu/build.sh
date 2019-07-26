#!/bin/bash
# Copyright (c) 2018, NVIDIA CORPORATION.
######################################
# cuDF CPU conda build script for CI #
######################################
set -e

# Logger function for build status output
function logger() {
  echo -e "\n>>>> $@\n"
}

cd ${WORKSPACE}

LOCAL_BUILD_DIR=${WORKSPACE}/build

. ci/common/build.sh ${LOCAL_BUILD_DIR}

if [ "${TEST_ON_CPU}" == '1' ]; then
  logger "Running CPU-based test..."
  cd ${LOCAL_BUILD_DIR}/bin

  logger "Test results..."
  ./racon_test
fi

