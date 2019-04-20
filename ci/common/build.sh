#!/bin/bash
# Copyright (c) 2018, NVIDIA CORPORATION.
######################################
# racon-gpu CPU/GPU conda build script for CI #
######################################
set -e

# Logger function for build status output
function logger() {
  echo -e "\n>>>> $@\n"
}

# Set path and build parallel level
export PATH=/conda/bin:/usr/local/cuda/bin:$PATH
export PARALLEL_LEVEL=4

export GW_REPO="ssh://git@gitlab-master.nvidia.com:12051/genomics/GenomeWorks.git"
export GW_NAME=GenomeWorks
export GW_DIR=$WORKSPACE/${GW_NAME}

# Set home to the job's workspace
export HOME=$WORKSPACE

# Switch to project root; also root of repo checkout
cd $WORKSPACE

################################################################################
# SETUP - Check environment
################################################################################

logger "Get env..."
env

logger "Activate conda env..."
source activate gdf

logger "Check versions..."
gcc --version
g++ --version

# FIX Added to deal with Anancoda SSL verification issues during conda builds
conda config --set ssl_verify False

CUDA_REL=${CUDA:0:3}
if [ "${CUDA:0:2}" == '10' ]; then
  # CUDA 10 release
  CUDA_REL=${CUDA:0:4}
fi

git clean -xdf

export LOCAL_BUILD_ROOT=${WORKSPACE}

CMAKE_COMMON_VARIABLES="-DCMAKE_BUILD_TYPE=Release"

if [ "${BUILD_FOR_GPU}" == '1' ]; then
  cd ${WORKSPACE}

  logger "Build racon-gpu for CUDA..."

  if [ ! -d "GenomeWorks" ]; then
    git clone ssh://git@gitlab-master.nvidia.com:12051/genomics/GenomeWorks.git
  fi
  cd GenomeWorks
  logger "Pull racon-gpu..."

  # pull from scratch each time
  cd ${WORKSPACE}
  rm -rf GenomeWorks
  mkdir ${GW_NAME}

  # is this is a merge request and there a branch with a name matching the MR branch in the
  # other repo, pull that
  export BRANCH_FOUND=""
  if [ "${BUILD_CAUSE_SCMTRIGGER}" == "true" ]; then
    logger "This is an SCM-caused build"
    if [ "${gitlabActionType}" == 'MERGE' ]; then
    logger "This is a merge-request-caused build"
      if [ "${gitlabSourceBranch}" != "" ]; then
        logger "The specified branch is: ${gitlabSourceBranch}"
        export BRANCH_FOUND=`git ls-remote -h ${GW_REPO} | grep "refs/heads/${gitlabSourceBranch}$"`
        logger "Branch found test ${BRANCH_FOUND}"
      fi
    fi
  fi

  if [ "${BRANCH_FOUND}" == "" ]; then
    logger "No specified branch - is there a target branch?: ${gitlabTargetBranch}"
    if [ "${gitlabTargetBranch}" != "" ]; then
      logger "A target branch is specified: ${gitlabTargetBranch}"
      export BRANCH_FOUND=`git ls-remote -h ${GW_REPO} | grep "refs/heads/${gitlabTargetBranch}$"`
      logger "Branch found test ${BRANCH_FOUND}"
      if [ "${BRANCH_FOUND}" != "" ]; then
        export MR_BRANCH=${gitlabTargetBranch}
      else
        export MR_BRANCH=master
      fi
    else
      export MR_BRANCH=master
    fi
  else
      export MR_BRANCH=${gitlabSourceBranch}
  fi

  git clone --branch ${MR_BRANCH} --single-branch --depth 1 ${GW_REPO}

  # Switch to project root; also root of repo checkout
  cd ${GW_DIR}

  git pull
  git submodule update --init --recursive


  CMAKE_BUILD_GPU="-Dracon_enable_cuda=ON -DGENOMEWORKS_SRC_PATH=${GW_DIR}"
else
  CMAKE_BUILD_GPU="-Dracon_enable_cuda=OFF -Dracon_build_tests=ON"
fi

export LOCAL_BUILD_ROOT=${WORKSPACE}

cd ${LOCAL_BUILD_ROOT}
export LOCAL_BUILD_DIR=${LOCAL_BUILD_ROOT}/build

# Use CMake-based build procedure
mkdir --parents ${LOCAL_BUILD_DIR}
cd ${LOCAL_BUILD_DIR}

# configure
cmake $CMAKE_COMMON_VARIABLES ${CMAKE_BUILD_GPU} ..
# build
make -j${PARALLEL_LEVEL} VERBOSE=1 all

if [ "${TEST_ON_CPU}" == '1' ]; then
  logger "Running CPU-based test..."
  logger "Test results..."
  cd ${LOCAL_BUILD_DIR}/bin
  ./racon_test
fi

if [ "${TEST_ON_GPU}" == '1' ]; then
  logger "Pulling GPU test data..."
  cd ${WORKSPACE}
  if [ ! -d "ont-racon-data" ]; then
    if [ ! -f "${ont-racon-data.tar.gz}" ]; then
      wget -q -L https://s3.us-east-2.amazonaws.com/racon-data/ont-racon-data.tar.gz
    fi
    tar xvzf ont-racon-data.tar.gz
  fi

  logger "Running GPU-based test..."
  logger "GPU config..."
  nvidia-smi

  logger "Test results..."
  cd ${LOCAL_BUILD_DIR}/bin
  ./cuda_test.sh
fi
