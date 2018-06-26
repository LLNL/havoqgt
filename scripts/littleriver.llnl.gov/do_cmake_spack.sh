#!/bin/bash

spack load cmake
spack load gcc
spack load openmpi
spack load boost

INSTALL_PREFIX=${PWD}

rm CMakeCache.txt
cmake ../../ \
  -DHAVOQGT_BUILD_TEST=FALSE \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DHAVOQGT_BUILD_TEST="ON" \
  -DCMAKE_CXX_FLAGS="-std=c++14 -O3 -lstdc++fs"

