#!/bin/bash

MACPORTS_PATH=/opt/local

INSTALL_PREFIX=${PWD}

rm CMakeCache.txt
cmake ../../ \
  -DHAVOQGT_BUILD_TEST=FALSE \
  -DCMAKE_CXX_COMPILER=${MACPORTS_PATH}/bin/g++-mp-4.7 \
  -DCMAKE_C_COMPILER=${MACPORTS_PATH}/bin/gcc-mp-4.7 \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI_C_COMPILER=${MACPORTS_PATH}/bin/mpicc-mpich-gcc47 \
  -DMPI_CXX_COMPILER=${MACPORTS_PATH}/bin/mpicxx-mpich-gcc47 \
  -DBOOST_ROOT=${MACPORTS_PATH}/include \
  -DHAVOQGT_BUILD_TEST="ON" \
  -DCMAKE_CXX_FLAGS="-std=c++11 -O3" \
  -DMPIEXEC_PREFLAGS="-hosts;localhost"

