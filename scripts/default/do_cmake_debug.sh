#!/bin/bash


INSTALL_PREFIX=${PWD}


rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Debug \
  -DMPI_C_COMPILER=mpicc \
  -DMPI_CXX_COMPILER=mpicxx \
  -DBOOST_ROOT=$HOME/boost_1_59_0/ \
  -DCMAKE_CXX_FLAGS="-std=c++11 -O0 -g3 -lrt " \
  -DCMAKE_C_FLAGS="-std=c++11 -O0 -g3  -lrt " \
