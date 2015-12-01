#!/bin/bash


INSTALL_PREFIX=${PWD}


rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-3/root/usr/bin/g++ \
  -DCMAKE_C_COMPILER=/opt/rh/devtoolset-3/root/usr/bin/gcc \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Debug \
  -DMPI_C_COMPILER=/usr/local/tools/mvapich2-gnu-1.9/bin/mpicc \
  -DMPI_CXX_COMPILER=/usr/local/tools/mvapich2-gnu-1.9/bin/mpicxx \
  -DBOOST_ROOT=$HOME/apps/boost_1_59_0/ \
  -DCMAKE_CXX_FLAGS="-std=c++11 -O0 -g3 -lrt " \
  -DCMAKE_C_FLAGS="-std=c++11 -O0 -g3  -lrt " \
