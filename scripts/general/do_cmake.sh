#!/bin/bash

rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=FALSE \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI_C_COMPILER=mpicc \
  -DMPI_CXX_COMPILER=mpicxx \
  -DBOOST_ROOT=$HOME/local/boost \
  -DMPIEXEC_NUMPROC_FLAG="-n" \
  -DCMAKE_CXX_FLAGS="-std=c++11 -lrt -lpthread -lm"