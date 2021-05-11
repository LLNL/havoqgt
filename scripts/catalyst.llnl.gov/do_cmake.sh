#!/bin/bash


INSTALL_PREFIX=${PWD}


rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/usr/tce/packages/gcc/gcc-8.3.1/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/gcc/gcc-8.3.1/bin/g++ \
  -DMPI_C_COMPILER=/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-8.3.1/bin/mpicc \
  -DMPI_CXX_COMPILER=/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-8.3.1/bin/mpicxx \
  -DBOOST_ROOT=/usr/tce/packages/boost/boost-1.66.0-mvapich2-2.2-gcc-4.9.3 \
  -DMPIEXEC_NUMPROC_FLAG="-n" \
  -DCMAKE_CXX_FLAGS="-std=c++17 -lrt -lstdc++fs -lpthread" \
