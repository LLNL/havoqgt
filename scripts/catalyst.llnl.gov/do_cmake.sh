#!/bin/bash


INSTALL_PREFIX=${PWD}


rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI_C_COMPILER=/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3/bin/mpicc \
  -DMPI_CXX_COMPILER=/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3/bin/mpicxx \
  -DBOOST_ROOT=/usr/tce/packages/boost/boost-1.58.0-mvapich2-2.2-intel-16.0.3 \
  -DMPIEXEC_NUMPROC_FLAG="-n" \
  -DCMAKE_CXX_FLAGS="-std=c++11 -lrt -lpthread -lm" \

