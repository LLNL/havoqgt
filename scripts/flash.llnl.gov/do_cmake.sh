#!/bin/bash


INSTALL_PREFIX=${PWD}


rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI_C_COMPILER=/usr/tce/packages/openmpi/openmpi-2.0.0-gcc-4.9.3/bin/mpicc \
  -DMPI_CXX_COMPILER=/usr/tce/packages/openmpi/openmpi-2.0.0-gcc-4.9.3/bin/mpicxx \
  -DBOOST_ROOT=/usr/tce/packages/boost/boost-1.58.0-mvapich2-2.2-intel-16.0.3 \
  -DCMAKE_CXX_FLAGS="-std=c++11 -lrt -lpthread -lm" \
