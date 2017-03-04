#!/bin/bash


INSTALL_PREFIX=${PWD}


rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI_C_COMPILER=/usr/tce/packages/openmpi/openmpi-2.0.0-gcc-4.9.3/bin/mpicc \
  -DMPI_CXX_COMPILER=/usr/tce/packages/openmpi/openmpi-2.0.0-gcc-4.9.3/bin/mpicxx \
  -DBOOST_ROOT=/g/g90/reza2/boost_1_57_0/ \
  -DCMAKE_CXX_FLAGS="-std=c++11 -lrt -lpthread -lm -lboost_system -lboost_filesystem -lboost_regex -L boost/lib -I /boost/include" \
