#!/bin/bash


INSTALL_PREFIX=${PWD}


rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_COMPILER=g++ \
  -DMPI_C_COMPILER=mpicc \
  -DMPI_CXX_COMPILER=mpicxx \
  -DMPIEXEC_NUMPROC_FLAG="-n" \
  -DCMAKE_CXX_FLAGS="-std=c++17 -Wredundant-decls -DDEBUG -DDEBUG_DPG -lrt -lpthread -lm" \
