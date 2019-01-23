#!/bin/bash

rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=FALSE \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Debug \
  -DMPI_C_COMPILER=mpicc \
  -DMPI_CXX_COMPILER=mpicxx \
  -DBOOST_ROOT=$HOME/local/boost \
  -DCMAKE_CXX_FLAGS="-std=c++17 -Wredundant-decls -DDEBUG -DDEBUG_DPG -lrt -lpthread -lm"

