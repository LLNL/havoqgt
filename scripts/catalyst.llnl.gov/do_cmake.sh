#!/bin/bash


INSTALL_PREFIX=${PWD}

if [[ ! -n $NUM_SOURCES ]]
then
    NUM_SOURCES=8
fi

rm CMakeCache.txt

cmake ../../ \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DMPI_C_COMPILER=mpicc \
  -DMPI_CXX_COMPILER=mpicxx \
  -DBOOST_ROOT=$HOME/local/boost \
  -DMPIEXEC_NUMPROC_FLAG="-n" \
  -DCMAKE_CXX_FLAGS="-std=c++11 -lrt -lpthread -lm -DNUM_SOURCES=${NUM_SOURCES}" \