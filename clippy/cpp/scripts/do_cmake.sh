#!/bin/bash

rm CMakeCache.txt

cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-8 -DCMAKE_CXX_FLAGS="-std=c++17 -O3 -g -fopenmp -lrt -lpthread -lstdc++fs"
