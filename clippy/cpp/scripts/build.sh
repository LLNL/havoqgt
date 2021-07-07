#!/bin/bash
cwd=$(pwd)
mkdir ~/clippy/cpp/build
mkdir ~/clippy/cpp/build/examples
cp ~/clippy/cpp/scripts/do_cmake.sh ~/clippy/cpp/build
cd ~/clippy/cpp/build
sh ~/clippy/cpp/build/do_cmake.sh
make st_path
make st_shortest_path
make sum
cd $cwd
