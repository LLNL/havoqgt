#!/usr/bin/env bash

for s in 17 18 19 20 21 22 23 24 25 26 27 28 29 30; do export DIR=/p/lustre1/iwabuchi/n2v/rmat/s${s}; mkdir -p ${DIR}; srun -t 2:00:00 -N8 --ntasks-per-node=24 ./src/generate_rmat_edge_list -s${s} -o ${DIR}/edge; done