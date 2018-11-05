#!/usr/bin/env bash


NODES=1
PROCS_PER_NODE=24

EDGE_LiST_PREFX="/dev/shm/edge"

SEED=0
GOAL=2
NUM_WK=$((2**20))
DIE_RATE=10
MAX_VISTS=$((2**16))
WARP_TO_SEED=""

while getopts "n:p:e:s:g:w:d:v:r" OPT
do
  case $OPT in
    n) NODES=$OPTARG;;
    p) PROCS_PER_NODE=$OPTARG;;
    e) EDGE_LiST_PREFX=$OPTARG;;
    s) SEED=$OPTARG;;
    g) GOAL=$OPTARG;;
    w) NUM_WK=$OPTARG;;
    d) DIE_RATE=$OPTARG;;
    v) MAX_VISTS=$OPTARG;;
    r) WARP_TO_SEED="-r";;
    :) echo  "[ERROR] Option argument is undefined.";;   #
    \?) echo "[ERROR] Undefined options.";;
  esac
done

#SCALE=20
#rm -f /dev/shm/edge*
#srun -N ${NODES} --ntasks-per-node=${PROCS_PER_NODE} ./src/generate_rmat_edge_list -s ${SCALE}  -o /dev/shm/edge

srun -N ${NODES} --ntasks-per-node=${PROCS_PER_NODE} --distribution=block ./src/ingest_edge_list -o /dev/shm/graph -d $((2**30)) -f 3 ${EDGE_LiST_PREFX}*

srun -N ${NODES} --ntasks-per-node=${PROCS_PER_NODE} --distribution=block --drop-caches=pagecache \
./src/run_rw_v0 -i /dev/shm/graph -s ${SEED} -g ${GOAL} -w ${NUM_WK} -d ${DIE_RATE} -v ${MAX_VISTS} ${WARP_TO_SEED}