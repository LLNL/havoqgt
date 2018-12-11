#!/usr/bin/env bash


NODES=1
PROCS_PER_NODE=24

EDGE_LIST_PREFX="/dev/shm/edge"
GRAPH="/dev/shm/graph"

PERSONALIZED=""
SEED=0
NUM_WK=$((2**20))
DIE_RATE=10
MAX_WK_LENGTH=$((2**16))
WARP_TO_SEED=""
CLOSING_RATE=""
SCORE_DUMP_FILE=""
NUM_TOP=10

while getopts "n:p:e:g:s:w:d:l:rc:o:t:" OPT
do
  case $OPT in
    n) NODES=$OPTARG;;
    p) PROCS_PER_NODE=$OPTARG;;
    e) EDGE_LIST_PREFX=$OPTARG;;
    g) GRAPH=$OPTARG;;
    s) SEED=$OPTARG PERSONALIZED="-p";;
    w) NUM_WK=$OPTARG;;
    d) DIE_RATE=$OPTARG;;
    l) MAX_WK_LENGTH=$OPTARG;;
    r) WARP_TO_SEED="-r";;
    c) CLOSING_RATE="-c "$OPTARG;;
    o) SCORE_DUMP_FILE="-o "$OPTARG;;
    t) NUM_TOP=$OPTARG;;
    :) echo  "[ERROR] Option argument is undefined.";;   #
    \?) echo "[ERROR] Undefined options.";;
  esac
done

#SCALE=20
#srun -N ${NODES} --ntasks-per-node=${PROCS_PER_NODE} ./src/generate_rmat_edge_list -s ${SCALE}  -o ${EDGE_LIST_PREFX}

#srun -N ${NODES} --ntasks-per-node=${PROCS_PER_NODE} --distribution=block ./src/ingest_edge_list -o ${GRAPH} -d $((2**30)) -f 3 ${EDGE_LIST_PREFX}*

srun -N ${NODES} --ntasks-per-node=${PROCS_PER_NODE} --distribution=block --drop-caches=pagecache \
./src/run_rw_v1 -i ${GRAPH} ${PERSONALIZED} -s ${SEED} -w ${NUM_WK} -d ${DIE_RATE} -l ${MAX_WK_LENGTH} ${WARP_TO_SEED} ${CLOSING_RATE} ${SCORE_DUMP_FILE} -t ${NUM_TOP}