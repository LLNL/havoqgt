#!/usr/bin/env bash

TAG=$1

NUM_P=20

unset USE_TAKE

if [${TAG} == "take"]
then
    export USE_TAKE="1" ## 1 is a dummy value
fi

function ingest_edge
{
  rm -rf /dev/shm/g*
  srun -N1 --ntasks-per-node=${NUM_P} ./src/ingest_edge_list -o /dev/shm/graph -f 3 -u 1 ${EDGE} |& tee ingest_${PREFIX}_${TAG}.log
}

function run
{
  if [${TAG} == "rnd"]
  then
    srun --drop-caches=pagecache -N1 --ntasks-per-node=${NUM_P} ./src/run_exact_ecc  -i /dev/shm/graph -e /p/lscratchf/iwabuchi/ecc/ecc_${PREFIX}_${TAG} -a 7 |& tee ecc_${PREFIX}_${TAG}.log
  else
    srun --drop-caches=pagecache -N1 --ntasks-per-node=${NUM_P} ./src/run_exact_ecc  -i /dev/shm/graph -e /p/lscratchf/iwabuchi/ecc/ecc_${PREFIX}_${TAG} -a 0:1:2:3:4:5:6:7:8 |& tee ecc_${PREFIX}_${TAG}.log
  fi
}

#function check
#{
#  ./src/exact_ecc_result_processor ${NUM_P} /p/lscratchf/iwabuchi/ecc/ecc_${PREFIX}_${TAG} ${EDGE} |& tee check_${PREFIX}_${TAG}.log
#}

# /p/lscratchf/iwabuchi/graph_dataset/as-skitter/as-skitter.txt
EDGE=/p/lscratchf/iwabuchi/graph_dataset/as-skitter/as-skitter.txt
PREFIX=sk
ingest_edge
run

# /p/lscratchf/iwabuchi/graph_dataset/flickr/flickr-links.txt
EDGE=/p/lscratchf/iwabuchi/graph_dataset/flickr/flickr-links.txt
PREFIX=fl
ingest_edge
run

# /p/lscratchf/iwabuchi/graph_dataset/wiki-talk/wiki-Talk.txt
EDGE=/p/lscratchf/iwabuchi/graph_dataset/wiki-talk/wiki-Talk.txt
PREFIX=wk
ingest_edge
run

# /p/lscratchf/iwabuchi/graph_dataset/com-youtube/com-youtube.ungraph.txt
EDGE=/p/lscratchf/iwabuchi/graph_dataset/com-youtube/com-youtube.ungraph.txt
PREFIX=yt
ingest_edge
run

# /p/lscratchf/iwabuchi/graph_dataset/roadca/roadNet-CA.txt
EDGE=/p/lscratchf/iwabuchi/graph_dataset/roadca/roadNet-CA.txt
PREFIX=rn
ingest_edge
run

