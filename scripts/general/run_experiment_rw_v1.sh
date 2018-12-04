#!/usr/bin/env bash

for src in {1..2}; do
    LOG=precision_src${src}.log
    echo -e cs '\t' ct '\t' dr '\t' rw '\t' pre > ${LOG}
    for cs in $(seq 0 20 100); do
        for ct in $(seq 0 20 100); do
            for dr in $(seq 15 15 90); do
                rw=10
                for i in {1..6}; do
                    rm -f /dev/shm/out_*
                    sh ../../scripts/general/run_rw_v1.sh -p24 -e /p/lustre1/iwabuchi/gc/edges -s 654 -w ${rw} -d ${dr} -l 9999999 -c ${cs}:${ct}:0 -o /dev/shm/out_
                    ret=$(python3 ../../scripts/general/calculate_precision.py 654 /p/lustre1/iwabuchi/gc/simulated_blockmodel_graph_5000000_nodes_truePartition.tsv /dev/shm/out__dead_score_${rw}_*)
                    echo -e ${cs} '\t' ${ct} '\t' ${dr} '\t' ${rw} '\t' ${ret} >> ${LOG}
                    rw=$((rw * 10))
                done
            done
        done
    done
done