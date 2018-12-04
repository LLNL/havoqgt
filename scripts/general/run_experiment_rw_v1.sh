#!/usr/bin/env bash

src_list=(1 2 3 4 5 67 8 9 10)
community_list=(8 52 10 106 118 219 191 132 165 106)
ref_file="/p/lustre1/iwabuchi/gc/simulated_blockmodel_graph_5000000_nodes_truePartition.tsv"

for (( s=0; s<${#src_list[@]}; s++ ))
do
    src=${src_list[$s]}
    cm=${community_list[$s]}

    num_v_in_cm=$(cat ${ref_file} | awk -v c="$cm" '{if ($2==c) count+=1} END{print count}')

    LOG=precision_src${src}.log
    echo -e cs '\t' ct '\t' dr '\t' rw '\t' pre > ${LOG}

    for cs in $(seq 0 20 100); do
        for ct in $(seq 0 20 100); do
            for dr in $(seq 15 15 90); do
                rw=10
                for i in {1..6}; do
                    rm -f /dev/shm/out_*
                    sh ../../scripts/general/run_rw_v1.sh -p24 -e /p/lustre1/iwabuchi/gc/edges -s ${src} -w ${rw} -d ${dr} -l 9999999 -c ${cs}:${ct}:0 -o /dev/shm/out -t ${num_v_in_cm}
                    ret=$(python3 ../../scripts/general/calculate_precision.py ${cm} ${ref_file} /dev/shm/out_dead_score)
                    echo -e ${cs} '\t' ${ct} '\t' ${dr} '\t' ${rw} '\t' ${ret} >> ${LOG}
                    rw=$((rw * 10))
                done
            done
        done
    done
done