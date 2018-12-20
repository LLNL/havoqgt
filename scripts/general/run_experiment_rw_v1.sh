#!/usr/bin/env bash


function read_src_community () {
    file_name=$1
    i=0
    while read line; do
        set -f
        set -- $line
        src_list[$i]=$1
        community_list[$i]=$2
        i=$((i+1))
    done < ${file_name}
}


# --- file name --- #
true_value_file="../../../tmp/true.txt"
graph_file="/dev/shm/graph"
# ------

# --- input --- #
# src_list=(1 2 3 4 5 67 8 9 10)
# community_list=(8 52 10 106 118 219 191 132 165 106)
# or
read_src_community ${true_value_file}
# ------

# --- main loop --- #
LOG_d=jordan/recall_dead.log
LOG_v=jordan/recall_visit.log

echo -e src '\t' cs '\t' ct '\t' dr '\t' rw '\t' recall > ${LOG_d}
echo -e src '\t' cs '\t' ct '\t' dr '\t' rw '\t' recall > ${LOG_v}

for (( s=0; s<${#src_list[@]}; s++ ))
do
    src=${src_list[$s]}
    cm=${community_list[$s]}

    num_v_in_cm=$(cat ${true_value_file} | awk -v c="$cm" '{if ($2==c) count+=1} END{print count}')

    echo $src " in " $cm " : " $num_v_in_cm


    for cs in $(seq 0 1 0); do
        for ct in $(seq 0 100 100); do
            for dr in $(seq 15 15 15); do
                rw=1000
                for i in {1..2}; do
                    rm -f /dev/shm/out_*
                    #LOG2=jordan/detail/rw_cs${cs}_ct${ct}_dr${dr}_rw${rw}_src${src}.log
                    sh ../../scripts/general/run_rw_v1.sh -p24 -g ${graph_file} -s ${src} -w ${rw} -d ${dr} -l 9999999 -c ${cs}:${ct}:0 -o /dev/shm/out -t ${num_v_in_cm} #| tee ${LOG2}

                    ret=$(python3 ../../scripts/general/calculate_precision.py ${cm} ${true_value_file} /dev/shm/out_dead_score)
                    echo -e ${src} '\t' ${cs} '\t' ${ct} '\t' ${dr} '\t' ${rw} '\t' ${ret} >> ${LOG_d}

                    ret=$(python3 ../../scripts/general/calculate_precision.py ${cm} ${true_value_file} /dev/shm/out_visit_score)
                    echo -e ${src} '\t' ${cs} '\t' ${ct} '\t' ${dr} '\t' ${rw} '\t' ${ret} >> ${LOG_v}

                    rw=$((rw * 10))
                    # cp /dev/shm/out_dead_score jordan/scores/dead_score_src${src}
                    # cp /dev/shm/out_visit_score jordan/scores/dead_score_src${src}
                done
            done
        done
    done
done