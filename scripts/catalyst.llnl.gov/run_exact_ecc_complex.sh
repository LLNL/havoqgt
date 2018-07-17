#!/usr/bin/env bash

# ---------------------------------------------------------------------------------------------------- #
# Usages
# ---------------------------------------------------------------------------------------------------- #
#sh scripts/run_exact_ecc_complex.sh 	  -g 	yt	 -n 	1	 -k 	64	 -t 	ds_fix	 -l 	20:00
#sh scripts/run_exact_ecc_complex.sh 	  -g 	rn	 -n 	8	 -k 	128	 -t 	ds_fix	 -l 	4:00:00
#sh scripts/run_exact_ecc_complex.sh 	  -g 	wk	 -n 	1	 -k 	64	 -t 	ds_fix	 -l 	20:00
#sh scripts/run_exact_ecc_complex.sh 	  -g 	ok	 -n 	8	 -k 	128	 -t 	ds_fix	 -l 	5:00:00
#sh scripts/run_exact_ecc_complex.sh 	  -g 	pt	 -n 	8	 -k 	128	 -t 	ds_fix	 -l 	2:00:00
#sh scripts/run_exact_ecc_complex.sh 	  -g 	wk2	 -n 	128	 -k 	128	 -t 	ds_fix	 -l 	9:00:00
#sh scripts/run_exact_ecc_complex.sh 	  -g 	fr	 -n 	64	 -k 	128	 -t 	ds_fix	 -l 	29:00
#sh scripts/run_exact_ecc_complex.sh 	  -g 	tw	 -n 	64	 -k 	128	 -t 	ds_fix	 -l 	4:00:00
#sh scripts/run_exact_ecc_complex.sh 	  -g 	wb	 -n 	128	 -k 	4	 -t 	tk	 -l 	1:00:00	 -b

# ---------------------------------------------------------------------------------------------------- #
# Command line options
# ---------------------------------------------------------------------------------------------------- #
num_nodes=1
k=128
time_limit="29:00"

while getopts "g:n:k:t:l:s:b" OPT
do
  case $OPT in
    g) graph_name=$OPTARG;;
    t) tag=$OPTARG;;
    n) num_nodes=$OPTARG;;
    k) k=$OPTARG;;
    l) time_limit=$OPTARG;;
    b) use_big_msg_size=true; ;;
    :) echo  "[ERROR] Option argument is undefined.";;   #
    \?) echo "[ERROR] Undefined options.";;
  esac
done

if [ "${graph_name}" == "" ]; then
    echo "Empty graph_name"
    exit
fi

if [ "${tag}" == "" ]; then
    echo "Empty tag"
    exit
fi

# ---------------------------------------------------------------------------------------------------- #
# Configurations
# ---------------------------------------------------------------------------------------------------- #f
storage_path="/p/lscratchh/iwabuchi/"
procs=36

# ---------------------------------------------------------------------------------------------------- #
# Functions
# ---------------------------------------------------------------------------------------------------- #
function make_script
{
    rm -f ${sbatch_out}
    echo "#!/bin/bash" > ${sh_file}
    echo "#SBATCH -N${num_nodes}" >> ${sh_file}
    echo "#SBATCH -o ${sbatch_out}" >> ${sh_file}
    echo "#SBATCH --ntasks-per-node=${procs}" >> ${sh_file}
    echo "#SBATCH -t ${time_limit}" >> ${sh_file}

    echo "${option}" >> ${sh_file}

    echo "srun --ntasks-per-node=${procs} --distribution=block ./src/transfer_graph ${graph_path} /dev/shm/" >> ${sh_file}

    if [ "$use_big_msg_size" ]; then
        echo "export HAVOQGT_MAILBOX_SHM_SIZE=16384" >> ${sh_file}
        echo "export HAVOQGT_MAILBOX_MPI_SIZE=131072" >> ${sh_file}
    fi
    # removed  -c ${two_core_file_name}
    echo "srun --drop-caches=pagecache -N${num_nodes} --ntasks-per-node=${procs} --distribution=block ./src/${exe_file_name} -i /dev/shm/graph  -e ${ecc_file_name} ${source_selection_algorithms}" >> ${sh_file}
}

function compile()
{
        num_sources=$1
        export NUM_SOURCES=${num_sources}
        sh scripts/do_cmake.sh
        make transfer_graph run_exact_ecc #ingest_edge_list
        cp ./src/run_exact_ecc ./src/${exe_file_name}
}

# ---------------------------------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------------------------------- #
base_graph_path="${storage_path}/${graph_name}"
ecc_out_path="${storage_path}/${graph_name}/ecc/"

unset USE_TAKE
unset USE_TAKE_PRUN
unset USE_NEW_ADP
unset USE_NEW_MAX_U
unset USE_SOFT_CONT_SCORE
unset USE_SKIP_STRATEGY

if [ "${tag}" == "tk" ]; then
    option="export USE_TAKE=1" ## 1 is a dummy value
    source_selection_algorithms=""
elif [ "${tag}" == "ds_fix" ]; then
    option="export USE_DS_FIX=1"
    source_selection_algorithms=""
elif [ "${tag}" == "ds_fix_tree" ]; then
    option="export USE_DS_FIX=1; export USE_TREE=1"
    source_selection_algorithms=""
elif [ "${tag}" == "ds_fix_ms" ]; then
    option="export USE_DS_FIX=1; export USE_DS_FIX_MS=1"
    source_selection_algorithms=""
elif [ "${tag}" == "ds_adp_ms" ]; then
    option="export USE_DS_ADP=1"
    source_selection_algorithms="-a 0:1:2:3:4:5:6"
elif [ "${tag}" == "ds_adp_ms_tree" ]; then
    option="export USE_DS_ADP=1; export USE_TREE=1"
    source_selection_algorithms="-a 0:1:2:3:4:5:6"
else
    echo "Invalid tag ${tag}"
    exit
fi

sbatch_out="sbatch_${graph_name}_n${num_nodes}_k${k}_${tag}.out"
sh_file="run_${graph_name}_n${num_nodes}_k${k}_${tag}.sh"
exe_file_name="run_exact_ecc_k${k}_${tag}"
ecc_file_name="${ecc_out_path}/ecc_k${k}_${tag}"
# two_core_file_name="${base_graph_path}/2core_table"
graph_path="${base_graph_path}/n${num_nodes}/graph"

make_script
compile ${k}

sbatch ${sh_file}