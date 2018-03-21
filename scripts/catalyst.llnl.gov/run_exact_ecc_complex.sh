#!/usr/bin/env bash

# ---------------------------------------------------------------------------------------------------- #
# Command line options
# ---------------------------------------------------------------------------------------------------- #
graph_name=$1
tag=$2

# ---------------------------------------------------------------------------------------------------- #
# Configurations
# ---------------------------------------------------------------------------------------------------- #f
storage_path="/p/lscratchh/iwabuchi/"
procs=36

# ---------------------------------------------------------------------------------------------------- #
# Functions
# ---------------------------------------------------------------------------------------------------- #
function init
{
    rm -f sbatch_out
    echo "#!/bin/bash" > ${sh_file}
    echo "#SBATCH -N${num_nodes}" >> ${sh_file}
    echo "#SBATCH -o ${sbatch_out}" >> ${sh_file}
    echo "#SBATCH --ntasks-per-node=${procs}" >> ${sh_file}
    echo "#SBATCH -t 12:00:00" >> ${sh_file}

    echo "srun --ntasks-per-node=${procs} --distribution=block ./src/transfer_graph ${graph_path} /dev/shm/" >> ${sh_file}
}

function make_script()
{
    num_sources=$1

#    echo "export NUM_SOURCES=${num_sources}" >> ${sh_file}
#    echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> ${sh_file}
#    echo "srun -n1 -N1 make -j4" >> ${sh_file}

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
ecc_out_path="${storage_path}/ecc/${graph_name}/"

unset USE_TAKE
unset USE_TAKE_PRUNING

if [${tag} == "tk"]
then
    export USE_TAKE="1" ## 1 is a dummy value
    source_selection_algorithms=""
elif [${tag} == "tk_pr"]
then
    export USE_TAKE_PRUNING="1" ## 1 is a dummy value
    source_selection_algorithms="-a 0:1:2:3:4:5:6:7"
elif [${tag} == "adp"]
then
    source_selection_algorithms="-a 0:1:2:3:4:5:6:7"
elif [${tag} == "adp_at"]
then
    source_selection_algorithms="-a 0:1:2:3:4:5:6:7:8"
fi

for ((num_nodes=128; num_nodes<=128; num_nodes*=2))
do
    for ((k=64; k<=64; k*=2))
    do
        sbatch_out="sbatch_${graph_name}_n${num_nodes}_k${k}_${tag}.out"
        sh_file="run_${graph_name}_n${num_nodes}_k${k}_${tag}.sh"
        exe_file_name="run_exact_ecc_k${num_sources}_${tag}"
        ecc_file_name="${ecc_out_path}/ecc_k${num_sources}_${tag}"
        graph_path="${base_graph_path}/n${num_nodes}/graph"

        init
        make_script ${k}
        compile ${k}

        #sbatch ${sh_file}
    done
done