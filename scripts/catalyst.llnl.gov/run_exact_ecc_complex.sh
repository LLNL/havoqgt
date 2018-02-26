#!/usr/bin/env bash

function init
{
    rm -f sbatch_out
    echo "#!/bin/bash" > ${sh_file}
    echo "#SBATCH -N${num_nodes}" >> ${sh_file}
    echo "#SBATCH -o ${sbatch_out}" >> ${sh_file}
    echo "#SBATCH --ntasks-per-node=24" >> ${sh_file}
    echo "#SBATCH -t 7:00:00" >> ${sh_file}

    echo "srun --ntasks-per-node=24 --distribution=block ./src/transfer_graph ${graph_path}${num_nodes}/graph /dev/shm/" >> ${sh_file}
}

function make_script()
{
    num_sources=$1

    echo "export NUM_SOURCES=${num_sources}" >> ${sh_file}
    echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> ${sh_file}
    echo "srun -n1 -N1 make -j4" >> ${sh_file}

    for i in $(seq 1 1)
    do
        echo "srun --drop-caches=pagecache -N${num_nodes} --ntasks-per-node=24 --distribution=block ./src/run_exact_ecc -i /dev/shm/graph  -e ${ecc_out_path}" >> ${sh_file}
    done
}


# ---------------------------------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------------------------------- #
prefix=$1
graph_path=$2
ecc_out_path="/p/lscratchf/iwabuchi/ecc/ecc_${prefix}"

for ((num_nodes=256; num_nodes<=256; num_nodes*=2))
do
    for ((k=64; k<=1024; k*=4))
    do
        sbatch_out="sbatch_${prefix}_n${num_nodes}_k${k}.out"
        sh_file="run_${prefix}_n${num_nodes}_k${k}.sh"

        init
        make_script ${k}

        #sbatch ${sh_file}
    done
done