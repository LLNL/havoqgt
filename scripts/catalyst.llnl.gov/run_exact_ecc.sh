#!/usr/bin/env bash

num_nodes=32
max_id=$((2**26))
graph_path="/p/lscfratchf/iwabuchi/havoqgt_graph/tw/n"

function init
{
    rm -f sbatch${num_nodes}.out
    echo "#!/bin/bash" > run_${num_nodes}.sh
    echo "#SBATCH -N${num_nodes}" >> run_${num_nodes}.sh
    echo "#SBATCH -o sbatch${num_nodes}.out" >> run_${num_nodes}.sh
    echo "#SBATCH --ntasks-per-node=24" >> run_${num_nodes}.sh
    echo "#SBATCH -t 4:00:00" >> run_${num_nodes}.sh

    echo "export HAVOQGT_MAILBOX_SHM_SIZE=16384" >> run_${num_nodes}.sh
    echo "export HAVOQGT_MAILBOX_MPI_SIZE=131072" >> run_${num_nodes}.sh

    echo "srun --ntasks-per-node=24 --distribution=block ./src/transfer_graph ${graph_path}${num_nodes}/graph /dev/shm/graph" >> run_${num_nodes}.sh
}

function make_script()
{
    num_sources=$1

    echo "export NUM_SOURCES=${num_sources}" >> run_${num_nodes}.sh
    echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> run_${num_nodes}.sh
    echo "srun -n1 -N1 make -j4" >> run_${num_nodes}.sh
    for i in $(seq 1 1)
    do
        echo "src_list=\$(./src/gen_src ${num_sources} ${max_id})" >> run_${num_nodes}.sh
        echo "srun --drop-caches=pagecache -N${num_nodes} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s \${src_list}" >> run_${num_nodes}.sh
    done
}


# ---------------------------------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------------------------------- #

init

for ((k=1; k<=256; k*=4))
do
    make_script ${k}
done

#sbatch run_${num_nodes}.sh