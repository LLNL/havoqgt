#!/usr/bin/env bash

num_nodes=$1
max_id=$((2**25))

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

    echo "srun --clear-ssd --ntasks-per-node=24 --distribution=block ./src/ingest_edge_list -o /dev/shm/graph -f 3 -u 1 /p/lscratchh/iwabuchi/tw/edgelist_cc_uniq/x*" >> run_${num_nodes}.sh
}

function make_script()
{
    num_sources=$1

    rm -f ./src/run_kbfs_sync
    NUM_SOURCES=${num_sources} sh scripts/do_cmake.sh
    make run_kbfs_sync
    cp ./src/run_kbfs_sync ./src/run_kbfs_sync_${num_nodes}_${num_sources}
    for i in $(seq 1 1)
    do
        src_list=$(./src/gen_src ${num_sources} ${max_id})
        echo "srun --drop-caches=pagecache -N${num_nodes} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync_${num_nodes}_${num_sources} -i /dev/shm/graph  -s ${src_list}" >> run_${num_nodes}.sh
    done
}


# ---------------------------------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------------------------------- #

init

for ((k=1; k<=256; k*=2))
do
    make_script ${k}
done

sbatch run_${num_nodes}.sh

