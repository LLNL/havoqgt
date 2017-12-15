#!/usr/bin/env bash

num_nodes=32
max_id=$((2**26))

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

    echo "srun --clear-ssd --ntasks-per-node=24 --distribution=block ./src/ingest_edge_list -o /dev/shm/graph -f 3 /p/lscratchf/havoqgtu/real/twitter/32x24parts/x0*" >> run_${num_nodes}.sh
}

function make_script()
{
    num_sources=$1

    echo "export NUM_SOURCES=${num_sources}" >> run_${num_nodes}.sh
    echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> run_${num_nodes}.sh
    echo "srun -n1 -N1 make -j4" >> run_${num_nodes}.sh
    for i in $(seq 1 10)
    do
        src_list=$(./src/gen_src ${num_sources} ${max_id})
        echo "srun --drop-caches=pagecache -N${num_nodes} --ntasks-per-node=24 --distribution=block ./src/run_kbfs_sync -i /dev/shm/graph  -s ${src_list}" >> run_${num_nodes}.sh
    done
}


# ---------------------------------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------------------------------- #

init

for ((k=1; k<=512; k*=2))
do
    make_script ${k}
done

sbatch run_${num_nodes}.sh

