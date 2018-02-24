#!/usr/bin/env bash

num_nodes=32
max_id=$((2**24))
graph_path="/p/lscratchf/iwabuchi/havoqgt_graph/tw/n"
ecc_out_path="/p/lscratchf/iwabuchi/havoqgt_graph/tw/n${num_nodes}/ecc"

sbatch_out="sbatch${num_nodes}.out"
sh_file="run_${num_nodes}.sh"

function init
{
    rm -f sbatch_out
    echo "#!/bin/bash" > ${sh_file}
    echo "#SBATCH -N${num_nodes}" >> ${sh_file}
    echo "#SBATCH -o ${sbatch_out}" >> ${sh_file}
    echo "#SBATCH --ntasks-per-node=24" >> ${sh_file}
    echo "#SBATCH -t 4:00:00" >> ${sh_file}

    echo "export HAVOQGT_MAILBOX_SHM_SIZE=16384" >> ${sh_file}
    echo "export HAVOQGT_MAILBOX_MPI_SIZE=131072" >> ${sh_file}

    echo "srun --ntasks-per-node=24 --distribution=block ./src/transfer_graph ${graph_path}${num_nodes}/graph /dev/shm/" >> ${sh_file}
}

function make_script()
{
    num_sources=$1

    echo "export NUM_SOURCES=${num_sources}" >> ${sh_file}
    echo "srun -n1 -N1 sh scripts/do_cmake.sh" >> ${sh_file}
    echo "srun -n1 -N1 make -j4" >> ${sh_file}
    #echo "srun -n1 -N1 g++ -std=c++11 ${HOME}/utility/src/gen_src.cpp -o ./src/gen_src" >> ${sh_file}
    for i in $(seq 1 1)
    do
        #echo "src_list=\$(./src/gen_src ${num_sources} ${max_id})" >> ${sh_file}
        echo "srun --drop-caches=pagecache -N${num_nodes} --ntasks-per-node=24 --distribution=block ./src/run_exact_ecc -i /dev/shm/graph  -e ${ecc_out_path}" >> ${sh_file}
    done
}


# ---------------------------------------------------------------------------------------------------- #
# Main
# ---------------------------------------------------------------------------------------------------- #

init

for ((k=64; k<=64; k*=4))
do
    make_script ${k}
done

#sbatch ${sh_file}