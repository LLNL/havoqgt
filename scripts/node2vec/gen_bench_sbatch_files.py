import os
import getopt
import sys

# python3 ../../scripts/node2vec_rw/gen_bench_sbatch_files.py -e edge_list -n 1

class Option:
    edge_list_path = ""
    sbatch_file_name = "sbatch.sh"
    out_file_name = "sbatch_out.txt"
    num_nodes = 1
    time_limit = "4:00:00"
    num_tasks_per_node = 24
    graph_path = "/dev/shm/graph"
    scratch_dir = "/dev/shm"
    k_bfs_sources = "4:900:183:296:550:493:67:366"


def parse_options(argv):
    option = Option()

    try:
        opts, args = getopt.getopt(argv, "e:s:o:n:")
    except getopt.GetoptError:
        print("Wrong arguments")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-e':
            option.edge_list_path = arg
        elif opt == '-s':
            option.sbatch_file_name = arg
        elif opt == '-o':
            option.out_file_name = arg
        elif opt == '-n':
            option.num_nodes = int(arg)

    return option


def init(sbatch_file, option):
    sbatch_file.write('#!/bin/bash\n')
    sbatch_file.write('#SBATCH -N %d\n' % option.num_nodes)
    sbatch_file.write('#SBATCH -o %s\n' % option.out_file_name)
    sbatch_file.write('#SBATCH --ntasks-per-node=%d\n' % option.num_tasks_per_node)
    sbatch_file.write('#SBATCH -t %s\n' % option.time_limit)
    # sbatch_file.write('#SBATCH --mail-type=all\n\n')


def add_ingest_edge_command(sbatch_file, option):
    command = 'srun --distribution=block -N%d --ntasks-per-node=%d ./src/ingest_edge_list -o%s -d%d -f%d' \
              % (option.num_nodes, option.num_tasks_per_node, option.graph_path, 2**30, 4)

    if not option.edge_list_path:
        sys.exit("Edge list path is not given")

    if os.path.isfile(option.edge_list_path):
        sbatch_file.write("%s %s\n" % (command, option.edge_list_path))

    # Edge list file
    if os.path.isdir(option.edge_list_path):
        edge_file_list = []
        with os.scandir(option.edge_list_path) as it:
            for entry in it:
                if not entry.name.startswith('.') and entry.is_file():
                    edge_file_list.extend(" %s " % entry.path)

        edge_files = ''.join(edge_file_list)
        sbatch_file.write("%s %s\n" % (command, edge_files))

    sbatch_file.write("\n")


def add_rw_command(sbatch_file, option):
    command = 'srun --distribution=block -N%d --ntasks-per-node=%d ./src/run_node2vec_rw -g%s -d%s -v%s\n' \
              % (option.num_nodes, option.num_tasks_per_node, option.graph_path, option.scratch_dir, option.k_bfs_sources)
    sbatch_file.write(command)
    sbatch_file.write("\n")


def main(argv):
    option = parse_options(argv)

    sbatch_file = open(option.sbatch_file_name, 'w')

    init(sbatch_file, option)
    sbatch_file.write('/bin/rm -rf /dev/shm/*\n\n')
    add_ingest_edge_command(sbatch_file, option)
    add_rw_command(sbatch_file, option)


if __name__ == '__main__':
    main(sys.argv[1:])

