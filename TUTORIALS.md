This is a tutorial to run a dynamic graph construction benchmark.

Program: run_dgs_bench

# Introduction
run_dgs_bench is a benchmark program to evaluate our graph store implementations, including two baseline implemenations for performance comparison, on distributed memory platform.
run_dgs_bench insert edges uniquely with generating edges or loading edges from files.
The total graph construction time will be printed out after finishing a graph construction with a prefix "only graph construction time (sec.)".

# Build
```
$ git clone git@github.com:LLNL/havoqgt.git
$ cd havoqgt
$ git checkout develop_keita
$ cd build/default/
# edit script/do_cmake.sh (cmake's configuretion) depends on your environment #
$ sh script/do_cmake.sh
$ make
```
To compile havoqgt, Boost C++ libraries are required. Please download from http://www.boost.org/ (at least version 1.59 is required).
You don't need to compile it; just specify its root directory in the script file (do_cmake.sh).

GCC 4.8 or more is also required.

# Options
+ `-o<string>`
	the base file name each process creates to store a graph
+ `-S<int>`
	the maximum size of the file each process creates in GB
+ `-g<string>`
	the name of graph store implementation (DegAwareRHH, Baseline or BaselieMap)
+ `-s<int>`
	the size of the RMAT graph to be constructed in SCALE (performs a test with generationg a RMAT graph)
+ `-E<string>`
	the name of a file that contains file names of edges (performs a test with loading a graph from the files)
+ `-d`
  load also delete flags when load edges from files

# Run Benchmark
## Example1 (with RMAT graph)
```
$ mpirun -np 2 ./run_dgs_bench -s20 -o/mnt/device/graph_out -S2 -gDegAwareRHH
```

* Run the dynamic graph construction program with 2 processes on a RMAT SCALE 20 graph.
* Each process creates 2 GB of a file to store the graph under '/mnt/device/' with the prefix 'graph_out'.
* Use DegAwareRHH

Note:

To create a SCALE 30 graph, 250 ~ 260 GB of space is required for DegAwareRHH and Baseline. BaselineMap requires 500 GB.
FYI, DegAwareRHH took 14 hours to create SCALE 30 graph with 40 processes with 128 GB of DRAM and enough size of NVRAM.

## Example 2 (loading edges from files)
```
$ mpirun -np 2 ./run_dgs_bench -E./edge_files -o/mnt/device/graph_out -S2 -gDegAwareRHH
```
* Run the dynamic graph construction program with 2 processes with loading a graph from files
```
$ cat ./edge_files
path/to/edge_file1
path/to/edge_file2

$ cat path/to/edge_file1
1 1
4 2
2 1
```
Each line of the file represents a pair of source vertex and target vertex.

Note:

The number of edge list files can be any number regardless of the number of MPI processes.
To achieve fair load balancing among MPI processes, it is recommended that the number of edge list files is multiples of the number of MPI processes.

## Example 3 (loading edges from files with delete requests)
```
$ mpirun -np 2 ./run_dgs_bench -E./edge_files -o/mnt/device/graph_out -S2 -gDegAwareRHH -d
```
* Run the dynamic graph construction program with 2 processes with loading a graph from files with delete flags
```
$ cat ./edge_files
path/to/edge_file1
path/to/edge_file2

$ cat path/to/edge_file1
1 1 1
4 2 0
2 1 0
```
Each line of the file represents a pair of source vertex and target vertex and a delete flag (0: add a edge, 1: delete a edge)

Note:

Please don't forget to specify -d option


# Reading our wikipedia graph (will be open soruced soon)
Please apply a patch file under the root directory, i.e., havoqgt, to read our wikipedia graph.
## apply patch
```
$ cd havoqgt
$ patch -p0 < ./scripts/default/edgelist_reader_wikipedia_support.patch
```
## undo patch
```
$ cd havoqgt
$ patch -R -p0 < ./scripts/default/edgelist_reader_wikipedia_support.patch
```
After applying the patch and building the program, execution commands are same as the example 2 (loading a graph from file).

