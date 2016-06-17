This is a instruction to run a dynamic graph construction test.

# Build
```
$ git clone git@bitbucket.org:PerMA/havoqgt.git
$ cd havoqgt
$ git checkout develop_keita
$ cd build/catalyst.llnl.gov/
# edit script/do_cmake.sh (cmake's configuretion) depends on your environment #
$ sh script/do_cmake.sh
$ make
```
To compile havoqgt, Boost libraries are required. Please download from http://www.boost.org/ (at least version 1.59 is required).
You don't need to compile it; just specify its root directory in the script file (do_cmake.sh).

GCC 4.8 or more is also required.

# Options
+ `-o<string>`
	the base file name each process creates to store the graph
+ `-S<int>`
	the size of the file each process create in the logarithm base two
+ `-g<string>`
	the name of  graph store (DegAwareRHH, Baseline or BaselieMap)
+ `-s<int>`
	the size of the RMAT graph to be constructed in SCALE (performs a test with generationg a RMAT graph)
+ `-E<string>`
	the name of a file that contains file names of edge lists (performs a test with reading a graph from files)

# Run
## Example1 (RMAT graph)
```
$ mpirun -np 2 ./dist_bench -s20 -o/mnt/device/graph_out -S30 -gDegAwareRHH
```

* Run the dynamic graph construction program with 2 processes on a RMAT SCALE 20 graph.
* Each process create 2^30 byte of a file to store the graph under '/mnt/device/' with the prefix 'graph_out'.
* Use DegAwareRHH

Note:

To create a SCALE 30 graph, 250 ~ 260 GB of space is required for DegAwareRHH and Baseline. BaselineMap requires 500 GB.
DegAwareRHH takes 14 hours to create SCALE 30 graph with 40 processes with 128 GB of DRAM.

## Example 2 (reading a graph from file)
```
$ mpirun -np 2 ./dist_bench -E./edgelists -o/mnt/device/graph_out -S30 -gDegAwareRHH
```
* Run the dynamic graph construction program with 2 processes with reading a graph from files
```
$ cat ./edgelists
path/to/edgelist1
path/to/edgelist2

$ cat path/to/edgelist1
1 1
4 2
2 1
```
Each line of the file represents a pair of source vertex and target vertex.

Note:

The number of edge list files can be any number regardless of the number of MPI processes.
To achieve fair load balancing among MPI processes, it is recommended that the number of edge list files is multiples of the number of MPI processes.

# Reading Suraj's wikipedia graph
Please apply a patch file under the root directory, i.e., havoqgt, to read the Suraj's wikipedia graph.
## apply patch
```
$ cd havoqgt
$ patch -p0 < ./scripts/catalyst.llnl.gov/edgelist_reader_wikipedia_support.patch
```
## undo patch
```
$ cd havoqgt
$ patch -R -p0 < ./scripts/catalyst.llnl.gov/edgelist_reader_wikipedia_support.patch
```
After applying the patch and building the program, execution commands are same as the example 2 (reading a graph from file).

