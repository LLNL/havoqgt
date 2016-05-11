This is a instruction to run a dynamic graph construction test.
Program: dist_bench

### Options:
+ `-s<int>`
	the size of the RMAT graph to be constructed in SCALE
+ `-o<string>`
	the base file name each process creates to store the graph
+ `-S<int>`
	the size of the file each process create in the logarithm base two
+ `-g<string>`
	the name of  graph store (DegAwareRHH, Baseline or BaselieMap)

### Example
``
$ mpirun -np 2 ./dist_bench -s20 -o/mnt/device/graph_out -S30 -gDegAwareRHH
``
1. Run the dynamic graph construction program with 2 MPI processes on a RMAT SCALE 20 graph.
2. Each process create 2^30 byte of a file to store the graph under /mnt/device directory with the suffix "graph_".
3. Use gDegAwareRHH.

To create a SCALE 30 graph, 250 ~ 260 GB of space is required for DegAwareRHH and Baseline. BaselineMap requires 500 GB.
DegAwareRHH takes 14 hours to create SCALE 30 graph with 40 processes with 128 GB of DRAM.


## Reading Suraj's wikipedia graph
To read Suraj's wikipedia graph dataset, applying a patch file is required
### apply patch
```
cd havoqgt
patch -p0 < scripts/catalyst.llnl.gov/edgelist_reader_wikipedia_support.patch
```
### undo patch
```
patch -R -p0 < scripts/catalyst.llnl.gov/edgelist_reader_wikipedia_support.patch
```
