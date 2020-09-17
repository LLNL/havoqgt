# CORAL2 Triangle Counting Benchmark

A triangle counting benchmark has been built using HavoqGT for CORAL2 evaluation.  
This benchmark evaluates the performance of computing triangles on large scale-free graphs
using the algorithm presented [here](https://doi.org/10.1109/HPEC.2017.8091051).

## Inputs
The input is designed to be similar to the Graph500, but is constructed by the kronecker product
of two previously generated graphs whose triangles have been previously counted.  After kroneckering
the two inputs, the triangle count of the input graph can be validated easily.

Inputs for the benchmark are in the /data/kron_inputs folder and are labeled as 'A' or 'B'.  Kroneckering two Scale 10 (i.e., A_G500_S10_tc.edges  B_G500_S10_tc.edges) will yield a Graph500-like Scale 20 graph.  Scales 20-36 can be created using the included data.  

## Example

Below is example of running the benchmark with a Graph500-like Scale20 graph.  Two Scale 10 graphs are kroneckered to form the input Scale 20.

```console
srun -n 16 src/CORAL2_benchmark ../../data/kron_inputs/A_G500_S10_tc.edges ../../data/kron_inputs/B_G500_S10_tc.edges
```

Example output:

```console
<snip graph construction details>

Time to construct DODGraph (seconds) = 1.70931
Largest original degree = 49504
Largest DODGraph out degree = 538
Count of nonzeros in original graph = 33547264
Count of directed edges in DODGraph = 16773624
Triangle Count = 139903248
Number of Wedge Checks = 1115054571
Total Triangle Count Time (seconds) = 18.4261
!!PASSED!!
```