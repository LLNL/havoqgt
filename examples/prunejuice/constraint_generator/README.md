
# Pattern Matching - Constraint Generator

This project generate the set of constraint used for graph pruning by [this asynchronous pattern matching code](https://github.com/hpcresearchanddevelopment/patternmatching) built on [HavoqGT](https://github.com/LLNL/havoqgt). 

## Usage
This script requires python 3.2 or a newer version and the module [NetworkX](https://networkx.github.io/) for graph processing.

The pattern are described through two different files :
 * An edge list file which contain all of the edges of the graph pattern. It can be a directed. An edge is described by two vertex identifier (int) : `5 6` is an edge from the vertex `5` to the vertex `6`.
``` 
# Example : data/test1/patterns/pattern_edge
# vertexFrom vertexTo
0 1
0 2
1 0
1 2
2 0
2 1
```
  * A label list which associate to each vertex of the graph a specific label.
```
# Example : data/test1/patterns/pattern_vertex_data
# vertex label
0 0
1 1
2 2
```

### Launch via the command line
The script can be launched through the python script `src/main.py`. It takes several input arguments :
 * `-ie` `--input_pattern_edges` The pattern edges file.
 * `-id` `--input_pattern_data` The pattern vertex label file.
 * `-o` `--output_directory` The output directory for the constraint generation

Example :
```shell
python src/main.py --input_pattern_edges data/test1/patterns/pattern_edge --input_pattern_data data/test1/patterns/pattern_vertex_data --output_directory data/test1/results/
```

### Launch via the test script
A script is available to test constraint generation on two pattern provided in the `data/` directory.
```shell
./run_test.sh
```

## Generated constraint 

The constraint are generated in different files with a csv format inside the _output directory_. 

### Local constraint checking
The file `local_constraint.txt` contains the output of the local constraint checking algorithm.  For each vertex of the pattern, this algorithms enumerates the labels of their neighbors.

### Leaf elimination
The leafs of the pattern whose label are unique are pruned from the pattern graph. The Local constraint checking steps will be enough to ensure that they exist.

### Cycle constraint
The file `cycle_constraint.txt` contains the output of the cycle generation algorithm.

A list of all of the cycles in the pattern is computed. From each vertex of the pattern, a Depth-First-Search is done. A list of the visited vertex is kept in memory to avoid backtracking. When the Depth-First-Search reaches the same vertex it originated from a cycle is created. 

To avoid duplicate cycle in the cycle list, an equality operator is defined, and the vertex list ![$(\mathcal{C}_i)$](http://latex.codecogs.com/gif.latex?%24%28%5Cmathcal%7BC%7D_i%29%24) of the constraint is ordered. The order is made such as  :
* ![$\mathcal{C}_0$](http://latex.codecogs.com/gif.latex?%24%5Cmathcal%7BC%7D_0%24) is the minimal index of the vertex list
* ![$\mathcal{C}_1<\mathcal{C}_{last}$](http://latex.codecogs.com/gif.latex?%24%5Cmathcal%7BC%7D_1%3C%5Cmathcal%7BC%7D_%7Blast%7D%24)

Finally the equality operator implies the equality of the size of the constraint and of the elements of the constraint: ![$\mathcal{C}^A=\mathcal{C}^B \iff |\mathcal{C}^A|=|\mathcal{C}^B|\land(\forall i, \mathcal{C}^A_i=\mathcal{C}^B_i)$ ](http://latex.codecogs.com/gif.latex?%24%5Cmathcal%7BC%7D%5EA%3D%5Cmathcal%7BC%7D%5EB%20%5Ciff%20%7C%5Cmathcal%7BC%7D%5EA%7C%3D%7C%5Cmathcal%7BC%7D%5EB%7C%5Cland%28%5Cforall%20i%2C%20%5Cmathcal%7BC%7D%5EA_i%3D%5Cmathcal%7BC%7D%5EB_i%29%24)


We could also check that ![$\mathcal{C}^A=\mathcal{C}^B \iff\mathcal{C}^A\oplus\mathcal{C}^B=\emptyset$](http://latex.codecogs.com/gif.latex?%24%5Cmathcal%7BC%7D%5EA%3D%5Cmathcal%7BC%7D%5EB%20%5Ciff%5Cmathcal%7BC%7D%5EA%5Coplus%5Cmathcal%7BC%7D%5EB%3D%5Cemptyset%24) .

### Path constraint
The file `path_constraint.txt` contains the output of the path generation algorithm.

A path constraint is defined as a path between two vertex with the same label. From each vertex of the pattern, a Depth-First-Search is done. A list of the visited vertex is kept in memory to avoid backtracking. A path is created when the Depth-First-Search reaches a vertex with the same label than the one it originated from. 

Duplicates are avoided through a similar method as the one presented in cycle constraint generation.

### TDS constraint
Template Driven Search constraint are more complicated as they involve searching for a subgraph of the pattern. The generated files contain a set of edges describing the subgraph. An edge is noted `(A B)` where `A` is the initial vertex or label, and `B` the target vertex or label.

### TDS cycle constraint
The file `tds_cycle_constraint.txt` contains the output of the tds cycle generation algorithm.

TDS Cycle constraint are made out of the aggregation of cycle constraint sharing at least an edge. The algorithm merges cycle constraint 2 by 2 after checking that they have a common edge. 

The specific implementation uses a doubly linked list for this part via `deque()`. When two subgraph of the list are found to share an edge, their union is computed and added at the beginning of the doubly linked list while the elements are removed. The algorithm stops when every element has been compare to the others. The elements that are not issued from a merge are removed because they are a duplicate of a cycle constraint.

### TDS path constraint
The file `tds_path_constraint.txt` contains the output of the tds path generation algorithm.

TDS Path constraint are made out of the aggregation of path constraint originating from the same label sharing at least a vertex. The algorithm merges path constraint 2 by 2 after checking that they have a common vertex. 

The specific implementation uses a doubly linked list for this part via `deque()`. When two subgraph of the list are found to share a vertex, their union is computed and added at the beginning of the doubly linked list while the elements are removed. The algorithm stops when every element has been compare to the others. The elements that are not issued from a merge are removed because they are a duplicate of a path constraint.

### TDS subtemplate constraint
The file `tds_subtemplate_constraint.txt` contains the output of the tds subtemplate generation algorithm.

TDS Subtemplate constraint are made out of the aggregation of TDS path constraint and TDS cycle constraint sharing at least a vertex. The algorithm merges the constraint 2 by 2 after checking that they have a common vertex. 

The specific implementation uses a doubly linked list for this part via `deque()`. When two subgraph of the list are found to share a vertex, their union is computed. If the difference between the union and each of the subgraph is not empty, the union is added at the beginning of the doubly linked list while the elements are removed. This avoids duplicate. The algorithm stops when every element has been compared to the others. The elements that are not issued from a merge are removed because they are a duplicate of a tds cyclic or path constraint.


