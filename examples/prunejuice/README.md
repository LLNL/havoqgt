<h2>PruneJuice : A Distribute Solution for Pattern Matching in Metadata Graphs</h2>
<p>We present an example of searching a pattern in an R-MAT generated graph using our program. The code is developed on top of HavoqGT.</p>
<p>Clone the code from https://github.com/LLNL/havoqgt/</p>
<p>git clone https://github.com/LLNL/havoqgt.git</p>
<!--<p>git checkout develop</p>-->
<p>Follow the <a href="https://github.com/LLNL/havoqgt/tree/develop">instructions</a> to setup HavoqGT.</p>
<p>You will require the latest releases of OpenMPI or MAVPICH2, and the Boost library (some Boost releases have bugs, e.g., 1.58) to run HavoqGT. The code has only been tested on latest generation of Linux distributions. <!--Once you have checked out the code, make sure you are on the develop branch.--></p>
<p>Go to the directory, havoqgt/build/catalyst.llnl.gov/:</p> 
<p>cd havoqgt/build/catalyst.llnl.gov/</p>
<p>Setup the CMake environment by running the following script:</p> 
<p>./scripts/catalyst.llnl.gov/do_cmake.sh</p>
<p>(Make necessary adjustments to the script for CMake to work within your environment.)</p>

<h4>Graph Generation</h4>
<p>The first step is to generate a graph in the HavoqGT format. Go to the directory, build/catalyst.llnl.gov/ and build the R-MAT generator:</p>
<p>make generate_rmat</p>
<p>Create a directory, e.g., /usr/graph/, to store the generated graph. Assuming you are working within the Slurm environment, run the following command to generate an R-MAT generated graph:
<p>srun -N1 --ntasks-per-node=4 --distribution=block ./src/generate_rmat -s 21 -p 1 -f 1 -o /dev/shm/rmat -b /urs/graph/rmat</p>
<p>This will create a graph with four partitions, to be run on four MPI processes. This is a Scale 21 graph (note the argument for the -s flag). The mmap/binary graph file will be stored in the directory /usr/graph/.</p>

<h4>Input Pattern</h4>
<p>We will search the following Tree pattern in the graph we just created. The numeric values on each vertex is the label of the respective vertex.</p>
<div align="center"><img src="doc/tree_0001.png" width="200" height="200"></div>
<p>We use degree information to create numeric vertex labels, computed using the formula ceil(log_2(d(v_i)+1)). Here, d(v_i) is the degree of a vertex v_i. The input pattern is available in the dircetory havoqgt/examples/prunejuice/patterns/tree_0001/.</p>

<h4>Searching a Pattern</h4>
<p>First, build the pattern matching executable:</p>
<p>make run_pattern_matching</p> 
<p>Next, use the following command to search the pattern stored in havoqgt/examples/prunejuice/patterns/tree_0001/.</p> 
<p>Note that we do not need to provide vertex labels for the Tree pattern as we will use labels based on vertex degree and the program will generate degree-based labels when no input label is provided, i.e., the -v flag is not set. The program requires a specific directory structure to output results. An example is available here: havoqgt/examples/prunejuice/output/.</p> 
<p>srun -N1 --ntasks-per-node=4 --distribution=block ./src/run_pattern_matching_beta_x.x -i /dev/shm/rmat -b /usr/graph/rmat -p ../../examples/prunejuice/patterns/tree_0001/ -o ../../examples/prunejuice/output/</p>
<p>The program logs status information to the standard output. <!--so you know the current state of the execution.--></p>

<!--<h4>Results</h4>
<p>Next, we discuss how to collect and interpret the outputs and retrieve the pruned graph. You will find (python) scripts in havoqgt/examples/prunejuice/scripts/ that should help you to parse the output files.</p>
<p>The files in the (output) directory examples/prunejuice/output/0/all_ranks_active_vertices_count/ contain the number of active vertices at the end of each iteration. However, the outputs are distributed among multiple files (one file for each MPI rank). The script examples/prunejuice/scripts/total_active_count.py produces the global statistics (to the standard output) from the distributed files:</p>
<p>python ../../examples/prunejuice/scripts/total_active_count.py ../../examples/prunejuice/output/0/all_ranks_active_vertices_count/ > /tmp/vertices_count</p>
<p>The last entry in the output file (e.g, /tmp/vertices_count) indicates the final number of active vertices. (You can use the same script to obtain edge statistics, output to the directory examples/prunejuice/output/0/all_ranks_active_edges_count/).</p>
<p>The file, examples/prunejuice/output/0/result_superstep, contains the global runtime, the time (in seconds) required to complete each LCC and NLCC iteration. Sum of time to complete all iterations is the time to complete a search.</p>
<p>The files in the (output) directory examples/prunejuice/output/0/all_ranks_active_vertices/ contain the number of active vertices after search/pruning has completed. The outputs can be easily merged into a single file:</p> 
<p>cat ../../examples/prunejuice/output/0/all_ranks_active_vertices/* > /tmp/vertices</p> 
<p>Following the same procedure, you can collect the list of final active edges (from the output in examples/prunejuice/output/0/all_ranks_active_edges/).</p>
!-->

<h4>Results</h4>
<p>The runtime, and vertex and edge paticipation information (i.e., match count) are dispalyed on the standard output. If requested, the program also produces full match enumeration and displays the result (i.e., match count) on the standard output. If the -o flag is set (i.e., it points to the output directory), the user can obtain the pruned solution subgraph, and collect additional data. The files in the directory, examples/prunejuice/output/all_ranks_active_vertices/, contain matching vertices, and, examples/prunejuice/output/all_ranks_active_edges/, contain matching edges, in the final solution subgraph.</p>
 
<p>For the Tree pattern used in this example, the last step in the execution enumerates the pattern in the pruned solution subgraph. The input for enumeration is the last entry in the file examples/prunejuice/patterns/tree_0001/pattern_nonlocal_constraint. The subgraphs are stored in the directory examples/prunejuice/output/all_ranks_subgraphs/. </p>

<h4>Constraint Generation</h4>
<p>examples/prunejuice/tools/query_preprocessor/cpp/</p>
