<h4>Distribute Pattern Matching on Large Metadata Graphs</h4>
<p>We present an example of searching a pattern in a R-MAT generated graph using our program. The code is developed on top of HavoqGT.</p>
<p>Clone (with SSH) the code from https://github.com/hpcresearchanddevelopment/patternmatching</p>
<p>git clone git@github.com:hpcresearchanddevelopment/patternmatching.git</p>
<p>You will require the latest releases of OpenMPI or MAVPICH2 and the Boost library (some Boost releases have bugs, e.g., 1.58; the code works fine with 1.57) to run HavoqGT. The code has only been tested on latest generation of Linux distributions. Once you have checked out the code, make sure you are on the master branch.</p>
<p>Go to the directory, patternmatching/build/quartz/:</p> 
<p>cd  patternmatching/build/quartz/</p>
<p>Setup the CMake environment by running the following script:</p> 
<p>./scripts/quartz/do_cmake.sh</p>
<p>(Make necessary adjustments to the script for CMake to work within your environment.)</p>

<h4>Graph Generation</h4>
<p>The first step is to generate a graph in HavoqGT format. Go to the directory, build/quartz/ and build the R-MAT generator:</p>
<p>make generate_rmat</p>
<p>Create a directory, e.g., /usr/graph/, to store the generated graph. Assuming you are in the Slurm environment, run the following command to generate a R-MAT graph:
<p>srun -N1 --ntasks-per-node=4 --distribution=block ./src/generate_rmat -s 21 -p 1 -f 1 -o /dev/shm/rmat -b /urs/graph/rmat</p>
<p>This will create a graph with four partitions, to be run on four MPI processes. This is a Scale 21 graph (notice the parameter for the -s flag). The mmap/binary graph file will be stored in the directory /usr/graph/.</p>

<h4>Input Pattern</h4>
<p>We will search the following Tree pattern on the graph we just created. The numeric values on each vertex is the label of the respective vertex.</p>
<div align="center"><img src="https://github.com/hpcresearchanddevelopment/patternmatching/blob/master/examples/doc/tree_0011.png" width="200" height="200"></div>
<p>We use degree information to create numeric vertex labels, computed using the formula ceil(log_2(d(v_i)+1)). Here, d(v_i) is the degree of a vertex v_i. The input pattern is available in the dircetory patternmatching/examples/rmat_log2_tree_pattern/.</p>

<h4>Searching a Pattern</h4>
<p>First, build the pattern matching executable:</p>
<p>make run_pattern_matching_beta</p> 
<p>Next, use the following command to search the pattern stored in patternmatching/examples/rmat_log2_tree_pattern/.</p> 
<p>(Note that we do not need to provide vertex labels for the Tree pattern as we will use labels based on vertex degree and the program will generate degree-based labels when no input label is provided, i.e., the -v flag is not set. The program requires a specific directory structure to output results. An example is available here, patternmatching/examples/results/.)</p> 
<p>srun -N1 --ntasks-per-node=4 --distribution=block ./src/run_pattern_matching_beta -i /dev/shm/rmat -b /usr/graph/rmat -p ../../examples/rmat_log2_tree_pattern/ -o ../../examples/results/</p>
<p>The program logs status information to the standard output so you know the current state of the execution.</p>

<h4>Results</h4>
<p>Next, we discuss how to collect and interpret the results and retrieve the pruned graph. You will find (python) scripts in patternmatching/examples/scripts/ that should help you to parse the result files.</p>
<p>The files in the (result) directory examples/results/0/all_ranks_active_vertices_count/ contain the number of active vertices at the end of each iteration. However, the results are distributed among multiple files (one file for each MPI rank). The script examples/scripts/total_active_count.py produces the global statistics (to the standard output) from the distributed files:</p>
<p>python ../../examples/scripts/total_active_count.py ../../examples/results/0/all_ranks_active_vertices_count/ > /tmp/vertices_count</p>
<p>The last entry in the output file (e.g, /tmp/vertices_count) indicates the final number of active vertices. (You can use the same script to obtain edge statistics, output to the directory examples/results/0/all_ranks_active_edges_count/).</p>
<p>The file, examples/results/0/result_superstep, contains the global runtime, the time (in seconds) required to complete each LCC and NLCC iteration. Sum of time to complete all iterations is the time to complete a search.</p>
<p>The files in the (result) directory examples/results/0/all_ranks_active_vertices/ contain the number of active vertices after search/pruning has completed. The results can be easily merged in to a single file:</p> 
<p>cat ../../examples/results/0/all_ranks_active_vertices/* > /tmp/vertices</p> 
<p>Following the same procedure, you can collect the list of final active edges (from the output in examples/results/0/all_ranks_active_edges/).</p>

<h4>Enumeration</h4>
<p>For the Tree pattern in this example, the last step in the execution enumerates the pattern in the pruned graph. The input for enumeration is the last entry in the file examples/rmat_log2_tree_pattern/0/pattern_non_local_constraint. Here, the same NLCC code walks the full template with work aggregation turned off. The subgraphs are output to the directory examples/results/0/all_ranks_subgraphs and are distributed among multiple files. </p>

<h4>Constraint Generation</h4>
<p>https://github.com/hpcresearchanddevelopment/patternmatching/tree/master/examples/constraint_generator</p>
