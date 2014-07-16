import sys
import atexit
import time
import subprocess
import os.path

RunTest = False
GenerateCSV = False

DEBUG_SCRIPT_TESTS = True
VERBOSE = True

log_dir = "/g/g17/mrdalek/havoqgt/build/catalyst.llnl.gov/logs/"
executable_dir = "/g/g17/mrdalek/havoqgt/build/catalyst.llnl.gov/src/"
executable = "run_bfs"

command_strings = []
high_node_threshold = 64
command_strings_high_nodes = []
test_count = 0

def log(s):
	if VERBOSE:
		print s
	with open(log_file_name, 'a') as f:
		f.write(s + "\n")

def init_test_dir():
	global log_dir
	global log_file_name
	global sbatch_file
	global executable

	time_stamp = str(time.time())

	while os.path.exists(log_dir+time_stamp):
		time_stamp = str(time.time())
	log_dir += time_stamp + "/"
	os.makedirs(log_dir)

	log_file_name = log_dir + "python_log.out"
	sbatch_file = log_dir + "batch.sh"


	cmd = ['cp', executable_dir+executable, log_dir+executable]
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	while (p.poll() == None):
		pass

	executable = log_dir+executable

def generate_shell_file():
	with open(sbatch_file, 'w') as f:
		f.write("#!/bin/bash\n")

		for s in command_strings:
			f.write(s + " &\n")

		f.write("#srun with more than " + str(high_node_threshold) + " nodes.\n")

		for s in command_strings_high_nodes:
			f.write(s + " \n")

		f.write("wait\n")

def execute_shell_file():
	pass

def add_command(nodes, processes, cmd):
	global test_count

	cmd_log_fname = log_dir+"test_"+str(test_count)+".out"

	cmd.append(">>")
	cmd.append(cmd_log_fname)


	temp = str(test_count) + ":\t" + " ".join(cmd)
	log(temp)

	with open(cmd_log_fname, 'w') as f:
		temp = "Test number: %d\n" %(test_count)
		f.write(temp)
		temp = "SRun Args: %s\n" %(" ".join(cmd))
		f.write(temp)
		temp = "Nodes: %d\n" %(nodes)
		f.write(temp)
		temp = "Processes: %d\n" %(processes)
		f.write(temp)
		f.write("\n")

	test_count += 1

	if nodes >= high_node_threshold:
		command_strings_high_nodes.append(" ".join(cmd))
	else:
		command_strings.append(" ".join(cmd))

def create_commands(initial_scale, scale_increments, max_scale,
	inital_nodes, node_multipler, max_nodes,
	intial_threshold, threshold_multiplier):

	graph_file = "out.graph"
	save_file = 0
	compare_files = 0
	test_type = "RMAT"


	scale = initial_scale
	nodes = inital_nodes
	degree_threshold = intial_threshold

	while (nodes <= max_nodes and (scale <= max_scale or max_scale == -1) ):

		processes = 24 * nodes
		str_processes = "-n %d" %(processes)
		str_nodes = "-N %d" %(nodes)

		cmd = ['srun', "--clear-ssd", "--di-mmap=" + str(96*1024*256) , str_processes, str_nodes, executable, test_type, str(scale), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files)]
		add_command(nodes, processes, cmd)

		nodes *= node_multipler
		scale += scale_increments
		degree_threshold *= threshold_multiplier


if RunTest:
	init_test_dir()

	if DEBUG_SCRIPT_TESTS:
		sleep_time = 5
		create_commands(17, 1, 20,		1, 1, 1,		1024, 1)
	else:
		sleep_time = 120
		#Data Scaling test spawning
		create_commands(17, 1, 30,		1, 1, 1,		1024, 1)

		#Weak Scaling test spawning
		create_commands(20, 2, -1,		1, 2, 64,		1024, 2)

		#make bash file and run it
		generate_shell_file()
		execute_shell_file()

		log("Generated %d Srun Tasks\n" %(test_count))


headers = ["Processes", "Nodes", "HAVOQGT_MAILBOX_NUM_IRECV", "HAVOQGT_MAILBOX_NUM_ISEND", "HAVOQGT_MAILBOX_AGGREGATION", "HAVOQGT_MAILBOX_TREE_AGGREGATION", "HAVOQGT_MAILBOX_PRINT_STATS", "Building graph type:", "Building graph Scale", "Hub threshold", "PA-beta", "File name ", "Load from disk", "Delete on Exit", "count_edge_degrees time", "partition_low_degree time", "calculate_overflow time", "partition_high_degree time", "delegate_partitioned_graph time", "Total MB Written:", "Total MB Read:" ,"Max Vertex Id", "Count of hub vertices", "Total percentage good hub edges", "total count del target", "Total percentage of localized edges", "Global number of edges", "Number of small degree", "Number of hubs", "oned imbalance", "hubs imbalance", "TOTAL imbalance ", "Max Degree ", "BFS Time", "Count BFS", "AVERAGE BFS", "Visited total", "Error"]

if GenerateCSV:
	if not RunTest:
		if len(sys.argv) == 2:
			pass_dir = str(sys.argv[1])
			if os.path.exists(pass_dir):
				log_dir = pass_dir
			elif os.path.exists(log_dir+pass_dir):
				log_dir += pass_dir
			else:
				log("Specified directory " + pass_dir +" not found"
					+"\n\tSearched:"+pass_dir +"\n\tSearched:"+log_dir+pass_dir+"\n")
				exit(-1)
		else:
			log("No directory specified and no jobs were spawned\n")
			exit(-1)


	if not os.path.exists(log_dir):
		log("Log directory("+log_dir+") not found\.\n")
		exit(-1)


	with open(log_dir+"results.csv") as fout:
		fout.write("\t ".join(headers)+"\n")
		for file in os.listdir(log_dir):
			if file.endswith(".out"):
				log("Parsing: " + file + "\n")
				with open(file, 'r') as fin:
					temp = []
					for h in headers:
						temp.append("")
					for line in fin:
						line = line.strip()
						for h in xrange(0, len(headers), 1):
							if headers[h] in l:
								words = l.split(' ')
								temp[h] = words[len(words)-1]
								break

					fout.write("\t".join(temp)+"\n")
