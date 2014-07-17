import sys
import atexit
import time
import subprocess
import os.path

RunTest = True
DEBUG_SCRIPT_TESTS = False
VERBOSE = True

log_dir = "/g/g17/mrdalek/havoqgt/build/catalyst.llnl.gov/logs/"
executable_dir = "/g/g17/mrdalek/havoqgt/build/catalyst.llnl.gov/src/"
executable = "run_bfs"

command_strings = []
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

	log_file_name = log_dir + "run_tests.log"
	log("Test Motivation:")
	if len(sys.argv) == 2:
		log(str(sys.argv[1]))
	else:
		var = raw_input("Please test motivation: ")
		log(var)

	if DEBUG_SCRIPT_TESTS:
		log("DEBUG_SCRIPT_TESTS = True")


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



def execute_shell_file():
	cmd = ['sh', sbatch_file]
	subprocess.call(cmd)

def add_command(nodes, processes, cmd):
	global test_count

	cmd_log_fname = log_dir+"test_"+str(test_count)+".out"

	cmd.append(">>")
	cmd.append(cmd_log_fname)
	cmd.append("2>&1")


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
		str_processes = "-n%d" %(processes)
		str_nodes = "-N%d" %(nodes)

		cmd = ['srun', "--clear-ssd", "--di-mmap=" + str(96*1024*256)]

		if (nodes >= 128):
			cmd.append("-pdit128")
		elif (nodes >= 64):
			cmd.append("-pdit64_1")
		elif (nodes >= 32):
			cmd.append("-pdit36")


		cmd.extend([str_processes, str_nodes, executable, test_type, str(scale), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files)])
		add_command(nodes, processes, cmd)

		nodes *= node_multipler
		scale += scale_increments
		degree_threshold *= threshold_multiplier


if RunTest:
	init_test_dir()

	if DEBUG_SCRIPT_TESTS:
		create_commands(17, 1, 20, 1, 1, 1, 1024, 1)
	else:
		#Data Scaling test spawning
		create_commands(17, 1, 28, 1, 1, 1, 1024, 1)

		#Weak Scaling test spawning
		#create_commands(20, 2, -1, 1, 2, 64, 1024, 2)

	#make bash file and run it
	generate_shell_file()
	execute_shell_file()

	log("Finished after generating %d Srun Tasks\n" %(test_count))


