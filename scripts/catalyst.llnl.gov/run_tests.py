import sys
import atexit
import time
import subprocess
import os.path

USE_PDEBUG = False
DEBUG = False

USE_DIMMAP = False
USE_DIMMAP_FOR_TUNE = True
NORUN = True

if USE_DIMMAP:
	graph_dir = "/dimmap/"
else:
	graph_dir = "/l/ssd/"

log_dir = "logs/"
executable_dir = "src/"
executable = "run_bfs" #"generate_graph"

command_strings = []
test_count = 0
test_motivation = ""


slurm_options = "-u --clear-ssd  -enable-hyperthread "
if USE_DIMMAP:
	slurm_options += " --di-mmap=" + str(96*1024*256) + " "
elif USE_DIMMAP_FOR_TUNE:
	slurm_options += " --di-mmap=" + str(10*256) + " "



def log(s):
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

	if DEBUG:
		log_dir += "debug/"
	else:
		log_dir += time_stamp + "/"
		os.makedirs(log_dir)

	log_file_name = log_dir + "run_tests.log"

	if len(sys.argv) > 1:
		test_motivation= str(sys.argv[1])
	elif DEBUG:
		test_motivation  = "Debuging..."
	else:
		test_motivation = raw_input("Please test motivation: ")
	log("Test Motivation:\n" + test_motivation)


	if len(sys.argv) > 2:
		line_name= str(sys.argv[2])
	elif DEBUG:
		line_name  = "Debuging..."
	else:
		line_name = raw_input("Please line name: ")

	log("Line Name:" + line_name)

	sbatch_file = log_dir + "batch.sh"

	if DEBUG:
		executable = executable_dir+executable
	else:
		cmd = ['cp', executable_dir+executable, log_dir+executable]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		while (p.poll() == None):
			pass
		executable = log_dir+executable

def generate_shell_file():
	global slurm_options
	with open(sbatch_file, 'w') as f:
		f.write("#!/bin/bash\n")
		i = 0
		for cmd in command_strings:
			nodes = cmd[0]
			processes = cmd[1]
			cmd_str = cmd[2]

			l_slurm_options = slurm_options
			if USE_PDEBUG:
				l_slurm_options += " -ppdebug " #-w catalyst324"
			# elif (nodes >= 128):
			# 	l_slurm_options += "-pdit128"
			# elif (nodes >= 64):
			# 	l_slurm_options += "-pdit64_1"
			# elif (nodes >= 32):
			# 	l_slurm_options += "-pdit36"

			if DEBUG:
				time.sleep(.01)
				cmd_log_fname = log_dir + "test_" + str(time.time()) + "_" + str(i) + ".out"
			else:
				cmd_log_fname = log_dir  +"test_" + str(i) + ".out"

			srun_cmd = "srun " + l_slurm_options \
				+ " -N" +str(nodes) + " -n" + str(processes) + " " \
				+ cmd_str  + " 2>&1 >> " + cmd_log_fname + "& \n"

			log(str(i) + ": " + srun_cmd)
			f.write(srun_cmd)

			i +=1

def execute_shell_file():
	if not NORUN:
		cmd = ['sh', sbatch_file]
		subprocess.call(cmd)

def add_command(nodes, processes, cmd):
	global test_count

	if not DEBUG:
		cmd_log_fname = log_dir+"test_"+str(test_count)+".out"

		log(str(test_count) + ":\t" + " ".join(cmd))

		with open(cmd_log_fname, 'w') as f:
			temp = "Test number: %d\n" %(test_count)
			f.write(temp)
			temp = "SRun Args: %s\n" %(" ".join(cmd))
			f.write(temp)
			temp = "Nodes: %d\n" %(nodes)
			f.write(temp)
			temp = "Processes: %d\n" %(processes)
			f.write(temp)

			temp = "Motivation: %s\n" %(test_motivation)
			f.write(temp)
			f.write("\n")

	test_count += 1

	command_strings.append([nodes, processes, " ".join(cmd)])

def create_commands(initial_scale, scale_increments, max_scale,
	inital_nodes, node_multipler, max_nodes,
	intial_threshold, threshold_multiplier):


	graph_file = graph_dir+"out.graph"

	save_file = 0
	compare_files = 0
	test_type = "RMAT"


	scale = initial_scale
	nodes = inital_nodes
	degree_threshold = intial_threshold

	while (nodes <= max_nodes and (scale <= max_scale or max_scale == -1) ):
		processes = 24 * nodes

		cmd = [executable, test_type, str(scale), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files)]
		add_command(nodes, processes, cmd)

		nodes *= node_multipler
		scale += scale_increments
		degree_threshold *= threshold_multiplier


init_test_dir()

if DEBUG:
	create_commands(17, 1, 21, 1, 1, 1, 1024, 1)
else:
	#create_commands(17, 1, 30, 1, 1, 1, 1024, 1)
	#create_commands(25, 1, 31, 1, 1, 1, 1024, 2)
	create_commands(31, 1, 37, 1, 2, 64, 65536*2, 2)


#Data Scaling test spawning
#create_commands(29, 1, 31, 1, 1, 1, 1024, 1)

#Weak Scaling test spawning
#create_commands(20, 2, -1, 1, 2, 64, 1024, 2)

#make bash file and run it
generate_shell_file()
execute_shell_file()

log("Finished after generating %d Srun Tasks\n" %(test_count))


