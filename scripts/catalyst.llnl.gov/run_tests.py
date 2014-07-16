import sys
import atexit
import time
import subprocess
import os.path


VERBOSE = False
EXECUTE = True

headers = ["Process", "Nodes", "HAVOQGT_MAILBOX_NUM_IRECV", "HAVOQGT_MAILBOX_NUM_ISEND", "HAVOQGT_MAILBOX_AGGREGATION", "HAVOQGT_MAILBOX_TREE_AGGREGATION", "HAVOQGT_MAILBOX_PRINT_STATS", "Building graph type:", "Building graph Scale", "Hub threshold", "PA-beta", "File name ", "Load from disk", "Delete on Exit", "count_edge_degrees time", "partition_low_degree time", "calculate_overflow time", "partition_high_degree time", "delegate_partitioned_graph time", "Total MB Written:", "Total MB Read:" ,"Max Vertex Id", "Count of hub vertices", "Total percentage good hub edges", "total count del target", "Total percentage of localized edges", "Global number of edges", "Number of small degree", "Number of hubs", "oned imbalance", "hubs imbalance", "TOTAL imbalance ", "Max Degree ", "BFS Time", "Count BFS", "AVERAGE BFS", "Visited total", "Error"]

csv_file = ""
out_file = ""

def log_output(values):
	global out_file
	try:
		p = values[0]
		out, err = p.communicate()
		out_file.write("=====================================================================\n")
		temp = "Test number: %d\n" %(values[2])
		out_file.write(temp)
		temp = "Exit Code: %s\n" %(p.returncode)
		out_file.write(temp)
		temp = "SRun Args: %s\n" %(values[1])
		out_file.write(temp)
		temp = "Output: %s\n" %(out)
		out_file.write(temp)
		temp = "Error: %s\n\n" %(err)
		out_file.write(temp)
		out_file.write("=====================================================================\n")
		out_file.flush()
	except:
		print "Exception 3.1: ", sys.exc_info()[0]
	try:
		log_output_csv(out, values[3], values[4])
	except:
		print "Exception 3.2: ", sys.exc_info()[0]



def log_output_csv(string, processes, nodes):
	global csv_file
	lines = string.split('\n')

	temp = []
	for h in headers:
		temp.append("")

	temp[0] = str(processes)
	temp[1] = str(nodes)

	for l in lines:
		l = l.strip()
		for h in xrange(2, len(headers), 1):
			if headers[h] in l:
				words = l.split(' ')
				temp[h] = words[len(words)-1]
				break

	csv_file.write("\t".join(temp)+"\n")
	csv_file.flush()

def generate_csv(ifile_name):
	global csv_file
	try:
		fin = open(ifile_name, 'r')
		csv_file = open(ifile_name + "_gen.csv", 'w')

		header = "\t, ".join(headers)
		csv_file.write(header+"\n")



		while True:
			fin.readline()
			fin.readline()

			line = fin.readline().strip()
			words = line.split(' ')
			exit_code = words[len(words)-1]

			line = fin.readline().strip()

			words = line.split(' ')
			nodes = words[len(words)-1]
			processes = words[len(words)-3]

			line = fin.readline().strip()
			string = ""
			while True:
				if "Error:" in line:
					break
				else:
					string += line + "\n"
					line = fin.readline().strip()

			while True:
				if "=====================================================================" in line:
					break
				else:
					line = fin.readline().strip()

			log_output_csv(exit_code, string, "", processes, nodes)
	except:
		print "Exception Gen: ", sys.exc_info()[0]

if False:
	generate_csv("gen_results_1404839527.99.out")
	exit(0)


log_dir = "/g/g17/mrdalek/havoqgt/build/catalyst.llnl.gov/logs/"
executable_dir = "/g/g17/mrdalek/havoqgt/build/catalyst.llnl.gov/src/"
executable = "run_bfs"


time_stamp = str(time.time())
log_file = executable + "_" + str(time_stamp) + ".out"
while os.path.isfile(log_dir+log_file):
	time_stamp = time.time()
	log_file = executable + "_" + str(time_stamp) + ".out"

log_file_csv = log_dir + executable + "_" + str(time_stamp) + ".csv"
log_file = log_dir + log_file

out_file = open(log_file, 'w')
csv_file = open(log_file_csv, 'w')


cmd = ['cp', executable_dir+executable, log_dir+executable+"_"+time_stamp]
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
while (p.poll() == None):
	pass

executable = log_dir+executable+"_"+time_stamp


test_count = 0
running_process = []
def kill_jobs():
	global running_process
	print ("Canceling Remaining jobs and Exiting")
	try:
		for values in running_process:
			values[0].kill()
			while (values[0].poll() == None):
				pass
			log_output(values)
	except:
		print "Exception 4: ", sys.exc_info()[0]
atexit.register(kill_jobs)


header = "\t ".join(headers)
csv_file.write(header+"\n")

graph_file = "out.graph"
save_file = 0
compare_files = 0
test_type = "RMAT"

def run_cmd(cmd):
	if VERBOSE:
		print " ".join(cmd)
	if EXECUTE:
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		return p
	return None

def spawn_data_scaling(initial_scale, scale_increments, max_scale, intial_threshold, threshold_multiplier, nodes):
	global test_count
	global running_process
	global executable
	global graph_file
	global save_file
	global compare_files
	global test_type

	processes = 24 * nodes
	scale = initial_scale
	degree_threshold = intial_threshold

	while (scale <= max_scale):

		str_processes = "-n %d" %(processes)
		str_nodes = "-N %d" %(nodes)
		cmd = ['srun', "--clear-ssd", "-di-mmap=" + str(96*1024*256) , str_processes, str_nodes, executable, test_type, str(scale), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files)]

		p = run_cmd(cmd)
		if p != None:
			running_process.append([p, str_processes + " " + str_nodes, test_count, processes, nodes])

		scale += scale_increments
		degree_threshold *= threshold_multiplier

		test_count += 1


def spawn_weak_scaling(initial_scale, scale_increments, inital_nodes, node_multipler, max_nodes, intial_threshold, threshold_multiplier):
	global test_count
	global running_process
	global executable
	global graph_file
	global save_file
	global compare_files
	global test_type

	nodes = scale_increments
	scale = initial_scale
	degree_threshold = intial_threshold

	while (nodes <= max_nodes):
		processes = 24 * nodes

		str_processes = "-n %d" %(processes)
		str_nodes = "-N %d" %(nodes)
		cmd = ['srun', "--clear-ssd", "--di-mmap=" + str(96*1024*256), str_processes, str_nodes, executable, test_type, str(scale), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files)]

		p = run_cmd(cmd)
		if p != None:
			running_process.append([p, str_processes + " " + str_nodes, test_count, processes, nodes])


		nodes *= node_multipler
		scale += scale_increments
		degree_threshold *= threshold_multiplier

		test_count += 1



#Data Scaling test spawning
test_script = False
if test_script:
	sleep_time = 5
	spawn_data_scaling(17, 1, 20, 1024, 1, 1)
else:
	sleep_time = 120
	#spawn_data_scaling(25, 1, 32, 1024, 1, 1)
	#spawn_data_scaling(26, 1, 30, 2048, 2, 1)
	spawn_data_scaling(17, 1, 30, 1024, 1, 1)

#Weak Scaling test spawning
if test_script:
	pass
	# spawn_weak_scaling(17, 1, 1, 2, 8, 1024, 1)
else:
	pass
	#spawn_weak_scaling(25, 1, 1, 2, 64, 1024, 1)
	#spawn_weak_scaling(26, 1, 2, 2, 64, 2048, 2) # skips first one, which is done above

print "Generated %d Srun Tasks" %(test_count)

i = 0
while (len(running_process) != 0):
	try:
		time.sleep(sleep_time)
	except:
		print "Exception 1: ", sys.exc_info()[0]

	if (i >= len(running_process)):
		i = 0

	values = running_process[i]
	if (values[0].poll() != None):
		try:
			log_output(values)
			running_process.remove(values)
		except ValueError as e:
			print "ValueError ({0}): {1}".format(e.errno, e.strerror)
		except:
			print "Exception 2: ", sys.exc_info()[0]

	i += 1
