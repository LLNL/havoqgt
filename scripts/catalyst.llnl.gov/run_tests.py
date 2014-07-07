import atexit
import time
import subprocess
import os.path

time_stamp = time.time()
base_string =  "gen_results_"
log_file = base_string + str(time_stamp) + ".out"
log_file_csv = base_string + str(time_stamp) + ".csv"
while os.path.isfile(log_file):
	log_file = base_string+ str(time_stamp) + ".out"
	log_file_csv = base_string + str(time_stamp) + ".csv"

out_file = open(log_file, 'w')
csv_file = open(log_file_csv, 'w')


test_count = 0
running_process = []
def kill_jobs():
	global running_process
	for values in running_process:
		values[0].kill()
		while (values[0].poll() == None):
			pass
		log_output(values)

	subprocess.Popen(["scancel", "-u mrdalek"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	print ("Canceling jobs, wait 60 seconds")
	time.sleep(60)
	print ("Exitd")

atexit.register(kill_jobs)

headers = ["Process", "Nodes", "HAVOQGT_MAILBOX_NUM_IRECV", "HAVOQGT_MAILBOX_NUM_ISEND", "HAVOQGT_MAILBOX_AGGREGATION", "HAVOQGT_MAILBOX_TREE_AGGREGATION", "HAVOQGT_MAILBOX_PRINT_STATS", "Building graph type:", "Building graph Scale", "Hub threshold", "PA-beta", "File name ", "Load from disk", "Delete on Exit", "count_edge_degrees time", "partition_low_degree time", "calculate_overflow time", "partition_high_degree time", "delegate_partitioned_graph time", "Max Vertex Id", "Count of hub vertices", "Total percentage good hub edges", "total count del target", "Total percentage of localized edges", "Global number of edges", "Number of small degree", "Number of hubs", "oned imbalance", "hubs imbalance", "TOTAL imbalance ", "Max Degree ", "Error"]

header = ", ".join(headers)
csv_file.write(header+"\n")

def log_output(values):
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

	log_output_csv(p.returncode, out, err, values[3], values[4])



def log_output_csv(exit_code, string, string_err, processes, nodes):
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

	csv_file.write(",".join(temp)+"\n")
	csv_file.flush()

executable = "/g/g17/mrdalek/havoqgt/build/catalyst.llnl.gov/src/generate_graph"
graph_file = "out.graph"
save_file = 0
compare_files = 0
test_type = "RMAT"

def spawn_data_scaling(initial_scale, scale_increments, max_scale, intial_threshold, threshold_multiplier):
	global test_count
	global running_process
	global executable
	global graph_file
	global save_file
	global compare_files
	global test_type

	nodes = 1
	processes = 24
	scale = initial_scale
	degree_threshold = intial_threshold

	while (scale <= max_scale):

		str_processes = "-n %d" %(processes)
		str_nodes = "-N %d" %(nodes)
		cmd = ['srun', "--clear-ssd", str_processes, str_nodes, executable, test_type, str(scale), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files)]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
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
		cmd = ['srun', "--clear-ssd", str_processes, str_nodes, executable, test_type, str(scale), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files)]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		running_process.append([p, str_processes + " " + str_nodes, test_count, processes, nodes])

		nodes *= node_multipler
		scale += scale_increments
		degree_threshold *= threshold_multiplier

		test_count += 1



#Data Scaling test spawning
spawn_data_scaling(25, 1, 29, 1024, 1)
spawn_data_scaling(25, 1, 29, 1024, 2)

#Weak Scaling test spawning
spawn_weak_scaling(28, 1, 2, 2, 32, 1024, 1)
spawn_weak_scaling(28, 1, 2, 2, 16, 1024, 2)

print "Generated %d Srun Tasks" %(test_count)

i = 0
while (len(running_process) != 0):
	time.sleep(120)
	if (i >= len(running_process)):
		i = 0

	values = running_process[i]
	if (values[0].poll() != None):
		log_output(values)
		running_process.remove(values)
	i += 1
