import atexit
import time
import subprocess
import os.path

time_stamp = time.time()
log_file = "GraphTest_" + str(time_stamp) + ".out"
while os.path.isfile(log_file):
	log_file = "GraphTest_" + str(time_stamp) + ".out"

out_file = open(log_file, 'w')


executable = "/g/g17/mrdalek/havoqgt/build/catalyst.llnl.gov/src/generate_graph"
graph_file = "out.graph"
save_file = 0
compare_files = 0
test_type = "RMAT"

max_nodes = 64
nodes = 1

running_process = []
test_count = 0

def kill_jobs():
	for pair in running_process:
		pair[0].kill()


atexit.register(kill_jobs)

#Data Scaling test spawning
if True:
	nodes = 1
	processes = 24
	scale = 28
	max_scale = 32
	degree_threshold = 1024

	while (scale <= max_scale):

		str_processes = "-n %d" %(processes)
		str_nodes = "-N %d" %(nodes)
		cmd = ['srun', "--clear-ssd", str_processes, str_nodes, executable, test_type, str(scale), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files)]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		running_process.append([p, str_processes + " " + str_nodes, test_count])

		scale += 1
		test_count += 1
		degree_threshold *= 2

#Weak Scaling test spawning
if True:
	max_nodes = 64
	nodes = 2
	scale = 28
	degree_threshold = 1024

	while (nodes <= max_nodes):
		processes = 24 * nodes

		str_processes = "-n %d" %(processes)
		str_nodes = "-N %d" %(nodes)
		cmd = ['srun', "--clear-ssd", str_processes, str_nodes, executable, test_type, str(scale), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files)]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		running_process.append([p, str_processes + " " + str_nodes, test_count])

		nodes *= 2
		scale += 1
		test_count += 1
		degree_threshold *= 2


print "Generated %d Srun Tasks" %(test_count)

i = 0
while (len(running_process) != 0):
	time.sleep(60)
	if (i >= len(running_process)):
		i = 0
	pair = running_process[i]
	p = pair[0]
	if (p.poll() != None):
		out, err = p.communicate()
		out_file.write("=====================================================================\n")
		temp = "Test number: %d\n" %(pair[2])
		out_file.write(temp)
		temp = "Exit Code: %s\n" %(p.returncode)
		out_file.write(temp)
		temp = "SRun Args: %s\n" %(pair[1])
		out_file.write(temp)
		temp = "Output: %s\n" %(out)
		out_file.write(temp)
		temp = "Error: %s\n\n" %(err)
		out_file.write(temp)
		out_file.write("=====================================================================\n")

		running_process.remove(pair)

	i += 1


