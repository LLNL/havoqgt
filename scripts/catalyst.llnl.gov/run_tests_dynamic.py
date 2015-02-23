import sys
import atexit
import time
import subprocess
import os.path
import datetime

USE_PDEBUG = False
DEBUG = False
if DEBUG:
	USE_PDEBUG = True
USE_DIMMAP = False
USE_DIMMAP_FOR_TUNE = False
NORUN = False
VERBOSE = True
USE_CATALYST = True
DELETE_WORK_FILES = False

MONITOR_IO = False

if USE_DIMMAP:
	graph_dir = "/dimmap/"
else:
	if USE_CATALYST:
		graph_dir = "/l/ssd/"
	else:
		graph_dir = "/usr/localdisk/fusion/"

log_dir = "logs/"
executable_dir = "src/"
executable = "generate_graph_dynamic"
#executable = "generate_graph" #"run_bfs"

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
	global io_monitoring_report_file

	if DEBUG:
		log_dir += "debug/"
		if not os.path.exists(log_dir):
			os.makedirs(log_dir)

	time_stamp = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
	while os.path.exists(log_dir+time_stamp):
		time_stamp = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

	log_dir += time_stamp + "/"
	os.makedirs(log_dir)

	log_file_name = log_dir + "run_tests.log"
	log("Test Motivation:")
	if len(sys.argv) == 2:
		log(str(sys.argv[1]))
	elif DEBUG:
		log("Debuging...")
	else:
		var = raw_input("Please test motivation: ")
		log(var)

	sbatch_file = log_dir + "batch.sh"

	if DEBUG:
		executable = executable_dir+executable
	else:
		cmd = ['cp', executable_dir+executable, log_dir+executable]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		while (p.poll() == None):
			pass
		executable = log_dir+executable


	if MONITOR_IO:
		io_monitoring_report_file = log_dir + "io_monitering_report.log"

def generate_shell_file():
	block_start = "echo -e \"\\n\\n------------------------------------\"\n"
	block_end = "echo -e \"------------------------------------\"\n"
	i = 0

	if USE_CATALYST:
		slurm_options = " --clear-ssd "
	else:
		slurm_options = ""

	if USE_DIMMAP:
		slurm_options += "--di-mmap=" + str(1024*256*12) + ",ver=1.1.20d --enable-hyperthreads "
	elif USE_DIMMAP_FOR_TUNE:
		slurm_options += "--di-mmap=" + str(15545139) + " "

	if USE_PDEBUG:
		slurm_options += "-ppdebug"

	with open(sbatch_file, 'w') as f:
		f.write("#!/bin/bash\n")
		for cmd in command_strings:
			nodes = cmd[0]
			processes = cmd[1]
			cmd_str = cmd[2]

			if DEBUG:
				cmd_log_fname = log_dir+"test_"+str(i)+".out"
				cmd_error_log_fname = log_dir+"test_"+str(i)+".out"
			else:
				cmd_log_fname = log_dir+"test_"+str(i)+".out"
				cmd_error_log_fname = log_dir+"test_"+str(i)+".out"

			sbatch = "sbatch " + slurm_options + " -N" +str(nodes) + " -o" + cmd_log_fname + " -e" + cmd_error_log_fname + " << EOF \n"

			s = "#!/bin/sh\n"

			s += block_start + "echo Nodes: \n" + block_end
			s += "echo \"SLURM_NODELIST = \$SLURM_NODELIST \"\n"

			s += block_start + "echo Tuned Info: \n" + block_end
			s += "echo \"/proc/sys/vm/dirty_ratio = \$(cat /proc/sys/vm/dirty_ratio)\" \n"
			s += "echo \"/proc/sys/vm/dirty_background_ratio = \$(cat /proc/sys/vm/dirty_background_ratio)\" \n"
			s += "echo \"/proc/sys/vm/dirty_expire_centisecs = \$(cat /proc/sys/vm/dirty_expire_centisecs)\" \n"

			s += block_start + "echo free -m \n" + block_end
			s += "free -m \n"

			s += block_start + "echo Top 10 for memory using process \n" + block_end
			s += "ps alx  | awk '{printf (\"%d\\t%s\\n\", \\$8, \\$13)}' | sort -nr | head -10 \n"
			if USE_CATALYST:
				s += block_start + "echo df -h /l/ssd \n" + block_end
				s += "df -h /l/ssd  \n"
			else:
				s += block_start + "echo df -h /usr/localdisk/fusion \n" + block_end
				s += "df -h /usr/localdisk/fusion \n"

			s += block_start + "echo io-stat -m | grep md0 2>&1\n" + block_end
			s += "iostat -m | grep Device 2>&1 \n"
			s += "iostat -m | grep md0 2>&1 \n"

			if MONITOR_IO:
				s += block_start + "echo start I/O monitoring \n" + block_end
				s += "iostat -d -m -t -x -p md0 10 > " + io_monitoring_report_file + " 2>&1 & \n"

			s += block_start + "echo Executable Log \n" + block_end
			s += "date \n"
			s += "srun -N" +str(nodes) + " -n" + str(processes) + " " + cmd_str  + " \n"
			s += "date \n"

			if MONITOR_IO:
				s += block_start + "echo stop I/O monitoring \n" + block_end
				s += "pkill iostat \n"
				s += "ps -a \n"

			s += block_start + "echo free -m \n" + block_end
			s += "free -m \n"

			if USE_CATALYST:
				s += block_start + "echo df -h /l/ssd \n" + block_end
				s += "df -h /l/ssd  \n"
			else:
				s += block_start + "echo df -h /usr/localdisk/fusion \n" + block_end
				s += "df -h /usr/localdisk/fusion \n"

			if USE_CATALYST:
				s += block_start + "echo du -sh /l/ssd/out.graph* \n" + block_end
				s += "du -sh /l/ssd/out.graph* \n"
			else:
				s += block_start + "echo du -sh /usr/localdisk/fusion/out.graph* \n" + block_end
				s += "du -sh /usr/localdisk/fusion/out.graph* \n"

			if USE_DIMMAP:
				s += block_start + "echo du -sh /dimmap/* \n" + block_end
				s += "du -sh /dimmap/*\n"

			if USE_CATALYST:
				s += block_start + "echo ls -lsth /l/ssd/ \n" + block_end
				s += "ls -lsth /l/ssd/\n"
			else:
				s += block_start + "echo ls -lsth /usr/localdisk/fusion/ \n" + block_end
				s += "ls -lsth /usr/localdisk/fusion/ \n"

			if USE_DIMMAP:
				s += block_start + "echo ls -lsth /dimmap/ \n" + block_end
				s += "ls -lsth /dimmap/\n"

			if USE_DIMMAP:
				s += block_start + "echo cat /proc/di-mmap-runtimeA-stats \n" + block_end
				s += "cat /proc/di-mmap-runtimeA-stats \n"

#			s += block_start + "echo dmesg \n" + block_end
#			s += "dmesg\n"

			s += block_start + "echo io-stat -m | grep md0 2>&1\n" + block_end
			s += "iostat -m | grep Device 2>&1 \n"
			s += "iostat -m | grep md0 2>&1 \n"

			if DELETE_WORK_FILES:
				if USE_CATALYST:
					s += "rm /l/ssd/out.graph*\n"
					s += "rm /dimmap/out.graph*\n"
				else:
					s += "rm /usr/localdisk/fusion/out.graph*\n"
					s += "rm /dimmap/out.graph*\n"

			s += "EOF\n\n"



			f.write(sbatch + s+ "\n\n")

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
			f.write("\n")

	test_count += 1

	command_strings.append([nodes, processes, " ".join(cmd)])

def create_commands(initial_scale, scale_increments, max_scale,
	inital_nodes, node_multipler, max_nodes,
	intial_threshold, threshold_multiplier, data_type,
	low_deg_tlh_s, low_deg_tlh_e, delete_ratio_list):

	graph_file = graph_dir+"out.graph"

	for k in delete_ratio_list:
		for i in range(low_deg_tlh_s, low_deg_tlh_e+1) :
			save_file = 0
			compare_files = 0
			test_type = "RMAT"
			chunk_size = 20
			edges_factor = 16
			scale = initial_scale
			nodes = inital_nodes
			degree_threshold = intial_threshold

			while (nodes <= max_nodes and (scale <= max_scale or max_scale == -1) ):
				processes = 1 * nodes

				cmd = [executable, test_type, str(scale), str(edges_factor), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files), str(chunk_size), data_type, str(i), str(k)]
				add_command(nodes, processes, cmd)

				nodes *= node_multipler
				scale += scale_increments
				degree_threshold *= threshold_multiplier

init_test_dir()

delete_ratio_list = [0]

create_commands(27, 1, 27, 1, 1, 1, 1024, 1, "HY_DA", 1, 1, delete_ratio_list)

#make bash file and run it
generate_shell_file()
execute_shell_file()

log("Finished after generating %d Srun Tasks\n" %(test_count))
