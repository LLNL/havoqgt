import sys
import atexit
import time
import subprocess
import os.path
import datetime

USE_PDEBUG = False

# ---- Configuration ----- #
N_NODES = 1
N_PROCS = 1


DEBUG = True
if DEBUG:
	USE_PDEBUG = True
USE_DIMMAP = True
USE_DIMMAP_FOR_TUNE = False
MONITOR_IO = False
MEMSIZE_DIMMAP = 1024*256*4
#MEMSIZE_DIMMAP = 1024*256*2*N_NODES*N_PROCS
GLOBAL_LOG_FILE = "/g/g90/iwabuchi/logs/sbatch_experiments_graphpartionning_test_0713.log"

NORUN = False
VERBOSE = True
USE_CATALYST = True
DELETE_WORK_FILES = False
SEGMENT_SIZE = 39
# --------------------------- #

if USE_DIMMAP:
	graph_dir = "/dimmap/"
else:
	if USE_CATALYST:
		graph_dir = "/l/ssd/"
	else:
		graph_dir = "/usr/localdisk/fusion/"

log_dir = "logs/"
executable_dir = "src/"
executable = "dynamicdegreeawaregraphstore_bench"


command_strings = []
test_count = 0

def log(s):
	if VERBOSE:
		print s
	with open(log_file_name, 'a') as f:
		f.write(s + "\n")

def log_global():
	fl = open(GLOBAL_LOG_FILE, 'a')
	if not NORUN:
	 fl.write("Job ID: " + job_id + "\n")
	fl.write("Log dir: " + log_dir + "\n")
	fl.write("graph_dir: " + graph_dir + "\n")
	fl.write("motivation: " + motivation + "\n")
	fl.write("DEBUG: " + str(DEBUG) + ", ")
	fl.write("USE_DIMMAP: " + str(USE_DIMMAP) + ", ")
	fl.write("USE_DIMMAP_FOR_TUNE: " + str(USE_DIMMAP_FOR_TUNE) + ", ")
	if USE_DIMMAP or USE_DIMMAP_FOR_TUNE:
		fl.write("MEMSIZE_DIMMAP: " + str(MEMSIZE_DIMMAP) + ", ")
	fl.write("USE_CATALYST: " + str(USE_CATALYST) + ", ")
	fl.write("DELETE_WORK_FILES: " + str(DELETE_WORK_FILES) + ", ")
	fl.write("MONITOR_IO: " + str(MONITOR_IO) + "\n")
	fl.write("---------------------------------------------\n\n")
	fl.close()

def init_test_dir():
	global log_dir
	global log_file_name
	global sbatch_file
	global executable
	global io_monitoring_report_file
	global motivation

	if DEBUG:
		log_dir += "debug/"
		if not os.path.exists(log_dir):
			os.makedirs(log_dir)

	time_stamp = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
	while os.path.exists(log_dir+time_stamp):
		time_stamp = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

	if USE_PDEBUG:
		log_dir += "pdebug_"

	log_dir += time_stamp + "/"
	os.makedirs(log_dir)

	log_file_name = log_dir + "run_tests.log"
	log("Test Motivation:")
	if len(sys.argv) == 2:
		log(str(sys.argv[1]))
		motivation = str(sys.argv[1])
	elif DEBUG:
		log("Debuging...")
		log(log_dir)
		motivation = "Debuging..."
	else:
		var = raw_input("Please input test motivation: ")
		log(log_dir + ": " + var)
		motivation = var

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
		slurm_options += "--di-mmap=npages=" + str(MEMSIZE_DIMMAP) + ",ver=1.1.21d,ra_tune=0 --enable-hyperthreads "
	elif USE_DIMMAP_FOR_TUNE:
		slurm_options += "--di-mmap=npages=" + str(MEMSIZE_DIMMAP) + ",ver=1.1.21d,ra_tune=0 "
	else:
	 	slurm_options += "--di-mmap=npages=" + str(1024*256*2) + ",ver=none,ra_tune=0 "

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

			# s += "export NUM_EDGES=157286400000\n"

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

			s += block_start + "dmesg | tee -n 200 \n" + block_end
			s += "dmesg | tee -n 200 \n"

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
	global job_id

	if not NORUN:
		cmd = ['sh', sbatch_file]
		job_id = subprocess.check_output(cmd)

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

def create_commands(initial_scale, scale_increments, max_scale, delete_ratio_list):

	graph_file = graph_dir+"out.graph"

	for k in delete_ratio_list:
		delete_segment_file = 0
		chunk_size = 20
		edges_factor = 16
		scale = initial_scale

		while (scale <= max_scale):
			cmd = [executable, str(scale), str(edges_factor), graph_file, str(SEGMENT_SIZE), str(delete_segment_file), str(chunk_size), str(k),
				"/g/g90/iwabuchi/tmp/file_list"]
			add_command(N_NODES, N_PROCS, cmd)
			scale = scale + scale_increments

init_test_dir()

low_deg_tlh_list = [0]
delete_ratio_list = [0]

create_commands(27, 1, 27, delete_ratio_list)

#make bash file and run it
generate_shell_file()
execute_shell_file()

log("Finished after generating %d Srun Tasks\n" %(test_count))
log_global()


