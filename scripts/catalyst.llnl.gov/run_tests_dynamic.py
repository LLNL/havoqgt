import sys
import atexit
import time
import subprocess
import os
import os.path
import datetime
import argparse

PWD = os.getcwd()
LOG_DIR = PWD + "/log"
EXEC_DIR = PWD + "/src"
EXEC_NAME = "graphconst_dynamicgstore_bench"
SEG_FILE_DIR = "/l/ssd"
SEG_FILE_NAME = "out.graph"
SEG_FILE_DIR_DIMMAP = "/dimmap"

NO_RUN = False

def parse_args():
	global args

	parser = argparse.ArgumentParser()


	parser.add_argument("gstore_name", help="graphstore name : DegAwareRHH or Baseline")

	# options for srun
	parser.add_argument("--num_procs", "-n", type=int,
							default=1,
							help="num processes")
	parser.add_argument("--num_nodes", "-N", type=int,
							default=1,
							help="num nodes")
	parser.add_argument("--time_limit", "-t", type=int,
							default=23*60+59,
							help="time limit in min.")
	parser.add_argument("--ppdebug", action="store_true",
							default=False,
							help="run on ppdebug")

	# options for dynamic
	parser.add_argument("--scale", "-s", type=int,
							default=18,
							help="time limit in min.")
	parser.add_argument("--edge_factor", "-e", type=int,
							default=16,
							help="time limit in min.")
	parser.add_argument("--segment_size", "-S", type=int,
							default=39,
							help="segment_size in GB")
	parser.add_argument("--chunk_size", "-c", type=int,
							default=6,
							help="chunk size in log10")
	parser.add_argument("--edge_files", "-E", default="",
							help="edge filies")

	# options for DI-MMAP
	parser.add_argument("--dimmap", type=int,
							default=0,
							help="DI-MMAP cache size (4KB)")
	parser.add_argument("--dimmap_tune", action="store_true",
							default=False,
							help="use DI-MMAP for tune")

	# options for log
	parser.add_argument("--motivation", "-m", default="motivation is not specified",
	                                                            help="motivation")
	parser.add_argument("--log_dir", "-l", default=LOG_DIR,
							help="path log directory")
	parser.add_argument("--global_log_file", default=PWD + "tests_dynamic_global_log.log",
							help="monitor I/O statistics using iostat")

	# options for this program
	parser.add_argument("--product_mode", "-p", action="store_true",
							default=False,
							help="production mode")

	parser.add_argument("--iostat", action="store_true",
							default=False,
							help="run iostat")
							
	parser.add_argument("--verbose", "-v", action="store_true",
							default=False,
							help="verbose")

	args = parser.parse_args()


def init_test_files():
 	global log_file_base
 	global log_file_name
 	global batch_file
 	global executable
 	global io_monitoring_report_file
	global motivation

	# init log directory
	log_dir = args.log_dir
	if not args.product_mode:
		log_dir += "/debug"
	if not os.path.exists(log_dir):
		os.makedirs(log_dir)
	
	# init log file base
	log_file_base = log_dir + "/"
	time_stamp = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
	while os.path.exists(log_file_base + time_stamp):
		time_stamp = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
	log_file_base += time_stamp
	
	# init log file	
	log_file_name = log_file_base + ".log"
	
	# init batch.sh	
	batch_file = log_file_base + "_batch.sh"
	
	# init exec file
	cmd = ['cp', EXEC_DIR + "/" + EXEC_NAME, log_file_base + "_" + EXEC_NAME]
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	while (p.poll() == None):
		pass
	executable = log_file_base + "_" + EXEC_NAME
	
	# iostat
	if args.iostat:
		io_monitoring_report_file = log_file_base + "_iostat.log"

	# log
	log_dsc("log_file_base", log_file_base)

	log_dsc("nodes", str(args.num_nodes))
	
	log_dsc("procs", str(args.num_procs))
	
	log_dsc("Test Motivation", args.motivation)

	if args.iostat:
		log("iostat running")	

# --- write to log file ---- #
def log_dsc(description, val):
	print description + ": " + val
	with open(log_file_name, 'a') as f:
		f.write(description + ": " + val + "\n")

def log(s):
	print s
	with open(log_file_name, 'a') as f:
		f.write(s + "\n")


# --- write to batch file --- #
def add_cmd(command):
	with open(batch_file, 'a') as f:
		f.write(command + "\n")

def add_dsc(description):
	block_start = "echo -e \"\\n\\n------------------------------------\""
	block_end = "echo -e \"------------------------------------\""

	with open(batch_file, 'a') as f:
		f.write("\n" + block_start + "\necho " + description + "\n" + block_end + "\n")
	
def add_dsc_cmd(command):
	add_dsc(command)
	add_cmd(command)


# --- generate batch file --- #
def generate_batch_file(exec_cmd):

	# - -- slurm options --- #
	slurm_options = " --clear-ssd "
	if args.dimmap > 0 and not args.dimmap_tune:
		slurm_options += " --di-mmap=npages=" + str(args.dimmap) + ",ver=1.1.21d,ra_tune=0 --enable-hyperthreads "
	elif args.dimmap > 0 and args.dimmap_tune:
		slurm_options += " --di-mmap=npages=" + str(args.dimmap) + ",ver=1.1.21d,ra_tune=0 "
	else:
		slurm_options += " --di-mmap=npages=" + str(1024*256*2) + ",ver=none,ra_tune=0 "

	if args.ppdebug:
		slurm_options += " -ppdebug "
		
	slurm_options += " -t" + str(args.time_limit) + " "
	slurm_options += " --msr-safe "
	log_dsc("slurm_options", slurm_options)
	
	# --- make batch file ---  #
	add_cmd("#!/bin/bash")
	add_cmd("sbatch " + slurm_options + \
			" -N" + str(args.num_nodes) + \
			" -o" + log_file_name + \
			" -e" + log_file_name + \
			" << EOF")
    
    	add_cmd("#!/bin/sh")
	add_dsc("Nodes: ")
	add_cmd("echo \"SLURM_NODELIST = \$SLURM_NODELIST\"")
	
	add_dsc("Tuned Info: ")
	add_cmd("echo \"/proc/sys/vm/dirty_ratio = \$(cat /proc/sys/vm/dirty_ratio)\"")
	add_cmd("echo \"/proc/sys/vm/dirty_background_ratio = \$(cat /proc/sys/vm/dirty_background_ratio)\"")
	add_cmd("echo \"/proc/sys/vm/dirty_expire_centisecs = \$(cat /proc/sys/vm/dirty_expire_centisecs)\"")

	add_dsc_cmd("free -m")

	add_dsc("Top 10 for memory using process")
	add_cmd("ps alx  | awk '{printf (\"%d\\t%s\\n\", \\$8, \\$13)}' | sort -nr | head -10")

	add_dsc_cmd("df -h " + SEG_FILE_DIR)

	add_dsc("io-stat -m | grep md0 2>&1")
	add_cmd("iostat -m | grep Device 2>&1")
	add_cmd("iostat -m | grep md0 2>&1")

	if args.iostat:
		add_dsc_cmd("iostat -d -m -t -x -p md0 10 > " + io_monitoring_report_file + "_" + str(i)+".log" + " 2>&1 &")

	add_dsc_cmd("date")
	add_dsc_cmd("srun -N" +str(args.num_nodes) + " -n" + str(args.num_procs) + " -W" + str(args.time_limit * 60) + " " + exec_cmd)
	add_dsc_cmd("date")

	if args.iostat:
		add_dsc("stop I/O monitoring")
		add_cmd("pkill iostat")
		add_cmd("ps -a")

	add_dsc_cmd("free -m")
	add_dsc_cmd("df -h " + SEG_FILE_DIR)
	add_dsc_cmd("du -sh " + SEG_FILE_DIR + "/*")

	if args.dimmap > 0 and not args.dimmap_tune:
		add_dsc_cmd("du -sh " +SEG_FILE_DIR_DIMMAP+ "/*")

	add_dsc_cmd("ls -lsth " + SEG_FILE_DIR)
	
	if args.dimmap > 0 and not args.dimmap_tune:
		add_dsc_cmd("ls -lsth " + SEG_FILE_DIR_DIMMAP)

	if args.dimmap > 0 and not args.dimmap_tune:
		add_dsc_cmd("cat /proc/di-mmap-runtimeA-stats")

	if args.verbose:
		add_dsc_cmd("dmesg | tail -n 500")
		
	add_dsc("io-stat -m | grep md0 2>&1")
	add_cmd("iostat -m | grep Device 2>&1")
	add_cmd("iostat -m | grep md0 2>&1")
	
	add_cmd("EOF\n")


def execute_shell_file():
	global job_id

	cmd = ['sh', batch_file]
	job_id = subprocess.check_output(cmd)

	log_dsc("job_id", job_id)

#def create_command(nodes, processes, cmd):

#    with open(log_file_name, 'w') as f:
 #       temp = "Test number: %d\n" %(test_count)
  #      f.write(temp)
   #     temp = "SRun Args: %s\n" %(" ".join(cmd))
    #    f.write(temp)
     #   temp = "Nodes: %d\n" %(nodes)
      #  f.write(temp)
       # temp = "Processes: %d\n" %(processes)
        #f.write(temp)
        #f.write("\n")

    #command_strings.append([nodes, processes, " ".join(cmd)])

	
def log_global():
	with open(args.global_log_file, 'a') as f:
		if not NO_RUN:
			f.write("Job ID: " + job_id + "\n")
		f.write("Log dir: " + log_dir + "\n")
		f.write("Log base: " + log_file_base + "\n")
		f.write("graph_dir: " + graph_dir + "\n")
		f.write("nodes " + args.num_nodes + "\n")
		f.write("procs " + args.num_procs + "\n")
		f.write("motivation: " + motivation + "\n")
		f.write("production mode: " + str(args.product_mode) + ", ")
		f.write("dimmap: " + str(args.dimmap) + ", ")
		f.write("dimmap_tune: " + str(args.dimmap_tune) + ", ")
		f.write("iostat: " + str(args.iostat) + "\n")
		f.write("---------------------------------------------\n\n")
		f.close()

if __name__ == '__main__':
	
	parse_args()
	init_test_files()
    
	if args.dimmap > 0 and not args.dimmap_tune:
		graph_path = SEG_FILE_DIR_DIMMAP + "/" + SEG_FILE_NAME
	else:
		graph_path = SEG_FILE_DIR + "/" + SEG_FILE_NAME

	exec_cmd = ""
	
	if args.edge_files:
		exec_cmd = executable + \
			" -g" + args.gstore_name + \
			" -o" + graph_path + \
			" -S" + str(args.segment_size) + \
			" -c" + str(args.chunk_size) + \
			" -E" + args.edge_files
	else:
		exec_cmd = executable + \
			" -g" + args.gstore_name + \
			" -s" + str(args.scale) + \
			" -e" + str(args.edge_factor) +\
			" -o" + graph_path + \
			" -S" + str(args.segment_size) + \
			" -c" + str(args.chunk_size)
	log_dsc("exec_cmd:", exec_cmd)

	generate_batch_file(exec_cmd)
	
	if not NO_RUN:
		execute_shell_file()

	if args.product_mode:
		log_global()

#EDGES_FILELIST="./work/file_list_rmat_s24"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_s24_64bit"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_rnd_1d_fine"
#create_commands(27, 1, 27, delete_ratio_list)

# -------------------------------------------- #
#EDGES_FILELIST="./work/file_list_srt_1d"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_srt_2d"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_srt_rnd"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_bfs_1d"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_bfs_2d"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_bfs_rnd"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_rnd_1d"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_rnd_2d"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_rnd_rnd"
#create_commands(27, 1, 27, delete_ratio_list)
# --------------------------------------------- #

