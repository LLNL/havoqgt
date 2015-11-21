import sys
import atexit
import time
import subprocess
import os
import os.path
import datetime
import argparse

GRAPH_PATH="/l/ssd/"
GRAPH_NAME="out.graph"
GRAPH_PATH_DIMMAP="/dimmap/"

def parse_options:
  global args

  parser = argparse.ArgumentParser()
  parser.add_argument("--non_debug", "-D", action="store_false",
                                  default=True,
                                  help="non debug mode")
  parser.add_argument("--dimmap", "-m", type=int,
                                  default=0,
                                  help="dimmap cache size (4KB)")
  parser.add_argument("--dimmap_tune", "-T", action="store_true",
                                       default=False,
                                       help="use dimmap for tune")
  parser.add_argument("--global_log_file", "-L", default="./tests_dynamic.log"
                                           help="monitor I/O statistics using iostat")
  parser.add_argument("--log_dir", "-l", default="./log"
                                           help="path log directory")
  parser.add_argument("--ppdebug", "-P", action="store_true",
                                   default=False,
                                   help="run on ppdebug")
  parser.add_argument("--segment_size", "-S", type=int,
                                        default="39"
                                        help="segment_size in GB")
  parser.add_argument("--chunk_size", "-c", type=int,
                                        default="6"
                                        help="chunk size in log10")
  parser.add_argument("--time_limit", "-t", type=int,
                                      default=23*60+59
                                      help="time limit in min.")
  parser.add_argument("--edge_files", "-E", default=""
                                      help="edge filies")
  parser.add_argument("--scale", "-s", type=int,
                                      default=18
                                      help="time limit in min.")
  parser.add_argument("--edge_factor", "-e", type=int,
                                      default=16
                                      help="time limit in min.")
  parser.add_argument("--verbose", "-v", action="store_true",
                                   default=False,
                                   help="verbose")
  parser.add_argument("--num_procs", "-n", type=int,
                                      default=1
                                      help="num processes")
  parser.add_argument("--num_nodes", "-N", type=int,
                                      default=1
                                      help="num nodes")
  parser.add_argument("executable", help="path for executable")
  parser.add_argument("motivation", default="debug",
                                    help="motivation")

  args = parser.parse_args()


command_strings = []
test_count = 0

def log(s):
	print s
	with open(log_file_name, 'a') as f:
		f.write(s + "\n")

def init_test_dir():
	global log_dir
	global log_file_name
	global sbatch_file
	global executable
	global io_monitoring_report_file
	global motivation

  log_dir += args.log_dir

	if !args.non_debug:
		log_dir += "/debug/"
		if not os.path.exists(log_dir):
			os.makedirs(log_dir)

	time_stamp = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
	while os.path.exists(log_dir+time_stamp):
		time_stamp = str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

	if args.ppdebug:
		log_dir += "/pdebug_"

	log_dir += time_stamp + "/"
	os.makedirs(log_dir)

	log_file_name = log_dir + "run_tests.log"
	log("Test Motivation:")
	motivation = str(arg.motivation)
  log(motivation)

	sbatch_file = log_dir + "batch.sh"

  exe_splt=args.executable.split("/")
  exe=exe_splt[-1]
	cmd = ['cp', args.executable, log_dir+exe]
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	while (p.poll() == None):
		pass
	executable = log_dir+exe

	if args.iostat:
		io_monitoring_report_file = log_dir + "io_monitering_report"


def generate_shell_file():
	block_start = "echo -e \"\\n\\n------------------------------------\"\n"
	block_end = "echo -e \"------------------------------------\"\n"
	i = 0

		slurm_options = " --clear-ssd "

	if args.dimmap > 0 && !args.dimmap_tune:
		slurm_options += " --di-mmap=npages=" + str(args.dimmap) + ",ver=1.1.21d,ra_tune=0 --enable-hyperthreads "
	elif args.dimmap > 0 && args.dimmap_tune:
		slurm_options += " --di-mmap=npages=" + str(args.dimmap) + ",ver=1.1.21d,ra_tune=0 "
	else:
	 	slurm_options += " --di-mmap=npages=" + str(1024*256*2) + ",ver=none,ra_tune=0 "

	if args.ppdebug:
		slurm_options += " -ppdebug "

	slurm_options += " -t" + str(args.time_limit) + " "
  slurm_options += " --msr-safe "

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
				s += "iostat -d -m -t -x -p md0 10 > " + io_monitoring_report_file + "_" + str(i)+".log" + " 2>&1 & \n"

			# s += "export NUM_EDGES=157286400000\n"

			s += block_start + "echo Executable Log \n" + block_end
			s += "date \n"
			s += "srun -N" +str(nodes) + " -n" + str(processes) + " -W" + str(args.time_limit * 60) + " " + cmd_str + " \n"
			s += "date \n"

			if args.iostat:
				s += block_start + "echo stop I/O monitoring \n" + block_end
				s += "pkill iostat \n"
				s += "ps -a \n"

			s += block_start + "echo free -m \n" + block_end
			s += "free -m \n"

			s += block_start + "echo df -h " + GRAPH_PATH + " \n" + block_end
			s += "df -h " + GRAPH_PATH + "  \n"

			s += block_start + "echo du -sh " + GRAPH_PATH + "/" + GRAPH_NAME + "* \n" + block_end
			s += "du -sh " + GRAPH_PATH + "/" + GRAPH_NAME + "* \n"

			if args.dimmap > 0 && !args.dimmap_tune:
				s += block_start + "echo du -sh " + GRAPH_PATH_DIMMAP + "* \n" + block_end
				s += "du -sh " +GRAPH_PATH_DIMMAP+ "/*\n"

			s += block_start + "echo ls -lsth " + GRAPH_PATH + " \n" + block_end
			s += "ls -lsth " + GRAPH_PATH + "\n"

			if args.dimmap > 0 && !args.dimmap_tune:
				s += block_start + "echo ls -lsth " + GRAPH_PATH_DIMMAP + " \n" + block_end
				s += "ls -lsth " + GRAPH_PATH_DIMMAP + "\n"

			if args.dimmap > 0 && !args.dimmap_tune:
				s += block_start + "echo cat /proc/di-mmap-runtimeA-stats \n" + block_end
				s += "cat /proc/di-mmap-runtimeA-stats \n"

      if args.verbose:
  			s += block_start + "dmesg | tail -n 500 \n" + block_end
  			s += "dmesg | tail -n 500 \n"

			s += block_start + "echo io-stat -m | grep md0 2>&1\n" + block_end
			s += "iostat -m | grep Device 2>&1 \n"
			s += "iostat -m | grep md0 2>&1 \n"

			s += "EOF\n\n"

			f.write(sbatch + s+ "\n\n")

			i +=1


def execute_shell_file():
	global job_id

	cmd = ['sh', sbatch_file]
	job_id = subprocess.check_output(cmd)


def add_command(nodes, processes, cmd):
	global test_count

	if args.non_debug:
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


def log_global():
  fl = open(args.global_log_file, 'a')
  fl.write("Job ID: " + job_id + "\n")
  fl.write("Log dir: " + log_dir + "\n")
  fl.write("graph_dir: " + graph_dir + "\n")
  fl.write("motivation: " + motivation + "\n")
  fl.write("non debug: " + str(args.non_debug) + ", ")
  fl.write("dimmap: " + str(args.dimmap) + ", ")
  fl.write("dimmap_tune: " + str(args.dimmap_tune) + ", ")
  fl.write("iostat: " + str(args.iostat) + "\n")
  fl.write("---------------------------------------------\n\n")
  fl.close()




parse_args()
init_test_dir()
cmd = [executable, "-s" + str(args.scale), "-e" + str(args.edge_factor),
       "-o" + GRAPH_PATH + GRAPH_NAME, "-f" + str(args.segment_size),
       "-c" + str(args.chunk_size), "-i" + args.edge_files]
add_command(args.num_nodes, args.num_procs, cmd)

#make bash file and run it
generate_shell_file()
execute_shell_file()

log("Finished after generating %d Srun Tasks\n" %(test_count))

if args.non_debug:
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

EDGES_FILELIST="./work/file_list_rnd_1d"
create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_rnd_2d"
#create_commands(27, 1, 27, delete_ratio_list)

#EDGES_FILELIST="./work/file_list_rnd_rnd"
#create_commands(27, 1, 27, delete_ratio_list)
# --------------------------------------------- #

