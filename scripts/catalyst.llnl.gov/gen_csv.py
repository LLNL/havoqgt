import sys
import os

headers = ["Processes", "Nodes", "HAVOQGT_MAILBOX_NUM_IRECV", "HAVOQGT_MAILBOX_NUM_ISEND", "HAVOQGT_MAILBOX_AGGREGATION", "HAVOQGT_MAILBOX_TREE_AGGREGATION", "HAVOQGT_MAILBOX_PRINT_STATS", "Building graph type:", "Building graph Scale", "Hub threshold", "PA-beta", "File name ", "Load from disk", "Delete on Exit", "count_edge_degrees time", "partition_low_degree time", "calculate_overflow time", "partition_high_degree time", "delegate_partitioned_graph time", "Total MB Written:", "Total MB Read:" ,"Max Vertex Id", "Count of hub vertices", "Total percentage good hub edges", "total count del target", "Total percentage of localized edges", "Global number of edges", "Number of small degree", "Number of hubs", "oned imbalance", "hubs imbalance", "TOTAL imbalance ", "Max Degree ", "BFS Time", "Count BFS", "AVERAGE BFS", "Visited total", "Error"]



if len(sys.argv) == 1:
	print "No directory specified and no jobs were spawned"
	exit(-1)

pass_dir = str(sys.argv[1])
if os.path.exists(pass_dir):
	log_dir = pass_dir
else:
	print "Specified directory " + pass_dir +" not found"+"\n\tSearched:"+pass_dir +"\n\tSearched:"+log_dir+pass_dir
	exit(-1)


log_dir += "/"

if not os.path.exists(log_dir):
	print "Log directory("+log_dir+") not found."
	exit(-1)

print "Searching " + log_dir

with open(log_dir+"results.csv", "w") as fout:
	fout.write("\t ".join(headers)+"\n")
	for fname in os.listdir(log_dir):
		if fname.endswith(".out"):
			print "Parsing: " + fname
			with open(log_dir+fname, 'r') as fin:
				temp = []
				for h in headers:
					temp.append("")
				for line in fin:
					line = line.strip()
					for h in xrange(0, len(headers), 1):
						if headers[h] in line:
							words = line.split(' ')
							temp[h] = words[len(words)-1]
							break

				fout.write("\t".join(temp)+"\n")
