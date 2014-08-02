import sys
import os

headers = ["Processes", "Nodes", "HAVOQGT_MAILBOX_NUM_IRECV", "HAVOQGT_MAILBOX_NUM_ISEND", "HAVOQGT_MAILBOX_AGGREGATION", "HAVOQGT_MAILBOX_TREE_AGGREGATION", "HAVOQGT_MAILBOX_PRINT_STATS", "Building graph type:", "Building graph Scale", "Hub threshold", "PA-beta", "File name ", "Load from disk", "Delete on Exit", "count_edge_degrees time", "partition_low_degree time", "calculate_overflow time", "partition_high_degree time", "delegate_partitioned_graph time", "Total MB Written:", "Total MB Read:" ,"Max Vertex Id", "Count of hub vertices", "Total percentage good hub edges", "total count del target", "Total percentage of localized edges", "Global number of edges", "Number of small degree", "Number of hubs", "oned imbalance", "hubs imbalance", "TOTAL imbalance ", "Max Degree ", "BFS Time", "Count BFS", "AVERAGE BFS", "Visited total", "Error"]


def generate_csv_file(dir):
	with open(dir+"results.csv", "w") as fout:
		fout.write("\t ".join(headers)+"\n")
		for fname in os.listdir(dir):
			if fname.endswith(".out"):
				print "\tParsing: " + fname
				with open(dir+fname, 'r') as fin:
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


def generate_graph_files(filename, title, x_axis, y_axis, lines):

	tick_set = {}
	for line in lines:
		line[2].sort()
		for p in line[2]:
			tick_set[p[0]] = p[0]

	ticks = []
	ticks.extend(tick_set)
	ticks.sort()
	tick_str = "{" + ",".join(ticks) + "}"





	with open(filename, "w") as fout:
		fout.write("\\begin{tikzpicture}[font=\\bfseries]\n")
		fout.write("\\begin{semilogxaxis}[")
		fout.write("title=" + title + ", \n")
		fout.write("xlabel="+ x_axis + ", \n")
		fout.write("xtick="+ tick_str + ", \n")
		fout.write("xticklabels="+ tick_str + ", \n")
		fout.write("ylabel style={align=center},\n")
		fout.write("ylabel="+ y_axis + ", \n")
		fout.write("legend style={fill=none, draw=none, at={(0.5,-0.17)},anchor=north,legend cell align=left}\n")
		fout.write("]\n")

		for line in lines:
			points = ""
			for p in line[2]:
				points += "("+p[0]+","+p[1]+")\n"

			fout.write("%" + line[1] + "\n")
			fout.write("\\addplot plot coordinates {\n")
			fout.write(points + "};\n")
			fout.write("\\addlegendentry{"+line[0]+"}\n")

		fout.write("\\end{semilogxaxis}\n\\end{tikzpicture}\n")


def get_line(dir, y_str, x_str):
	comment = ""
	line_name = "not in file"
	with open(dir+"run_tests.log", "r") as fin:
		next_line = False
		for line in fin:
			line = line.strip()
			if next_line:
				comment = "%" + line + "\n"
				next_line = False
			if "Test Motivation:" in line:
				next_line = True
			if "Line name:" in line:
				words = line.split(' ')
				line_name = words[len(words)-1]

	data = []
	for fname in os.listdir(dir):
		if fname.endswith(".out"):
			count = 0
			y_val = ""
			x_val = ""
			print "\tParsing: " + fname
			with open(dir+fname, 'r') as fin:
				for line in fin:
					line = line.strip()
					if y_str in line:
						words = line.split(' ')
						y_val = words[len(words)-1]
						count += 1

					if x_str in line:
						words = line.split(' ')
						x_val = words[len(words)-1]
						count += 1

					if count == 2:
						break
				if count == 2:
					data.append([y_val, x_val])

	return [line_name, comment, data]

rmat_scaling_time = []
def get_rmat_time_line(dir):
	print "Get Rmat Time Scaling Line:"
	line = get_line(dir, "Building graph Scale:", "delegate_partitioned_graph time =")
	rmat_scaling_time.append(line)

def gen_rmat_time_tex(dir):
	generate_graph_files(dir+"rmat_time.tex", "RMAT Graph Scaling", "RMAT", "Time(s)", rmat_scaling_time)

rmat_scaling_written = []
def get_rmat_written_line(dir):
	print "Get Rmat Written Scaling Line:"
	line = get_line(dir, "Building graph Scale:", "Total MB Written:")
	rmat_scaling_written.append(line)

def gen_rmat_written_tex(dir):
	generate_graph_files(dir+"rmat_written.tex", "RMAT Graph Scaling", "RMAT", "Data Written (MB)", rmat_scaling_written)


rmat_scaling_read = []
def get_rmat_read_line(dir):
	print "Get Rmat Rad Scaling Line:"
	line = get_line(dir, "Building graph Scale:", "Total MB Read:")
	rmat_scaling_read.append(line)

def gen_rmat_read_tex(dir):
	generate_graph_files(dir+"rmat_read.tex", "RMAT Graph Scaling", "RMAT", "Data Rad (MB)", rmat_scaling_read)


if len(sys.argv) == 1:
	print "No directory specified and no jobs were spawned"
	exit(-1)

for i in xrange(1, len(sys.argv)):
	pass_dir = str(sys.argv[i]) + "/"
	if not os.path.exists(pass_dir):
		"Directory " + pass_dir + " not found, skipping..."
		continue

	print "Processing: " + pass_dir

	generate_csv_file(pass_dir)
	get_rmat_time_line(pass_dir)
	get_rmat_written_line(pass_dir)
	get_rmat_read_line(pass_dir)

tex_dir = "tex_graphs/"
if not os.path.exists(tex_dir):
	os.makedirs(tex_dir)

gen_rmat_time_tex(tex_dir)
gen_rmat_written_tex(tex_dir)
gen_rmat_read_tex(tex_dir)
