import os
import glob
import sys
import thread
from urlparse import urlparse
from multiprocessing.dummy import Pool as ThreadPool 
from multiprocessing import Pool

def build_file_list(file_list, skip_filename) : 
	for filename in glob.glob('*'):
		if filename == skip_filename:
			continue
		file_list.append(filename)
"""			
def to_undir_edgelist(filename, output_filepath) :
	# filename example part-r-00063
	file_number = int(filename.split("-")[2])
	out_file_number = file_number + len(dir_file_list)
	output_filename = out_filepath + "/part-r-000" + str(out_file_number)
 	print(" > Reading from " + filename + ", Writing to " + output_filename + " ... ")

	line_count = 0

	with open(filename, 'r') as readfile, open(output_filename, 'w') as writefile:
		for line in readfile:
			line_count = line_count + 1
			source = line.split(" ")[0].strip()
			target = line.split(" ")[1].strip()
			writefile.write(target + " " + source + "\n")
			#if line_count == 10:
			#	break
		print(" >> wrote " + str(line_count) + " lines to " + output_filename + ".")

def thread_function(in_filename) :
	to_undir_edgelist(in_filename, out_filepath)
"""

def active_vertices_count(dir_file_list, out_filename) :
	print ("Total active vertices count")

	# read one file and count number of iterations
	print ("Counting total number of iterations ... ")
	max_itr_count = 0
	with open (dir_file_list[0], 'r') as readfile:
		#print ("> Reading file: " + dir_file_list[0] + " ... ") 
		for line in readfile:
			#print (line)
			#line = line.strip()	
			max_itr_count = max_itr_count + 1;
		print ("Total number of iterations: " + str(max_itr_count))
	
	active_vertices_count_string_list = max_itr_count*[""]
	active_vertices_count_list = max_itr_count*[0]
	#print ("active_vertices_count_list size: " + str(len(active_vertices_count_list)))

	line_count = 0
	with open (dir_file_list[0], 'r') as readfile:
		for line in readfile:
			line = line.strip()
			tokens = line.split(",")
			tokens.pop(-1) # remove the last element
			for t in tokens:
				active_vertices_count_string_list[line_count] = active_vertices_count_string_list[line_count] + t.strip() + "," 
			line_count = line_count + 1 

	# read all files and accumulate results	
	file_count = 0
	for f in dir_file_list:			
		with open (f, 'r') as readfile:
			#print ("> [" + str(file_count) + "] Reading file: " + f + " ... ")
			line_count = 0
			for line in readfile:

				line = line.strip()
				tokens = line.split(",")							

				active_vertices_count_list[line_count] = active_vertices_count_list[line_count] + int(tokens[len(tokens) - 1]); 

				#print (">> Iteration: " + str(line_count) + " - " + line)

				line_count = line_count + 1
		file_count = file_count + 1
		#if file_count == 100:
		#	break

	# write output to file
	# TODO: write output to file
	# out_filename
	for i in range(0, len(active_vertices_count_list)) :
		print (active_vertices_count_string_list[i] + str(active_vertices_count_list[i]))	
	print ("Done.")

#in_filename = sys.argv[1]
#out_filename = sys.argv[2]

filepath = sys.argv[1]
skip_filename = "" #sys.argv[2]
#out_filepath = sys.argv[2] # gloabal variable 
out_filename = "" #sys.argv[2]
#thread_count = int(sys.argv[4])

dir_file_list = []

os.chdir(filepath)

build_file_list(dir_file_list, skip_filename)

#for f in dir_file_list:
#	print(f) 	
print(str(len(dir_file_list)) + " files to process ... ")

active_vertices_count(dir_file_list, out_filename)

"""
#pool = ThreadPool(thread_count)
pool = Pool(processes=thread_count)
print("Created pool with " + str(thread_count) + " processes.")

print("Creating undirected edgelist file(s) ... ")

results = pool.map(thread_function, dir_file_list)
pool.close()
pool.join()

print("Done creating undirected edgelist file(s).")	
"""
