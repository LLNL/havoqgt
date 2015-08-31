#!/usr/bin/python

from operator import itemgetter
import sys
import os
import re

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.pyplot as plt

#------------------------
num_files = 1
#keyword = "segment size"
keyword = "Execution time"
col_num = 5 # NOTE: [ awk -1 ]
#is_cumulate = False
is_cumulate = True
#------------------------


if len(sys.argv) != 2:
	print "# of arguments is wrong"
	sys.exit(1)

target_dir = sys.argv[1]
# num_files = int(sys.argv[2])
# keyword = sys.argv[3]
# col_num = int(sys.argv[4]) # NOTE: [ awk -1 ]
# is_cumulate = bool(int(sys.argv[5]))

flist = []

for i in range(num_files):
	fname = target_dir + "/test_" + str(i) + ".out"
	flist.append(open(fname))

results = []
counts = []

for i in range(num_files):
	result = []
	for line in flist[i]:
		if re.search(keyword, line):
			items =	line.split()
			#print items
			if len(items) > col_num:
				result.append(float(items[-1]))
	results.append(result)
	counts.append(len(result))

max_length = max(counts)

if is_cumulate:
	for i in range (1, max_length):
		for j in range (num_files):
			if i < len(results[j]):
				results[j][i] = results[j][i-1] + results[j][i]


# --- print result --- #
fnameresults = target_dir + "/value_" + str(num_files) + "_" + keyword.replace(" ", "_") + "_" + str(col_num) + "_" + str(is_cumulate) + ".txt"
fresult = open(fnameresults, "w")
for i in range (max_length):
	for j in range (num_files):
		if i < len(results[j]):
			fresult.write(str(results[j][i]))
		fresult.write("\t")
	fresult.write("\n")


# --- plot result --- #
x_value = range(1, max_length+1)
for i in range (num_files):
	x_value = range(1, len(results[i])+1)
	plt.plot(x_value, results[i], label="test_" + str(i) + ".out")
plt.legend()

x1,x2,y1,y2 = plt.axis()
plt.axis((0, x2, 0, y2))

plt.ylabel(keyword)
plt.xlabel('inserted chunks')

plt.grid()

outname = target_dir + "/chart_" + str(num_files) + "_" + keyword.replace(" ", "_") + "_" + str(col_num) + "_" + str(is_cumulate) + ".pdf"
plt.savefig(outname ,format = 'pdf')
