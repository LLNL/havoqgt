#!/usr/bin/python

from operator import itemgetter
import sys

is_dupication = 1

f = open(sys.argv[1])
edges = []

for line in f:
	items = line.split()
	edges.append((int(items[0]), int(items[1])))
f.close()
print "Reading file1 done."

edges.sort(key=lambda x: x[1])
print "Sort1-1 done."
edges.sort(key=lambda x: x[0])
print "Sort1-2 done."

#print edges

f = open(sys.argv[2])
edges2 = []

for line in f:
	items = line.split()
	edges2.append((int(items[0]), int(items[1])))
f.close()
print "Reading file2 done"

edges2.sort(key=lambda x: x[1])
print "Sort2-1 done."
edges2.sort(key=lambda x: x[0])
print "Sort2-2 done."

old_edge_src = -1;
old_edge_trg = -1;
i2 = 0;

diff_count = 0;
for i, item in enumerate(edges):
	if is_dupication and old_edge_src == item[0] and old_edge_trg == item[1]:
		continue

	if item[0] != edges2[i2][0] or item[1] != edges2[i2][1]:
		print "(", item[0], item[1], ") != (", edges2[i2][0], edges2[i2][1], ")"
		diff_count = diff_count + 1

	old_edge_src = item[0]
	old_edge_trg = item[1]
	
	#print item[0], item[1], edges2[i2][0], edges2[i2][1]

	i2 = i2 + 1

if diff_count == 0:
	print "Success !!"
else:
	print "# of diff =", diff_count
