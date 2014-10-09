#!/usr/bin/python

from operator import itemgetter
import sys
import bisect


f = open(sys.argv[1])
elems = {}

for line in f:
  items = line.split()
  if int(items[1]) == 1:
    if elems.get(int(items[0])) == 1:
      print "!!!!"
    else:
      elems[int(items[0])] = 1
  elif int(items[1]) == 2:
      elems[int(items[0])] = 0

f.close()