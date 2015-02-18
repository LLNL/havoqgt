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

make_figure_separately = True
target_dir_list = ["./20150212_181159/", "./20150212_181356/", "./20150212_182833/", "./20150212_182937/"]
legend_labels_list = ["normal-mmap, without-sort", "normal-mmap, sort-chunk", "di-mmap, without-sort", "di-mmap, sort-chunk"]

target_file = "io_monitering_report.log"
column_labels = ["rrqm/s", "wrqm/s", "r/s", "w/s", "rMB/s", "wMB/s", "avgrq-sz", "avgqu-sz", "await", "svctm", "%%util"]
device_name = "md0"

def read_data(_target_dir_list):
  global target_files_list
  global values_root

  target_files_list = []
  for tdir in _target_dir_list:
    fname = tdir + target_file
    target_files_list.append(fname)

  # etract results from all targert files
  values_root =  {}
  for label in column_labels:
    columnwise_values_all_filse = {}
    for f in target_files_list:
      columnwise_values_all_filse[f] = []
    values_root[label] = columnwise_values_all_filse

  for f in target_files_list: # for all files

    columnwise_values_for_a_file = {}
    for label in column_labels:
      columnwise_values_for_a_file[label] = []

    for line in open(f, "r"):
      if re.search(device_name, line):
        rowwise_values = line.split()
        rowwise_values = rowwise_values[1:] # Delete first element since it is a device name
        for i, val in enumerate(rowwise_values):
          # if val == device_name:
          #   continue
          columun = columnwise_values_for_a_file[column_labels[i]]
          columun.append(float(val))

    for colum_label in columnwise_values_for_a_file:
      columnwise_values_all_filse = values_root[colum_label]
      columnwise_values_all_filse[f] = columnwise_values_for_a_file[colum_label]
      #values_root[column_labels[i]] = columnwise_values_all_filse


# --- plot result --- #
def plot_result(_target_dir_list, _index):

  if make_figure_separately:
    f_log = open(_target_dir_list[0]+"iostat_result.log", "w")
    for ylabel in column_labels:
      column_for_a_file = values_root[ylabel][target_files_list[0]]
      f_log.write(ylabel + " ")
      for d in column_for_a_file:
        f_log.write(str(d) + " ")
      f_log.write("\n")
    f_log.close()

  for ylabel in column_labels:
    columnwise_values_all_filse = values_root[ylabel]
    f_pos = _index
    for f in target_files_list:
      column_for_a_file = columnwise_values_all_filse[f]
      x_value = range(0, len(column_for_a_file)*10, 10)
      #legend_label = _target_dir_list[i].split('/')[-2]
      legend_label = legend_labels_list[_index]
      column_for_a_file[0] = 0
      plt.plot(x_value, column_for_a_file, label=legend_label)
      plt.legend()

      # ------------------------------ #
      x1,x2,y1,y2 = plt.axis()
      if ylabel == "r/s" or ylabel == "w/s":
        y2 = 80000
      elif  ylabel == "rMB/s" or ylabel == "wMB/s":
        y2 = 1200
      x2 = 250000
      # ------------------------------ #

      plt.axis((0, x2, 0, y2))

      plt.title("iostat-" + ylabel)
      plt.ylabel(ylabel)
      plt.xlabel("Time (sec.)")

      plt.grid()

      tmp = ylabel.replace("/", "")
      tmp = tmp.replace("%%", "")
      outname = _target_dir_list[0] + tmp + ".png"
      plt.savefig(outname, format = 'png')

      f_pos = f_pos + 1
    plt.clf()


def main():
  if make_figure_separately:
    for i, t_d in enumerate(target_dir_list):
      target_dir = []
      target_dir.append(t_d)
      read_data(target_dir)
      plot_result(target_dir, i)
  else:
    read_data(target_dir_list)
    plot_result(target_dir_list, 0)

if __name__ == "__main__":
    main()