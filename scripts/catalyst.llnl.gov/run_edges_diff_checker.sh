#!/bin/sh

# $1 -> scale
# $2 -> DATA_type
# $3 -> low_deg_threshold
# $4 -> delete ratio

mkdir /g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/work 2> /dev/null
mkdir /g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/work/debuglog 2> /dev/null
# /g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/src/generate_graph_dynamic RMAT $1 16 /l/ssd/out.graph 0 20 $2 $3 $4 2>&1 | tee  /l/ssd/edge_diff_check.log
#337641472
#168820736
export NUM_EDGES=16882073600
/g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/src/generate_graph_dynamic 20 16 /l/ssd/out.graph 0 20 HY_DA 0 0 logs/debug/webgraph_file_list.txt  2>&1 | tee  /l/ssd/edge_diff_check.log
#/g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/src/generate_graph_dynamic 20 16 /dimmap/out.graph 0 20 HY_DA 0 0 logs/debug/webgraph_file_list.txt  2>&1 | tee  /l/ssd/edge_diff_check.log

ls -lsth /l/ssd/*  2>&1 | tee -a  /l/ssd/edge_diff_check.log

echo "making checkfile..."  2>&1 | tee -a  /l/ssd/edge_diff_check.log
/opt/rh/devtoolset-1.1/root/usr/bin/g++ -std=c++11 -O3 /g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/scripts/delete_uniqe_processer.cpp -o /g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/work/delete_uniqe_processer  2>&1 | tee -a  /l/ssd/edge_diff_check.log
time /g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/work/delete_uniqe_processer < /l/ssd/graph_out.debug_edges_raw > /l/ssd/tmp_del_uniq.dat  2>&1 | tee -a  /l/ssd/edge_diff_check.log

echo "Sorting 1 ..."  2>&1 | tee -a  /l/ssd/edge_diff_check.log
sort -k 1,1n -k 2,2n /l/ssd/tmp_del_uniq.dat > /l/ssd/tmp_sortr1.dat
echo "Sorting 2 ..."  2>&1 | tee -a  /l/ssd/edge_diff_check.log
sort -k 1,1n -k 2,2n /l/ssd/graph_out.debug_edges_graph > /l/ssd/tmp_sortr2.dat

echo "Checking diff..." 2>&1 | tee -a  /l/ssd/edge_diff_check.log
diff /l/ssd/tmp_sortr1.dat /l/ssd/tmp_sortr2.dat > /l/ssd/edge_diff_out.log
echo "# of different edges"  2>&1 | tee -a  /l/ssd/edge_diff_check.log
wc /l/ssd/edge_diff_out.log  2>&1 | tee -a  /l/ssd/edge_diff_check.log

echo ""  2>&1 | tee -a  /l/ssd/edge_diff_check.log
echo "Check done!!"  2>&1 | tee -a  /l/ssd/edge_diff_check.log

echo "-- wc --"  2>&1 | tee -a  /l/ssd/edge_diff_check.log
wc /l/ssd/tmp_del_uniq.dat  2>&1 | tee -a  /l/ssd/edge_diff_check.log
wc /l/ssd/graph_out.debug_edges_graph  2>&1 | tee -a  /l/ssd/edge_diff_check.log
