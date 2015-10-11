#!/bin/sh

mkdir ./work 2> /dev/null
mkdir ./work/debuglog 2> /dev/null
srun --clear-ssd echo "Clear the ssd"
./src/rhhda_bench 21 16 /l/ssd/graph_out 38 1 6 10  2>&1 | tee  /l/ssd/edge_diff_check.log


ls -lsth /l/ssd/*  2>&1 | tee -a  /l/ssd/edge_diff_check.log

echo "making checkfile..."  2>&1 | tee -a  /l/ssd/edge_diff_check.log
/opt/rh/devtoolset-1.1/root/usr/bin/g++ -std=c++11 -O3 ./scripts/delete_uniqe_processer.cpp -o ./work/delete_uniqe_processer  2>&1 | tee -a  /l/ssd/edge_diff_check.log
time ./work/delete_uniqe_processer < /l/ssd/graph_out.debug_edges_raw > /l/ssd/tmp_del_uniq.dat  2>&1 | tee -a  /l/ssd/edge_diff_check.log

echo "Sorting 1 ..."  2>&1 | tee -a  /l/ssd/edge_diff_check.log
sort -k 1,1n -k 2,2n /l/ssd/tmp_del_uniq.dat > /l/ssd/tmp_sortr1.dat
echo "Sorting 2 ..."  2>&1 | tee -a  /l/ssd/edge_diff_check.log
sort -k 1,1n -k 2,2n /l/ssd/graph_out.debug_edges_graph > /l/ssd/tmp_sortr2.dat

echo "Checking diff..." 2>&1 | tee -a  /l/ssd/edge_diff_check.log
diff -Eb /l/ssd/tmp_sortr1.dat /l/ssd/tmp_sortr2.dat > /l/ssd/edge_diff_out.log
echo "# of different edges"  2>&1 | tee -a  /l/ssd/edge_diff_check.log
wc /l/ssd/edge_diff_out.log  2>&1 | tee -a  /l/ssd/edge_diff_check.log

echo ""  2>&1 | tee -a  /l/ssd/edge_diff_check.log
echo "Check done!!"  2>&1 | tee -a  /l/ssd/edge_diff_check.log

echo "-- wc --"  2>&1 | tee -a  /l/ssd/edge_diff_check.log
wc /l/ssd/tmp_del_uniq.dat  2>&1 | tee -a  /l/ssd/edge_diff_check.log
wc /l/ssd/graph_out.debug_edges_graph  2>&1 | tee -a  /l/ssd/edge_diff_check.log
