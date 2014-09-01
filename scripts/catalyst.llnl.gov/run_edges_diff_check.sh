#!/bin/sh

/g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/src/generate_graph_dynamic RMAT $1 16 0 1024 /l/ssd/out.graph 0 0 20 RB_MX 1 $2 | tee  /l/ssd/edge_diff_check.log

ls -lsth /l/ssd/*  | tee -a  /l/ssd/edge_diff_check.log

echo "making checkfile..."  | tee -a  /l/ssd/edge_diff_check.log
g++ -std=c++0x /g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/scripts/delete_uniqe_processer.cpp -o /g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/work/delete_uniqe_processer
/g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/work/delete_uniqe_processer < /l/ssd/graph_out.debug_edges_raw > /l/ssd/tmp_del_uniq.dat  | tee -a  /l/ssd/edge_diff_check.log

echo "Sorting 1 ..."  | tee -a  /l/ssd/edge_diff_check.log
sort -k 1,1n -k 2,2n /l/ssd/tmp_del_uniq.dat > /l/ssd/tmp_sortr1.dat
echo "Sorting 2 ..."  | tee -a  /l/ssd/edge_diff_check.log
sort -k 1,1n -k 2,2n /l/ssd/graph_out.debug_edges_graph > /l/ssd/tmp_sortr2.dat

echo "Checking diff..." | tee -a  /l/ssd/edge_diff_check.log
diff /l/ssd/tmp_sortr1.dat /l/ssd/tmp_sortr2.dat | tee -a /l/ssd/edge_diff_check.log

echo ""  | tee -a  /l/ssd/edge_diff_check.log
echo "Check done!!"  | tee -a  /l/ssd/edge_diff_check.log

echo "-- wc --"  | tee -a  /l/ssd/edge_diff_check.log
wc /l/ssd/tmp_del_uniq.dat  | tee -a  /l/ssd/edge_diff_check.log
wc /l/ssd/graph_out.debug_edges_graph  | tee -a  /l/ssd/edge_diff_check.log
echo "Fin."  | tee -a  /l/ssd/edge_diff_check.log