#!/bin/sh

# $1 -> scale
# $2 -> DATA_type
# $3 -> low_deg_threshold
# $4 -> delete ratio

# /g/g90/iwabuchi/havoqgt/build/catalyst.llnl.gov/src/generate_graph_dynamic RMAT $1 16 0 1024 /home/iwabuchi/tmp/test_havoc/tmp/out.graph 0 0 20 $2 $3 $4 | tee  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log

# ls -lsth /home/iwabuchi/tmp/test_havoc/tmp/*  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log

g++ -std=c++11 -I$HOME/apps/include -L$HOME/apps/lib -lrt -pthread -g3 -O0 /home/iwabuchi/tmp/test_havoc/RHH/RHHMain_test.cpp -o /home/iwabuchi/tmp/test_havoc/tmp/RHHMain_test  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log
/home/iwabuchi/tmp/test_havoc/tmp/RHHMain_test  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log

echo "making checkfile..."  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log
g++ -std=c++0x /home/iwabuchi/tmp/test_havoc/work/delete_uniqe_processer.cpp -o /home/iwabuchi/tmp/test_havoc/tmp/delete_uniqe_processer
/home/iwabuchi/tmp/test_havoc/tmp/delete_uniqe_processer < /home/iwabuchi/tmp/test_havoc/tmp/input.txt > /home/iwabuchi/tmp/test_havoc/tmp/tmp_del_uniq.dat  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log

echo "Sorting 1 ..."  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log
sort -k 1,1n -k 2,2n /home/iwabuchi/tmp/test_havoc/tmp/tmp_del_uniq.dat > /home/iwabuchi/tmp/test_havoc/tmp/tmp_sortr1.dat
echo "Sorting 2 ..."  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log
sort -k 1,1n -k 2,2n /home/iwabuchi/tmp/test_havoc/tmp/output.txt > /home/iwabuchi/tmp/test_havoc/tmp/tmp_sortr2.dat

echo "Checking diff..." | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log
diff /home/iwabuchi/tmp/test_havoc/tmp/tmp_sortr1.dat /home/iwabuchi/tmp/test_havoc/tmp/tmp_sortr2.dat | tee -a /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log

echo ""  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log
echo "Check done!!"  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log

echo "-- wc --"  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log
wc /home/iwabuchi/tmp/test_havoc/tmp/tmp_del_uniq.dat  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log
wc /home/iwabuchi/tmp/test_havoc/tmp/output.txt  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log
echo "Fin."  | tee -a  /home/iwabuchi/tmp/test_havoc/tmp/edge_diff_check.log