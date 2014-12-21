#!/bin/bash

# $1 scale
# $2  Data Type
# $3 highdeg_thshold
# $4 delete_ratio

FNAME=./work/debuglog/`date +%Y%m%d`_`date +%H%M%S`_scale"$1"_"$2"_highdeg"$3"_deleteratio"$4".log
echo $FNAME

sbatch  --clear-ssd -N1 -o $FNAME -e $FNAME << EOF
#!/bin/sh
echo -e "\n\n------------------------------------"
echo Nodes:
echo -e "------------------------------------\n\n"
echo "SLURM_NODELIST = \$SLURM_NODELIST "
echo -e "\n\n------------------------------------"
echo Tuned Info:
echo -e "------------------------------------\n\n"
echo "/proc/sys/vm/dirty_ratio = \$(cat /proc/sys/vm/dirty_ratio)"
echo "/proc/sys/vm/dirty_background_ratio = \$(cat /proc/sys/vm/dirty_background_ratio)"
echo "/proc/sys/vm/dirty_expire_centisecs = \$(cat /proc/sys/vm/dirty_expire_centisecs)"
echo -e "\n\n------------------------------------"
echo free -m
echo -e "------------------------------------\n\n"
free -m
echo -e "\n\n------------------------------------"
echo Top 10 for memory using process
echo -e "------------------------------------\n\n"
ps alx  | awk '{printf ("%d\t%s\n", \$8, \$13)}' | sort -nr | head -10
echo -e "\n\n------------------------------------"
echo df -h /l/ssd
echo -e "------------------------------------\n\n"
df -h -h /l/ssd
echo -e "\n\n------------------------------------"
echo io-stat -m | grep md0 2>&1
echo -e "------------------------------------\n\n"
iostat -m | grep Device 2>&1
iostat -m | grep md0 2>&1
date
echo -e "\n\n------------------------------------"
echo Executable Log
echo -e "------------------------------------\n\n"
srun -N1 -n1 ./scripts/run_edges_diff_check.sh $1 $2 $3 $4
date
echo -e "\n\n------------------------------------"
echo free -m
echo -e "------------------------------------\n\n"
free -m
echo -e "\n\n------------------------------------"
echo df -h /l/ssd
echo -e "------------------------------------\n\n"
df -h /l/ssd
echo -e "\n\n------------------------------------"
echo du -sh /l/ssd/out.graph*
echo -e "------------------------------------\n\n"
du -sh /l/ssd/out.graph*
echo -e "\n\n------------------------------------"
echo io-stat -m | grep md0 2>&1
echo -e "------------------------------------\n\n"
iostat -m | grep Device 2>&1
iostat -m | grep md0 2>&1
echo -e "\n\n------------------------------------"
echo ls -lst /l/ssd/
echo -e "------------------------------------\n\n"
ls -lst /l/ssd/
rm /l/ssd/out.graph*
EOF
