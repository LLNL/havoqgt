#!/bin/sh


# ---- Configuration ----- #
export EXECUTABLE="rhhda_bench"
export N_NODES=1
export N_PROCS = 1
export DEBUG=True
export USE_DIMMAP=True
export USE_DIMMAP_FOR_TUNE=False
export MONITOR_IO=False
export MEMSIZE_DIMMAP=1024*256*4
export GLOBAL_LOG_FILE="/g/g90/iwabuchi/logs/sc_ipdps_etc.log"
export NORUN=False
export VERBOSE=True
export USE_CATALYST=True
export DELETE_WORK_FILES=False
export SEGMENT_SIZE=39
export EDGES_FILELIST=os.getenv('EFILE_LIST', "./work/file_list")
export TIME_LIMIT=60 * 23
# ---------------------
python run_tests_dynamic.py
