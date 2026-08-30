#!/bin/bash
cd "$(dirname "$0")"
source /disk/rh_usr_local/intel/oneapi/mkl/2025.3/env/vars.sh > setvars.log 2>&1 || true
export LD_LIBRARY_PATH="/disk/rh_usr_local/intel/oneapi/2025.3/lib:/disk/rh_usr_local/intel/oneapi/2024.0/lib:$LD_LIBRARY_PATH"
set -eo pipefail
ulimit -s unlimited
export OMP_STACKSIZE=2G
export OMP_NUM_THREADS=16
mpirun -np 1 ./lhbiem_O perf.in > job.log 2>&1
echo "exit=$?" > exit_code
