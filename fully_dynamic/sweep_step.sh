#!/bin/sh
# Stepover sweep: both senses, several widths, QD then FD. Each configuration
# writes step_<mode>_W<km>_ov<km>_dx<m>.npz after every event, so a killed run
# keeps what it has. Usage: ./sweep_step.sh [events] [dx] [overlap_km]
events=${1-6}
dx=${2-100}
ov=${3-2.0}
for W in -0.5 0.5 -1.0 1.0 -2.0 2.0 -3.0 3.0; do
  python3 run_step.py --qd --W $W --overlap $ov --dx $dx --events $events > step_qd_W${W}_ov${ov}.log 2>&1
done
for W in -0.5 0.5 -1.0 1.0 -2.0 2.0 -3.0 3.0; do
  python3 run_step.py --W $W --overlap $ov --dx $dx --events $events > step_fd_W${W}_ov${ov}.log 2>&1
done
python3 summarize_step.py "step_*_dx${dx}.npz"
