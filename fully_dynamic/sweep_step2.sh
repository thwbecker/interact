#!/bin/sh
# Transition-region sweep: W = 1, 1.5, 2 km, both senses, margins 0.5 and 1 km,
# 2 km overlap. Usage: ./sweep_step2.sh [events] [dx] [mode: qd|fd|both]
events=${1-8}
dx=${2-100}
mode=${3-both}
for vs in 0.5 1.0; do
  for W in -1.0 1.0 -1.5 1.5 -2.0 2.0; do
    if [ "$mode" != "fd" ]; then
      python3 run_step.py --qd --W $W --overlap 2.0 --vs-end $vs --dx $dx --events $events \
        --out step2_qd_W${W}_vs${vs}_dx${dx}.npz > step2_qd_W${W}_vs${vs}_dx${dx}.log 2>&1
    fi
    if [ "$mode" != "qd" ]; then
      python3 run_step.py --W $W --overlap 2.0 --vs-end $vs --dx $dx --events $events \
        --out step2_fd_W${W}_vs${vs}_dx${dx}.npz > step2_fd_W${W}_vs${vs}_dx${dx}.log 2>&1
    fi
  done
done
python3 summarize_step.py "step2_*_dx${dx}.npz"
