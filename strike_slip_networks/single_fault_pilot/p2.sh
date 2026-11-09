#!/bin/bash
# P2 pilot: single strip, aging law, a 0.01 b 0.02, 10 % VS ends, L 60 km,
# L/h* set by dc, far-field against margin loading
S=$(cd "$(dirname "$0")/.." && pwd)
tmax=${1-10000}
for X in 5 10 20 30 40; do
  dc=$(python3 -c "print(0.5515/$X)")
  for load in backslip farfield; do
    if grep -q "^wall" X${X}_$load/rsf_solve.log 2>/dev/null; then continue; fi
    $S/run_model X${X}_$load $tmax 1 --a 0.01 --b 0.02 --dc $dc --vs_margin 0.1 --a_vs 0.03 --load $load --cells_per_lb 4 -- -rsf_series_stride 8
  done
done
echo P2_DONE
