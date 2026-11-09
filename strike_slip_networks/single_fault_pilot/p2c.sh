#!/bin/bash
while pgrep -f "p2.sh" > /dev/null; do sleep 10; done
S=$(cd "$(dirname "$0")/.." && pwd)
c="--a 0.01 --b 0.02 --vs_margin 0.1 --a_vs 0.03 --cells_per_lb 4"
# twin of X20_backslip with a 1 Pa kick on the middle cell
$S/run_model X20_backslip_kick 10000 1 $c --dc $(python3 -c "print(0.5515/20)") --load backslip --ic_kick 1 -- -rsf_series_stride 8
# long far-field runs (about 100 cycles)
for X in 20 30; do
  $S/run_model X${X}_farfield_long 40000 1 $c --dc $(python3 -c "print(0.5515/$X)") --load farfield -- -rsf_series_stride 16
done
echo P2C_DONE
