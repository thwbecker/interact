#!/usr/bin/env python3
"""
bifurcation_scan.py rsf_solve_binary outfile [vpert] [ttotal] [tsettle]:
rsf_solve two-state single-patch slider (see slider2sv.py) over a list of
kappa' values; writes, per kappa', the stress extrema y = (tau-tau0)/(a
sigma) after tsettle (units of dc/v0) for a bifurcation diagram, and a
summary line (number of distinct extrema at a tolerance of 3e-3 in y,
i.e. 3e-5 in tau/sigma for a = 0.01, permutation entropy of the maxima)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "strike_slip_networks"))
import slider2sv as s  # noqa: E402
from complexity import perm_entropy, n_distinct  # noqa: E402

KAPPAS = [0.95, 0.92, 0.90, 0.88, 0.87, 0.865, 0.86, 0.858, 0.8577, 0.8574, 0.857, 0.856, 0.8555,
          0.8552, 0.855, 0.8545, 0.854, 0.8536, 0.8534, 0.8533, 0.853, 0.852, 0.851, 0.85]

if __name__ == "__main__":
    binary, out = sys.argv[1], sys.argv[2]
    s.VPERT = float(sys.argv[3]) if len(sys.argv) > 3 else 1.001
    ttot = float(sys.argv[4]) if len(sys.argv) > 4 else 1e5
    tset = float(sys.argv[5]) if len(sys.argv) > 5 else 5e4
    with open(out, "w") as f:
        f.write("# kappa' kind(+1 max,-1 min) y   [vpert %g, total %g, settle %g dc/v0]\n" % (s.VPERT, ttot, tset))
    for kp in KAPPAS:
        d = "k%.5f" % kp
        t, y, nf = s.run_rsf(binary, kp, ttot, d)
        w = t >= tset
        tw, yw = t[w], y[w]
        mx, mn = [], []
        for i in range(1, len(yw) - 1):
            if (yw[i] - yw[i - 1]) * (yw[i + 1] - yw[i]) < 0:
                c = np.polyfit(tw[i - 1:i + 2] - tw[i], yw[i - 1:i + 2], 2)
                v = c[2] - c[1] ** 2 / (4 * c[0]) if c[0] != 0 else yw[i]
                (mx if c[0] < 0 else mn).append(v)
        with open(out, "a") as f:
            for v in mx:
                f.write("%.6f 1 %.6f\n" % (kp, v))
            for v in mn:
                f.write("%.6f -1 %.6f\n" % (kp, v))
        allv = np.array(mx + mn)
        print("kappa' %.5f: %i extrema, distinct %i, y %.3f to %.3f, PE(maxima) %.2f" %
              (kp, len(allv), n_distinct(allv, 3e-3), allv.min(), allv.max(), perm_entropy(mx, 6, 3e-3)),
              flush=True)
        os.remove(os.path.join(d, "s.bin"))
