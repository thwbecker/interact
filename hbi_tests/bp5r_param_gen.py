#!/usr/bin/env python3
"""
bp5r_param_gen.py -- generate HBI's per-element BP5 (rectangular) parameter file
at any resolution. Reproduces the exact frictional distribution of HBI's
examples/bp5r_param.py (a, dc, tau, vel; VW core + VS borders + nucleation patch),
but parametrized by resolution and without the matplotlib/stray-f0 issues so it
runs headless.

Element ordering is i (along-strike) outer, j (down-dip) inner -- matching HBI's
internal numbering and the shipped bp5r_param.dat. Physical coordinates (x, dep in
km) are resolution-independent, so finer ds just samples the same 100 x 40 km fault
more densely (a proper convergence refinement).

Usage:
    bp5r_param_gen.py IMAX JMAX DS  OUTFILE
    e.g.  bp5r_param_gen.py 100 40 1.0  bp5r_1km_param.dat     # N=4000  (1 km)
          bp5r_param_gen.py 200 80 0.5  bp5r_0.5km_param.dat   # N=16000 (0.5 km)
"""
import sys, numpy as np

if len(sys.argv) != 5:
    sys.exit(__doc__)
imax, jmax = int(sys.argv[1]), int(sys.argv[2])
ds0 = float(sys.argv[3])
outfile = sys.argv[4]

# --- BP5 frictional constants (identical to HBI examples/bp5r_param.py) ---
a0, b0, dc0 = 0.004, 0.03, 0.14
a_max = 0.04
mu0, vref = 0.6, 1e-6
sigma0, rigid, cs = 25.0, 32.04, 3.464
vel0 = 1e-9

a   = np.empty(imax * jmax)
dc  = np.empty(imax * jmax)
tau = np.empty(imax * jmax)
vel = np.empty(imax * jmax)

k = -1
for i in range(imax):
    for j in range(jmax):
        k += 1
        x   = (i + 0.5 - imax / 2) * ds0      # along-strike, km  (-50..+50)
        dep = (j + 0.5) * ds0                 # depth, km          (0..40)
        dc[k] = dc0
        r = max(abs(dep - 10) - 6, abs(x) - 30) / 2
        a[k] = min(a0 + r * (a_max - a0), a_max)
        a[k] = max(a[k], a0)
        vel[k] = vel0
        if abs(x + 24) < 6 and abs(dep - 10) < 6:   # nucleation patch
            vel[k] = 3e-2
            dc[k]  = 0.13
        tau[k] = sigma0 * (a[k] * np.arcsinh(0.5 * vel[k] / vref *
                 np.exp((mu0 + b0 * np.log(vref / vel0)) / a[k]))
                 + rigid / (2 * cs) * vel[k])

with open(outfile, "w") as f:
    f.write("a\tdc\ttau\tvel\n")
    for k in range(imax * jmax):
        f.write(f"{a[k]:.18g}\t{dc[k]:.18g}\t{tau[k]:.18g}\t{vel[k]:.18g}\n")

print(f"wrote {outfile}: {imax*jmax} rows (imax={imax} jmax={jmax} ds={ds0} km)")
