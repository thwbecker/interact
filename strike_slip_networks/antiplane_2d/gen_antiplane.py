#!/usr/bin/env python3
"""
gen_antiplane.py: 2-D antiplane (mode III, out-of-plane slip) single-fault
setups for rsf_solve that follow two published parameter studies, as a
check of the 2.5-D strip results against the 1-D antiplane literature.

  wuchen    Wu & Chen (2014, JGR 119, 4890), model A: uniform fault of
            length W, aging law, a 0.008, b 0.016, dc 0.1 mm, sigma 50 MPa,
            G 30 GPa, cs 3464 m/s, f0 0.6, v0 = vpl = 0.035 m/yr; the fault
            is "surrounded by a steady-slipping region", which in the
            backslip formulation means backslip at vpl on the fault cells
            only. Their whole space is approximated by burying the fault
            (top at --depth_W fault lengths below the free surface).
            X = W/hc with hc = (2/pi) G dc/((b-a) sigma) = 4.77 m; their
            regime boundaries are at W/h* = 1.2, 5.5, 33, 105 with
            h* = hc/(1 - a/b) = 2 hc here.
  cattania  Cattania (2019, GRL 46, 10384): vertical fault reaching the
            free surface, velocity weakening from the surface to depth W
            (a - b = -0.005), velocity strengthening below (a - b =
            +0.005) to depth (1 + --vs_W) W, b 0.02, sigma 50 MPa, dc 0.1
            mm, aging law, loaded by creep below (backslip at vpl on all
            cells). X = W/L_inf with L_inf = b G dc/(pi (b-a)^2 sigma).
            G, cs, vpl, v0 and the extent of the creeping region are not
            stated in the paper; here G 30 GPa, cs 3464 m/s, vpl = v0 =
            1e-9 m/s, f0 0.6, and a creeping region of --vs_W W.

The --law option switches to the slip law (rsf_solve -state_law 2) with
otherwise identical parameters. Cells: --cpl cells per L_b = G dc/(b
sigma) (default 4).

usage: gen_antiplane.py model X [--law aging|slip] [--cpl n] [--vs_W f]
       [--depth_W d] [--out dir]
writes geom.in, rsf.dat, ic.in, run.opts, params.txt in --out
"""
import argparse
import math
import os

import numpy as np

YR = 365.25 * 86400.0

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("model", choices=("wuchen", "cattania"))
p.add_argument("X", type=float, help="W/hc (wuchen) or W/L_inf (cattania)")
p.add_argument("--law", choices=("aging", "slip"), default="aging")
p.add_argument("--cpl", type=float, default=4.0, help="cells per L_b")
p.add_argument("--vs_W", type=float, default=1.0, help="cattania: creeping-region extent / W")
p.add_argument("--depth_W", type=float, default=50.0, help="wuchen: burial depth of the fault top / W")
p.add_argument("--seed", type=int, default=1)
p.add_argument("--out", default=".")
a = p.parse_args()

G, cs, f0, sig = 30e9, 3464.0, 0.6, 50e6
if a.model == "wuchen":
    av, bv, dc = 0.008, 0.016, 1e-4
    vpl = 0.035 / YR
    v0 = vpl
    hc = 2.0 / math.pi * G * dc / ((bv - av) * sig)
    W = a.X * hc
    ztop = a.depth_W * W
    lscale = hc
else:
    bv, dc = 0.02, 1e-4
    av = bv - 0.005
    avs = bv + 0.005
    vpl = v0 = 1e-9
    linf = bv * G * dc / (math.pi * (bv - av) ** 2 * sig)
    W = a.X * linf
    ztop = 0.0
    lscale = linf
Lb = G * dc / (bv * sig)
ds = Lb / a.cpl
nvw = max(2, int(math.ceil(W / ds)))
ds = W / nvw
nvs = 0 if a.model == "wuchen" else int(math.ceil(a.vs_W * W / ds))
eta = G / (2.0 * cs)
rng = np.random.default_rng(a.seed)
os.makedirs(a.out, exist_ok=True)
with open(os.path.join(a.out, "geom.in"), "w") as fg, \
     open(os.path.join(a.out, "rsf.dat"), "w") as fr, \
     open(os.path.join(a.out, "ic.in"), "w") as fi:
    for j in range(nvw + nvs):
        z = ztop + (j + 0.5) * ds
        aa = av if j < nvw else avs
        # steady state at a slightly perturbed plate rate (breaks the
        # unstable uniform equilibrium; 0.1 percent random)
        v = vpl * (1.0 + 1e-3 * rng.uniform(-1, 1))
        tau = sig * aa * math.asinh(v / (2 * v0) * math.exp((f0 + bv * math.log(v0 / v)) / aa)) + eta * v
        fg.write("0 %.8e 0 0 -90 %.8e 0 0\n" % (-z, ds / 2.0))
        fr.write("%.6f %.6f\n" % (aa, bv))
        fi.write("%.10e %.6e\n" % (tau, v))
# rough recurrence estimate for the run length: stress drop of a few
# (b-a) sigma over the fault, slip ~ dtau W/G, at vpl
t_est = 20.0 * (bv - av) * sig * W / G / vpl / YR
opts = ("-geom_file geom.in -rsf_file rsf.dat -rsf_ic_file ic.in -dc %.6e -sigma_init %.6e -f0 %g "
        "-v0 %.6e -vpl %.6e -shear_modulus %.6e -s_wave_speed %g"
        % (dc, sig, f0, v0, vpl, G, cs))
if a.law == "slip":
    opts += " -state_law 2"
with open(os.path.join(a.out, "run.opts"), "w") as f:
    f.write(opts + "\n")
with open(os.path.join(a.out, "params.txt"), "w") as f:
    f.write("model %s law %s X %g W %.6e ds %.6e nvw %i nvs %i Lb %.6e lscale %.6e a %g b %g dc %g "
            "vpl %.6e t_est_yr %.6e\n" % (a.model, a.law, a.X, W, ds, nvw, nvs, Lb, lscale, av, bv, dc, vpl, t_est))
print("%s %s X %g: W %.1f m, %i VW + %i VS cells of %.3f m (L_b/%.1f), recurrence estimate %.3g yr"
      % (a.model, a.law, a.X, W, nvw, nvs, ds, Lb / ds, t_est))
