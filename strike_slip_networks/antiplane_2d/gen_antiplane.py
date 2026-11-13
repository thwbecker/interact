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

The --law option switches to the slip law (rsf_solve -state_law 2) or the
PRZ law (-state_law 3) with otherwise identical parameters, or (two) to
two state variables with the slip law: b1 + b2 = b, b2/b1 = 0.84,
dc2 = dc1/0.048, dc1 such that the two-state spring-slider critical
stiffness equals the one-state one (X then refers to the one-state hc or
L_inf); cells G dc2/(b2 sigma)/(8 cpl), which does not resolve the
first variable's G dc1/(b1 sigma). Cells: --cpl cells per L_b = G dc/(b
sigma) (default 4).

usage: gen_antiplane.py model X [--law aging|slip|prz] [--a_over_b r] [--cpl n]
       [--vs_W f] [--depth_W d] [--out dir]

--a_over_b changes a of the velocity-weakening part at fixed b (Wu &
Chen varied a and b over their models A to E; Cattania varied a/b at two
values of b), for the parameter-space maps.
writes geom.in, rsf.dat, ic.in, run.opts, stations.in, params.txt in --out
"""
import argparse
import math
import os

import numpy as np

YR = 365.25 * 86400.0

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("model", choices=("wuchen", "cattania"))
p.add_argument("X", type=float, help="W/hc (wuchen) or W/L_inf (cattania)")
p.add_argument("--law", choices=("aging", "slip", "prz", "two"), default="aging",
               help="state evolution law: aging, slip or PRZ (rsf_solve -state_law 1, 2, 3), or two state "
                    "variables with the slip law (see below)")
p.add_argument("--a_over_b", type=float, default=0.0,
               help="a/b of the velocity-weakening part (b fixed); 0: paper value (0.5 Wu & Chen model A, "
                    "0.75 Cattania)")
p.add_argument("--cpl", type=float, default=4.0, help="cells per L_b")
p.add_argument("--vs_W", type=float, default=1.0, help="cattania: creeping-region extent / W")
p.add_argument("--depth_W", type=float, default=50.0, help="wuchen: burial depth of the fault top / W")
p.add_argument("--seed", type=int, default=1)
p.add_argument("--out", default=".")
a = p.parse_args()

G, cs, f0, sig = 30e9, 3464.0, 0.6, 50e6
if a.model == "wuchen":
    bv, dc = 0.016, 1e-4
    av = (a.a_over_b if a.a_over_b > 0 else 0.5) * bv
    vpl = 0.035 / YR
    v0 = vpl
    hc = 2.0 / math.pi * G * dc / ((bv - av) * sig)
    W = a.X * hc
    ztop = a.depth_W * W
    lscale = hc
else:
    bv, dc = 0.02, 1e-4
    av = (a.a_over_b if a.a_over_b > 0 else 0.75) * bv
    avs = bv + 0.005
    vpl = v0 = 1e-9
    linf = bv * G * dc / (math.pi * (bv - av) ** 2 * sig)
    W = a.X * linf
    ztop = 0.0
    lscale = linf
Lb = G * dc / (bv * sig)
ds = Lb / a.cpl
if a.law == "two":
    # two state variables (slip law) with b1 + b2 = b and b2/b1 = 0.84,
    # dc2 = dc1/0.048 (the ratios of Gu et al. 1984 and Becker 2000);
    # dc1 is chosen so that the critical stiffness of the two-state spring
    # slider, kappa_c a sigma/dc1, equals the one-state (b - a) sigma/dc,
    # so that X keeps the meaning of the one-state W/hc or W/L_inf
    rho = 0.048
    b1, b2 = bv / 1.84, 0.84 * bv / 1.84
    be1, be2 = b1 / av, b2 / av
    k0 = be1 - 1.0
    kc = (k0 + rho * (2 * be1 + (be2 - 1) * (2 + rho)) +
          math.sqrt(4 * rho ** 2 * (k0 + be2) + (k0 + rho ** 2 * (be2 - 1)) ** 2)) / (2 * (1 + rho))
    dc1 = kc * av * dc / (bv - av)
    dc2 = dc1 / rho
    # cells: G dc2/(b2 sigma)/(8 cpl), i.e. /32 for the default cpl 4
    ds = G * dc2 / (b2 * sig) / (8.0 * a.cpl)
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
        fr.write("%.6f %.6f\n" % (aa, b1 if a.law == "two" else bv))
        fi.write("%.10e %.6e\n" % (tau, v))
# stations for rsf_solve -rsf_stations (0-based cell indices): for
# Cattania the loading point (first creeping cell below the
# velocity-weakening part), the deepest velocity-weakening cell and the
# top cell; for Wu & Chen one end and the centre
with open(os.path.join(a.out, "stations.in"), "w") as f:
    if a.model == "cattania":
        f.write("load %i\nvwbot %i\ntop 0\n" % (nvw, nvw - 1))
    else:
        f.write("end 0\nmid %i\n" % (nvw // 2))
# rough recurrence estimate for the run length: stress drop of a few
# (b-a) sigma over the fault, slip ~ dtau W/G, at vpl
t_est = 20.0 * (bv - av) * sig * W / G / vpl / YR
opts = ("-geom_file geom.in -rsf_file rsf.dat -rsf_ic_file ic.in -dc %.6e -sigma_init %.6e -f0 %g "
        "-v0 %.6e -vpl %.6e -shear_modulus %.6e -s_wave_speed %g"
        % (dc1 if a.law == "two" else dc, sig, f0, v0, vpl, G, cs))
if a.law == "two":
    opts += " -state_law 2 -rsf_two_state 1 -rsf_b2 %.8e -rsf_dc2 %.8e" % (b2, dc2)
if a.law == "slip":
    opts += " -state_law 2"
elif a.law == "prz":
    opts += " -state_law 3"
with open(os.path.join(a.out, "run.opts"), "w") as f:
    f.write(opts + "\n")
with open(os.path.join(a.out, "params.txt"), "w") as f:
    f.write("model %s law %s X %g W %.6e ds %.6e nvw %i nvs %i Lb %.6e lscale %.6e a %g b %g dc %g "
            "vpl %.6e t_est_yr %.6e\n" % (a.model, a.law, a.X, W, ds, nvw, nvs, Lb, lscale, av, bv, dc, vpl, t_est))
print("%s %s X %g: W %.1f m, %i VW + %i VS cells of %.3f m (L_b/%.1f), recurrence estimate %.3g yr"
      % (a.model, a.law, a.X, W, nvw, nvs, ds, Lb / ds, t_est))
