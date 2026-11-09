#!/usr/bin/env python3
"""
single-patch rsf_solve with two state variables (slip law) against the
nondimensional spring-slider of Gu et al. (1984) / Becker (2000), eqs.
(4)-(6) of Sun, Gilpin & Becker: beta1 = 1, beta2 = 0.84, rho = 0.048.
usage: slider2sv.py rsf_solve_binary kappa_prime [kappa_prime ...]

runs rsf_solve on the 10 km test_slider patch (self-stiffness from
-mat_view, 1.0327845339680668 G/1e4 m) with the parameters mapped onto
the nondimensional system, and integrates the nondimensional ODE with
scipy from the same initial velocity perturbation (VPERT V0); prints
the number of distinct stress extrema (Poincare section x = 0) of both
in the second half of the run
"""
import os
import subprocess
import sys

import numpy as np
from scipy.integrate import solve_ivp

B1, B2, RHO = 1.0, 0.84, 0.048
K_SELF = 1.0327845339680668 * 3.204e10 / 1e4      # Pa/m, 10 km slider patch (-mat_view)
A, SIG, F0, V0 = 0.01, 50e6, 0.6, 1e-9
YR = 365.25 * 86400.0
VPERT = 1.05


def kappa_c():
    k0 = B1 - 1.0
    return (k0 + RHO * (2 * B1 + (B2 - 1) * (2 + RHO)) +
            np.sqrt(4 * RHO**2 * (k0 + B2) + (k0 + RHO**2 * (B2 - 1))**2)) / (2 * (1 + RHO))


def extrema(t, y, tmin):
    """values of local maxima and minima of y after tmin, parabolic refinement"""
    out = []
    for i in range(1, len(y) - 1):
        if t[i] < tmin:
            continue
        if (y[i] - y[i - 1]) * (y[i + 1] - y[i]) < 0:
            # parabola through three points
            t3, y3 = t[i - 1:i + 2], y[i - 1:i + 2]
            c = np.polyfit(t3 - t3[1], y3, 2)
            out.append(c[2] - c[1]**2 / (4 * c[0]) if c[0] != 0 else y[i])
    return np.array(out)


def distinct(v, tol):
    v = np.sort(v)
    if len(v) == 0:
        return 0
    return 1 + int(np.sum(np.diff(v) > tol))


def run_rsf(binary, kp, tmax_tstar, tdir):
    kap = kp * kappa_c()
    dc1 = kap * A * SIG / K_SELF
    tstar = dc1 / V0
    os.makedirs(tdir, exist_ok=True)
    with open(os.path.join(tdir, "rsf.dat"), "w") as f:
        f.write("%g %g\n" % (A, B1 * A))
    with open(os.path.join(tdir, "geom.in"), "w") as f:
        f.write("0 0 -10000 0 90 5000 5000 0\n")
    cmd = [binary, "-geom_file", "geom.in", "-rsf_file", "rsf.dat", "-shear_modulus", "3.204e10",
           "-f0", str(F0), "-dc", "%.10e" % dc1, "-vpl", str(V0), "-v0", str(V0), "-sigma_init", str(SIG),
           "-vel_init", "%.6e" % (VPERT * V0), "-state_law", "2", "-rsf_two_state", "1",
           "-rsf_b2", "%.10e" % (B2 * A), "-rsf_dc2", "%.10e" % (dc1 / RHO),
           "-stop_time_yr", "%.6f" % (tmax_tstar * tstar / YR), "-rtol", "1e-8", "-ts_rk_type", "5dp",
           "-rsf_series_file", "s.bin", "-rsf_series_vdense", "0", "-dt_monitor_yr", "1e6"]
    with open(os.path.join(tdir, "log"), "w") as lg:
        subprocess.run(cmd, cwd=tdir, stdout=lg, stderr=subprocess.STDOUT, check=True)
    with open(os.path.join(tdir, "s.bin"), "rb") as f:
        h = np.fromfile(f, "<i4", 3)
        d = np.fromfile(f, np.dtype([("t", "<f8"), ("f", "<f4", (h[2], h[1]))]))
    t = d["t"] / tstar
    # y = (tau - tau0)/A; tau stored in MPa as float32, so recompute from
    # log10 v, theta1, theta2 would be more precise; float32 tau of ~30 MPa
    # gives about 2 Pa, i.e. 4e-6 in y, sufficient here
    y = (d["f"][:, 1, 0].astype(float) * 1e6 - F0 * SIG) / (A * SIG)
    return t, y, h[2]


def run_ode(kp, tmax):
    kap = kp * kappa_c()

    def rhs(t, u):
        x, y, z = u
        e = np.exp(x)
        f1 = kap * (1 - e)
        f2 = -e * RHO * (B2 * x + z)
        f0 = e * ((B1 - 1) * x + y - z) + f1 - f2
        return [f0, f1, f2]
    s = solve_ivp(rhs, (0, tmax), [np.log(1.05), 0.0, 0.0], method="LSODA", rtol=1e-10, atol=1e-12,
                  max_step=0.5)
    return s.t, s.y[1]


if __name__ == "__main__":
    binary = sys.argv[1]
    tmax = 2e4
    print("kappa_c %.6f" % kappa_c())
    for kp in [float(v) for v in sys.argv[2:]]:
        t, y, nf = run_rsf(binary, kp, tmax, "k%.5f" % kp)
        to, yo = run_ode(kp, tmax)
        er = extrema(t, y, 0.5 * tmax)
        eo = extrema(to, yo, 0.5 * tmax)
        tol = 2e-3
        print("kappa' %.5f: rsf_solve %i records (nfield %i), distinct extrema %i, y range %.4f %.4f | "
              "ODE distinct extrema %i, y range %.4f %.4f" %
              (kp, len(t), nf, distinct(er, tol), er.min(), er.max(), distinct(eo, tol), eo.min(), eo.max()))
        print("   rsf_solve extrema:", np.unique(np.round(er, 3))[:12])
        print("   ODE       extrema:", np.unique(np.round(eo, 3))[:12])
