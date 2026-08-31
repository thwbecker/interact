#!/usr/bin/env python3
"""T6: Maxwell-time scaling symmetry of the relaxation kernel.

The substrate enters the wavenumber-domain system only through
mu2(s) = mu2 s/(s + 1/tau) and lam2(s) = K2 - 2 mu2(s)/3, and the
buoyancy rows carry no s, so substituting s = sigma/tau leaves the
system dependent on s only through s*tau.  The step response is then
a function of t/tau alone,

    dK(t; tau) = dKtilde(t/tau),      dK(t) = K(t) - K(0+)

which is what bp3_ve_kernels.py exploits to sample once and emit a
kernel file per Maxwell time by scaling the ladder (checked at file
level by run_tau_scaling_test).  This gate tests the underlying
statement directly on the sampled arrays, for both traction families
and with gravity on and off, at integer and non-integer Maxwell-time
ratios.  Note that the relaxation part dK is the scale-invariant
object: K(0+) itself is elastic and Maxwell-time independent, so the
same holds for K(t) - K(0+) however K(0+) is carried."""
import numpy as np
from inplane2d import fields_pairs, ve_fields_pairs, segment_sources

yr = 3.15576e7
mu = lam = 32.04e9                     # BP3 elastic parameters
K2 = lam + 2*mu/3
rho1, rho2, gon = 2800.0, 3300.0, 9.81
H = 60e3
kw = dict(kmin=1e-7, kmax=8.0/H, npanel=24, nquad=4)
mur = 1e-10*mu
tol = 1e-9                             # relative to the relaxed limit
ok = True

# a dipping (60 deg) contained fault, as in the BP3 cross section:
# receivers at the segment centers, tractions resolved the way
# bp3_ve_kernels.py resolves them
ds, dip = 5e3, np.deg2rad(60.0)
d = (np.arange(4) + 3.5)*ds                    # 4 receiver segments
xc, zc = d*np.cos(dip), d*np.sin(dip)          # z down
tvec = np.array([np.cos(dip), np.sin(dip)])    # down-dip tangent
nvec = np.array([-tvec[1], tvec[0]])

def tractions(F):
    """(shear, normal) tractions on the receiver segments"""
    tx_ = F["sxx"]*nvec[0] + F["sxz"]*nvec[1]
    tz_ = F["sxz"]*nvec[0] + F["szz"]*nvec[1]
    return (tx_*tvec[0] + tz_*tvec[1], tx_*nvec[0] + tz_*nvec[1])

j = 1                                          # the source segment
src = segment_sources(xc[j] - ds/2*tvec[0], zc[j] - ds/2*tvec[1],
                      xc[j] + ds/2*tvec[0], zc[j] + ds/2*tvec[1],
                      1.0, ns=2)
that = np.array([0.1, 1.0, 10.0])              # sampling times t/tau

def sampled(tM, gg):
    """relaxation part of the two traction families at t = that*tM"""
    el = fields_pairs(xc, zc, src, lam, mu, lam, mu, H,
                      rho1=rho1, rho2=rho2, g=gg, **kw)
    s0, n0 = tractions(el)
    elr = fields_pairs(xc, zc, src, lam, mu, K2 - 2*mur/3, mur, H,
                       rho1=rho1, rho2=rho2, g=gg, **kw)
    sinf, ninf = tractions(elr)
    dS = np.empty((that.size, xc.size))
    dN = np.empty((that.size, xc.size))
    for k, th in enumerate(that):
        ve = ve_fields_pairs(xc, zc, src, lam, mu, lam, mu, H, tM,
                             th*tM, M=10, rho1=rho1, rho2=rho2,
                             g=gg, **kw)
        vr = {c: np.real(ve[c]) for c in ("sxx", "sxz", "szz")}
        s, nn = tractions(vr)
        dS[k] = s - s0
        dN[k] = nn - n0
    return dS, dN, sinf - s0, ninf - n0

for gg, gl in ((gon, "g on "), (0.0, "g off")):
    ratios = (3.0, 7.3) if gg else (3.0,)
    dS0, dN0, dSi, dNi = sampled(15.0*yr, gg)
    sS = np.max(np.abs(dSi)); sN = np.max(np.abs(dNi))
    for r in ratios:
        dS, dN, dSi2, dNi2 = sampled(r*15.0*yr, gg)
        for lab, a, b, s in (("shear ", dS, dS0, sS),
                             ("normal", dN, dN0, sN)):
            dev = float(np.max(np.abs(a - b)))/s
            print(f"T6 {gl} tau x {r:4.1f} {lab}: max |dK(t/tau) dev| "
                  f"= {dev:8.2e} of {s:.3e} Pa/m")
            ok = ok and dev < tol
        # the relaxed limit is elastic and must not move at all
        dev = float(np.max(np.abs(dSi2 - dSi)))/sS
        print(f"T6 {gl} tau x {r:4.1f} relaxed limit: dev "
              f"= {dev:8.2e}")
        ok = ok and dev < 1e-13
print("T6 " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
