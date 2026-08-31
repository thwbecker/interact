#!/usr/bin/env python3
"""T7: the batched matched-pair evaluation path against the reference.

fields_pairs() sums the wavenumber quadrature as a matrix product over
a separable form of the region bases (see the BATCHED EVALUATION PATH
comment in inplane2d.py).  Correctness rests on three things, each
gated here:

 (a) the separation itself: B(z) = sum_m f_m(z) C_m must reproduce
     basis() and basis_cs() exactly, for real and for Laplace-domain
     (complex) moduli;
 (b) the batched matcher solve_k_b must solve the same 10x10 systems as
     solve_k.  The equilibrated systems reach condition numbers of
     1e7 to 1e9 here, so two backward-stable solutions of the same
     (differently rounded) matrix may differ by much more than
     machine epsilon; the gate is therefore on the RESIDUAL of both
     solutions, with the coefficient difference reported;
 (c) the assembled fields must agree with fields_pairs_ref, which is
     the original per-wavenumber loop, kept for exactly this purpose;
     and with the diagonal of fields_xz, which the batched path does
     not touch and which is therefore an independent reference."""
import numpy as np
from inplane2d import (basis, basis_cs, region_basis, _split_cs,
                       _split_exp, jumps,
                       jumps_b, solve_k, solve_k_b, kgrid, fields_xz,
                       fields_pairs, fields_pairs_ref, maxwell_moduli,
                       segment_sources)

yr = 3.15576e7
mu = lam = 32.04e9
H = 40e3
rho1, rho2, gon = 2800.0, 3300.0, 9.81
ok = True
rng = np.random.default_rng(11)

# ---------------------------------------------------------------- (a)
wcs = wex = 0.0
for trial in range(200):
    k = 10**rng.uniform(-7, -2)
    lm, mm = lam*(1 + 0.3*rng.normal()), mu*(1 + 0.3*rng.normal())
    if trial % 2:
        lm, mm = lm*(1 - 0.2j), mm*(1 + 0.4j)
    C = _split_cs(k, lm, mm)
    zl = rng.uniform(0.0, 2.0/k)
    kz = k*zl
    ch, sh = np.cosh(kz), np.sinh(kz)
    B = ch*C[0][0] + sh*C[1][0] + kz*ch*C[2][0] + kz*sh*C[3][0]
    Br = basis_cs(zl, k, lm, mm)
    wcs = max(wcs, np.max(np.abs(B - Br))/np.max(np.abs(Br)))
    D = _split_exp(k, lm, mm)
    za = rng.uniform(0.0, 40e3)
    zb = za + rng.uniform(1e3, 40e3)
    z = rng.uniform(za, zb)
    ed, eg, kz = np.exp(-k*(z - za)), np.exp(k*(z - zb)), k*z
    B = ed*D[0][0] + eg*D[1][0] + kz*ed*D[2][0] + kz*eg*D[3][0]
    Br = basis(z, k, lm, mm, za, zb)
    wex = max(wex, np.max(np.abs(B - Br))/np.max(np.abs(Br)))
print(f"T7a separated basis vs basis_cs: {wcs:.2e}, vs basis: {wex:.2e}")
ok = ok and wcs < 1e-13 and wex < 1e-13

# ---------------------------------------------------------------- (b)
wres = wcoef = wcond = 0.0
for trial in range(4):
    l2, m2 = lam, mu
    if trial % 2:
        l2, m2 = maxwell_moduli((0.3 + 1.7j)/1e10, mu, lam, 45*yr)
    gg = gon if trial < 2 else 0.0
    d = rng.uniform(2e3, 35e3)
    ks, ws = kgrid(1e-7, 8.0/H, 24, 4)
    sh, nh, x0 = (0.5, -0.866), (0.866, 0.5), rng.uniform(-2e4, 2e4)
    Jb = jumps_b(ks, lam, mu, 1.0, sh, nh, x0)
    cb = solve_k_b(ks, lam, mu, l2, m2, H, d, Jb, rho1, rho2, gg)
    for i in range(0, ks.size, 5):
        J = jumps(ks[i], lam, mu, 1.0, sh, nh, x0)
        wj = np.max(np.abs(J - Jb[i]))/np.max(np.abs(J))
        ok = ok and wj < 1e-14
        c = solve_k(ks[i], lam, mu, l2, m2, H, d, J, rho1, rho2, gg)
        wcoef = max(wcoef, np.max(np.abs(c - cb[i]))/np.max(np.abs(c)))
        # both must solve the reference system (assembled the way
        # solve_k assembles it) to machine precision
        k = ks[i]
        M = np.zeros((10, 10), dtype=complex)
        rhs = np.zeros(10, dtype=complex)
        B1 = region_basis(k, lam, mu, 0.0, d)
        B2 = region_basis(k, lam, mu, d, H)
        B3 = basis(H, k, l2, m2, H, H)[:, :2]
        M[0, 0:4] = B1(0.0)[2, :]
        M[1, 0:4] = B1(0.0)[3, :] - rho1*gg*B1(0.0)[1, :]
        M[2:6, 4:8] = B2(d)
        M[2:6, 0:4] = -B1(d)
        rhs[2:6] = J
        M[6:10, 8:10] = B3
        M[6:10, 4:8] = -B2(H)
        M[9, 8:10] -= (rho2 - rho1)*gg*B3[1, :]
        sc = np.max(np.abs(M), axis=0)
        sc[sc == 0.0] = 1.0
        wcond = max(wcond, np.linalg.cond(M/sc))
        rn = np.max(np.abs(rhs))
        for cc in (c, cb[i]):
            wres = max(wres, np.max(np.abs(M @ cc - rhs))/rn)
print(f"T7b jumps_b and solve_k_b vs scalars: worst coefficient "
      f"deviation {wcoef:.2e} at condition number up to {wcond:.1e};")
print(f"    worst scaled residual of either solution {wres:.2e}")
ok = ok and wcoef < 1e-7 and wres < 1e-12

# ---------------------------------------------------------------- (c)
ds, dip = 1e3, np.deg2rad(60.0)
dd = (np.arange(40) + 0.5)*ds
tx, tz = np.cos(dip), np.sin(dip)
xc = np.concatenate([dd*tx, [5e3, 20e3]])       # two substrate receivers
zc = np.concatenate([dd*tz, [45e3, 60e3]])
j = 17
src = segment_sources(xc[j] - ds/2*tx, zc[j] - ds/2*tz,
                      xc[j] + ds/2*tx, zc[j] + ds/2*tz, 1.0, ns=2)
kw = dict(kmin=1e-7, kmax=8.0/H, npanel=24, nquad=4)
worst = 0.0
for gg in (gon, 0.0):
    for ve in (False, True):
        l2, m2 = (maxwell_moduli((0.4 + 2.1j)/1e10, mu, lam, 45*yr)
                  if ve else (lam, mu))
        a = fields_pairs_ref(xc, zc, src, lam, mu, l2, m2, H,
                             rho1=rho1, rho2=rho2, g=gg, **kw)
        b = fields_pairs(xc, zc, src, lam, mu, l2, m2, H,
                         rho1=rho1, rho2=rho2, g=gg, **kw)
        for cmp_ in ("ux", "uz", "sxx", "sxz", "szz"):
            A, B = np.asarray(a[cmp_]), np.asarray(b[cmp_])
            if np.iscomplexobj(A) != np.iscomplexobj(B):
                print(f"  dtype mismatch in {cmp_}  FAIL")
                ok = False
            worst = max(worst, np.max(np.abs(B - A))/np.max(np.abs(A)))
print(f"T7c batched vs reference matched-pair fields: {worst:.2e}")
ok = ok and worst < 1e-10

# independent: the matched pairs must be the diagonal of fields_xz
xs = xc[[3, 9, 21, 33, 40]]
zs = zc[[3, 9, 21, 33, 40]]
l2, m2 = maxwell_moduli((0.4 + 2.1j)/1e10, mu, lam, 45*yr)
F = fields_xz(xs, zs, src, lam, mu, l2, m2, H,
              rho1=rho1, rho2=rho2, g=gon, **kw)
P = fields_pairs(xs, zs, src, lam, mu, l2, m2, H,
                 rho1=rho1, rho2=rho2, g=gon, **kw)
wd = 0.0
for cmp_ in ("ux", "uz", "sxx", "sxz", "szz"):
    dg = np.array([F[cmp_][i, i] for i in range(xs.size)])
    wd = max(wd, np.max(np.abs(P[cmp_] - dg))/np.max(np.abs(dg)))
print(f"T7d batched pairs vs diagonal of fields_xz: {wd:.2e}")
ok = ok and wd < 1e-10
print("T7 " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
