#!/usr/bin/env python3
"""T2: uniform-Maxwell correspondence gate.  In an unbounded Maxwell
medium (elastic bulk, Maxwell shear) every stress component of a held
edge dislocation evolves with ONE exact scalar time factor,
  f(t) = L^-1[ mu(s)(lam(s)+mu(s)) / ((lam(s)+2mu(s)) s) ] / D0,
  D0 = mu(lam+mu)/(lam+2mu),
because the full-space field is proportional to mu/(1-nu).  We place
the fault 100 km deep (free-surface correction < 1 percent), run the
LAYERED solver in uniform-Maxwell mode through the Talbot inversion,
and compare each stress component at several points and times against
sigma_elastic * f(t) with f(t) from sympy (exact partial fractions).
This gates the s-plane moduli plumbing AND the Talbot inversion."""
import numpy as np
import sympy as sp
from inplane2d import fields_xz, ve_fields_xz, segment_sources

mu_, lam_ = 30e9, 30e9
tau = 10.0*3.15576e7                     # 10 yr Maxwell time
d1, d2, H = 98e3, 102e3, 150e3
src = segment_sources(0.0, d1, 0.0, d2, 1.0, ns=16)
xs = np.array([8e3, 14e3])
zs = np.array([100e3 - 16e3, 100e3 + 16e3])
kw = dict(kmin=1e-7, kmax=2e-3, npanel=40, nquad=6)

# exact f(t) by sympy
s, t_ = sp.symbols("s t", positive=True)
a = 1.0/tau
K = lam_ + 2.0*mu_/3.0
mus = mu_*s/(s + a)
lams = K - 2*mus/3
D0 = mu_*(lam_ + mu_)/(lam_ + 2*mu_)
Fs = mus*(lams + mus)/((lams + 2*mus)*s)/D0
f_exact = sp.inverse_laplace_transform(sp.nsimplify(sp.simplify(Fs)),
                                       s, t_)
f = sp.lambdify(t_, f_exact, "numpy")

el = fields_xz(xs, zs, src, lam_, mu_, lam_, mu_, H, **kw)
ok = True
for tfac in (0.3, 1.0, 3.0):
    t = tfac*tau
    ve = ve_fields_xz(xs, zs, src, lam_, mu_, lam_, mu_, H, tau,
                      t, tau1=tau, M=14, **kw)
    fx = float(f(t))
    for comp in ("sxx", "sxz", "szz"):
        r = ve[comp]/(el[comp]*fx)
        err = np.max(np.abs(r - 1.0))
        print(f"t = {tfac:3.1f} tau  {comp}: max |ratio-1| = {err:8.2e}"
              f"   (f_exact = {fx:.5f})")
        ok = ok and err < 0.02
print("T2 " + ("PASS" if ok else "FAIL"))
raise SystemExit(0 if ok else 1)
