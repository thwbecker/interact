"""Dump a PairTable of the prototype as a text reference for the C port
(src/testing/rsf_solve/test_fd_tables.c).

usage: python3 dump_tables.py <case> <outfile>
cases: inplane_self, inplane_offset, antiplane_self, antiplane_offset
Header line: mode ns nr dx dt et alpha beta rho tol lmax M L X0s_x X0s_y ss_x ss_y X0r_x X0r_y sr_x sr_y
then nm lines: m ks kn, then L lines of nm values for W shear (per lag,
offsets in order), then (inplane) L lines for W normal."""
import sys, numpy as np
import shells as sh

AL, BE, RHO = 6000.0, 3464.0, 2700.0
case, out = sys.argv[1], sys.argv[2]
mode = "inplane" if case.startswith("inplane") else "antiplane"
dx = 100.0; dt = 0.5 * dx / BE; et = 1.0; tol = 1e-2; lmax = 600
fs = sh.Fault([-1000.0, 0.0], [1.0, 0.0], 2000.0, 20)
if case.endswith("self"):
    fr = fs
else:
    fr = sh.Fault([-300.0, 600.0], [1.0, 0.0], 1600.0, 16)   # parallel, offset along and across
M = int(2 ** np.ceil(np.log2(fr.n + fs.n - 1)))
p = sh.PairTable(fs, fr, dt, AL, BE, RHO, tol=tol, lmax=lmax, M=M, mode=mode, cap_ok=True, et=et)
with open(out, "w") as f:
    f.write(f"{0 if mode == 'inplane' else 1} {fs.n} {fr.n} {dx:.17g} {dt:.17g} {et:.17g} {AL:.17g} {BE:.17g} {RHO:.17g} {tol:.17g} {lmax} {M} {p.L} "
            f"{fs.X0[0]:.17g} {fs.X0[1]:.17g} {fs.s[0]:.17g} {fs.s[1]:.17g} {fr.X0[0]:.17g} {fr.X0[1]:.17g} {fr.s[0]:.17g} {fr.s[1]:.17g}\n")
    for mm, ks, kn in zip(p.m, p.Ks_row, p.Kn_row):
        f.write(f"{mm} {ks:.17g} {kn:.17g}\n")
    for comp in ("s", "n"):
        if comp == "n" and not p.has_normal:
            break
        hat = p.Ws_hat if comp == "s" else p.Wn_hat
        for l in range(p.L):
            c = np.fft.ifft(hat[l]).real
            vals = [c[mm] if mm >= 0 else c[p.M + mm] for mm in p.m]
            f.write(" ".join(f"{v:.12e}" for v in vals) + "\n")
print(case, "L =", p.L, "nm =", len(p.m), "M =", p.M)
