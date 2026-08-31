#!/usr/bin/env python3
"""
tau_scale_compare.py: compare two -ve_prony_file kernel files that
differ only in the substrate Maxwell time, one written by scaling the
tau ladder of a single sampling pass and one generated directly.

The claim under test (see bp3_ve_kernels.py, "MAXWELL-TIME SCALING"):
the substrate enters only through mu2(s) = mu2 s/(s + 1/tau), so the
relaxation kernel is a function of t/tau alone and the fitted
amplitude blocks C_p do not depend on the Maxwell time; only tau_p
scales.  Because the sampling times and the ladder are held in
dimensionless form, agreement is expected at roundoff level rather
than at fit-noise level.

usage: tau_scale_compare.py ref_file scaled_file ratio [tol]
       ratio  expected tau_scaled/tau_ref (e.g. 3 for 15 -> 45 yr)
       tol    tolerance on the amplitude blocks, relative to the
              largest amplitude in the family (default 1e-7, against
              a held-out fit residual of 1e-3 to 1e-2 in the same
              units, so this is a roundoff-level gate).  The
              blocks are written with 8 decimal digits, so this
              comparison cannot resolve better than about 1e-9; the
              arithmetic itself is checked at array level by
              t6_tau_scaling.py.
exit 0 if every check passes
"""
import sys
import numpy as np

def parse(path):
    """returns (n, npr, nfam, taus, blocks); blocks are (n, n) in
    file order: K0 gate, C_1..C_np [, N0 gate, Cn_1..Cn_np]"""
    with open(path) as f:
        rec = [l for l in f if not l.lstrip().startswith("#")
               and l.strip()]
    n, npr, nfam = (int(x) for x in rec[0].split())
    taus = np.array([float(x) for x in rec[1].split()])
    if taus.size != npr:
        sys.exit(f"{path}: ladder has {taus.size} entries, np = {npr}")
    nb = nfam*(npr + 1)
    need = 2 + nb*n
    if len(rec) < need:
        sys.exit(f"{path}: expected {need} records, found {len(rec)}")
    vals = np.array([[float(x) for x in rec[2 + i].split()]
                     for i in range(nb*n)])
    if vals.shape[1] != n:
        sys.exit(f"{path}: rows have {vals.shape[1]} entries, n = {n}")
    return n, npr, nfam, taus, [vals[b*n:(b + 1)*n] for b in range(nb)]

def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    fref, fscl = sys.argv[1], sys.argv[2]
    ratio = float(sys.argv[3])
    tol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-7
    nr, pr, fr, tr, Br = parse(fref)
    ns, ps, fs, tsc, Bs = parse(fscl)
    fail = 0
    if (nr, pr, fr) != (ns, ps, fs):
        print(f"  header mismatch: ({nr},{pr},{fr}) vs ({ns},{ps},{fs})"
              "  FAIL")
        return 1
    print(f"  n = {nr}, np = {pr}, nfam = {fr}")

    # 1. the ladder, and only the ladder, scales.  The files carry
    # tau_p with 8 decimal digits, so agreement is bounded by the
    # output format rather than by the arithmetic.
    dev = float(np.max(np.abs(tsc/(ratio*tr) - 1.0)))
    if dev < 1e-7:
        print(f"  ladder ratio: max |tau_s/({ratio:g} tau_r) - 1| = "
              f"{dev:.2e}  PASS")
    else:
        print(f"  ladder ratio: max |tau_s/({ratio:g} tau_r) - 1| = "
              f"{dev:.2e}  FAIL")
        fail = 1

    # 2. the elastic gate blocks carry no Maxwell time at all, so
    # they must come out of the two runs identically
    for fam in range(fr):
        b = fam*(pr + 1)
        scl = np.max(np.abs(Br[b])) + 1e-30
        dev = float(np.max(np.abs(Bs[b] - Br[b])))/scl
        lab = "K0" if fam == 0 else "N0"
        ok = dev < 1e-14
        print(f"  {lab} gate block: max dev {dev:.2e} of "
              f"{scl:.3e} Pa/m  {'PASS' if ok else 'FAIL'}")
        fail = fail or (not ok)

    # 3. the amplitude blocks must be Maxwell-time independent
    for fam in range(fr):
        b0 = fam*(pr + 1) + 1
        amp = [Br[b0 + p] for p in range(pr)]
        scl = max(float(np.max(np.abs(a))) for a in amp) + 1e-30
        dev = max(float(np.max(np.abs(Bs[b0 + p] - Br[b0 + p])))
                  for p in range(pr))/scl
        lab = "shear" if fam == 0 else "normal"
        ok = dev < tol
        print(f"  {lab} amplitudes C_p: worst dev {dev:.2e} "
              f"(of {scl:.3e} Pa/m, tol {tol:.0e})  "
              f"{'PASS' if ok else 'FAIL'}")
        fail = fail or (not ok)
    return 1 if fail else 0

if __name__ == "__main__":
    sys.exit(main())
