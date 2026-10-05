#!/usr/bin/env python3
"""check_restart.py: compare the catalogs of a continuous run and a
split (checkpoint-extended) run of the same case.

usage: check_restart.py cont_dir split_dir [tM] [tol_yr]
       tM      Maxwell time tag to compare besides el   (50)
       tol_yr  pass threshold on |onset difference|     (0.05)

For each tag (el, tm<tM>) it reports the event counts, the largest
onset difference, the largest relative differences of mean slip and
mean stress drop, and whether any event in the split catalog is
duplicated across the restart (onsets closer than 1e-6 yr).  Exit
status 0 if every tag passes, 1 otherwise.
"""
import sys, os
import numpy as np

cont = sys.argv[1]
split = sys.argv[2]
tm = sys.argv[3] if len(sys.argv) > 3 else '50'
tol = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
ONSET, SLIP, DROP = 1, 6, 8


def load(fn):
    d = np.loadtxt(fn, comments='#', ndmin=2)
    return d[np.argsort(d[:, ONSET])]


ok = True
for tag in ['el', f'tm{tm}']:
    fa, fb = os.path.join(cont, f'cat_{tag}.dat'), os.path.join(split, f'cat_{tag}.dat')
    if not (os.path.isfile(fa) and os.path.isfile(fb)):
        print(f'{tag}: missing catalog'); ok = False
        continue
    a, b = load(fa), load(fb)
    dup = int((np.diff(b[:, ONSET]) < 1e-6).sum())
    n = min(len(a), len(b))
    don = np.abs(a[:n, ONSET] - b[:n, ONSET])
    rsl = np.abs(a[:n, SLIP] / b[:n, SLIP] - 1)
    rdr = np.abs(a[:n, DROP] / b[:n, DROP] - 1)
    worst = int(np.argmax(don)) if n else -1
    passed = (len(a) == len(b)) and (dup == 0) and (n > 0) and (don.max() < tol)
    ok &= passed
    print(f'{tag}: events {len(a)} vs {len(b)}, duplicated across restart {dup}, '
          f'max |d onset| {don.max() if n else np.nan:.3e} yr (event {worst + 1}, '
          f't = {a[worst, ONSET] if n else np.nan:.1f} yr), '
          f'max rel d slip {rsl.max() if n else np.nan:.2e}, '
          f'max rel d drop {rdr.max() if n else np.nan:.2e}  -> {"PASS" if passed else "FAIL"}')
    if n and don.max() >= tol:
        # where does the divergence start
        first = int(np.argmax(don >= tol))
        print(f'      first onset difference >= {tol} yr at event {first + 1}, '
              f't = {a[first, ONSET]:.1f} yr')
print('restart test', 'PASSED' if ok else 'FAILED')
sys.exit(0 if ok else 1)
