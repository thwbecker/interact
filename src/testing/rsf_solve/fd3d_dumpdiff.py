#!/usr/bin/env python3
"""Compare two FD3D_DUMP files (fd3d_test, hexadecimal values).

Usage: fd3d_dumpdiff.py a.txt b.txt [tolerance]
Prints the largest difference of each quantity relative to its largest
absolute value (over all faults, or over the vector); with a tolerance the
exit status is 1 if any relative difference exceeds it.
"""
import sys


def read(fn):
    d = {}
    for line in open(fn):
        w = line.split()
        if w[0] == 'f':
            for c, name in enumerate(('V', 'slip', 'tau', 'sigma')):
                d.setdefault(name, []).append(float.fromhex(w[3 + c]))
        else:
            d.setdefault(w[0], []).append(float.fromhex(w[2]))
    return d


a = read(sys.argv[1])
b = read(sys.argv[2])
tol = float(sys.argv[3]) if len(sys.argv) > 3 else None
worst = 0.0
out = []
for k in a:
    if k not in b or len(a[k]) != len(b[k]):
        print(f"{k}: sizes differ")
        sys.exit(1)
    m = max(abs(x) for x in a[k])
    e = max(abs(x - y) for x, y in zip(a[k], b[k]))
    r = e/m if m > 0 else e
    worst = max(worst, r)
    out.append(f"{k} {r:.1e}")
print(' '.join(out))
if tol is not None and worst > tol:
    sys.exit(1)
