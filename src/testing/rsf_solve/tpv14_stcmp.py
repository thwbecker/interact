#!/usr/bin/env python3
"""Compare an on-fault station file of the SCEC benchmark (time, h-slip, h-slip-rate,
h-shear-stress, v-slip, v-slip-rate, v-shear-stress, n-stress, as saved from the SCEC
web site with the HTML removed or as plain text) with station files of fd3d_test
tpv14/tpv14s (FD3D_TPV14_STATIONS; time, slip, slip rate, shear, normal).

Usage: tpv14_stcmp.py scec_file our_file [our_file ...]
Prints rupture time (slip rate above 1 mm/s), peak slip rate, slip at 6, 8, 10 and
12 s, shear and normal stress at 8 and 12 s for each file.
"""
import re
import sys
import numpy as np


def load_scec(fn):
    rows = []
    for line in open(fn):
        line = re.sub(r'<[^>]+>', ' ', line).strip()
        if re.match(r'^-?\d\.\d+E[+-]\d+(\s+-?\d\.\d+E[+-]\d+){7}$', line):
            rows.append([float(x) for x in line.split()])
    return np.array(rows)


def summ(t, sl, sr, ts, ns):
    i = np.argmax(np.abs(sr) > 1e-3) if np.any(np.abs(sr) > 1e-3) else -1
    out = dict(rt=t[i] if i >= 0 else np.nan, pk=np.abs(sr).max(), tpk=t[np.abs(sr).argmax()], sl=sl[-1])
    for tt in (6, 8, 10, 12):
        k = np.argmin(abs(t - tt))
        out['s%d' % tt] = sl[k]; out['ts%d' % tt] = ts[k]; out['ns%d' % tt] = ns[k]
    return out
ref = load_scec(sys.argv[1])
r = summ(ref[:, 0], ref[:, 1], ref[:, 2], ref[:, 3], ref[:, 7])
rows = [('SCEC ' + sys.argv[1].split('/')[-1].replace('.dat', ''), r)]
for f in sys.argv[2:]:
    d = np.loadtxt(f, comments='#')
    rows.append((f, summ(d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4])))
for name, v in rows:
    print(f"{name:42s} rupture {v['rt']:6.2f} s, peak rate {v['pk']:5.2f} m/s at {v['tpk']:5.2f} s, slip at 6/8/10/12 s {v['s6']:.2f}/{v['s8']:.2f}/{v['s10']:.2f}/{v['s12']:.2f} m, "
          f"shear 8/12 s {v['ts8']:.1f}/{v['ts12']:.1f}, normal 8/12 s {v['ns8']:.1f}/{v['ns12']:.1f} MPa")
