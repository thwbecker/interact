#!/usr/bin/env python3
"""plot_p2.py rundir [rundir ...] -o out.png: per run, fault-averaged
tau/sigma against time (last 4000 yr) and the recurrence return map
T(n+1) against T(n) of the first fault group (after a quarter of the run)"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from rsf_series import Run, cat_columns  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("rundirs", nargs="+")
p.add_argument("-o", "--out", default="p2.png")
a = p.parse_args()
n = len(a.rundirs)
fig, ax = plt.subplots(n, 2, figsize=(12, 1.9 * n), gridspec_kw=dict(width_ratios=[4, 1]), squeeze=False)
for i, d in enumerate(a.rundirs):
    r = Run(d)
    sel = r.group == 0
    mu = (r.tau[:, sel] / r.sigma[:, sel]).mean(axis=1)
    w = r.t >= r.t[-1] - 4000
    ax[i, 0].plot(r.t[w], mu[w], "k-", lw=0.6)
    ax[i, 0].set_ylabel("tau/sigma")
    ax[i, 0].set_title(os.path.basename(os.path.abspath(d)), fontsize=9, loc="left")
    c = r.gcat[0] if r.gcat[0] is not None else r.cat
    gc = cat_columns(None, group=r.gcat[0] is not None)
    if c is not None:
        on = c[:, gc["onset"]]
        on = on[on > r.t[0] + 0.25 * (r.t[-1] - r.t[0])]
        dt = np.diff(on)
        if len(dt) > 2:
            ax[i, 1].plot(dt[:-1], dt[1:], "k.", ms=3)
            m = max(dt.max(), 1.0)
            ax[i, 1].plot([0, m], [0, m], ":", color="0.6")
    ax[i, 1].set_xlabel("T(n) [yr]", fontsize=8)
    ax[i, 1].set_ylabel("T(n+1)", fontsize=8)
ax[-1, 0].set_xlabel("time [yr]")
fig.tight_layout()
fig.savefig(a.out, dpi=120)
print("wrote", a.out)
