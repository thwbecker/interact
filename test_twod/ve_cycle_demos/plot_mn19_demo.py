#!/usr/bin/env python3
"""
plot_mn19_demo.py: recurrence vs t_c/t_load from a run_mn19_demo
directory, with the Miyake & Noda (2019) divergence law
t_rec = alpha/(x - x_cr) + beta fitted to our quasi-dynamic points
and their fully-dynamic periodic-fault reference curve
(alpha = 2.63 yr, beta = 7.28 yr, x_cr = 1.5182) for regime
comparison (QD isolated strips diverge near x_cr ~ 0.86; the offset
from their 1.52 is inertia plus periodicity, see noda_mn_tests.md).

usage: plot_mn19_demo.py rundir t_load_yr
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rundir = sys.argv[1]
tload = float(sys.argv[2])
os.chdir(rundir)

def trec(f, nskip=3):
    t = [float(l.split()[1]) for l in open(f) if not l.startswith("#")]
    if len(t) < nskip + 3:
        return None
    iv = np.diff(t[nskip:])
    return float(np.mean(iv))

xs, tr = [], []
for f in sorted(glob.glob("cat_x*.dat")):
    x = float(f.replace("cat_x", "").replace(".dat", ""))
    r = trec(f)
    if r:
        xs.append(x); tr.append(r)
xs, tr = np.array(xs), np.array(tr)
te = trec("cat_el.dat")
print("elastic t_rec:", te, "yr")
for x, r in zip(xs, tr):
    print(f"  x = {x:5.2f}: t_rec = {r:7.2f} yr")

fig, ax = plt.subplots(figsize=(7.5, 5.5))
ax.plot(xs, tr, "ks", ms=7, mfc="none", label="rsf_solve (QD, this demo)")
if te:
    ax.axhline(te, color="0.6", lw=0.8, ls=":",
               label=f"elastic ({te:.1f} yr)")
# fit our points
if xs.size >= 4:
    from scipy.optimize import curve_fit
    def law(x, a, b, xc):
        return a/(x - xc) + b
    try:
        p, _ = curve_fit(law, xs, tr, p0=[2.0, te or 8.0, 0.8],
                         maxfev=20000)
        xf = np.linspace(max(p[2] + 0.02, xs.min() - 0.3), xs.max(), 200)
        ax.plot(xf, law(xf, *p), "k-", lw=1,
                label=f"fit: {p[0]:.2f}/(x - {p[2]:.2f}) + {p[1]:.2f}")
        print(f"fit: alpha {p[0]:.3f} beta {p[1]:.3f} x_cr {p[2]:.3f}")
    except Exception as e:
        print("fit failed:", e)
# their reference law (fully dynamic, periodic)
xf = np.linspace(1.60, max(3.2, xs.max()), 200)
ax.plot(xf, 2.63/(xf - 1.5182) + 7.28, "r--", lw=1.2,
        label="Miyake & Noda 2019 (FD, periodic):\n2.63/(x-1.52)+7.28")
ax.set_xlabel("t_c / t_load")
ax.set_ylabel("recurrence interval [yr]")
ax.set_title("Viscoelastic recurrence divergence "
             "(Miyake & Noda 2019 systematics)", fontsize=11)
ax.legend(fontsize=8)
ax.set_ylim(0, max(np.max(tr)*1.3, 40))
fig.tight_layout()
fig.savefig("mn19_demo.png", dpi=150, bbox_inches="tight")
print("wrote mn19_demo.png")
