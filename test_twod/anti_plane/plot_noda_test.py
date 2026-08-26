#!/usr/bin/env python3
"""
plot_noda_test.py: render the two summary figures for run_noda_test from
the runs in the current directory (tmp_noda), together with precomputed
full-sweep values from this study (labeled as such) and the published
reference values of Kato (2002) Table 1 and Miyake & Noda (2019) Fig. 5.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def iv_cat(f, nmin=30):
    t = []
    for l in open(f):
        if l.startswith("#"): continue
        c = l.split()
        if int(c[5]) > nmin: t.append(float(c[1]))
    return np.diff(np.array(t))

def iv_chain(f):
    on = []
    for l in open(f):
        if l.startswith(" onsets"):
            on = sorted(set(float(x) for x in l.split(":")[1].split()))
    return np.diff(np.array(on))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

# ---------------- Kato panel -------------------------------------------------
Tel = 137.87
tr_sweep = np.array([1, 3, 10, 30, 100, 300])
# full-resolution sweep, digitized Fig. 2, virgin start, cycles 2-6
# (precomputed in this study; run gen_kato02.py + rsf_solve to reproduce)
Tcy_sweep = np.array([105.5, 106.9, 110.3, 116.9, 125.6, 131.9])
ax1.axhline(1.0, color="0.6", lw=1)
ax1.plot(tr_sweep, np.ones(6), "s--", color="0.35",
         label="Kato (2002) Table 1: 138 yr for all tr")
ax1.plot([300], [137/138.], "s", color="0.35")
ax1.plot(tr_sweep, Tcy_sweep/Tel, "o-", color="crimson",
         label="this study, full sweep (n=200, cycles 2-6)")
# fresh gating runs (reduced n=80 problem, own elastic reference)
el = iv_chain("chain_el.log"); ex = iv_chain("chain_ex.log")
sf = iv_chain("chain_sf.log"); rs = iv_cat("cat_kato_ve.dat")
ax1.plot([1], [rs[:2].mean()/el[0]], "D", color="darkred", ms=8,
         label=f"gate run: rsf_solve (n=80) {rs[:2].mean():.1f} yr")
ax1.plot([1], [ex[:2].mean()/el[0]], "x", color="k", ms=10, mew=2,
         label=f"gate run: exact Erlang chain {ex[:2].mean():.1f} yr")
ax1.plot([1], [sf[:2].mean()/el[0]], "*", color="royalblue", ms=14,
         label=f"chain, image family A only: {sf[:2].mean():.1f} yr\n"
               "(the single-family bug reproduces his null)")
ax1.set_xscale("log")
ax1.set_xlabel("tr = 2$\\eta$/G [yr]")
ax1.set_ylabel("Tcy / Tcy(elastic)")
ax1.set_ylim(0.6, 1.1)
ax1.set_title("Kato (2002) replica, digitized Fig. 2, virgin start")
ax1.legend(fontsize=7.2, loc="lower right")

# ---------------- M&N panel --------------------------------------------------
tload = 2.113386
# full scan at Rc/R = 0.35 (this study, single patch, late intervals)
xs = np.array([2.37, 1.7, 1.4, 1.2, 1.0, 0.9])
ts = np.array([9.35, 9.99, 11.35, 14.23, 28.2, 91.0])
# multi-patch (4 patches every 4.096R, approximating their periodicity)
xp = np.array([1.4, 1.2])
tp = np.array([13.6, 19.0])
ax2.axhline(12.47, color="0.6", lw=1)
ax2.annotate("elastic (single patch)", xy=(2.0, 12.9), fontsize=7, color="0.4")
ax2.plot(xs, ts, "o-", color="crimson", label="this study, QD single patch")
ax2.plot(xp, tp, "^-", color="royalblue",
         label="QD, 4 patches every 4.096R (stuck at 1.0)")
# divergence fit to the QD single-patch branch
from scipy.optimize import curve_fit
def f(x, a, b, xc): return a/(x-xc) + b
try:
    po, _ = curve_fit(f, xs[2:], ts[2:], p0=(2.0, 9.0, 0.8))
    xx = np.linspace(po[2]+0.02, 2.4, 200)
    ax2.plot(xx, f(xx, *po), ":", color="crimson", lw=1,
             label=f"fit: t_cr(QD, strip) = {po[2]:.2f}")
except Exception:
    pass
# their FD fit (Fig. 5): alpha 2.63 yr, beta 7.28 yr, tcr 1.5182
xx = np.linspace(1.55, 2.4, 100)
ax2.plot(xx, 2.63/(xx-1.5182)+7.28, "--", color="0.3",
         label="M&N (2019) Fig. 5 fit (FD, periodic):\n"
               "t_rec = 2.63/(x-1.518) + 7.28, elastic 7.15 yr")
# fresh gating runs
ive = iv_cat("cat_mn_ve.dat"); ich = iv_chain("chain_mn.log")
ax2.plot([1.0], [ive[:3].mean()], "D", color="darkred", ms=8,
         label=f"gate: rsf_solve mode 1 {ive[:3].mean():.1f} yr")
ax2.plot([1.0], [ich[:3].mean()], "x", color="k", ms=10, mew=2,
         label=f"gate: exact chain {ich[:3].mean():.1f} yr")
ax2.annotate("stuck (ST):\nQD strip 0.8\nQD periodic 1.0\nM&N FD 1.52",
             xy=(0.62, 55), fontsize=8)
ax2.set_yscale("log")
ax2.set_xlabel("$t_c/t_{load}$")
ax2.set_ylabel("$t_{rec}$ [yr]")
ax2.set_title("Miyake & Noda (2019) replica, $R_c^{EL}/R$ = 0.35, QD")
ax2.legend(fontsize=7.2, loc="upper right")

plt.tight_layout()
plt.savefig("noda_kato02_test.png", dpi=160,
            bbox_inches="tight") if False else None
fig.savefig("noda_tests.png", dpi=160)
# split into the two named figures for convenience: save each axis region
extent1 = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
extent2 = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
fig.savefig("noda_kato02_test.png", dpi=160,
            bbox_inches=extent1.expanded(1.35, 1.32))
fig.savefig("noda_mn19_test.png", dpi=160,
            bbox_inches=extent2.expanded(1.35, 1.32))
print("plot_noda_test: wrote noda_tests.png, noda_kato02_test.png, noda_mn19_test.png")
