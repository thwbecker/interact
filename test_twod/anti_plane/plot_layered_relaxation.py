#!/usr/bin/env python3
"""
plot_layered_relaxation.py: render the Nur-Mavko single-event
relaxation comparison from the ve_layered_check profile files
({prefix}_elastic.dat, {prefix}_ve_disp.dat, {prefix}_ve_vel.dat):
analytic series as lines, the element-based machinery / np = 6 Prony
reconstruction as symbols, plus residual panels. Exits nonzero if
machinery and analytic disagree beyond tolerance.

usage: plot_layered_relaxation.py prefix [outdir]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    prefix = sys.argv[1] if len(sys.argv) > 1 else "prof"
    outdir = sys.argv[2] if len(sys.argv) > 2 else "."
    gams = [1, 0.75, 0.5, 0, -0.5, -1]
    tps = [0, 0.5, 1, 5, 10]
    de = np.loadtxt(prefix + "_elastic.dat")
    dd = np.loadtxt(prefix + "_ve_disp.dat")
    dv = np.loadtxt(prefix + "_ve_vel.dat")
    xh = de[:, 0]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    worst = [0.0, 0.0, 0.0]
    for i, g in enumerate(gams):
        c = plt.cm.coolwarm(i/(len(gams)-1.0))
        ana = de[:, 1+2*i]; num = de[:, 2+2*i]
        axes[0].plot(xh, ana, "-", color=c, lw=1.2, label=f"Gamma={g}")
        axes[0].plot(xh[::2], num[::2], "o", color=c, ms=3.5, mfc="none")
        worst[0] = max(worst[0], np.max(np.abs(num-ana)))
    axes[0].set_title("elastic layered profiles (Rybicki)\nlines analytic, symbols elements")
    axes[0].set_xlabel("x / H"); axes[0].set_ylabel("u / slip")
    axes[0].legend(fontsize=8)
    for i, tp in enumerate(tps):
        c = plt.cm.viridis(i/(len(tps)-1.0))
        ana = dd[:, 1+2*i]; num = dd[:, 2+2*i]
        axes[1].plot(xh, ana, "-", color=c, lw=1.2, label=f"t = {tp} tau_r")
        axes[1].plot(xh[::2], num[::2], "o", color=c, ms=3.5, mfc="none")
        worst[1] = max(worst[1], np.max(np.abs(num-ana)))
        ana = dv[:, 1+2*i]; num = dv[:, 2+2*i]
        axes[2].plot(xh, ana, "-", color=c, lw=1.2)
        axes[2].plot(xh[::2], num[::2], "o", color=c, ms=3.5, mfc="none")
        worst[2] = max(worst[2], np.max(np.abs(num-ana)))
    axes[1].set_title(f"postseismic displacement (Nur-Mavko)\nworst |fit-analytic| = {worst[1]:.1e}")
    axes[1].set_xlabel("x / H"); axes[1].set_ylabel("u / slip")
    axes[1].legend(fontsize=8)
    axes[2].set_title(f"postseismic velocity, scaled by tau_r\nworst |fit-analytic| = {worst[2]:.1e}")
    axes[2].set_xlabel("x / H"); axes[2].set_ylabel("v tau_r / slip")
    fig.tight_layout()
    fig.savefig(f"{outdir}/layered_relaxation.png", dpi=160)
    print(f"plot_layered_relaxation: elastic {worst[0]:.2e}, disp {worst[1]:.2e}, vel {worst[2]:.2e}")
    ok = (worst[0] < 1e-3) and (worst[1] < 3e-3) and (worst[2] < 2e-2)
    print("plot_layered_relaxation:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
