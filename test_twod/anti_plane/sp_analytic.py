#!/usr/bin/env python3
"""
sp_analytic.py: independent analytic Savage-Prescott reference (route
C) for the ve_sp_cycle testbed, plus comparison plots.

The analytic solution is built from scratch here (no interact code):
surface velocity and displacement of the spun-up earthquake cycle on
a through-plate antiplane fault in an elastic plate (thickness H,
rigidity G) over a Maxwell half space of equal rigidity, from the
classic image expansion (Savage and Prescott, 1978; Segall, 2010,
ch. 12): a slip event du at time 0 has step response

  u(x,t)/du = (1/pi) { pi/2 sgn(x) - atan(x/H)
              + sum_{n>=1} W_n(t) [atan(x/((2n-1)H)) - atan(x/((2n+1)H))] }

with W_n(t) = P(n, b t) = 1 - exp(-bt) sum_{k<n} (bt)^k / k! the
Erlang ramp of order n and b = 1/(2 t_M); its rate uses the Erlang
density Wdot_n = b exp(-bt)(bt)^{n-1}/(n-1)!. The spun-up cycle with
events every T is the sum of these responses over all past events;
no other loading is required, since the relaxed response of each
event is rigid block motion +/- du/2.

Fault-plane stress (x = 0, depth z, G = 1) uses the telescoped image
families: n = 0 gives 1/(z+H) - 1/(z-H), family n >= 1 gives
1/(z-(2n-1)H) - 1/(z-(2n+1)H) + 1/(z+(2n+1)H) - 1/(z+(2n-1)H), each
weighted W_n(t), all times (1/2 pi) du.

Reads the ve_sp_cycle output files {prefix}_vel_profiles.dat,
{prefix}_stress_hist.dat, {prefix}_surf_hist.dat and
{prefix}_checks.dat, overlays the machinery (route A), the exact
element-based route B, and this analytic route C, writes PNG figures
and a deviation summary, and exits nonzero if C disagrees with B
beyond tolerance (B and C share the exact time dependence, so B-C
measures spatial-part and truncation consistency only).

usage: sp_analytic.py prefix H t_M T_over_taur [outdir]
"""
import sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def erlang_all(nmax, bt):
    """W_n = P(n, bt) and Wdot_n/b, n = 0..nmax, by upward recurrence;
    bt may be a numpy array (vectorized over ages)"""
    bt = np.atleast_1d(np.asarray(bt, dtype=float))
    w = np.zeros((nmax+1, bt.size)); wd = np.zeros((nmax+1, bt.size))
    w[0] = 1.0
    eterm = np.exp(-bt); psum = np.zeros_like(bt)
    for n in range(1, nmax+1):
        psum = psum + eterm
        w[n] = np.maximum(1.0 - psum, 0.0)
        wd[n] = eterm          # exp(-bt)(bt)^(n-1)/(n-1)!, times b later
        eterm = eterm*bt/n
    return w, wd

def families_disp(x, H, nmax):
    """elastic image-family surface displacement brackets, unit slip"""
    X = np.zeros((nmax+1, x.size))
    X[0] = (0.5*np.pi*np.sign(x) - np.arctan(x/H))/np.pi
    for n in range(1, nmax+1):
        X[n] = (np.arctan(x/((2*n-1)*H)) - np.arctan(x/((2*n+1)*H)))/np.pi
    return X

def families_stress(z, H, nmax):
    """elastic image-family fault-plane stress at depth z, unit slip, G=1"""
    S = np.zeros((nmax+1, z.size))
    S[0] = (1.0/(z+H) - 1.0/(z-H))/(2*np.pi)
    for n in range(1, nmax+1):
        S[n] = (1.0/(z-(2*n-1)*H) - 1.0/(z-(2*n+1)*H)
                + 1.0/(z+(2*n+1)*H) - 1.0/(z+(2*n-1)*H))/(2*np.pi)
    return S

def spun_weights(nmax, b, t, T, K, rate=False):
    """sum over K past events of W_n (or b-scaled Wdot_n) at ages t + j T"""
    ages = t + T*np.arange(K)
    w, wd = erlang_all(nmax, b*ages)
    return (b*wd).sum(axis=1) if rate else w.sum(axis=1)

def main():
    if len(sys.argv) < 5:
        print(__doc__); sys.exit(1)
    prefix = sys.argv[1]
    H = float(sys.argv[2]); t_M = float(sys.argv[3]); Ttr = float(sys.argv[4])
    outdir = sys.argv[5] if len(sys.argv) > 5 else "."
    b = 1.0/(2.0*t_M); taur = 1.0/b; T = Ttr*taur
    Vpl = 1.0; du = Vpl*T

    vel = np.loadtxt(prefix + "_vel_profiles.dat")
    sth = np.loadtxt(prefix + "_stress_hist.dat")
    sfh = np.loadtxt(prefix + "_surf_hist.dat")
    checks = dict()
    with open(prefix + "_checks.dat") as f:
        for line in f:
            if line.startswith("#"): continue
            k, v = line.split(); checks[k] = float(v)
    # match the testbed's truncations exactly: same image count and
    # the same number of superposed events
    nmax = int(checks["n_img"]); K = int(checks["spinup_cycles"])
    tfrac = [0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 0.98]
    xh = vel[:, 0]; x = xh*H
    Xd = families_disp(x, H, nmax)

    # ---- route C velocity profiles and comparison -------------------
    vC = np.zeros((len(tfrac), x.size))
    for k, tf in enumerate(tfrac):
        swd = spun_weights(nmax, b, tf*T, T, K, rate=True)
        vC[k] = du*(swd[1:, None]*Xd[1:]).sum(axis=0)
    vA = vel[:, 1::3].T; vA2 = vel[:, 2::3].T; vB = vel[:, 3::3].T
    devAC = np.max(np.abs(vA - vC))/Vpl
    devBC = np.max(np.abs(vB - vC))/Vpl

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharex=True)
    cmap = plt.cm.viridis
    for k, tf in enumerate(tfrac):
        c = cmap(k/(len(tfrac)-1.0))
        axes[0].plot(xh, vC[k], "-", color=c, lw=1.2,
                     label=f"t/T={tf}" if k % 2 == 0 else None)
        axes[0].plot(xh[::2], vA[k][::2], "o", color=c, ms=3.5, mfc="none")
        axes[1].semilogy(xh, np.abs(vA[k]-vC[k])/Vpl + 1e-12, "-", color=c, lw=1)
    axes[0].set_xlabel("x / H"); axes[0].set_ylabel("v / V_pl")
    axes[0].set_title(f"spun-up cycle velocity, T/tau_r = {Ttr:g}\n"
                      "lines: analytic (Savage-Prescott); symbols: h-state machinery")
    axes[0].legend(fontsize=8, loc="upper right")
    axes[1].set_xlabel("x / H"); axes[1].set_ylabel("|machinery - analytic| / V_pl")
    axes[1].set_title(f"worst {devAC:.1e} (exact-time route B - C: {devBC:.1e})")
    fig.tight_layout()
    fig.savefig(f"{outdir}/sp_cycle_velocity_T{Ttr:g}.png", dpi=160)

    # ---- route C stress histories ------------------------------------
    # actual receiver patch depths from the testbed; TRANSIENT image
    # weights (W_n - 1): the exact relaxed stress of the through-plate
    # rupture is zero, so the event count drops out
    z = np.array([checks["depth_0"], checks["depth_1"], checks["depth_2"]])
    zh = z/H
    Sd = families_stress(z, H, nmax)
    tt = sth[:, 0]
    sC = np.zeros((tt.size, z.size))
    for i, t in enumerate(tt):
        nev = K if t >= 0 else K-1
        j0 = 0 if t >= 0 else 1
        t0 = t*T if t >= 0 else (t+1)*T
        ages = t0 + T*np.arange(j0, K)
        w, _ = erlang_all(nmax, b*ages)
        sw = w.sum(axis=1)
        sC[i] = du*((sw[:, None] - nev)*Sd).sum(axis=0)
    sA = sth[:, 1::2]; sB = sth[:, 2::2]
    # element/analytic stress can differ by an overall sign convention of
    # the resolved slip mode; calibrate once
    sgn = 1.0 if np.sum(sB*sC) >= 0 else -1.0
    sC *= sgn
    sswing = np.max(np.abs(sB))
    devS_AC = np.max(np.abs(sA - sC))/sswing
    devS_BC = np.max(np.abs(sB - sC))/sswing
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for j, zf in enumerate(zh):
        c = plt.cm.plasma(j/2.0)
        axes[0].plot(tt, sC[:, j], "-", color=c, lw=1.2, label=f"z/H={zf}")
        axes[0].plot(tt[::4], sA[::4, j], "o", color=c, ms=3.5, mfc="none")
        axes[1].semilogy(tt, np.abs(sA[:, j]-sC[:, j])/sswing + 1e-12, "-", color=c)
    axes[0].set_xlabel("t / T (negative: previous cycle)")
    axes[0].set_ylabel("resolved fault stress (G = 1)")
    axes[0].set_title(f"fault stress through the spun-up cycle, T/tau_r = {Ttr:g}\n"
                      "lines: analytic; symbols: h-state machinery")
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel("t / T"); axes[1].set_ylabel("|machinery - analytic| / swing")
    axes[1].set_title(f"worst {devS_AC:.1e} (B - C: {devS_BC:.1e})")
    fig.tight_layout()
    fig.savefig(f"{outdir}/sp_cycle_stress_T{Ttr:g}.png", dpi=160)

    # ---- surface displacement histories ------------------------------
    ttd = sfh[:, 0]
    xsel = [checks["xobs_1"]/H, checks["xobs_3"]/H]
    iA = [np.argmin(np.abs(xh - xs)) for xs in xsel]
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    dev_u = 0.0
    for j, xs in enumerate(xsel):
        c = plt.cm.viridis(j/1.0*0.7)
        xg = np.array([xh[iA[j]]*H])
        Xg = families_disp(xg, H, nmax)
        uC = np.zeros(ttd.size)
        sw0 = spun_weights(nmax, b, 1e-12*T, T, K)
        u0 = du*(sw0*Xg[:, 0]).sum()
        for i, t in enumerate(ttd):
            sw = spun_weights(nmax, b, t*T, T, K)
            uC[i] = du*(sw*Xg[:, 0]).sum() - u0
        uC /= du
        uA = sfh[:, 1+2*j]
        ax.plot(ttd, uC, "-", color=c, lw=1.2, label=f"x/H={xh[iA[j]]:.2f} analytic")
        ax.plot(ttd[::4], uA[::4], "o", color=c, ms=3.5, mfc="none")
        dev_u = max(dev_u, np.max(np.abs(uA-uC)))
    ax.set_xlabel("t / T"); ax.set_ylabel("u / du")
    ax.set_title(f"surface displacement through the cycle, T/tau_r = {Ttr:g}\n"
                 f"worst |machinery - analytic| = {dev_u:.1e} du")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{outdir}/sp_cycle_surfdisp_T{Ttr:g}.png", dpi=160)

    print(f"sp_analytic T/taur={Ttr:g}: velocity |A-C| {devAC:.2e} |B-C| {devBC:.2e} V_pl; "
          f"stress |A-C| {devS_AC:.2e} |B-C| {devS_BC:.2e} of swing; surf disp |A-C| {dev_u:.2e} du")
    # gates: B and C share exact time dependence; B-C tests spatial parts
    # and truncations, A-C the full machinery
    ok = (devBC < 3e-3) and (devAC < 1.2e-2) and (devS_BC < 3e-3) and \
         (devS_AC < 2e-2) and (dev_u < 6e-3)
    print("sp_analytic:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
