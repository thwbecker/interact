"""
Run a quasi-dynamic and a fully dynamic earthquake sequence on the same
single fault (run_cycle1.setup) and plot the comparison.

Usage
    python3 compare_qd_fd.py                 # 4 events each, dx = 100 m, then plot
    python3 compare_qd_fd.py --events 2      # fewer events
    python3 compare_qd_fd.py --plot-only     # redraw from the saved results

Results are written to <prefix>_fd.npz and <prefix>_qd.npz (one file per
sequence, written after every event so an interrupted run keeps what it has)
and the figure to <prefix>.png. The second event is recorded in detail
(slip rate and stress snapshots) for the x-t and time-series panels; if only
one event is requested, the first is used.
"""
import argparse, os, time, warnings
import numpy as np
# rejected trial steps of the adaptive QD solver can evaluate the closure at
# wild states; the step is then rejected, so the resulting warnings are noise
warnings.filterwarnings("ignore", category=RuntimeWarning, module="scipy.integrate")

YR = 3.15576e7


def run_sequence(dynamic, n_events, dx, V_sw, prefix, record_event=1, verbose=True, friction="aging"):
    import run_cycle1
    tag = "fd" if dynamic else "qd"
    model, y0, x = run_cycle1.setup(V_sw=V_sw, dx=dx, dynamic=dynamic, friction=friction)
    model.record = 20
    model.verbose = False
    out = dict(x=x, dx=dx, dynamic=dynamic, t_start=[], Vmax=[], slip=[], nsteps=[])
    y, t = y0, 0.0
    t0 = time.time()
    for k in range(n_events):
        y, t = model.run(y, t_end=t + 400 * YR, max_events=1, t0=t)
        if len(model.events) <= k:
            print(f"[{tag}] no further event before t_end"); break
        e = model.events[k]
        out["t_start"].append(e["t_start"]); out["Vmax"].append(e["Vmax"])
        out["slip"].append(e["slip"]); out["nsteps"].append(e["nsteps"])
        if k == min(record_event, n_events - 1):
            out["snap_t"] = np.array([s[0] for s in e["snaps"]])
            out["snap_V"] = np.array([s[1] for s in e["snaps"]])
            out["snap_tau"] = np.array([s[2] for s in e["snaps"]])
        e["snaps"] = []   # free memory
        if verbose:
            print(f"[{tag}] event {k+1}: t = {e['t_start']/YR:8.3f} yr, Vmax {e['Vmax']:.2f} m/s, "
                  f"mean slip {e['slip'].mean():.3f} m, extent {np.sum(e['slip'] > 1e-2) * dx / 1e3:.1f} km, "
                  f"{e['nsteps']} steps ({time.time() - t0:.0f} s)")
        np.savez(f"{prefix}_{tag}.npz", **{kk: np.array(v) for kk, v in out.items()})
    return out


def load(prefix, tag):
    d = np.load(f"{prefix}_{tag}.npz")
    return {k: d[k] for k in d.files}


def plot(prefix, out_png=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fd, qd = load(prefix, "fd"), load(prefix, "qd")
    x = fd["x"] / 1e3; dx = float(fd["dx"])
    fig, ax = plt.subplots(3, 2, figsize=(12, 12))
    # slip per event
    a = ax[0, 0]
    for k in range(len(fd["slip"])):
        a.plot(x, fd["slip"][k], color="C3", lw=1.2, label="fully dynamic" if k == 0 else None)
    for k in range(len(qd["slip"])):
        a.plot(x, qd["slip"][k], color="C0", lw=1.2, ls="--", label="quasi-dynamic" if k == 0 else None)
    a.axvspan(-8, 8, color="0.9", zorder=0); a.set_xlabel("x (km)"); a.set_ylabel("coseismic slip per event (m)")
    a.set_title("slip per event (grey: velocity-weakening patch)"); a.legend()
    # recurrence
    a = ax[0, 1]
    for d, c, m, lab in ((fd, "C3", "o-", "FD"), (qd, "C0", "s--", "QD")):
        ts = d["t_start"] / YR
        rec = np.diff(ts).mean() if len(ts) > 2 else np.nan
        a.plot(np.arange(1, len(ts) + 1), ts, m, color=c, label=f"{lab}: recurrence {rec:.1f} yr" if np.isfinite(rec) else lab)
    a.set_xlabel("event number"); a.set_ylabel("event time (yr)"); a.set_title("recurrence"); a.legend()
    # x-t slip rate
    def aligned(d):
        t = d["snap_t"]; V = d["snap_V"]
        i0 = np.argmax(V.max(axis=1) > 1e-2); return t, V, t[i0]
    for j, (d, name) in enumerate(((fd, "fully dynamic"), (qd, "quasi-dynamic"))):
        t, V, t0 = aligned(d)
        m = (t >= t0 - 5) & (t <= t0 + 25)
        a = ax[1, j]
        im = a.pcolormesh(x, t[m] - t0, np.log10(np.maximum(V[m], 1e-6)), cmap="magma", vmin=-4, vmax=0.5, shading="nearest")
        a.set_xlabel("x (km)"); a.set_ylabel("time after max V exceeds 1e-2 m/s (s)")
        a.set_title(f"{name}: log10 slip rate, recorded event"); plt.colorbar(im, ax=a, label="log10 V (m/s)")
    # slip-rate time series
    a = ax[2, 0]
    for d, c, ls, lab in ((fd, "C3", "-", "FD"), (qd, "C0", "--", "QD")):
        t, V, t0 = aligned(d)
        a.plot(t - t0, V.max(axis=1), color=c, ls=ls, label=f"{lab} max over fault")
        ix = np.argmin(np.abs(x - 6.0))
        a.plot(t - t0, V[:, ix], color=c, ls=ls, lw=0.8, alpha=0.7, label=f"{lab} at x = 6 km")
    a.set_xlim(-2, 20); a.set_yscale("log"); a.set_ylim(1e-4, 5); a.set_xlabel("time (s)"); a.set_ylabel("slip rate (m/s)")
    a.set_title("slip rate during the recorded event"); a.legend(fontsize=8)
    # stress change
    a = ax[2, 1]
    for d, c, ls, lab in ((fd, "C3", "-", "FD"), (qd, "C0", "--", "QD")):
        tau = d["snap_tau"]; a.plot(x, (tau[-1] - tau[0]) / 1e6, color=c, ls=ls, label=lab)
    a.axvspan(-8, 8, color="0.9", zorder=0); a.set_xlabel("x (km)"); a.set_ylabel("shear stress change (MPa)")
    a.set_title("coseismic stress change, recorded event"); a.legend()
    plt.tight_layout()
    out_png = out_png or f"{prefix}.png"
    plt.savefig(out_png, dpi=110)
    # summary numbers
    for d, lab in ((fd, "FD"), (qd, "QD")):
        t, V, t0 = aligned(d)
        front = np.array([t[np.argmax(V[:, i] > 1e-2)] if (V[:, i] > 1e-2).any() else np.nan for i in range(len(x))])
        sel = (np.abs(x) > 2) & (np.abs(x) < 7) & np.isfinite(front)
        speed = 1 / np.polyfit(np.abs(x[sel]), front[sel], 1)[0] if sel.sum() > 3 else np.nan
        k = min(1, len(d["slip"]) - 1)
        drop = -(d["snap_tau"][-1] - d["snap_tau"][0])[np.abs(x) < 8].mean() / 1e6
        print(f"{lab}: Vmax {V.max():.2f} m/s, mean slip {d['slip'][k].mean():.3f} m, max slip {d['slip'][k].max():.3f} m, "
              f"extent {np.sum(d['slip'][k] > 1e-2) * dx / 1e3:.1f} km, patch stress drop {drop:.2f} MPa, "
              f"front speed 2-7 km {speed:.2f} km/s")
    print("figure:", out_png)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--events", type=int, default=4)
    p.add_argument("--dx", type=float, default=100.0)
    p.add_argument("--vsw", type=float, default=1e-3, help="QD to FD switch slip rate (m/s)")
    p.add_argument("--prefix", default="qd_fd")
    p.add_argument("--plot-only", action="store_true")
    p.add_argument("--only", choices=["fd", "qd"], default=None, help="run just one of the two sequences")
    p.add_argument("--friction", choices=["aging", "flash"], default="aging", help="aging law, or strong rate weakening (flash heating, slip law)")
    args = p.parse_args()
    if not args.plot_only:
        if args.only in (None, "qd"):
            run_sequence(False, args.events, args.dx, args.vsw, args.prefix, friction=args.friction)
        if args.only in (None, "fd"):
            run_sequence(True, args.events, args.dx, args.vsw, args.prefix, friction=args.friction)
    plot(args.prefix)
