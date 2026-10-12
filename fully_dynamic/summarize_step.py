"""Table of the stepover sequence results (step_*.npz written by run_step.py)."""
import glob, sys, numpy as np
YR = 3.15576e7
files = sorted(glob.glob(sys.argv[1] if len(sys.argv) > 1 else "step_*.npz"))
print(f"{'file':40s} {'mode':4s} {'W km':>6s} {'ov km':>6s}  events: kind (lag after source arrest / delay)")
for f in files:
    d = np.load(f, allow_pickle=True)
    mode = "QD" if "_qd_" in f else "FD"
    kinds = [str(k) for k in d["kinds"]]; lags = d["lags"]
    desc = []
    for k, l, t in zip(kinds, lags, d["t_start"]):
        if k in ("dynamic", "postseismic", "simultaneous"):
            desc.append(f"{k} ({l:+.1f} s)")
        else:
            desc.append(k)
    print(f"{f:40s} {mode:4s} {float(d['W'])/1e3:+6.1f} {float(d['overlap'])/1e3:6.1f}  " + "; ".join(desc))
    rows = d["rows"]
    if rows.size:
        for j in (0, 1):
            r = rows[rows[:, 1] == j]
            if r.size:
                print(f"{'':40s}      seg{j}: n={len(r)} events, slip {r[:,4].mean()/1e3:.1f} m km, extent {r[:,5].mean()/1e3:.1f} km, Vmax {r[:,6].mean():.2f} m/s (means)")
