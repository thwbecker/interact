"""Table of the stepover sequence results (step_*.npz written by run_step.py)."""
import glob, sys, numpy as np
YR = 3.15576e7
files = sorted(glob.glob(sys.argv[1] if len(sys.argv) > 1 else "step_*.npz"))
print(f"{'file':44s} {'mode':4s} {'W km':>6s} {'ov km':>6s}  events: kind (receiver onset after source front at step / delay; single: other segment's peak V)")
for f in files:
    d = np.load(f, allow_pickle=True)
    import os; mode = "QD" if "qd" in os.path.basename(f) else "FD"
    kinds = [str(k) for k in d["kinds"]]; lags = d["lags"]
    rcv = d["rcv_Vpeak"] if "rcv_Vpeak" in d.files else np.full(len(kinds), np.nan)
    desc = []
    for k, l, v in zip(kinds, lags, rcv):
        if k in ("dynamic", "postseismic", "joint", "simultaneous"):
            desc.append(f"{k} ({l:+.1f} s)")
        elif k == "single" and np.isfinite(v):
            desc.append(f"single ({v:.0e})")
        else:
            desc.append(k)
    print(f"{f:44s} {mode:4s} {float(d['W'])/1e3:+6.1f} {float(d['overlap'])/1e3:6.1f}  " + "; ".join(desc))
    rows = d["rows"]
    if rows.size:
        for j in (0, 1):
            r = rows[rows[:, 1] == j]
            if r.size:
                print(f"{'':44s}      seg{j}: n={len(r)} events, slip {r[:,4].mean()/1e3:.1f} m km, extent {r[:,5].mean()/1e3:.1f} km, Vmax {r[:,6].mean():.2f} m/s (means)")
