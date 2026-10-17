"""Table of the stepover sequence results (npz files written by run_step.py).

For each event: source segment (earliest onset), receiver response within the
same dynamic phase (full rupture or minor re-slip, with the receiver onset
relative to the source front's arrival at the step), and the delay since the
previous event whose source was the other segment. Jumps are then binned:
immediate (< 60 s after the source front reached the step), minutes (< 1 h),
hours (< 2 d), days; computed from the saved rows, so older files work."""
import glob, sys, os, numpy as np
YR = 3.15576e7
JUMP_WINDOW = YR          # a rupture of the other segment within a year counts as a jump


def events_from_rows(d):
    rows = d["rows"]            # t_start, seg, t_on, t_off, slip, extent, Vmax[, t_end]
    out = []
    for ts in np.unique(rows[:, 0]):
        r = rows[rows[:, 0] == ts]
        segs = {int(x[1]): x for x in r}
        src = min(segs, key=lambda k: segs[k][2])
        e = dict(t=ts, src=src, slip_src=segs[src][4], rcv=None)
        if len(segs) == 2:
            rcv = 1 - src; x = segs[rcv]
            t_end_src = segs[src][7] if r.shape[1] > 7 else segs[src][3]
            e.update(rcv=rcv, slip_rcv=x[4], lag=x[2] - t_end_src, frac=x[4] / max(segs[src][4], 1e-9))
        out.append(e)
    return out


def fmt(sec):
    if sec < 60: return f"{sec:.1f} s"
    if sec < 3600: return f"{sec/60:.0f} min"
    if sec < 2 * 86400: return f"{sec/3600:.0f} h"
    if sec < 2 * 365.25 * 86400: return f"{sec/86400:.0f} d"
    return f"{sec/YR:.0f} yr"


def binname(sec):
    return "immediate" if sec < 60 else "minutes" if sec < 3600 else "hours" if sec < 2 * 86400 else "days" if sec < YR else "none"


files = sorted(glob.glob(sys.argv[1] if len(sys.argv) > 1 else "step*.npz"))
print(f"{'file':44s} {'mode':4s} {'W':>5s} {'margin':>6s}  per event: single, or delay of the jump to the other segment (seconds after the source front reached the step if within one dynamic phase)")
for f in files:
    d = np.load(f, allow_pickle=True)
    mode = "QD" if "qd" in os.path.basename(f) else "FD"
    import json
    try:
        vs = float(json.loads(str(d["params"])).get("vs_end", 4000.0)) / 1e3
    except Exception:
        vs = np.nan
    ev = events_from_rows(d)
    jumps, desc, cycles = [], [], []
    for i, e in enumerate(ev):
        if e["rcv"] is not None and e["frac"] > 0.1:              # both segments rupture in one phase
            jumps.append(e["lag"]); desc.append(f"{fmt(max(e['lag'], 0.0))}" + (" [joint]" if e["lag"] < 0 else ""))
        else:
            tag = f" (+re-slip {e['frac']*100:.0f}%)" if e["rcv"] is not None else ""
            # delay since the previous event with the other source
            if i > 0 and ev[i - 1]["src"] != e["src"] and e["t"] - ev[i - 1]["t"] < JUMP_WINDOW:
                dt = e["t"] - ev[i - 1]["t"]
                jumps.append(dt); desc.append(f"{fmt(dt)}{tag}")
            else:
                cycles.append(e["t"] - ev[i - 1]["t"] if i > 0 else np.nan)
                desc.append(f"single{tag}")
    bins = {b: sum(1 for j in jumps if binname(max(j, 0.0)) == b) for b in ("immediate", "minutes", "hours", "days")}
    print(f"{f:44s} {mode:4s} {float(d['W'])/1e3:+5.1f} {vs:6.2f}  " + "; ".join(desc))
    src_events = len(ev) - len(jumps)
    rec = np.nanmedian(cycles) / YR if len(cycles) > 1 else np.nan
    print(f"{'':44s}      {src_events} source events, {len(jumps)} followed by a jump within a year "
          f"({100*len(jumps)/max(src_events,1):.0f}%): " + ", ".join(f"{k} {v}" for k, v in bins.items() if v) +
          (f"; median delay {fmt(float(np.median(np.maximum(jumps, 0))))}" if jumps else "") +
          (f"; recurrence {rec:.0f} yr" if np.isfinite(rec) else ""))
