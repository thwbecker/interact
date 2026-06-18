#!/usr/bin/env python3
"""
Score one rsf_solve run against a dense reference run, using the accuracy
metrics already adopted in rsf_solve.md: the log10(max|V|) time trace and the
event onset times.

Both runs write rsf_monitor.dat (columns: step, time[s], time[yr], dt[s],
log10(max|v|), mean_slip, mean_mu, max_sigma, min_sigma) and, with
-track_events, rsf_events.dat (columns: time[s], time[yr], onset(1)/arrest(-1),
log10(max|v|), mean_slip, mean_mu).

We compare on the trace by interpolating the run's log10(max|V|) onto the
reference run's time samples over their overlapping time span, then reporting
the maximum and RMS absolute deviation in log10 units.  We compare on timing by
matching event onsets (the +1 rows) and reporting the first-onset-time
difference and, when at least two onsets are present, the mean-recurrence
difference.  All differences are signed run-minus-reference.

This is a difference relative to the dense solve, not an absolute error, and it
is specific to the run configuration (resolution, rtol, dt_monitor, friction
setup).  A small deviation here means the H-matrix run reproduces the dense
sequence for this problem; it does not by itself certify the underlying physics.

Usage:
  python3 rsf_accuracy.py REF_DIR RUN_DIR
prints a single space-separated key=value line (parsed by the sweep harness)
and a short human-readable summary on stderr.
"""
import sys, os
import numpy as np

def load_trace(d):
    p = os.path.join(d, "rsf_monitor.dat")
    if not os.path.isfile(p):
        return None, None
    try:
        arr = np.loadtxt(p, comments="#")
    except Exception:
        return None, None
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[0] < 2:
        return None, None
    t = arr[:, 2]          # time [yr]
    lv = arr[:, 4]         # log10(max|v|)
    # keep strictly increasing time for interpolation
    keep = np.concatenate(([True], np.diff(t) > 0))
    return t[keep], lv[keep]

def load_onsets(d):
    p = os.path.join(d, "rsf_events.dat")
    if not os.path.isfile(p):
        return np.array([])
    try:
        arr = np.loadtxt(p, comments="#")
    except Exception:
        return np.array([])
    if arr.size == 0:
        return np.array([])
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    onset = arr[arr[:, 2] == 1]      # +1 rows are onsets
    return np.sort(onset[:, 1])      # onset times [yr]

def trace_dev(tr, lr, tx, lx):
    # interpolate run (tx,lx) onto reference times tr over the overlap
    lo = max(tr[0], tx[0]); hi = min(tr[-1], tx[-1])
    m = (tr >= lo) & (tr <= hi)
    if m.sum() < 2:
        return float("nan"), float("nan")
    li = np.interp(tr[m], tx, lx)
    dv = np.abs(li - lr[m])
    return float(dv.max()), float(np.sqrt(np.mean(dv**2)))

def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: rsf_accuracy.py REF_DIR RUN_DIR\n"); sys.exit(2)
    ref, run = sys.argv[1], sys.argv[2]

    tr, lr = load_trace(ref)
    tx, lx = load_trace(run)
    if tr is None or tx is None:
        print("maxdev=NA rmsdev=NA dt1=NA drec=NA nev_ref=NA nev_run=NA rec_ref=NA")
        sys.stderr.write("  accuracy: missing/short monitor trace\n"); return
    maxdev, rmsdev = trace_dev(tr, lr, tx, lx)

    er = load_onsets(ref); ex = load_onsets(run)
    nev_ref, nev_run = len(er), len(ex)
    dt1 = (ex[0] - er[0]) if (nev_ref >= 1 and nev_run >= 1) else float("nan")
    rec_ref = float(np.mean(np.diff(er))) if nev_ref >= 2 else float("nan")
    rec_run = float(np.mean(np.diff(ex))) if nev_run >= 2 else float("nan")
    drec = (rec_run - rec_ref) if (nev_ref >= 2 and nev_run >= 2) else float("nan")

    def f(x): return "NA" if (x != x) else f"{x:.6g}"   # NaN -> NA
    print(f"maxdev={f(maxdev)} rmsdev={f(rmsdev)} dt1={f(dt1)} drec={f(drec)} "
          f"nev_ref={nev_ref} nev_run={nev_run} rec_ref={f(rec_ref)}")
    sys.stderr.write(
        f"  accuracy vs dense: max d(log10|V|)={f(maxdev)}, rms={f(rmsdev)}, "
        f"first-event d={f(dt1)} yr, recurrence d={f(drec)} yr "
        f"(events ref={nev_ref} run={nev_run})\n")

if __name__ == "__main__":
    main()
