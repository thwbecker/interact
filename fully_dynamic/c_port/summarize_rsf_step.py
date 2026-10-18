"""Jump table of a stepover sequence run with rsf_solve (run_step_fd.sh):
events per segment from the station files fltst_s1.dat and fltst_s2.dat
(centre of each segment; an event is log10|v| crossing -2 upward or a slip
step above 0.1 m between samples, merged within 300 s), then for each event the delay since the previous event on the
other segment if within a year. usage: summarize_rsf_step.py <run_dir> [...]"""
import sys, os, numpy as np
YR = 3.15576e7
def events(fn, thr=-2.0, dslip=None):
    """an event at the station: log10|v| crosses thr upward (QD, dense output),
    or the slip jumps by more than dslip between consecutive samples (FD: the
    coseismic slip happens inside the dynamic phase, invisible to the QD
    station writer, and appears as a step)"""
    a = np.loadtxt(fn); t = a[:, 0]; lv = a[:, 2]; sl = a[:, 1]
    on = (lv[1:] > thr) & (lv[:-1] <= thr)
    if dslip is not None: on = on | (np.diff(sl) > dslip)
    ts = t[1:][on]; out = []
    for s in ts:
        if not out or s - out[-1] > 300.0: out.append(s)
    return np.array(out)
def fmt(sec):
    if sec < 60: return f"{sec:.1f} s"
    if sec < 3600: return f"{sec/60:.0f} min"
    if sec < 2 * 86400: return f"{sec/3600:.0f} h"
    return f"{sec/86400:.0f} d"
def events_fd(d, n1=None, dslip=0.1):
    """FD runs: one event per segment per dynamic phase, from the phase lines of
    run.log (time to 1e-8 yr) and fd_phase_<n>_slip.dat (which segment slipped);
    both segments in one phase are listed as the source first (the one with the
    larger slip) and the other at the same time (joint rupture)"""
    import re, glob
    if n1 is None:
        f = os.path.join(d, "segments.txt")
        n1 = int(open(f).read().split()[0]) if os.path.exists(f) else 240
    ev = []
    for line in open(os.path.join(d, "run.log")):
        m = re.match(r"rsf_fd: phase (\d+): t ([0-9.]+) yr", line)
        if not m: continue
        n = int(m.group(1)); t = float(m.group(2)) * YR
        fn = os.path.join(d, f"fd_phase_{n}_slip.dat")
        if not os.path.exists(fn): continue
        a = np.loadtxt(fn); g = a[:, 1].astype(int); sl = a[:, 2]
        s0 = sl[g < n1].mean(); s1 = sl[g >= n1].mean()
        segs = [(s0, 0), (s1, 1)]
        segs.sort(reverse=True)
        for k, (m_, s) in enumerate(segs):
            if m_ > dslip: ev.append((t + 1e-3 * k, s))
    return ev
for d in sys.argv[1:]:
    if os.path.exists(os.path.join(d, "fd_phase_1_slip.dat")):
        ev = events_fd(d)
    else:
        e1 = events(os.path.join(d, "fltst_s1.dat")); e2 = events(os.path.join(d, "fltst_s2.dat"))
        ev = sorted([(t, 0) for t in e1] + [(t, 1) for t in e2])
    desc, jumps, cyc = [], [], []
    for i, (t, s) in enumerate(ev):
        if i > 0 and ev[i - 1][1] != s and t - ev[i - 1][0] < YR:
            jumps.append(t - ev[i - 1][0]); desc.append(f"seg{s} {fmt(t - ev[i-1][0])} after seg{1-s}")
        else:
            desc.append(f"seg{s} at {t/YR:.2f} yr")
            if i > 0: cyc.append(t - ev[i - 1][0])
    nsrc = len(ev) - len(jumps)
    bins = {"immediate (<60 s)": sum(j < 60 for j in jumps), "minutes": sum(60 <= j < 3600 for j in jumps),
            "hours": sum(3600 <= j < 2 * 86400 for j in jumps), "days": sum(j >= 2 * 86400 for j in jumps)}
    print(f"{d}: {len(ev)} events; " + "; ".join(desc))
    p1d = sum(j < 86400 for j in jumps) / max(nsrc, 1); p30d = sum(j < 30 * 86400 for j in jumps) / max(nsrc, 1)
    print(f"    P(jump within 1 d) = {p1d:.2f}, within 30 d = {p30d:.2f} ({nsrc} source events)")
    print(f"    {nsrc} source events, {len(jumps)} followed by a jump within a year: " +
          ", ".join(f"{k} {v}" for k, v in bins.items() if v) +
          (f"; median delay {fmt(float(np.median(jumps)))}" if jumps else "") +
          (f"; recurrence {np.median(cyc)/YR:.1f} yr" if cyc else ""))
