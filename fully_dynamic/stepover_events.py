"""Event pairs on a two-segment stepover from rsf_solve station files.

usage: stepover_events.py <run directory> [station 1 name] [station 2 name] [jump window s]

Reads <dir>/fltst_<s1>.dat and <dir>/fltst_<s2>.dat (columns: t [s], slip [m], log10 |v| [m/s],
...; rsf_solve -rsf_stations), takes an event onset where log10 |v| first rises through -3
(1 mm/s) after having been below -6, and its coseismic slip as the slip increase over the
following 0.01 yr. Onsets on the two segments are paired in time order; a pair whose onsets
differ by less than the jump window (default 60 s) is reported as a jump (one rupture
across the step), otherwise as delayed triggering with its delay. Works for QD and FD runs
(the station series are dense while slip rates exceed 1e-6 m/s). For FD runs the station
series do not resolve the dynamic phases; there the phases in <dir>/run.log are used
instead (a plane counts as ruptured in a phase when its mean slip exceeds 0.1 m)."""
import sys, numpy as np

def load(f):
    rows = [l.split() for l in open(f) if l.strip() and not l.startswith("#")]
    return np.array([r for r in rows if len(r) >= 3], float)

def onsets(a):
    t, sl, lv = a[:, 0], a[:, 1], a[:, 2]
    out, armed = [], True
    for k in range(1, len(t)):
        if lv[k] < -6: armed = True
        if armed and lv[k - 1] < -3 <= lv[k]:
            tc = t[k - 1] + (t[k] - t[k - 1])*(-3 - lv[k - 1])/(lv[k] - lv[k - 1])
            s1 = np.interp(tc + 0.01*3.15576e7, t, sl)
            out.append((tc, s1 - np.interp(tc, t, sl)))
            armed = False
    return out

def fd_phases(d):
    import re
    ev = []
    for l in open(f"{d}/run.log"):
        m = re.search(r"3-D phase \d+ \(two planes\): t ([\d.]+) yr, \d+ steps over ([\d.]+) s, peak slip rate ([\d.]+) m/s, mean slip ([\d.]+) m \(plane 1\) and ([\d.]+) m", l)
        if m: ev.append(tuple(float(x) for x in m.groups()))
    return ev

if __name__ == "__main__":
    import os
    d = sys.argv[1]
    if os.path.exists(f"{d}/run.log") and any("3-D phase" in l for l in open(f"{d}/run.log")):
        yr = 3.15576e7; ph = [e for e in fd_phases(d) if max(e[3], e[4]) > 0.1]
        print(f"{d} (FD phases): {len(ph)} rupture phases")
        k = 0
        while k < len(ph):
            t, dur, vp, m1, m2 = ph[k]
            both = m1 > 0.1 and m2 > 0.1
            if both: print(f"  {t:12.6f} yr  both segments in one phase (jump): mean slip {m1:.2f} / {m2:.2f} m, peak {vp:.2f} m/s"); k += 1; continue
            if k + 1 < len(ph) and (ph[k + 1][0] - t) < 0.5 and ((m1 > 0.1) != (ph[k + 1][3] > 0.1)):
                dt = (ph[k + 1][0] - t)*yr
                print(f"  {t:12.6f} yr  segment {1 if m1 > 0.1 else 2} first ({max(m1, m2):.2f} m), the other {dt/3600:.2f} h = {dt/86400:.2f} days later ({max(ph[k + 1][3], ph[k + 1][4]):.2f} m): delayed")
                k += 2
            else:
                print(f"  {t:12.6f} yr  segment {1 if m1 > 0.1 else 2} alone ({max(m1, m2):.2f} m)"); k += 1
        sys.exit(0)
    d = sys.argv[1]; n1 = sys.argv[2] if len(sys.argv) > 2 else "s0"; n2 = sys.argv[3] if len(sys.argv) > 3 else "s1"
    win = float(sys.argv[4]) if len(sys.argv) > 4 else 60.0
    e1, e2 = onsets(load(f"{d}/fltst_{n1}.dat")), onsets(load(f"{d}/fltst_{n2}.dat"))
    yr = 3.15576e7
    print(f"{d}: {len(e1)} onsets on segment 1, {len(e2)} on segment 2")
    j = 0
    for (t1, s1) in e1:
        while j < len(e2) and e2[j][0] < t1 - 0.5*yr: j += 1
        if j < len(e2) and abs(e2[j][0] - t1) < 0.5*yr:
            t2, s2 = e2[j]; dt = t2 - t1
            kind = "jump" if abs(dt) < win else "delayed"
            print(f"  {t1/yr:10.4f} yr  segment 1 {s1:5.2f} m  segment 2 {s2:5.2f} m  {kind:7s} "
                  + (f"(second segment after {dt/3600:.2f} h = {dt/86400:.2f} days)" if kind == "delayed" else f"({dt:+.1f} s)"))
            j += 1
        else:
            print(f"  {t1/yr:10.4f} yr  segment 1 {s1:5.2f} m  segment 2 alone or not at all within 0.5 yr")
