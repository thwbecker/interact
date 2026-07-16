#!/usr/bin/env python3
"""
analyze the single-patch slider benchmark produced by run_slider.

outputs, to stdout and to two whitespace tables:

  slider_physics.dat  : per law variant, from the reference runs:
                        recurrence, stress drop, peak slip rate, duration
  slider_wp.dat       : work-precision data, per (variant, solver, rtol):
                        cost measures and the event-time (phase) error
                        against that variant's reference run

the phase error is the RMS of |t_i - t_i_ref| over events matched by
index; runs whose event count differs from the reference are flagged
(the error is then computed over the common prefix and marked with n<).
cost is reported as accepted steps, rejected steps, RHS evaluations,
and (for imex) implicit function evaluations; total_rhs_equiv sums RHS
and implicit-function evaluations as a crude common currency, but note
an implicit evaluation is local (no matvec), so for LARGE systems the
matvec count (= RHS evaluations) is the better cost proxy, while for
this single-patch problem wall time is dominated by per-step overhead
and neither proxy maps to wall time directly.
"""
import numpy as np, os, re, glob

rdir = "results"

def events(d):
    p = os.path.join(d,"rsf_catalog.dat")
    if not os.path.exists(p): return None
    try:
        c = np.atleast_2d(np.loadtxt(p))
    except Exception:
        return np.zeros((0,13))
    if c.size == 0: return np.zeros((0,13))
    return c

def logstats(d):
    s = dict(steps=np.nan,rej=np.nan,rhs=np.nan,ifunc=np.nan,nli=np.nan,lin=np.nan)
    p = os.path.join(d,"log.dat")
    if not os.path.exists(p): return s
    txt = open(p, errors="ignore").read()
    for key,pat in [("rej",r"rejected steps=(\d+)"),
                    ("rhs",r"RHS function evaluations=(\d+)"),
                    ("ifunc",r"I function evaluations=(\d+)"),
                    ("nli",r"nonlinear solver iterations=(\d+)"),
                    ("lin",r"linear solver iterations=(\d+)")]:
        m = re.search(pat,txt)
        if m: s[key] = int(m.group(1))
    mp = os.path.join(d,"rsf_monitor.dat")
    if os.path.exists(mp):
        try:
            with open(mp) as f:
                last = f.readlines()[-1].split()
            s["steps"] = int(float(last[0])); s["t_end"] = float(last[2])
        except Exception: pass
    return s

# physics table from the references
variants = sorted(set(os.path.basename(d).rsplit("_",1)[0]
                      for d in glob.glob(f"{rdir}/*_ref")))
print("== physics (reference runs, 5dp, tight rtol; single 10 km patch,")
print("   one configuration; values are specific to these parameters) ==")
hdr = f"{'variant':8s} {'n_ev':>4s} {'recur[yr]':>10s} {'drop[MPa]':>10s} {'peakv[m/s]':>10s} {'dur[s]':>8s}"
print(hdr)
fp = open("slider_physics.dat","w"); fp.write("# "+hdr+"\n")
ref_ev = {}
for v in variants:
    c = events(f"{rdir}/{v}_ref")
    if c is None or len(c) < 2:
        print(f"{v:8s}  reference missing or fewer than 2 events"); continue
    on = c[:,1]; ref_ev[v] = on
    rec = np.diff(on)[1:] if len(on) > 2 else np.diff(on)  # skip first interval (IC transient)
    line = (f"{v:8s} {len(on):4d} {rec.mean():10.2f} {c[1:,8].mean():10.2f} "
            f"{c[1:,10].mean():10.3f} {c[1:,3].mean():8.1f}")
    print(line); fp.write(line+"\n")
fp.close()

# work-precision
print("\n== work-precision (phase error of matched events vs reference) ==")
hdr = (f"{'variant':8s} {'solver':6s} {'rtol':>6s} {'steps':>8s} {'rej':>7s} {'rhs':>8s} "
       f"{'ifunc':>8s} {'err_rms[yr]':>12s} {'err_max[yr]':>12s} {'n':>4s}")
print(hdr)
fw = open("slider_wp.dat","w"); fw.write("# "+hdr+"\n")
for d in sorted(glob.glob(f"{rdir}/*_*_1e-*")):
    b = os.path.basename(d)
    v,slab,rtol = b.rsplit("_",2)
    if v not in ref_ev: continue
    c = events(d); s = logstats(d)
    flag = ""
    if c is None or len(c) == 0:
        err = errm = np.nan; n = 0
    else:
        on = c[:,1]; ref = ref_ev[v]
        n = min(len(on),len(ref))
        if len(on) != len(ref): flag = "n<" if len(on) < len(ref) else "n>"
        dt = np.abs(on[:n]-ref[:n])
        err, errm = np.sqrt((dt**2).mean()), dt.max()
    line = (f"{v:8s} {slab:6s} {rtol:>6s} {s['steps']:8.0f} {s['rej']:7.0f} {s['rhs']:8.0f} "
            f"{s['ifunc']:8.0f} {err:12.4e} {errm:12.4e} {n:3d}{flag}")
    print(line); fw.write(line+"\n")
fw.close()
print("\nwrote slider_physics.dat and slider_wp.dat")
