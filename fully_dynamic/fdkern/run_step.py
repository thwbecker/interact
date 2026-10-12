"""
WP6: earthquake sequences on two parallel segments (stepover) with
rate-and-state friction, QD or FD, with per-event bookkeeping and a jump
classification.

Geometry (right-lateral convention of shells.py): segment 1 along +x from
x = -L1 to 0 at y = 0; segment 2 from x = -overlap to L2 - overlap at
y = W (W < 0: right step, dilational for right-lateral slip; W > 0: left
step, compressional). Each segment is velocity weakening except for
velocity-strengthening margins of length `vs_end` at both ends.

Event bookkeeping (from the FD phase snapshots): for each event, per
segment, time of first slip above V_thr, time of arrest (last slip above
V_thr), moment-like slip integral, peak slip rate. Classification for
events that involve both segments:
    dynamic     : receiver starts within t_dyn of the source arrest
                  (t_dyn = W/beta + 3 s by default)
    postseismic : receiver starts later in the same FD phase
Events with only one segment are "single".
"""
import argparse, time
import numpy as np
import shells as sh, stepper as st, cycle as cy

YR = 3.15576e7


def build(W=-1000.0, overlap=0.0, L1=24e3, L2=24e3, dx=100.0, vs_end=4e3, a_vw=0.010, a_vs=0.020,
          b=0.015, dc=0.016, sigma0=50e6, Vpl=1e-9, friction="aging", Vw=0.1, fw=0.2, tol=1e-2,
          V_sw=1e-3, V_back=1e-4, dynamic=True, perturb_x=None):
    AL, BE, RHO = 6000.0, 3464.0, 2700.0
    f1 = sh.Fault([-L1, 0.0], [1.0, 0.0], L1, int(L1 / dx), "f1")
    f2 = sh.Fault([-overlap, W], [1.0, 0.0], L2, int(L2 / dx), "f2")
    dt = 0.5 * dx / BE
    sys_ = sh.System([f1, f2], dt, AL, BE, RHO, tol=tol)
    x = np.concatenate([f1.centers[:, 0], f2.centers[:, 0]])
    n1 = f1.n
    seg = np.concatenate([np.zeros(n1, int), np.ones(f2.n, int)])
    # along-segment coordinate for the VS margins
    s_along = np.concatenate([f1.centers[:, 0] + L1, f2.centers[:, 0] + overlap])
    Lseg = np.where(seg == 0, L1, L2)
    vs = (s_along < vs_end) | (s_along > Lseg - vs_end)
    a = np.where(vs, a_vs, a_vw); bb = np.full(sys_.N, b)
    V0 = 1e-6; f0 = 0.6
    fr = st.RateState(a, bb, dc, V0, f0) if friction == "aging" else st.FlashRateState(a, bb, dc, V0, f0, Vw=Vw, fw=fw)
    psi_ss = f0 + bb * np.log(V0 / Vpl)
    Vinit = np.full(sys_.N, Vpl)
    if perturb_x is None:
        perturb_x = -L1 / 2          # centre of segment 1
    Vinit[(seg == 0) & (np.abs(x - perturb_x) < 2e3)] *= 1.5
    tau0 = fr.tau_of_V(Vinit, psi_ss, sigma0)
    model = cy.CycleModel(sys_, fr, tau0, sigma0, Vpl, normal_coupling=True, V_sw=V_sw, V_back=V_back,
                          dynamic=dynamic, verbose=False)
    model.record = 10
    y0 = np.concatenate([np.zeros(sys_.N), psi_ss])
    geom = dict(W=W, overlap=overlap, L1=L1, L2=L2, dx=dx, x=x, seg=seg, n1=n1, vs=vs, BE=BE)
    return model, y0, geom


def classify(ev, geom, V_thr=1e-2, t_extra=3.0):
    t = np.array([s[0] for s in ev["snaps"]]); V = np.array([s[1] for s in ev["snaps"]])
    seg = geom["seg"]; dx = geom["dx"]
    out = dict(t_start=ev["t_start"], Vmax=ev["Vmax"])
    per = {}
    for k in (0, 1):
        m = seg == k
        active = (V[:, m] > V_thr).any(axis=1)
        if active.any():
            i0 = np.argmax(active); i1 = len(active) - 1 - np.argmax(active[::-1])
            per[k] = dict(t_on=t[i0], t_off=t[i1], Vmax=V[:, m].max(), slip=ev["slip"][m].sum() * dx,
                          extent=np.sum(ev["slip"][m] > 1e-2) * dx)
    out["segments"] = per
    if len(per) == 2:
        src = 0 if per[0]["t_on"] <= per[1]["t_on"] else 1
        rcv = 1 - src
        lag_after_arrest = per[rcv]["t_on"] - per[src]["t_off"]
        lag_after_onset = per[rcv]["t_on"] - per[src]["t_on"]
        t_dyn = abs(geom["W"]) / geom["BE"] + t_extra
        kind = "dynamic" if lag_after_arrest < t_dyn else "postseismic"
        if lag_after_onset < 0.5:
            kind = "simultaneous"
        out.update(kind=kind, source=src, receiver=rcv, lag_after_arrest=lag_after_arrest, lag_after_onset=lag_after_onset)
    else:
        out.update(kind="single", source=list(per.keys())[0] if per else None)
    return out


def run(n_events=4, t_end_yr=600.0, verbose=True, delayed_window_yr=1.0, **kw):
    model, y0, geom = build(**kw)
    t0 = time.time()
    y, t = y0, 0.0
    results = []
    for k in range(n_events):
        y, t = model.run(y, t_end=min(t + 400 * YR, t_end_yr * YR), max_events=1, t0=t)
        if len(model.events) <= k:
            break
        ev = model.events[k]
        c = classify(ev, geom)
        ev["snaps"] = []
        # a single-segment event shortly after one on the other segment: delayed triggering
        if c["kind"] == "single" and results:
            prev = results[-1]
            if prev["segments"] and c["source"] is not None and c["source"] not in prev["segments"] \
                    and c["t_start"] - prev["t_start"] < delayed_window_yr * YR:
                c["kind"] = f"delayed ({(c['t_start'] - prev['t_start']) / 86400:.0f} d)"
        results.append(c)
        if verbose:
            segs = ", ".join(f"seg{j}: on {v['t_on']-c['t_start']:.1f}s off {v['t_off']-c['t_start']:.1f}s slip {v['slip']/1e3:.1f} m km extent {v['extent']/1e3:.1f} km Vmax {v['Vmax']:.2f}" for j, v in c["segments"].items())
            lag = f", lag after arrest {c['lag_after_arrest']:+.2f} s" if "lag_after_arrest" in c else ""
            print(f"event {k+1} at {c['t_start']/YR:8.2f} yr [{c['kind']}{lag}]: {segs}  ({time.time()-t0:.0f}s)")
        if t >= t_end_yr * YR:
            break
    return results, geom


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--W", type=float, default=-1.0, help="step width in km, negative = right step (dilational)")
    p.add_argument("--overlap", type=float, default=0.0, help="overlap in km (negative = gap)")
    p.add_argument("--events", type=int, default=4)
    p.add_argument("--dx", type=float, default=100.0)
    p.add_argument("--qd", action="store_true")
    p.add_argument("--friction", choices=["aging", "flash"], default="aging")
    a = p.parse_args()
    run(n_events=a.events, W=a.W * 1e3, overlap=a.overlap * 1e3, dx=a.dx, dynamic=not a.qd, friction=a.friction)
