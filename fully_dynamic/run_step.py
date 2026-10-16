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
import argparse, time, warnings
import numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning, module="scipy.integrate")
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
    """Per segment: onset (first V > V_thr anywhere on it), arrival of the
    rupture front at the end nearest the step (first V > V_thr on that end
    element), last slip above V_thr, slip, extent, peak slip rate. For events
    involving both segments the receiver onset is referred to the source
    front's arrival at the step:
        joint       : receiver starts before the source front reaches the step
        dynamic     : within t_dyn = |W|/beta + t_extra after that arrival
        postseismic : later in the same FD phase
    For single-segment events the other segment's peak slip rate during the
    rupture (onset to arrival + t_dyn) is recorded as rcv_Vpeak."""
    t = np.array([s_[0] for s_ in ev["snaps"]]); V = np.array([s_[1] for s_ in ev["snaps"]])
    seg = geom["seg"]; dx = geom["dx"]; x = geom["x"]
    # element of each segment nearest the step: seg0 ends at x = 0, seg1 starts at x = -overlap
    end_el = {0: np.where(seg == 0)[0][np.argmax(x[seg == 0])], 1: np.where(seg == 1)[0][np.argmin(x[seg == 1])]}
    out = dict(t_start=ev["t_start"], Vmax=ev["Vmax"])
    per = {}
    for k in (0, 1):
        m = seg == k
        active = (V[:, m] > V_thr).any(axis=1)
        if active.any():
            i0 = np.argmax(active); i1 = len(active) - 1 - np.argmax(active[::-1])
            at_end = V[:, end_el[k]] > V_thr
            t_end = t[np.argmax(at_end)] if at_end.any() else t[i1]
            per[k] = dict(t_on=t[i0], t_end=t_end, t_off=t[i1], Vmax=V[:, m].max(),
                          slip=ev["slip"][m].sum() * dx, extent=np.sum(ev["slip"][m] > 1e-2) * dx)
    out["segments"] = per
    t_dyn = abs(geom["W"]) / geom["BE"] + t_extra
    if len(per) == 2:
        src = 0 if per[0]["t_on"] <= per[1]["t_on"] else 1
        rcv = 1 - src
        lag = per[rcv]["t_on"] - per[src]["t_end"]
        frac = per[rcv]["slip"] / max(per[src]["slip"], 1e-9)
        if frac < 0.1:
            kind = f"single, minor re-slip of the other segment ({frac*100:.0f}% of the source slip)"
        else:
            kind = "joint" if lag < 0 else ("immediate" if lag < 60.0 else "minutes" if lag < 3600.0 else "postseismic")
        out.update(kind=kind, source=src, receiver=rcv, lag_after_arrival=lag,
                   lag_after_onset=per[rcv]["t_on"] - per[src]["t_on"])
    elif len(per) == 1:
        src = list(per.keys())[0]; rcv = 1 - src
        win = (t >= per[src]["t_on"]) & (t <= per[src]["t_end"] + t_dyn)
        out.update(kind="single", source=src, rcv_Vpeak=float(V[win][:, seg == rcv].max()) if win.any() else 0.0)
    else:
        out.update(kind="none")
    return out


def _fmt_delay(seconds):
    return f"{seconds/3600:.0f} h" if seconds < 2 * 86400 else f"{seconds/86400:.0f} d"


def run(n_events=4, t_end_yr=600.0, verbose=True, delayed_window_yr=1.0, out=None, **kw):
    model, y0, geom = build(**kw)
    model.record = 5
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
                c["kind"] = f"delayed ({_fmt_delay(c['t_start'] - prev['t_start'])})"
        results.append(c)
        if verbose:
            segs = ", ".join(f"seg{j}: on {v['t_on']-c['t_start']:.1f}s off {v['t_off']-c['t_start']:.1f}s slip {v['slip']/1e3:.1f} m km extent {v['extent']/1e3:.1f} km Vmax {v['Vmax']:.2f}" for j, v in c["segments"].items())
            lag = f", receiver onset {c['lag_after_arrival']:+.1f} s after the source front reached the step" if "lag_after_arrival" in c else (f", other segment peak V during rupture {c['rcv_Vpeak']:.1e} m/s" if "rcv_Vpeak" in c else "")
            print(f"event {k+1} at {c['t_start']/YR:8.2f} yr [{c['kind']}{lag}]: {segs}  ({time.time()-t0:.0f}s)")
        if out is not None:
            save(out, results, geom, kw)
        if t >= t_end_yr * YR:
            break
    return results, geom


def save(out, results, geom, kw):
    """One npz per configuration: event table and per-event, per-segment numbers."""
    import json
    rows = []
    for c in results:
        for j, v in c["segments"].items():
            rows.append((c["t_start"], j, v["t_on"], v["t_off"], v["slip"], v["extent"], v["Vmax"], v["t_end"]))
    kinds = [c["kind"] for c in results]
    lags = [c.get("lag_after_arrival", np.nan) for c in results]
    rcvV = [c.get("rcv_Vpeak", np.nan) for c in results]
    np.savez(out, rows=np.array(rows), kinds=np.array(kinds), lags=np.array(lags), rcv_Vpeak=np.array(rcvV),
             t_start=np.array([c["t_start"] for c in results]), Vmax=np.array([c["Vmax"] for c in results]),
             W=geom["W"], overlap=geom["overlap"], dx=geom["dx"], params=json.dumps({k: (v if isinstance(v, (int, float, str, bool)) else str(v)) for k, v in kw.items()}))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--W", type=float, default=-1.0, help="step width in km, negative = right step (dilational)")
    p.add_argument("--overlap", type=float, default=0.0, help="overlap in km (negative = gap)")
    p.add_argument("--events", type=int, default=4)
    p.add_argument("--dx", type=float, default=100.0)
    p.add_argument("--qd", action="store_true")
    p.add_argument("--friction", choices=["aging", "flash"], default="aging")
    p.add_argument("--vs-end", type=float, default=4.0, help="velocity-strengthening margin at each segment end (km)")
    p.add_argument("--tol", type=float, default=1e-2, help="window tolerance")
    p.add_argument("--out", default=None, help="npz file for the results (default: step_<mode>_W<km>_ov<km>_dx<m>.npz)")
    a = p.parse_args()
    out = a.out or f"step_{'qd' if a.qd else 'fd'}_W{a.W:+.1f}_ov{a.overlap:.1f}_dx{int(a.dx)}.npz"
    run(n_events=a.events, W=a.W * 1e3, overlap=a.overlap * 1e3, dx=a.dx, dynamic=not a.qd, friction=a.friction,
        vs_end=a.vs_end * 1e3, tol=a.tol, out=out)
