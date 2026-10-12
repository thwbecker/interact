"""Run SEAS BP1 (QD or FD) and record station series at 2.5, 7.5, 12.5 km."""
import numpy as np, sys, time, pickle
import bp1
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="scipy.integrate")

YR = 3.15576e7

def run(dynamic, dz=50.0, t_end_yr=520.0, lmax=6000, V_sw=1e-3, tag=None, resume=False):
    model, y0, z, fr = bp1.build(dz=dz, dynamic=dynamic, lmax=lmax, V_sw=V_sw)
    model.record = 5
    idx = {d: int(np.argmin(np.abs(z - d))) for d in (2.5e3, 7.5e3, 12.5e3)}
    tag = tag or ("fd" if dynamic else "qd")
    y, t = y0, 0.0
    events = []
    import os
    if resume and os.path.exists(f"bp1_{tag}.pkl"):
        d = pickle.load(open(f"bp1_{tag}.pkl", "rb"))
        events = d["events"]; y = d["y"]; t = d["t"]
        print(f"[{tag}] resuming after event {len(events)} at t = {t/YR:.2f} yr", flush=True)
    t0 = time.time()
    n_before = len(events)
    while t < t_end_yr * YR:
        y, t = model.run(y, t_end=t_end_yr * YR, max_events=1, t0=t)
        if len(model.events) <= len(events) - n_before:
            break
        e = model.events[-1]
        snaps = e["snaps"]; e["snaps"] = []
        st = np.array([s[0] for s in snaps]); Vs = np.array([s[1] for s in snaps]); taus = np.array([s[2] for s in snaps])
        rec = {d: dict(t=st, V=Vs[:, i], tau=taus[:, i]) for d, i in idx.items()}
        events.append(dict(t_start=e["t_start"], t_end=e["t_end"], Vmax=e["Vmax"], slip=e["slip"], stations=rec, nsteps=e["nsteps"]))
        print(f"[{tag}] event {len(events)} at {e['t_start']/YR:8.2f} yr: Vmax {e['Vmax']:.2f} m/s, max slip {e['slip'].max():.2f} m, "
              f"slip at 7.5 km {e['slip'][idx[7.5e3]]:.2f} m, {e['nsteps']} steps ({time.time()-t0:.0f} s)", flush=True)
        pickle.dump(dict(events=events, z=z, log=model.log, y=y, t=t), open(f"bp1_{tag}.pkl", "wb"))
    return events

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="SEAS BP1 with the space-time BIEM prototype")
    p.add_argument("mode", choices=["qd", "fd"])
    p.add_argument("--dz", type=float, default=50.0, help="element size (m); the benchmark suggests 25")
    p.add_argument("--tend", type=float, default=520.0, help="end time (yr)")
    p.add_argument("--lmax", type=int, default=6000, help="cap on the history window (fine-grid lags); 0 = no cap")
    p.add_argument("--vsw", type=float, default=1e-3, help="QD to FD switch slip rate (m/s)")
    p.add_argument("--resume", action="store_true", help="continue from bp1_<tag>.pkl")
    p.add_argument("--tag", default=None, help="output tag (default: mode_dz<dz>)")
    a = p.parse_args()
    tag = a.tag or f"{a.mode}_dz{int(a.dz)}"
    run(a.mode == "fd", dz=a.dz, t_end_yr=a.tend, lmax=(a.lmax if a.lmax > 0 else None), V_sw=a.vsw, tag=tag, resume=a.resume)
