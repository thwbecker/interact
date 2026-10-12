"""Compare our BP1 QD/FD sequences (bp1_qd.pkl, bp1_fd.pkl) with the uploaded
reference series (ref/bp1qd_*.npy: groups 3 and 4; ref/bp1fd_*.npy: groups 1
and 2) at 2.5, 7.5 and 12.5 km depth."""
import numpy as np, pickle, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

YR = 3.15576e7
DEPTHS = (2.5e3, 7.5e3, 12.5e3)
REFS = {"QD": {"Jiang SBIM": "bp1qd_3", "group 4": "bp1qd_4"},
        "FD": {"Jiang 25 m (file 1)": "bp1fd_1", "Yang SCycle 25 m": "bp1fd_2",
               "sem2dpack 50 m": "bp1fd_5", "Yang SCycle 50 m": "bp1fd_6"}}


def ref_events(a, V_thr=1e-3):
    t = a[:, 0]; lv = a[:, 2]; seis = lv > np.log10(V_thr)
    starts = t[1:][seis[1:] & ~seis[:-1]]
    merged = []
    for s_ in starts:
        if not merged or s_ - merged[-1] > 300.0:
            merged.append(s_)
    starts = np.array(merged)
    out = []
    for s in starts:
        m = (t >= s - 2) & (t <= s + 120)
        i0 = np.argmax(t >= s - 2); i1 = min(len(t) - 1, np.argmax(t >= s + 120) if (t >= s + 120).any() else len(t) - 1)
        out.append(dict(t=s, Vmax=10 ** lv[m].max(), slip=a[i1, 1] - a[i0, 1], taumin=a[m, 3].min(), taumax=a[m, 3].max()))
    return out


def summarize():
    rows = []
    for mode in ("QD", "FD"):
        fn = f"bp1_{mode.lower()}.pkl"
        if os.path.exists(fn):
            d = pickle.load(open(fn, "rb")); z = d["z"]
            ev = d["events"]
            times = np.array([e["t_start"] for e in ev]) / YR
            rows.append((f"ours {mode} (50 m)", times, {dep: [e["stations"][dep]["V"].max() for e in ev] for dep in DEPTHS},
                         {dep: [e["slip"][int(np.argmin(np.abs(z - dep)))] for e in ev] for dep in DEPTHS}))
        for name, pref in REFS[mode].items():
            try:
                arrs = {dep: np.load(f"ref/{pref}_{k}.npy") for dep, k in zip(DEPTHS, (1, 2, 3))}
            except FileNotFoundError:
                continue
            evs = {dep: ref_events(a) for dep, a in arrs.items()}
            times = np.array([e["t"] for e in evs[2.5e3]]) / YR
            rows.append((f"{name} {mode}", times, {dep: [e["Vmax"] for e in evs[dep]] for dep in DEPTHS},
                         {dep: [e["slip"] for e in evs[dep]] for dep in DEPTHS}))
    for name, times, pk, sl in rows:
        n = min(5, len(times))
        print(f"{name:22s} events (yr): {np.round(times[:n], 2)}  recurrence {np.round(np.diff(times[:n]), 1)}")
        for dep in DEPTHS:
            print(f"    {dep/1e3:4.1f} km: peak V {np.round(pk[dep][:n], 2)} m/s, coseismic slip {np.round(sl[dep][:n], 2)} m")
    return rows


def plot(event_index=1, out_png="bp1_compare.png"):
    fig, ax = plt.subplots(2, 3, figsize=(14, 8), sharex="col")
    for r, mode in enumerate(("QD", "FD")):
        fn = f"bp1_{mode.lower()}.pkl"
        ours = pickle.load(open(fn, "rb")) if os.path.exists(fn) else None
        for c, dep in enumerate(DEPTHS):
            a = ax[r, c]
            if ours is not None and len(ours["events"]) > event_index:
                e = ours["events"][event_index]; s = e["stations"][dep]
                i0 = np.argmax(s["V"] > 1e-3)
                a.plot(s["t"] - s["t"][i0], np.log10(np.maximum(s["V"], 1e-12)), "k-", lw=1.8, label=f"ours {mode} 50 m")
            for name, pref in REFS[mode].items():
                try:
                    arr = np.load(f"ref/{pref}_{c + 1}.npy")
                except FileNotFoundError:
                    continue
                evs = ref_events(arr)
                if len(evs) > event_index:
                    t0 = evs[event_index]["t"]
                    m = (arr[:, 0] > t0 - 5) & (arr[:, 0] < t0 + 40)
                    a.plot(arr[m, 0] - t0, arr[m, 2], lw=1.0, label=name)
            a.set_ylim(-4, 1.3); a.set_xlim(-5, 40)
            a.set_title(f"{mode}: {dep/1e3:.1f} km depth, event {event_index + 1}")
            a.set_xlabel("time after slip rate exceeds 1e-3 m/s at the station (s)")
            if c == 0:
                a.set_ylabel("log10 slip rate (m/s)")
            a.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(out_png, dpi=110); print("figure:", out_png)


if __name__ == "__main__":
    summarize()
    plot()
