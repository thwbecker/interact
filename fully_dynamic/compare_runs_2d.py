#!/usr/bin/env python3
"""Compare quasi-dynamic (QD) and fully dynamic (FD) cycle runs of the same section: the events
of each run (rsf_events.dat: onsets and arrests at |v| = vel_event), their event slip (mean slip
at arrest minus at onset), recurrence intervals, and the FD phases (run.log), and the QD against
the FD events in sequence: onset time differences (in years and in percent of the QD
recurrence) and slip differences. Optionally a figure of the mean slip against time
(rsf_monitor.dat) with the events marked.

Usage: compare_runs_2d.py outdir [qd_run] [fd_run ...] [--plot file.png] [--merge seconds]
  onsets that follow the previous arrest within --merge seconds (default one day: the slip rate
  flickering about the threshold, or short dynamic phases after the main one) count as one event
  outdir: a run_2d_cascadia.sh output directory; runs default to qd and fd_omp (and fd_mpi if
  present). Writes the comparison to stdout."""
import sys
import os
import numpy as np


MERGE_S = 86400.0          # sub-events (re-crossings of the threshold) closer than this are merged


def merge(ev, gap_s):
    """events whose onset follows the previous arrest within gap_s seconds are one event"""
    out = []
    for e in ev:
        if out and not np.isnan(out[-1][2]) and (e[0] - out[-1][2])*3.15576e7 < gap_s:
            out[-1][2] = e[2]
            out[-1][3] = e[3]
        else:
            out.append(list(e))
    return out


def events(d):
    """[(t_onset_yr, slip_onset, t_arrest_yr or nan, slip_arrest or nan, log10 vmax at onset)],
    sub-events within MERGE_S merged"""
    f = os.path.join(d, "rsf_events.dat")
    if not os.path.exists(f):
        return None
    a = np.loadtxt(f, comments="#", ndmin=2)
    ev = []
    cur = None
    for r in a:
        if r[2] == 1:
            if cur is not None:
                ev.append(cur)
            cur = [r[1], r[4], np.nan, np.nan, r[3]]
        elif r[2] == -1 and cur is not None:
            cur[2] = r[1]
            cur[3] = r[4]
            ev.append(cur)
            cur = None
    if cur is not None:
        ev.append(cur)
    return merge(ev, MERGE_S)


def phases(d):
    """FD phases from run.log: [(t_yr, steps, duration_s, peak, mean_slip)]"""
    f = os.path.join(d, "run.log")
    ph = []
    if not os.path.exists(f):
        return ph
    for line in open(f):
        if line.startswith("rsf_fd: phase"):
            w = line.replace(",", " ").split()
            try:
                t = float(w[w.index("t") + 1])
                st = int(w[w.index("steps") - 1])
                du = float(w[w.index("over") + 1])
                pk = float(w[w.index("rate") + 1])
                ms = float(w[w.index("mean") + 2])
                ph.append((t, st, du, pk, ms))
            except (ValueError, IndexError):
                pass
    return ph


def describe(name, ev, ph):
    if ev is None:
        print(f"{name}: no rsf_events.dat")
        return
    print(f"{name}: {len(ev)} events")
    for k, e in enumerate(ev):
        rec = (e[0] - ev[k - 1][0]) if k > 0 else np.nan
        line = f"  {k + 1:3d}: onset {e[0]:10.3f} yr, slip {e[3] - e[1]:8.3f} m, duration {(e[2] - e[0])*3.15576e7:10.1f} s"
        if k > 0:
            line += f", after {rec:8.2f} yr"
        print(line)
    for k, p in enumerate(ph):
        print(f"  phase {k + 1}: t {p[0]:10.3f} yr, {p[1]} steps over {p[2]:.1f} s, peak {p[3]:.3f} m/s, phase slip {p[4]:.3f} m")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    png = None
    if "--plot" in sys.argv:
        png = sys.argv[sys.argv.index("--plot") + 1]
        args = [a for a in args if a != png]
    global MERGE_S
    if "--merge" in sys.argv:
        MERGE_S = float(sys.argv[sys.argv.index("--merge") + 1])
        args = [a for a in args if a != sys.argv[sys.argv.index("--merge") + 1]]
    out = args[0]
    qd = args[1] if len(args) > 1 else "qd"
    fds = args[2:] if len(args) > 2 else [r for r in ("fd_omp", "fd_mpi") if os.path.isdir(os.path.join(out, r))]
    eq = events(os.path.join(out, qd))
    describe(qd, eq, [])
    for fd in fds:
        ef = events(os.path.join(out, fd))
        describe(fd, ef, phases(os.path.join(out, fd)))
        if not eq or not ef:
            continue
        n = min(len(eq), len(ef))
        print(f"{fd} against {qd}, events in sequence (FD minus QD):")
        for k in range(n):
            dt = ef[k][0] - eq[k][0]
            rq = (eq[k][0] - eq[k - 1][0]) if k > 0 else np.nan
            sq = eq[k][3] - eq[k][1]
            sf = ef[k][3] - ef[k][1]
            pct = f"{100.0*dt/rq:+6.2f} % of the QD recurrence" if k > 0 else "(first event)"
            print(f"  {k + 1:3d}: onset {dt:+9.3f} yr {pct}, slip {sf - sq:+7.3f} m ({100.0*(sf - sq)/sq:+6.2f} %)")
        if len(eq) != len(ef):
            print(f"  ({len(eq)} QD and {len(ef)} FD events: matched in sequence up to {n})")
    if png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 4))
        for r, c in [(qd, "k")] + [(f, col) for f, col in zip(fds, ("r", "b"))]:
            f = os.path.join(out, r, "rsf_monitor.dat")
            if os.path.exists(f):
                m = np.loadtxt(f, comments="#", ndmin=2)
                ax.plot(m[:, 2], m[:, 5], c, lw=1, label=r)
            ev = events(os.path.join(out, r)) or []
            for e in ev:
                ax.axvline(e[0], color=c, lw=0.5, ls=":")
        ax.set_xlabel("time [yr]")
        ax.set_ylabel("mean slip [m]")
        ax.legend()
        fig.tight_layout()
        fig.savefig(png, dpi=150)
        print(f"figure: {png}")


if __name__ == "__main__":
    main()
