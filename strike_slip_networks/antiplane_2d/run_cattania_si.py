#!/usr/bin/env python3
"""
run_cattania_si.py [rundir] [maxcells] [ncycles] [laws ...]: the Cattania
(2019) setup with the settings of her supporting text S1: earthquakes
above 0.1 m/s (-vel_event 0.1), cells with min(L_inf, L_b)/dx = 7.7
(gen_antiplane.py --cmin 7.7; two state variables keep their cell rule),
b = 0.02 and 0.01, a/b = 0.3, 0.5, 0.625, 0.75 and W/L_inf = 1 ... 32,
for the aging, slip and two-state slip laws (default all three). Runs
are made in order of estimated cost (cells, times 2 for the slip law and
4 for two state variables), skipping runs above maxcells cells (default
1300) and completed runs; no all-cell series. Directory names
si_<law>_b<b>_ab<a/b>_X<X>. rsf_solve is taken from ../config.sh if
present, else from the PATH.
"""
import os
import subprocess
import sys

rundir = sys.argv[1] if len(sys.argv) > 1 else "."
maxcells = int(sys.argv[2]) if len(sys.argv) > 2 else 1300
ncyc = float(sys.argv[3]) if len(sys.argv) > 3 else 30.0
laws = sys.argv[4:] if len(sys.argv) > 4 else ["aging", "slip", "two"]
sdir = os.path.dirname(os.path.abspath(__file__))
BS = (0.02, 0.01)
ABS = (0.3, 0.5, 0.625, 0.75)
XS = (1, 2, 4, 8, 16, 32)
FAC = {"aging": 1.0, "slip": 2.0, "two": 4.0}
rsf_solve = "rsf_solve"
cfg = os.path.join(sdir, "..", "config.sh")
if os.path.exists(cfg):
    for line in open(cfg):
        if line.startswith("rsf_solve="):
            rsf_solve = line.split("=", 1)[1].strip().strip('"')


def gen(law, b, ab, x, out):
    cmd = ["python3", os.path.join(sdir, "gen_antiplane.py"), "cattania", str(x), "--law", law, "--b", str(b),
           "--a_over_b", str(ab), "--out", out]
    if law != "two":
        cmd += ["--cmin", "7.7"]
    w = subprocess.run(cmd, capture_output=True, text=True).stdout.split()
    return int(w[w.index("VW") - 1]) + int(w[w.index("VS") - 1])


def done(d):
    """finished and exited normally (a failed or killed run also ends with
    the "wall" line of /usr/bin/time, after a non-zero status or a signal)"""
    try:
        tail = open(os.path.join(d, "rsf_solve.log")).read().splitlines()[-3:]
    except OSError:
        return False
    if not tail or not tail[-1].startswith("wall"):
        return False
    return not any("non-zero status" in t or "terminated by signal" in t for t in tail)


os.makedirs(rundir, exist_ok=True)
os.chdir(rundir)
jobs = []
for law in laws:
    for b in BS:
        for ab in ABS:
            for x in XS:
                n = gen(law, b, ab, x, "/tmp/si_ncells")
                if n <= maxcells:
                    jobs.append((n * FAC[law] * x ** 0.5, law, b, ab, x, n))
for cost, law, b, ab, x, n in sorted(jobs):
    d = "si_%s_b%g_ab%g_X%g" % (law, b, ab, x)
    if done(d):
        continue
    gen(law, b, ab, x, d)
    t = float(open(os.path.join(d, "params.txt")).read().split("t_est_yr")[1].split()[0])
    opts = open(os.path.join(d, "run.opts")).read().split()
    opts += ["-stop_time_yr", repr(ncyc * t), "-rtol", "1e-5", "-ts_rk_type", "3bs", "-dt_monitor_yr", repr(t / 20),
             "-rdx_monitor", "0.1", "-rsf_catalog", "-vel_event", "0.1", "-rsf_stations", "stations.in",
             "-rsf_station_dt_yr", repr(t / 500), "-rsf_station_vdense", "1e-3", "-rsf_checkpoint_wall", "120"]
    with open(os.path.join(d, "rsf_solve.log"), "w") as log:
        subprocess.run(["/usr/bin/time", "-f", "wall %e s", rsf_solve] + opts, cwd=d, stdout=log,
                       stderr=subprocess.STDOUT)
    tail = open(os.path.join(d, "rsf_solve.log")).read().splitlines()
    last = " | ".join(tail[-2:]) if tail else "(empty log)"
    print("%s (%i cells): %s" % (d, n, last), flush=True)
print("CATTANIA_SI_DONE", flush=True)
