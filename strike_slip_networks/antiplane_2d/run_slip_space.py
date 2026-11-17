#!/usr/bin/env python3
"""
run_slip_space.py [options]: Cattania-setup runs with her S1 settings (as
run_cattania_si.py: events above 0.1 m/s, min(L_inf, L_b)/dx = 7.7 for one
state variable, creeping region of depth W, 30 recurrence estimates, no
series) over user-chosen ranges, meant to extend the slip-law parameter
space; all combinations of the given values are run in order of
estimated cost, completed runs are skipped.

  --laws aging slip two   evolution laws
  --b 0.02                b of the velocity-weakening part
  --ab 0.1 0.2 0.3        a/b values
  --X 16 32 64            W/L_inf values
  --rho 0.048             two state variables: dc1/dc2 (one or more)
  --b2b1 0.84             two state variables: b2/b1 (one or more)
  --maxcells 4000         skip larger runs
  --cycles 30             run length in recurrence estimates
  --np 1                  MPI ranks (mpirun from ../config.sh)
  --dir .                 run directory

Directory names si_<law>_b<b>_ab<a/b>_X<X>, with the law written as
two-r<rho>-q<b2b1> when the two-state ratios differ from 0.048 and 0.84,
so that plot_cattania_si.py and plot_slip_space.py read them.
"""
import argparse
import os
import subprocess

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--laws", nargs="+", default=["slip"])
p.add_argument("--b", nargs="+", type=float, default=[0.02])
p.add_argument("--ab", nargs="+", type=float, default=[0.1, 0.2, 0.3])
p.add_argument("--X", nargs="+", type=float, default=[16, 32, 64])
p.add_argument("--rho", nargs="+", type=float, default=[0.048])
p.add_argument("--b2b1", nargs="+", type=float, default=[0.84])
p.add_argument("--maxcells", type=int, default=4000)
p.add_argument("--cycles", type=float, default=30.0)
p.add_argument("--np", type=int, default=1)
p.add_argument("--dir", default=".")
a = p.parse_args()
sdir = os.path.dirname(os.path.abspath(__file__))
FAC = {"aging": 1.0, "slip": 2.0, "two": 4.0}
rsf_solve, mpirun = "rsf_solve", "mpirun"
cfg = os.path.join(sdir, "..", "config.sh")
if os.path.exists(cfg):
    for line in open(cfg):
        if line.startswith("rsf_solve="):
            rsf_solve = line.split("=", 1)[1].strip().strip('"')
        if line.startswith("mpirun="):
            mpirun = line.split("=", 1)[1].strip().strip('"')


def gen(law, b, ab, x, rho, q, out):
    cmd = ["python3", os.path.join(sdir, "gen_antiplane.py"), "cattania", "%g" % x, "--law", law, "--b", "%g" % b,
           "--a_over_b", "%g" % ab, "--rho", "%g" % rho, "--b2_b1", "%g" % q, "--out", out]
    if law != "two":
        cmd += ["--cmin", "7.7"]
    w = subprocess.run(cmd, capture_output=True, text=True).stdout.split()
    return int(w[w.index("VW") - 1]) + int(w[w.index("VS") - 1])


def done(d):
    try:
        return open(os.path.join(d, "rsf_solve.log")).read().splitlines()[-1].startswith("wall")
    except (OSError, IndexError):
        return False


os.makedirs(a.dir, exist_ok=True)
os.chdir(a.dir)
jobs = []
for law in a.laws:
    combos = [(r, q) for r in a.rho for q in a.b2b1] if law == "two" else [(0.048, 0.84)]
    for b in a.b:
        for ab in a.ab:
            for x in a.X:
                for rho, q in combos:
                    tag = law if (rho, q) == (0.048, 0.84) else "two-r%g-q%g" % (rho, q)
                    d = "si_%s_b%g_ab%g_X%g" % (tag, b, ab, x)
                    n = gen(law, b, ab, x, rho, q, "/tmp/slip_space_ncells")
                    if n <= a.maxcells:
                        jobs.append((n * FAC[law] * x ** 0.5, d, law, b, ab, x, rho, q, n))
                    else:
                        print("skipping %s (%i cells)" % (d, n), flush=True)
for cost, d, law, b, ab, x, rho, q, n in sorted(jobs):
    if done(d):
        continue
    gen(law, b, ab, x, rho, q, d)
    t = float(open(os.path.join(d, "params.txt")).read().split("t_est_yr")[1].split()[0])
    opts = open(os.path.join(d, "run.opts")).read().split()
    opts += ["-stop_time_yr", repr(a.cycles * t), "-rtol", "1e-5", "-ts_rk_type", "3bs", "-dt_monitor_yr",
             repr(t / 20), "-rdx_monitor", "0.1", "-rsf_catalog", "-vel_event", "0.1", "-rsf_stations", "stations.in",
             "-rsf_station_dt_yr", repr(t / 500), "-rsf_station_vdense", "1e-3", "-rsf_checkpoint_wall", "120"]
    cmd = ["/usr/bin/time", "-f", "wall %e s"]
    if a.np > 1:
        cmd += mpirun.split() + ["-np", str(a.np)]
    with open(os.path.join(d, "rsf_solve.log"), "w") as log:
        subprocess.run(cmd + [rsf_solve] + opts, cwd=d, stdout=log, stderr=subprocess.STDOUT)
    last = open(os.path.join(d, "rsf_solve.log")).read().splitlines()[-1]
    print("%s (%i cells): %s" % (d, n, last), flush=True)
print("SLIP_SPACE_DONE", flush=True)
