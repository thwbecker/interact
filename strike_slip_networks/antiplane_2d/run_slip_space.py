#!/usr/bin/env python3
"""
run_slip_space.py [options]: Cattania-setup runs with her S1 settings (as
run_cattania_si.py: events above 0.1 m/s, min(L_inf, L_b)/dx = 7.7 for one
state variable, creeping region of depth W, no all-cell series), either
over user-chosen ranges or for the predefined final study (--plan final);
all runs are made in order of estimated cost.

  --plan grid|final       grid (default): all combinations of the values
                          below; final: the set of runs of the slip-law
                          study (see below)
  --laws aging slip two   evolution laws
  --b 0.02                b of the velocity-weakening part
  --ab 0.1 0.2 0.3        a/b values
  --X 16 32 64            W/L_inf values
  --rho 0.048             two state variables: dc1/dc2 (one or more)
  --b2b1 0.84             two state variables: b2/b1 (one or more)
  --maxcells 6000         skip larger runs
  --cycles 30             run length in recurrence estimates
  --np 1                  MPI ranks (mpirun from ../config.sh)
  --dir .                 run directory
  --dry                   only list what would be done
  --checkpoint_wall 120   seconds of wall time between checkpoints

Before running, each existing run directory is classified (runs whose
directory is the working directory of a running rsf_solve are left
alone):

  finished  log ends with the /usr/bin/time "wall" line of a normal exit:
            skipped
  killed    stopped by a signal or without a "wall" line, with a
            checkpoint: continued from rsf_checkpoint.bin (the earlier
            log is kept as rsf_solve.log.<n>)
  failed    non-zero exit status, empty log, or killed without a
            checkpoint: deleted and run again from the start

A run that fails again in this invocation is reported and not retried.

The final plan (about 110 runs; the largest have about 5000 cells):
  slip law,   b 0.02, a/b 0.1 0.2 0.3 0.4 0.5, W/L_inf 16 ... 256
  slip law,   b 0.02, a/b 0.625, W/L_inf 16 ... 128; a/b 0.75, 16 ... 64
  slip law,   b 0.01, a/b 0.3 0.5, W/L_inf 16 ... 256
  aging law,  b 0.02, a/b 0.1 0.2 0.3 0.4 0.5, W/L_inf 16 ... 128
              (reference at the same settings)
  two states, b 0.02, a/b 0.3 0.5, W/L_inf 16 ... 128 (dc1/dc2 0.048,
              b2/b1 0.84)
  two states, b 0.02, a/b 0.5, W/L_inf 4 ... 32, dc1/dc2 0.048 0.2 0.5,
              b2/b1 0.5 0.84 2

Directory names si_<law>_b<b>_ab<a/b>_X<X>, with the law written as
two-r<rho>-q<b2b1> when the two-state ratios differ from 0.048 and 0.84,
so that plot_cattania_si.py and plot_slip_space.py read them.
"""
import argparse
import glob
import os
import shutil
import subprocess

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--plan", choices=("grid", "final"), default="grid")
p.add_argument("--laws", nargs="+", default=["slip"])
p.add_argument("--b", nargs="+", type=float, default=[0.02])
p.add_argument("--ab", nargs="+", type=float, default=[0.1, 0.2, 0.3])
p.add_argument("--X", nargs="+", type=float, default=[16, 32, 64])
p.add_argument("--rho", nargs="+", type=float, default=[0.048])
p.add_argument("--b2b1", nargs="+", type=float, default=[0.84])
p.add_argument("--maxcells", type=int, default=6000)
p.add_argument("--cycles", type=float, default=30.0)
p.add_argument("--np", type=int, default=1)
p.add_argument("--dir", default=".")
p.add_argument("--dry", action="store_true")
p.add_argument("--checkpoint_wall", type=float, default=120.0, help="seconds between checkpoints")
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


def combos(laws, bs, abs_, xs, rhos=(0.048,), qs=(0.84,)):
    out = []
    for law in laws:
        rq = [(r, q) for r in rhos for q in qs] if law == "two" else [(0.048, 0.84)]
        for b in bs:
            for ab in abs_:
                for x in xs:
                    for rho, q in rq:
                        out.append((law, b, ab, x, rho, q))
    return out


def plan():
    if a.plan == "grid":
        return combos(a.laws, a.b, a.ab, a.X, a.rho, a.b2b1)
    big = (16, 32, 64, 128, 256)
    jobs = combos(["slip"], [0.02], [0.1, 0.2, 0.3, 0.4, 0.5], big)
    jobs += combos(["slip"], [0.02], [0.625], big[:4]) + combos(["slip"], [0.02], [0.75], big[:3])
    jobs += combos(["slip"], [0.01], [0.3, 0.5], big)
    jobs += combos(["aging"], [0.02], [0.1, 0.2, 0.3, 0.4, 0.5], big[:4])
    jobs += combos(["two"], [0.02], [0.3, 0.5], big[:4])
    jobs += combos(["two"], [0.02], [0.5], (4, 8, 16, 32), (0.048, 0.2, 0.5), (0.5, 0.84, 2.0))
    return sorted(set(jobs))


def name(law, b, ab, x, rho, q):
    tag = law if (law != "two" or (rho, q) == (0.048, 0.84)) else "two-r%g-q%g" % (rho, q)
    return "si_%s_b%g_ab%g_X%g" % (tag, b, ab, x)


def gen(law, b, ab, x, rho, q, out):
    cmd = ["python3", os.path.join(sdir, "gen_antiplane.py"), "cattania", "%g" % x, "--law", law, "--b", "%g" % b,
           "--a_over_b", "%g" % ab, "--rho", "%g" % rho, "--b2_b1", "%g" % q, "--out", out]
    if law != "two":
        cmd += ["--cmin", "7.7"]
    w = subprocess.run(cmd, capture_output=True, text=True).stdout.split()
    return int(w[w.index("VW") - 1]) + int(w[w.index("VS") - 1])


def running_dirs():
    out = set()
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            cmd = open("/proc/%s/cmdline" % pid, "rb").read().split(b"\0")
            if cmd and os.path.basename(cmd[0].decode(errors="ignore")) == "rsf_solve":
                out.add(os.path.realpath(os.readlink("/proc/%s/cwd" % pid)))
        except OSError:
            continue
    return out


def status(d):
    """finished, killed, failed or new"""
    log = os.path.join(d, "rsf_solve.log")
    if not os.path.isdir(d):
        return "new"
    try:
        lines = open(log).read().splitlines()
    except OSError:
        return "failed"
    if not lines:
        return "failed"
    tail = lines[-3:]
    if any("non-zero status" in t for t in tail):
        return "failed"
    if tail[-1].startswith("wall") and not any("terminated by signal" in t for t in tail):
        return "finished"
    return "killed" if os.path.exists(os.path.join(d, "rsf_checkpoint.bin")) else "failed"


def run(d, job, restart):
    law, b, ab, x, rho, q = job
    gen(law, b, ab, x, rho, q, d)
    t = float(open(os.path.join(d, "params.txt")).read().split("t_est_yr")[1].split()[0])
    opts = open(os.path.join(d, "run.opts")).read().split()
    opts += ["-stop_time_yr", repr(a.cycles * t), "-rtol", "1e-5", "-ts_rk_type", "3bs", "-dt_monitor_yr",
             repr(t / 20), "-rdx_monitor", "0.1", "-rsf_catalog", "-vel_event", "0.1", "-rsf_stations", "stations.in",
             "-rsf_station_dt_yr", repr(t / 500), "-rsf_station_vdense", "1e-3", "-rsf_checkpoint_wall", repr(a.checkpoint_wall)]
    if restart:
        opts += ["-rsf_restart", "rsf_checkpoint.bin"]
    cmd = ["/usr/bin/time", "-f", "wall %e s"]
    if a.np > 1:
        cmd += mpirun.split() + ["-np", str(a.np)]
    with open(os.path.join(d, "rsf_solve.log"), "w") as log:
        subprocess.run(cmd + [rsf_solve] + opts, cwd=d, stdout=log, stderr=subprocess.STDOUT)
    tail = open(os.path.join(d, "rsf_solve.log")).read().splitlines()
    return " | ".join(tail[-2:]) if tail else "(empty log)"


os.makedirs(a.dir, exist_ok=True)
os.chdir(a.dir)
busy = running_dirs()
todo = []
for job in plan():
    d = name(*job)
    if os.path.realpath(d) in busy:
        print("%-36s running elsewhere, left alone" % d, flush=True)
        continue
    st = status(d)
    if st == "finished":
        continue
    n = gen(*job, "/tmp/slip_space_ncells_%i" % os.getpid())
    if n > a.maxcells:
        print("%-36s skipped (%i cells)" % (d, n), flush=True)
        continue
    todo.append((n * FAC[job[0]] * job[3] ** 0.5, d, job, st, n))
todo.sort()
for cost, d, job, st, n in todo:
    print("%-36s %-8s %5i cells -> %s" % (d, st, n, "continue from checkpoint" if st == "killed" else "run"),
          flush=True)
if a.dry:
    raise SystemExit
for cost, d, job, st, n in todo:
    if st == "killed":
        k = len(glob.glob(os.path.join(d, "rsf_solve.log.*")))
        shutil.copy(os.path.join(d, "rsf_solve.log"), os.path.join(d, "rsf_solve.log.%i" % (k + 1)))
        msg = run(d, job, True)
    else:
        if st == "failed":
            shutil.rmtree(d)
        msg = run(d, job, False)
    st2 = status(d)
    print("%-36s %s: %s" % (d, "done" if st2 == "finished" else "FAILED", msg), flush=True)
print("SLIP_SPACE_DONE", flush=True)
