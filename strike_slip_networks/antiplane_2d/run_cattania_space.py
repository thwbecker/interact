#!/usr/bin/env python3
"""
run_cattania_space.py [rundir] [maxcells] [ncycles] [laws ...] [--series]: the
Cattania (2019) setup over W/L_inf = 2 ... 64 and a/b = 0.5, 0.625, 0.75
for the given laws (default aging slip two), run in order of estimated
cost (cells, times 2 for the slip law and 4 for two state variables);
runs with more than maxcells cells (default 800) are skipped, completed
runs are skipped, so the script can be restarted. Directory names as in
run_parameter_maps (cattania_<law>_ab<a/b>_X<X>). The all-cell series is
not written unless --series is given (plot_cattania_space.py needs only
the catalogs).
"""
import os
import subprocess
import sys

series = "1" if "--series" in sys.argv else "0"
sys.argv = [v for v in sys.argv if v != "--series"]
rundir = sys.argv[1] if len(sys.argv) > 1 else "."
maxcells = int(sys.argv[2]) if len(sys.argv) > 2 else 800
ncyc = sys.argv[3] if len(sys.argv) > 3 else "30"
laws = sys.argv[4:] if len(sys.argv) > 4 else ["aging", "slip", "two"]
sdir = os.path.dirname(os.path.abspath(__file__))
XS = (2, 3, 4, 6, 8, 11, 16, 23, 32, 45, 64)
ABS = (0.5, 0.625, 0.75)
FAC = {"aging": 1.0, "slip": 2.0, "two": 4.0, "prz": 1.5}


def ncells(law, x, ab):
    out = subprocess.run(["python3", os.path.join(sdir, "gen_antiplane.py"), "cattania", str(x), "--law", law,
                          "--a_over_b", str(ab), "--out", "/tmp/ncells_cs"], capture_output=True, text=True).stdout
    w = out.split()
    return int(w[w.index("VW") - 1]) + int(w[w.index("VS") - 1])


def done(d):
    try:
        return open(os.path.join(d, "rsf_solve.log")).read().splitlines()[-1].startswith("wall")
    except (OSError, IndexError):
        return False


os.makedirs(rundir, exist_ok=True)
os.chdir(rundir)
jobs = []
for law in laws:
    for x in XS:
        for ab in ABS:
            n = ncells(law, x, ab)
            if n <= maxcells:
                jobs.append((n * FAC[law] * (x ** 0.5), law, x, ab, n))
for cost, law, x, ab, n in sorted(jobs):
    d = "cattania_%s_ab%g_X%g" % (law, ab, x)
    if done(d):
        continue
    print("running %s (%i cells)" % (d, n), flush=True)
    subprocess.run([os.path.join(sdir, "run_antiplane"), "cattania", str(x), law, ncyc, "4", str(ab), d, series])
print("CATTANIA_SPACE_DONE", flush=True)
