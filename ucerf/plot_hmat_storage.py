#!/usr/bin/env python3
#
# plot_hmat_storage.py
#
# Analyze the hmat_storage.dat table written by sweep_hmat_storage.sh and
# plot, per backend and relative to the dense baseline, compression,
# matvec speedup, assembly speedup, and the Ax = b GMRES timing against
# the backend tolerance. The tolerances are backend-internal accuracy
# parameters, not verified operator errors, so curves are comparable
# across backends only to the extent the settings were calibrated
# against dense at smaller N.
#
# Solve-column honesty: the sweep runs GMRES(30) without a
# preconditioner, capped at 10000 iterations. If a run (or the dense
# baseline) hits the cap, its time measures iteration throughput, not
# time to solution; the script detects this and relabels the solve
# panel accordingly, since comparing non-converged "solve times" as if
# they were solutions would overstate what was measured.
#
# Parameters are plain assignments below.

infile = "hmat_storage.dat"
outfile = "hmat_storage.pdf"
maxit = 10000            # the tool's GMRES iteration cap

import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def fget(tokens, key):
    try:
        v = tokens[tokens.index(key) + 1]
        return float(v) if v != "NA" else None
    except (ValueError, IndexError):
        return None

dense_ms = dense_as = dense_ss = dense_si = None
dense_tag = ""
npatch = None
data = {}   # backend -> list of (eps, mb, matvec_ms, assembly_s, solve_s, solve_its)

for line in open(infile):
    t = line.split()
    if line.startswith("#"):
        if "dense baseline:" in line:
            dense_ms = fget(t, "per_matvec")
            dense_as = fget(t, "assembly")
            dense_ss = fget(t, "per_solve")
            dense_si = fget(t, "solve_its")
            if line.rstrip().endswith("]"):
                dense_tag = line[line.rindex("[") + 1:line.rindex("]")]
        for f in t:
            if f.startswith("npatch="):
                npatch = int(f.split("=")[1])
        continue
    if len(t) < 12:
        continue
    try:
        eps = float(t[1])
        mb = float(t[3]) if t[3] != "NA" else None
        ms = float(t[5]) if t[5] != "NA" else None
        asm = float(t[7]) if t[7] != "NA" else None
        ss = float(t[9]) if t[9] != "NA" else None
        si = float(t[10]) if t[10] != "NA" else None
    except ValueError:
        continue
    data.setdefault(t[0], []).append((eps, mb, ms, asm, ss, si))

if not data:
    sys.exit(f"no data rows found in {infile}")
if npatch is None:
    sys.exit(f"no npatch entry found in the header of {infile}")

dense_mb = npatch * float(npatch) * 8.0 / 1048576.0
dense_capped = (dense_si is not None and dense_si >= maxit)
any_capped = dense_capped or any(r[5] is not None and r[5] >= maxit
                                 for rows in data.values() for r in rows)

def human_mem(mb):
    if mb >= 1048576.0: return f"{mb/1048576.0:.1f} TB"
    if mb >= 1024.0:    return f"{mb/1024.0:.1f} GB"
    return f"{mb:.0f} MB"

def human_time_ms(ms):
    if ms >= 1000.0: return f"{ms/1000.0:.2f} s"
    return f"{ms:.1f} ms"

def human_time_s(s):
    if s >= 3600.0: return f"{s/3600.0:.1f} h"
    if s >= 60.0:   return f"{s/60.0:.1f} min"
    return f"{s:.2f} s"

if any_capped:
    solve_title = (f"GMRES(30) throughput at the {maxit}-iteration cap\n"
                   f"(NO run converged; dense: {human_time_s(dense_ss)}/{maxit} its)"
                   if dense_ss else "GMRES (iteration-capped, not converged)")
    solve_ylab = "iteration-capped GMRES time ratio vs dense (x)"
else:
    solve_title = ("Ax=b GMRES solve speedup"
                   + (f" (dense: {human_time_s(dense_ss)}, {dense_si:.0f} its)"
                      if dense_ss else ""))
    solve_ylab = "Ax=b solve speedup vs dense (x)"

fig, axs = plt.subplots(2, 2, figsize=(10.5, 8.4))
(axc, axm), (axa, axsv) = axs

panels = [
    (axc, 1, dense_mb, "storage reduction vs dense (x)",
     f"compression (dense: {human_mem(dense_mb)})"),
    (axm, 2, dense_ms, "matvec speedup vs dense (x)",
     "matvec speedup" + (f" (dense: {human_time_ms(dense_ms)}/apply)"
                         if dense_ms else "")),
    (axa, 3, dense_as, "assembly speedup vs dense (x)",
     "assembly speedup" + (f" (dense: {human_time_s(dense_as)})"
                           if dense_as else "")),
    (axsv, 4, dense_ss, solve_ylab, solve_title),
]

for ax, col, ref, ylab, title in panels:
    for b in sorted(data):
        pts = [(r[0], ref / r[col]) for r in sorted(data[b])
               if ref is not None and r[col] is not None and r[col] > 0]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", label=b)
    ax.axhline(1.0, color="k", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("backend tolerance")
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=9.5)
    ax.grid(alpha=0.3, which="both")
    if ax.get_lines():
        ax.legend(fontsize=8)

fig.suptitle(f"H-matrix backends relative to dense, N = {npatch} [{dense_tag}]",
             fontsize=11)
fig.subplots_adjust(top=0.90, hspace=0.38, wspace=0.26,
                    left=0.08, right=0.98, bottom=0.07)
fig.savefig(outfile)
print(f"wrote {outfile}")

if any_capped:
    print(f"\nNOTE: GMRES hit the {maxit}-iteration cap in the dense baseline"
          " and/or sweep runs; the solve panel shows iteration-capped"
          " throughput, not time to solution. Unpreconditioned GMRES does"
          " not converge on this operator at this size; a preconditioner"
          " would be needed for a true solve comparison.")

print(f"\n{'backend':12s} {'best eps':>8s} {'compress_x':>10s} {'matvec_x':>8s} "
      f"{'assembly_x':>10s}")
for b in sorted(data):
    rows = [r for r in sorted(data[b]) if r[2] is not None]
    if not rows:
        continue
    r = min(rows, key=lambda r: r[2])
    def x(ref, v): return f"{ref/v:.1f}" if (ref and v) else "NA"
    print(f"{b:12s} {r[0]:8.0e} {x(dense_mb,r[1]):>10s} {x(dense_ms,r[2]):>8s} "
          f"{x(dense_as,r[3]):>10s}")
