#!/usr/bin/env python3
#
# plot_hmat_storage.py
#
# Analyze the hmat_storage.dat table written by sweep_hmat_storage.sh and
# plot, per backend and relative to the dense baseline, compression,
# matvec speedup, assembly speedup, and Ax = b solve speedup against the
# backend tolerance. The tolerances are backend-internal accuracy
# parameters, not verified operator errors, so curves are comparable
# across backends only to the extent the settings were calibrated
# against dense at smaller N. Rows from the pre-solve table format (no
# assembly or solve columns) are handled; those panels then stay empty.
#
# Parameters are plain assignments below.

infile = "hmat_storage.dat"
outfile = "hmat_storage.pdf"

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
data = {}     # backend -> list of (eps, mb, ms, assembly_s, solve_s)

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
    if len(t) < 7:
        continue
    b = t[0]
    try:
        eps = float(t[1])
        mb = float(t[3]) if t[3] != "NA" else None
        ms = float(t[5]) if t[5] != "NA" else None
        asm = float(t[7]) if len(t) >= 13 and t[7] != "NA" else None
        ss = float(t[9]) if len(t) >= 13 and t[9] != "NA" else None
    except ValueError:
        continue
    data.setdefault(b, []).append((eps, mb, ms, asm, ss))

if not data:
    sys.exit(f"no data rows found in {infile}")
if npatch is None:
    sys.exit(f"no npatch entry found in the header of {infile}")

dense_mb = npatch * float(npatch) * 8.0 / 1048576.0

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

fig, axs = plt.subplots(2, 2, figsize=(10.5, 8.2))
(axc, axm), (axa, axsv) = axs

panels = [
    (axc,  1, dense_mb, "storage reduction vs dense (x)",
     f"compression (dense: {human_mem(dense_mb)})"),
    (axm,  2, dense_ms, "matvec speedup vs dense (x)",
     "matvec speedup" + (f" (dense: {human_time_ms(dense_ms)}/apply)"
                         if dense_ms else " (no dense baseline)")),
    (axa,  3, dense_as, "assembly speedup vs dense (x)",
     "assembly speedup" + (f" (dense: {human_time_s(dense_as)})"
                           if dense_as else " (no dense assembly time)"),),
    (axsv, 4, dense_ss, "Ax=b solve speedup vs dense (x)",
     "Ax=b GMRES solve speedup"
     + (f" (dense: {human_time_s(dense_ss)}, {dense_si:.0f} its)"
        if dense_ss else " (no dense solve time)")),
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
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3, which="both")
    if ax.get_lines():
        ax.legend(fontsize=8)

fig.suptitle(f"H-matrix backends relative to dense, N = {npatch} [{dense_tag}]",
             fontsize=11)
fig.subplots_adjust(top=0.92, hspace=0.32, wspace=0.26,
                    left=0.08, right=0.98, bottom=0.07)
fig.savefig(outfile)
print(f"wrote {outfile}")

# console summary: best point per backend by solve time (fall back to matvec)
print(f"{'backend':12s} {'best eps':>8s} {'matvec_x':>8s} {'assembly_x':>10s} "
      f"{'solve_x':>8s}")
for b in sorted(data):
    rows = [r for r in sorted(data[b]) if r[2] is not None]
    if not rows:
        continue
    key = (lambda r: r[4]) if all(r[4] is not None for r in rows) else (lambda r: r[2])
    r = min(rows, key=key)
    def x(ref, v): return f"{ref/v:.1f}" if (ref and v) else "NA"
    print(f"{b:12s} {r[0]:8.0e} {x(dense_ms,r[2]):>8s} {x(dense_as,r[3]):>10s} "
          f"{x(dense_ss,r[4]):>8s}")
