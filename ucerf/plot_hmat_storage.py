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
else:
    solve_title = ("inverse solve x = A\\b (GMRES)"
                   + (f" (dense: {human_time_s(dense_ss)}, {dense_si:.0f} its)"
                      if dense_ss else ""))

fig, axs = plt.subplots(2, 2, figsize=(10.5, 8.4))
(axc, axm), (axa, axsv) = axs

# right-column panels state the operation they time: the top one is the
# forward product b = A x (a single operator application), the bottom
# one the inverse solve x = A\b (iterative GMRES solution); their axis
# labels sit on the right edge to make the pairing visually explicit
panels = [
    (axc, 1, dense_mb, "storage reduction vs dense (x)",
     f"compression (dense: {human_mem(dense_mb)})", False),
    (axm, 2, dense_ms, "forward product  b = A x   speedup vs dense (x)",
     "matvec: forward product b = A x"
     + (f" (dense: {human_time_ms(dense_ms)}/apply)" if dense_ms else ""), True),
    (axa, 3, dense_as, "assembly speedup vs dense (x)",
     "assembly speedup" + (f" (dense: {human_time_s(dense_as)})"
                           if dense_as else ""), False),
    (axsv, 4, dense_ss, "inverse solve  x = A\\b   speedup vs dense (x)",
     solve_title, True),
]

for ax, col, ref, ylab, title, right in panels:
    vlo = 1e30
    vhi = 0.0
    for b in sorted(data):
        pts = [(r[0], ref / r[col]) for r in sorted(data[b])
               if ref is not None and r[col] is not None and r[col] > 0]
        if pts:
            ln, = ax.plot([p[0] for p in pts], [p[1] for p in pts],
                          "o-", label=b)
            for p in pts:
                if p[1] < vlo: vlo = p[1]
                if p[1] > vhi: vhi = p[1]
            # annotate each backend's best point with its factor
            pb = max(pts, key=lambda p: p[1])
            ax.annotate(f"{pb[1]:.0f}x", xy=pb, xytext=(0, 5),
                        textcoords="offset points", ha="center",
                        fontsize=7.5, color=ln.get_color())
    ax.set_xscale("log")
    ax.set_yscale("log")
    if vhi > 0.0:
        ax.set_ylim(vlo/1.35, vhi*1.45)
    # limited log range: place explicit 1-2-3-5-7 subdecade ticks within
    # the data range so the axis stays readable at roughly one order of
    # magnitude span
    if vhi > 0.0:
        tickv = []
        dec = -2
        while 10.0**dec < vhi*1.45:
            for sub in (1.0, 2.0, 3.0, 5.0, 7.0):
                v = sub*10.0**dec
                if (v >= vlo/1.35) and (v <= vhi*1.45):
                    tickv.append(v)
            dec += 1
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(tickv))
        ax.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(
            [f"{v:g}" for v in tickv]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("backend tolerance")
    ax.set_ylabel(ylab)
    if right:
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.yaxis.set_ticks_position("both")
    ax.set_title(title, fontsize=9.5)
    ax.grid(alpha=0.3, which="both")
    if ax.get_lines():
        ax.legend(fontsize=8, loc="best")

fig.suptitle(f"H-matrix backends relative to dense, N = {npatch} [{dense_tag}]",
             fontsize=11)
fig.subplots_adjust(top=0.90, hspace=0.38, wspace=0.30,
                    left=0.08, right=0.91, bottom=0.07)
fig.savefig(outfile, bbox_inches="tight")
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
