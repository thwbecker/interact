# rsf_solve vs. HBI — BP5-QD earthquake-cycle scaling (1 km & 0.5 km)

All numbers are specific to this configuration, platform, code versions, and tolerances,
and are not general statements about either code. Same single node throughout.

## Setup

Both codes integrate the SCEC SEAS **BP5-QD** benchmark (Jiang et al., 2022) to 1000 yr.

- **rsf_solve** (`interact`): PETSc `TS`, Dormand–Prince RK5(4); backends `dense`,
  `HTOOL`, `HACApK`, `hmmvp`.
- **HBI** (Ozawa, `github.com/sozawa94/hbi`): Cash–Karp RK5(4), **lattice H-matrix** on
  HACApK; built from a pristine clone with a one-line bounds guard. The lattice
  H-matrix targets very large rank counts (Ozawa's demonstrations reach >10⁴), so a
  test on ≤48 cores is well below its design regime.

The two integrators take different step counts (rsf ~4.6k/10k vs HBI ~4.0k/8.8k at
1 km/0.5 km), so **per-matvec time is the cleaner cross-code metric than total wall**.

## Matched conditions

ODE tolerance matched (`ts_rtol`/`eps_r` = 1e-4). H-matrix accuracy matched by **stored
bytes** (a proxy — neither cycle run computes a forward error against a dense reference):

| size | rsf-hacapk `ztol 1e-1` | HBI `eps_h 1d-1` |
|---|---|---|
| 1 km (N=4000) | ~12.1 MB | ~12.3 MB |
| 0.5 km (N=16000) | ~73 MB | ~68 MB (~7% apart) |

HBI's lattice needs a square MPI grid, so it runs only on perfect-square np (1, 4, 9, 16,
25, 36); rsf runs on 1, 2, 4, 8, 16, 24, 48. Shared points: np = 1, 4, 16.

## Notes / corrections (reproducibility)

- **HBI `eps_h` is a per-run `.in` keyword** (read before assembly, `main_LH.f90:171`),
  not compile-time. Earlier "1d-1" runs were actually near-dense `1d-4` (~25 MB at 1 km);
  the scaling script now injects `eps_h` from its `$3` argument.
- **HBI matvec time = `timeH`** (field 4 of the `time(s)` line); `sum(st_ctl%time)` is
  uninitialised in this build (~1e-309/NaN).
- **Same node required** for cross-code absolute timing — an early mixed-host run
  flattered HBI by 6–25%, growing with np.

## Results (matched storage, same node)

Best total wall, and HBI-vs-rsf-hacapk head-to-head (the like-for-like HACApK pairing):

| size | rsf best (backend@np) | HBI best (@np) | wall ratio HBI/rsf-hac @1 / @4 / @16 | per-matvec ratio @1 / @4 / @16 |
|---|---|---|---|---|
| 1 km   | ~7.5 s (htool@48), 8.1 s (hacapk@16) | 11.6 s (@36) | 1.00 / 1.90 / 1.75 | 1.70 / 3.28 / 2.78 |
| 0.5 km | 57 s (hacapk@16), 59 s (htool@48)    | 104 s (@36)  | 0.90 / 1.45 / 2.48 | 1.37 / 2.20 / 3.68 |

(HBI per-matvec is `matvec_s / (6·nsteps)` — six Cash–Karp stage applies/step, rejected
steps ignored — so an **upper bound**; rsf per-matvec is exact from `-log_view`.)
Figures: `rsf_vs_hbi_1km.png`, `rsf_vs_hbi_0p5km.png`.

## Interpretation (this configuration, ≤48 cores)

At both sizes, with matched accuracy, **rsf_solve (`hacapk`/`htool`/`hmmvp`) is faster
than HBI's lattice on ≤48 cores.** HBI is competitive — even slightly faster — at np=1,
where the larger N amortizes its serial overhead (0.5 km: 628 vs 696 s), but standard
HACApK's stronger parallel scaling wins by np=4–16.

**The 1 km data alone suggested an imminent crossover; the 0.5 km data does not support
that.** At 1 km, HBI's per-matvec kept descending while rsf-hacapk's bottomed near np=16,
hinting the lattice would overtake at larger N. Instead, the larger problem helped
*standard* HACApK far more: its np=16 parallel efficiency rose from ~0.46 (1 km) to ~0.76
(0.5 km, near-linear), while HBI's barely moved (~0.26 → ~0.28). So the gap *widens* with
core count at 0.5 km, and any crossover is pushed to **higher** np, not lower. This is
consistent with the lattice's advantage being asymptotic (its regime is ~10³–10⁴ ranks,
not reached here): a larger N gives standard HACApK more runway before its own
communication wall, rather than handing the lattice an early win.

**Hardware caveat.** rsf-hacapk's total wall bottoms near np=16 and then *rises* (np=24,
48), with `htool`/`hmmvp` — which keep scaling — overtaking it at high core counts. This
saturation-and-degradation, and the exact np at which backends cross, may substantially
reflect **this specific node's memory bandwidth and interconnect** rather than an
intrinsic property of standard HACApK's algorithm. On a machine with more bandwidth per
core or a faster network, the saturation point would likely move and the post-np=16
degradation could shrink or vanish. We therefore do **not** read the np≈16 roll-off as a
fixed property of the method; only the broad ordering (H-backends >> dense; lattice not
competitive at ≤48 cores at matched accuracy) is robust here.

## Caveats

- Matched **stored bytes** is a proxy for matched accuracy; no forward error is computed
  in the cycle runs.
- HBI per-matvec is an **upper bound** (÷6 stages, rejections ignored).
- Total wall mixes assembly + solve and a ~10–15% step-count difference; per-matvec is
  the cleaner engine metric.
- ≤48 cores says nothing about HBI at the ~10³–10⁴ ranks the lattice targets; this
  comparison is about the small/medium-node regime only.

## Reproduce

```
# HBI inputs (place in ./hbi/examples/): bp5r_param_gen.py 100 40 1.0 -> 1 km (N=4000),
#                                        bp5r_param_gen.py 200 80 0.5 -> 0.5 km (N=16000)
MPIF90=/path/to/mpif90 ./build_hbi_bp5.sh          # build once; eps_h is set per run
./rsf_solve_scaling_test.sh 0.5km 1000             # rsf (drop dense for N=16000: very slow)
./hbi_bp5_scaling_test.sh   0.5km 1000 1d-1        # HBI; $3 injects eps_h into the .in
```
