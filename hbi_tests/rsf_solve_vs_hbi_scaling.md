# rsf_solve vs. HBI — BP5-QD earthquake-cycle scaling (initial, 1 km)

**Status: interim.** This records the 1 km (N = 4000) cross-code comparison only. The
0.5 km (N = 16000) runs are in progress and will be added; the conclusions below are
expected to shift at the larger size (see *Interpretation* and *Open items*). All
numbers are specific to the configuration, platform, code versions, and tolerances
described here and should not be read as general statements about either code.

## What is being compared

Both codes integrate the SCEC SEAS **BP5-QD** quasi-dynamic benchmark (Jiang et al.,
2022) on a planar vertical strike-slip fault, rate-and-state friction, run to 1000 yr.

- **rsf_solve** (this repository, `interact`): PETSc `TS` time stepping with the
  Dormand–Prince RK5(4) integrator (`rk 5dp`), and a choice of four interaction-matrix
  backends — `dense`, `HTOOL`, `HACApK`, `hmmvp`.
- **HBI** (T. Ozawa, `github.com/sozawa94/hbi`): an independent boundary-element SEAS
  code using a Cash–Karp RK5(4) integrator and a **lattice H-matrix** built on HACApK.
  HBI is externally maintained; we build it from a pristine clone with a single one-line
  bounds guard (see `build_hbi_bp5.sh`). The lattice H-matrix is explicitly designed for
  very large process counts (Ozawa's demonstrations reach >10^4 ranks), so a small
  N = 4000 test on ≤ 36 cores is well outside its intended operating regime; the results
  here should be read with that in mind.

The two codes use different (but same-order) RK integrators with independent step
controllers, so they take somewhat different step counts for the same nominal ODE
tolerance. This is intrinsic and not tuned away; it is the main reason **per-matvec time
is treated as the cleaner cross-code metric than total wallclock**.

## Matched conditions

| Knob | rsf_solve | HBI | Notes |
|---|---|---|---|
| ODE rel. tolerance | `ts_rtol 1e-4` | `eps_r 1d-4` | matched |
| H-matrix tolerance | HACApK `ztol 1e-1` | lattice `eps_h 1d-1` | matched by *stored bytes*, see below |
| H-matrix storage @ N=4000 | ~12.1 MB | ~12.3 MB | the matching target |
| Integrator | Dormand–Prince 5(4) | Cash–Karp 5(4) | both order-5 adaptive |
| Steps to 1000 yr | ~4590 | ~3983 | differ by ~13% (integrator/controller) |
| Platform | (single node) | same node | see methodological note |

H-matrix accuracy is matched by **stored bytes per DOF**, not by the nominal tolerance
label, because the lattice and standard H-matrix partition the operator differently. At
N = 4000, HACApK `ztol 1e-1` and lattice `eps_h 1d-1` both land near 12 MB (≈ 10×
compression vs. the 122 MB dense operator), which is the closest proxy for matched
achieved accuracy available here — neither cycle run computes a forward error against a
dense reference, so this remains a proxy rather than a verified accuracy match.

## Methodological corrections made along the way

Three issues were found and fixed while setting this up; they materially changed the
comparison and are recorded so the earlier (now-retired) numbers are not reused.

1. **HBI `eps_h` is a per-run `.in` keyword, not a compile-time constant.**
   `read_inputfile()` (`main_LH.f90:171`) runs before matrix generation (`:493`), so an
   `eps_h <value>` line in the input overrides the hardcoded default with no recompile.
   Earlier runs that were *intended* to be `eps_h 1d-1` were in fact running the stock
   `1d-4` (a near-dense operator, ~25.3 MB at N=4000), because the compile-time change
   had not taken effect. `eps_h 1d-1` via the keyword gives ~12.3 MB. The scaling script
   now injects this from its `$3` argument. Storage vs. `eps_h` at N=4000:
   `1d-4 → 25.3 MB`, `1d-1 → 12.3 MB`, `3d-1 → 10.9 MB`, `1d0 → 9.8 MB`.

2. **HBI matvec timing.** The `time for matvec(s)` line prints `sum(st_ctl%time)`, which
   is uninitialized in this build (values ~1e-309 / NaN). The usable matvec wallclock is
   `timeH`, field 4 of the `time(s)` line; the scaling script reads that.

3. **Hardware must match for cross-code absolute timing.** An initial comparison
   inadvertently used two different nodes; the HBI host was ~6% faster at np=1 growing to
   ~25% at np=36 (better memory bandwidth/interconnect), which flattered HBI exactly where
   the matvec is communication-bound. All numbers below are same-node.

A fourth, structural constraint (not a bug): **HBI's lattice H-matrix requires a square
MPI process grid** (`npgl == npgt`), so only perfect-square rank counts run (1, 4, 9, 16,
25, 36, …); non-squares abort with "Process grid must have square shape!". HBI is
therefore sampled on squares and rsf_solve on its usual `1 2 4 8 16 24 48`; the two share
np = 1, 4, 16.

## Results — 1 km (N = 4000), matched ~12 MB, same node

rsf_solve, total wallclock [s] (includes assembly) and per-matvec [ms]:

| np | dense tot | htool tot | hacapk tot | hmmvp tot | hacapk mv | htool mv | hmmvp mv |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1  | 217.9 | 212.7 | 59.77 | 305.3 | 1.340 | 4.826 | 2.651 |
| 4  | 131.7 | 30.79 | 15.17 | 39.51 | 0.329 | 0.526 | 0.464 |
| 16 | 20.96 | 8.096 | 8.090 | 10.90 | 0.193 | 0.161 | 0.106 |
| 48 | 6.352 | 7.476 | 10.00 | 10.78 | 0.260 | 0.166 | 0.061 |

HBI (lattice HACApK, `eps_h 1d-1`):

| np | wall [s] (asm+loop) | matvec [s] | nsteps | per-matvec [ms]* | mem [MB] |
|---:|---:|---:|---:|---:|---:|
| 1  | 59.75 | 54.57 | 3988 | 2.281 | 12.27 |
| 4  | 28.78 | 25.90 | 3997 | 1.080 | 12.33 |
| 9  | 21.66 | 19.93 | 3992 | 0.832 | 12.54 |
| 16 | 14.14 | 12.81 | 3977 | 0.537 | 13.38 |
| 25 | 12.91 | 11.81 | 3981 | 0.494 | 13.23 |
| 36 | 11.59 | 10.38 | 3983 | 0.434 | 13.23 |

\* HBI per-matvec is estimated as `matvec_s / (6 · nsteps)` — six Cash–Karp stage
applies per step, ignoring rejected steps — so it is an **upper bound** (true per-matvec
is somewhat lower). rsf_solve per-matvec is exact (PETSc `-log_view` MatMult count).

Head-to-head against the like-for-like HACApK backend (`rsf hacapk`):

| np | wall HBI / rsf-hacapk | per-matvec HBI / rsf-hacapk |
|---:|---:|---:|
| 1  | 1.00× | 1.70× |
| 4  | 1.90× | 3.28× |
| 16 | 1.75× | 2.78× |

## Interpretation (specific to this configuration)

**Matching the operator removed the serial deficit and roughly halved the rest.** With
both operators near 12 MB, HBI ties `rsf hacapk` in total wall at np=1 (and is far ahead
of rsf's `htool`/`hmmvp`/`dense` there), and sits ~1.75–1.9× behind `rsf hacapk` at
np = 4–16 — versus the ~1.5–2.5× total / ~3–4.6× per-matvec seen when HBI was
inadvertently near-dense.

**The residual per-matvec gap is communication-bound, not accuracy-bound.** Reducing HBI
storage from ~25 to ~12 MB cut its *serial* per-matvec sharply (np=1: 3.77 → 2.28 ms,
×0.61) but barely moved its *parallel* per-matvec (np=16: 0.578 → 0.537 ms, ×0.93). The
serial matvec is flop-bound, so fewer stored entries help directly; the parallel matvec
is dominated by the lattice's communication pattern and its full near-field blocks, which
`eps_h` does not touch. The remaining ~2.8× at np=16 is therefore best read as **intrinsic
lattice communication overhead at this problem size**, the component expected to amortize
as N grows and each rank holds more work per unit of communication.

**Scaling slopes differ as expected.** In this test the standard HACApK matvec
parallelizes faster at low np but reaches its strong-scaling wall at N=4000: `rsf hacapk`
per-matvec bottoms at 0.193 ms (np=16) and then rises (0.216 at 24, 0.260 at 48). HBI's
lattice descends monotonically through np=36 (0.537 → 0.494 → 0.434 ms) with no floor yet
reached. The two have not crossed at this size — HBI's best total (11.59 s at np=36) is
still 1.43× `rsf hacapk`'s best (8.09 s at np=16) — but the slopes are consistent with a
crossover beyond the core counts reachable here.

**Net, for BP5 at 1 km on ≤ 36 cores and matched accuracy:** rsf_solve (`hacapk` or
`htool`, ~8 s at np=16) is the faster choice, while HBI is competitive rather than several-
fold slower, ties at np=1, and carries the better scaling slope. This is consistent with
the lattice H-matrix being built for much larger problems and core counts than this test;
a different problem size, core count, kernel, or HBI/HACApK version could change the
picture, and the 0.5 km runs below are the more relevant test of the lattice's design
regime.

## Open items / caveats

- **0.5 km (N = 16000) pending.** Prediction to test: the communication-bound residual
  shrinks and any crossover moves to lower relative np, so HBI's curve should sit closer
  to — possibly below — `rsf hacapk` at the higher squares. Re-check HBI's stored MB
  against rsf's at 0.5 km, since compression-per-DOF shifts with N (the matched `eps_h`
  may differ from `1d-1`).
- **Matched storage ≠ verified matched error.** A forward-error check against a dense
  matvec for both codes would make the accuracy match rigorous; neither cycle run does
  this today.
- **HBI per-matvec is an upper bound** (÷6 stages, rejections ignored).
- **Total wall mixes assembly + solve + a ~13% step-count difference**; per-matvec is the
  cleaner engine comparison.

## Reproduce

Inputs (place HBI files in `./hbi/examples/`):

- `bp5r_param_gen.py` — generates HBI per-element BP5 params at any resolution
  (`100 40 1.0` → 1 km / N=4000; `200 80 0.5` → 0.5 km / N=16000); reproduces HBI's
  shipped 1 km file bit-for-bit.
- `bp5r_1km.in` + `bp5r_1km_param.dat`, `bp5r_0.5km.in` + `bp5r_0.5km_param.dat`.

Build (once; `eps_h` is set per run, not here):

```
MPIF90=/path/to/mpif90 ./build_hbi_bp5.sh
```

Run scaling:

```
./rsf_solve_scaling_test.sh 1km 1000          # rsf_solve, 4 backends
./hbi_bp5_scaling_test.sh   1km 1000 1d-1      # HBI; $3 injects eps_h into the .in
```

Figure: `rsf_vs_hbi_1km.png` (total wall and per-matvec vs. np; HBI on square np only).
