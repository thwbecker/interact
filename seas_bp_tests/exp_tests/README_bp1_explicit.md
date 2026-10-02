# Explicit pseudo-transient SEAS BP1 prototype: status, results, open items

Updated: 2026-10-02, revision 2 (supersedes the 2026-09-29 version, kept as README_bp1_explicit_2026-09-29.md)

## 1. What this is

A test of whether an explicit, accelerated pseudo-transient (PT) solver in the JustRelax/ParallelStencil style can serve as the basis for quasi-dynamic seismic cycle simulation with rate-and-state friction, and whether a quasi-viscous (Herrendörfer et al. 2018 style) fault treatment reproduces the sharp-interface result. Test problem: SEAS BP1-type setup, 2D antiplane strain, vertical strike-slip fault at x = 0, free surface, regularized rate-and-state friction with the aging law, radiation damping, plate loading at Vp. The half-space is truncated to a box with Dirichlet loading on the far boundaries; only x >= 0 is computed (symmetry).

### Method

- Quasi-static elasticity solved at every physical time step by a damped-wave PT iteration (nu = 4), warm-started by linear extrapolation from the two previous steps. Convergence when max residual times dx < tol_pa (1 Pa in all production runs).
- Two fault treatments, selectable with `--fault`:
  - `robin`: displacement discontinuity at x = 0, traction from a second-order one-sided stencil, slip rate from a per-node safeguarded Newton solve inside every PT iteration.
  - `layer`: fault layer of width D = dx straddling x = 0 (nodes at (i + 1/2) dx) carrying plastic slip; layer stress mu (2 u0 - delta_p) / dx; the quasi-static form of the viscoplastic viscosity eta_f = tau dx / V. Node 0 is a PT unknown; slip rate from the same Newton solve.
- State variable held fixed within a step (explicit in state, implicit in slip rate); updated afterwards with the exact aging-law integral for constant V. First order in time.
- Adaptive time step dt = 0.1 Dc / Vmax.
- Dense logging (every step) while Vmax > 1e-4 m/s; every `out_every` steps otherwise.

### Files (`bp1/`)

| File | Content |
|---|---|
| `bp1_pt.py` | NumPy reference plus fused numba kernels (`--numba`, `--threads`) for both fault treatments; checkpoint/restart; dense logging; `--stop_Vmax`; `--wall_limit` |
| `bp1_pt.jl` | ParallelStencil translation of the robin scheme, CPU threads or CUDA via `USE_GPU`; CPU version verified against Python (identical iteration counts); no layer mode, no checkpointing yet |
| `run_compare` | Runs robin and layer at 250 and 125 m in the 30 x 50 km box and calls the comparison |
| `compare_runs.py` | Event catalogues, intervals, iteration statistics, figure; `--align` shifts each run so its first full event is at t = 0 |
| `run_*_{250,125,60}.csv` | Logs of the production runs described below |

## 2. What has been tested

### Solver

- Manufactured-solution test of the elastic PT solve: second-order spatial convergence; iterations to 1e-10 scale linearly with n (about 12.5 n at nu = 4).
- Fault Newton solve verified against brute-force bisection to 1e-12 relative; 4 iterations warm-started.
- Numba kernels reproduce the NumPy paths bit for bit (both fault modes); 1 vs 4 threads give identical iteration counts.
- Julia CPU version reproduces Python (robin) with identical iteration counts and Vmax within 5e-11 over 237 steps.
- Tolerance: 1 Pa vs 0.1 Pa from a common post-event state at 250 m changed event onsets by under 0.5 days over 137 years. 1 Pa is adequate.

### Two findings about the formulation

- The state variable must be explicit within the fault equation. An implicit state makes g(V) non-monotone under velocity weakening and root finding unsafe.
- The first event in BP1 grows out of numerical noise on an unstable initial steady state. Its timing (8 to 28 yr across runs) carries no information and must be removed before comparing runs (`--align`). Everything after it is reproducible.

### Resolution, 30 x 50 km box, 1 Pa

h* (Rubin-Ampuero) is about 2 km; the quasi-static cohesive zone about 300 m.

| dx | cells per h* | robin | layer |
|---|---|---|---|
| 1 km | 2 | checkerboard after first event; unusable | not run |
| 500 m | 4 | single-node fast-slip flicker between events; unusable | same |
| 250 m | 8 | partial deep events (10 to 16 km) alternating with full ruptures; 137 yr full-to-full | same pattern, weaker partials; 144 yr |
| 125 m | 16 | period-4 cycle: intervals 63.4 / 88.0 / 67.3 / 77.7 yr, slips 2.26 / 2.69 / 0.85 / 3.55 m, the third event arrested (peak 0.76 m/s); strictly periodic to 988 yr | period-2 cycle: 67.7 / 88.4 yr, 2.23 / 2.69 m, peaks 4.05 / 4.40 m/s; strictly periodic to 990 yr |
| 60 m | 33 | period-2: 67.8 / 88.7 yr, 2.25 / 2.71 m, 4.07 / 4.44 m/s (5 events) | period-2: 67.8 / 89.1 yr, 2.28 / 2.71 m, 4.10 / 4.46 m/s (5 events) |

Conclusions:

- The converged solution in this box is the period-2 cycle. Layer at 125 m, robin at 60 m and layer at 60 m agree to 0.3 yr in intervals, 0.03 m in slip, 1 percent in peak slip rate, and 0.01 MPa in post-event stress (27.41 MPa).
- The robin scheme at 125 m is the unconverged one: its period-4 cycle with an arrested event disappears at 60 m. The layer scheme is converged at 125 m where the sharp-interface scheme is not. This is the reverse of the expectation from Herrendörfer et al.'s first-order convergence for a D = dx layer.
- Post-event stress and event durations (125 to 155 s for full ruptures, quasi-dynamic) are consistent across all converged runs.
- The slip budget closes in every periodic cycle (sum of slips per period equals Vp times the period).

### Cost

| n along fault (dx) | interseismic PT iterations per step | coseismic mean / max | steps per 320 yr |
|---|---|---|---|
| 80 (1 km) | 50 to 100 | 82 / 220 | |
| 200 (250 m) | 300 (robin), 320 (layer) | 290 / 500, 330 / 545 | 129k, 117k |
| 400 (125 m) | 539, 618 | 569 / 1132, 662 / 1197 | 138k, 137k |
| 800 (60 m) | 1099, 1328 | 1248 / 2168, 1352 / 2417 | 133k, 128k |

Iterations per step scale linearly with n through n = 800. Step count is resolution-independent (dt set by Dc / V), so cost per resolution doubling is about 7x. Layer costs 15 to 20 percent more iterations than robin. Thread scaling of the numba kernel at 120 x 200 was 1.35x on 4 threads; larger grids scale better.

## 3. Open items

### Physics / benchmark

1. **Box size.** The period-2 cycle is resolution-converged but the BP1 reference is period-1 at about 78 yr. The remaining suspect is the 30 x 50 km box with Dirichlet loading. Run layer mode at 125 m (adequate resolution) in a 60 x 100 km box:

   ```
   python3 bp1_pt.py --numba --threads 4 --fault layer --nx 480 --nz 800 --Lx 60e3 --Lz 100e3 \
       --t_end_yr 500 --tol_pa 1.0 --out_every 500 \
       --log run_layer_125_bigbox.csv --checkpoint ck_layer_125_bigbox.npz
   ```

   Expected cost: per step about that of the 60 m runs (n along the fault doubles to 800), roughly 2 h per 100 yr on 4 threads. If period-1 near 78 yr appears, the prototype reproduces BP1 up to loading geometry and the comparison against the reference values (recurrence, peak slip rate, stress drop, slip per event) can be made quantitatively. If period-2 persists, the method is the suspect and item 2 is next.
2. **Time integration order.** State is first order in time. Test a predictor-corrector (Heun) step on slip and state at 125 m, layer mode, in the smaller box, where the period-2 result is known to 0.3 yr. Doubles the PT solves per step.
3. **Robin at 60 m beyond five events.** Run in progress (2026-10-02); the five completed events are on the period-2 cycle.
4. **Direct comparison to the BP1-QD reference** (Erickson et al. 2020 time series at the SEAS stations) once item 1 is settled.

### Numerics / cost

5. **Iteration scaling.** Linear in n is the main cost driver and is now measured to n = 800. Options to break it: better initial guess (second-order extrapolation in time), multigrid-type acceleration of the PT iteration, or nonuniform grid refinement toward the fault. Any of these should be tested on the 125 m layer case where the answer is known.
6. **Switch to physical inertia during rupture.** The one structural advantage of the PT framework (set pseudo-inertia and damping to physical values and the same loop becomes an explicit elastodynamic solver) has not been exercised. The 60 m runs give the quasi-dynamic reference to compare against.
7. **Healing and Maxwell time step limits** (Herrendörfer et al. eqs. 31 and 33) are not implemented; needed once a viscous off-fault medium is added, not for the elastic benchmark.

### Software

8. **Julia version**: add the layer mode and checkpoint/restart, then run the GPU equivalence test and the time-per-iteration scan versus grid size. The 60 m runs on CPU took hours; 30 m is the next resolution step and belongs on the GPU.
9. **JustRelax integration**: the layer treatment is a pointwise nonlinear viscosity update, the same structure as JustRelax's plasticity Picard step, so the natural next step is a JustRelax miniapp with the fault layer as a viscoplastic material in the velocity-stress Maxwell formulation, which also brings the viscous off-fault medium for free.

## 4. Assessment

The elliptic PT solve was never the limiting factor: it converged reliably with warm starts, and 1 Pa tolerance is sufficient. The resolution requirement (16 cells per h* for the layer scheme, 33 for the sharp interface) is set by the physics, not the method, and matches experience with other volume methods. The quasi-viscous fault treatment is validated against the sharp-interface one at the level of 0.3 yr in recurrence over several cycles and converges earlier, which makes it the preferred basis for the visco-elasto-plastic extension. What the prototype has not yet shown is agreement with the BP1 reference sequence; that question is now reduced to one run (the larger box). The cost of the method is well characterized: linear iteration growth with resolution on top of the physics-set step count.
