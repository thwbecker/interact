# Explicit pseudo-transient SEAS BP1 prototype: status, setup, next tests

Date: 2026-09-29

## 1. What this is

A test of whether an explicit, accelerated pseudo-transient (PT) solver in the JustRelax/ParallelStencil style is a workable basis for quasi-dynamic seismic cycle simulation with rate-and-state friction. The test problem is a SEAS BP1-type setup: 2D antiplane strain, vertical strike-slip fault at x = 0, free surface, regularized rate-and-state friction with the aging law, radiation damping, plate loading at Vp. The half-space is truncated to a box with Dirichlet loading on the far boundaries and the fault is modeled with symmetry, so only x >= 0 is computed.

Method summary:

- The quasi-static elastic problem is solved at each physical time step by a damped-wave PT iteration (damping parameter nu = 4, CFL 0.9), warm-started by linear extrapolation from the previous two steps.
- The fault is a nonlinear Robin condition. At every PT iteration the slip rate V(z) is found per node from the current traction by a safeguarded Newton iteration in log V, with the state variable held at its value from the previous time step. The fault displacement is then set to u0_old + V dt / 2.
- State is updated after the step with the exact integral of the aging law for constant V.
- Adaptive time step dt = xi Dc / Vmax with xi = 0.1.
- Convergence: maximum residual times dx below tol_pa (default 1 Pa).

Files (`bp1/`):

| File | Content |
|---|---|
| `bp1_pt.py` | NumPy reference implementation with optional fused numba kernel (`--numba`), checkpoint/restart (`--checkpoint`), incremental CSV log, `--threads`, `--stop_Vmax`, `--wall_limit` |
| `bp1_pt.jl` | ParallelStencil translation, CPU threads or CUDA backend via `USE_GPU`; verified equivalent to the Python version on CPU (see below) |
| `run80*.csv`, `run160*.csv`, `run250*.csv` | Logs from the 1 km, 500 m, and 250 m runs described below |
| `ck250*.npz` | Restartable 250 m states (1 Pa branch at 296 yr, 0.1 Pa branch, and the post-event-3 state at 159 yr) |

## 2. Results so far

All runs use BP1 parameters (mu = 32.04 GPa, sigma_n = 50 MPa, a = 0.010 to 0.025, b = 0.015, Dc = 8 mm, Vp = 1e-9 m/s, VW to 15 km, transition to 18 km, RSF to 40 km).

### Solver verification

- Manufactured-solution test of the elastic PT solve: second-order spatial convergence (L2 error 3.3e-4, 8.3e-5, 2.1e-5 for n = 32, 64, 128). Iterations to 1e-10 residual scale linearly with n, about 12.5 n at nu = 4.
- Fault Newton solve verified against brute-force bisection to 1e-12 relative. Warm-started it converges in about 4 iterations.
- Numba kernel reproduces the NumPy path bit for bit; 1 and 4 threads give identical iteration counts and Vmax.
- Julia CPU version (`Threads` backend, `check_every = 1`) reproduces the Python version over 237 steps at 120 x 200: identical iteration counts, Vmax within 5e-11 relative.

### Two design findings

- Holding the state implicit in the fault equation makes g(V) non-monotone under velocity weakening and root finding unsafe. State must be explicit within the step.
- The first event in BP1 starts from an unstable steady state, so its timing depends on numerical noise. First-event times at 1 km varied non-monotonically between 8.8 and 28 yr as the tolerance changed from 10 Pa to 0.01 Pa. First-event timing is therefore not a usable accuracy metric.

### Resolution

Box 40 x 60 km unless noted; h* (Rubin-Ampuero) is about 2 km.

| dx | cells per h* | Behaviour |
|---|---|---|
| 1 km (80 x 80, 80 km box) | 2 | First event 8.8 yr, then node-to-node checkerboard in slip and state; catalogue meaningless after the first event |
| 500 m (80 x 120) | 4 | First event smooth (23.6 yr, 4.2 m at 7.5 km). Afterwards single nodes in the shallow VW zone slide at 1e-4 to 1e-2 m/s for extended periods while neighbours are locked; 40 percent of interseismic samples show fast slip. Reducing xi from 0.1 to 0.02 did not change this, so it is spatial, not temporal |
| 250 m (120 x 200, 30 x 50 km box) | 8 | Quiet interseismic periods; smooth post-event profiles. Catalogue to 296 yr: full ruptures at 22.6, 159.1 and 296.2 yr (about 4 m at 7.5 km, stress 32 to 34 MPa dropping to 27.4 MPa), partial deep events (10 to 16 km only) at 82.3 and 218.7 yr |

The alternating partial/full pattern at 250 m with a 137 yr full-to-full period does not match the BP1 reference (full ruptures only, 78 yr). Two candidate causes, not yet separated: the small box with Dirichlet loading, and residual underresolution. A 125 m run is in progress (Python, numba, threaded).

### Tolerance

At 250 m, restarting from the identical state after the third event with 1 Pa and 0.1 Pa tolerance gave event onsets within 0.5 days over the following 137 years (218.6936 vs 218.6935 yr; 296.1376 vs 296.1388 yr). 1 Pa is adequate at this resolution. Iterations per step were about 325 (1 Pa) vs 545 (0.1 Pa) interseismic.

### Cost

| n along fault | PT iterations per step, interseismic (warm start) | coseismic mean / max |
|---|---|---|
| 80 | about 50 to 100 | 82 / 220 |
| 200 | about 350 | 285 / 531 |

Iterations scale roughly linearly with n. The 250 m run needed 77,000 steps and 2.5 million PT iterations for 159 years, about 35 minutes single-core with numba. Extrapolated to 25 m in this box: 3000 to 4000 iterations per step and about 1e5 coseismic steps per event.

Thread scaling of the numba kernel at 120 x 200: 1.35x on 4 threads (barrier overhead dominates at this size).

## 3. Julia installation

### Without GPU

```
curl -fsSL https://install.julialang.org | sh
```

Restart the shell, then:

```
julia -e 'using Pkg; Pkg.add("ParallelStencil")'
```

Run with N threads:

```
julia -t N -e 'include("bp1_pt.jl"); run_bp1(nx=120, nz=200, Lx=30e3, Lz=50e3, t_end_yr=6.0, out_every=1, check_every=1, logfile="jl6.csv")'
```

`bp1_pt.jl` must have `const USE_GPU = false`.

### With NVIDIA GPU

Requires a CUDA-capable driver on the machine; CUDA.jl downloads its own toolkit.

```
julia -e 'using Pkg; Pkg.add(["ParallelStencil", "CUDA"])'
julia -e 'using CUDA; CUDA.versioninfo()'
```

Set `const USE_GPU = true` in `bp1_pt.jl`. The number of Julia threads is irrelevant for the GPU backend.

For AMD GPUs, replace `CUDA` with `AMDGPU` in both the `Pkg.add` line and the `@static if USE_GPU` block of `bp1_pt.jl` (`using AMDGPU; @init_parallel_stencil(AMDGPU, Float64, 2)`). Not tried.

### Notes on the Julia file

- `Data.Array` takes the number type from `@init_parallel_stencil`; integer index arrays are stored as Float64 and converted in the kernels.
- `check_every` sets how often the residual reduction is evaluated. 1 reproduces the Python iteration counts; 10 saves reductions at a cost of up to 9 extra iterations per step.
- The Julia log has six columns (`step,t_yr,dt_s,it,resid_pa,Vmax`); the Python log has eleven. Comparisons use `it` and `Vmax`.
- No checkpoint/restart yet.

## 4. What to test next on a GPU

In order. Each item has a pass criterion.

1. **Correctness on the GPU backend.** Same 6 yr, 120 x 200 test as on CPU, compare against the CPU `jl6.csv`. Pass: identical `it` column, Vmax within about 1e-11 relative. Watch for `Performing scalar indexing` warnings; the boundary assignments in the driver are the first suspects if any appear.

2. **Time per PT iteration versus grid size.** Run about 100 steps at 120 x 200, 240 x 400, 480 x 800, 960 x 1600 with `check_every = 10` and record wall time divided by total iterations. Expected: nearly flat time per iteration until the GPU saturates, then linear in node count. The knee tells you which resolutions are effectively free. At small sizes the run is launch-bound (three kernels plus a reduction per iteration).

3. **Reduction cost.** Compare `check_every = 1` and `check_every = 10` at one large size. If the reduction is a significant fraction of the iteration time, it is worth fusing the residual and update kernels and keeping the reduction sparse.

4. **Checkpoint/restart.** Add saving of `u`, `u_prev`, `V`, `theta`, `t`, `step`, `dt_prev` (as host arrays) before any run longer than an hour. Pass: a restarted run continues the log with identical values.

5. **Resolution convergence in the 30 x 50 km box.** Run 125 m and 60 m to at least 320 yr with 1 Pa tolerance. Compare event catalogues (onset times, peak slip rate, slip at 7.5 km, stress drop) against the 250 m Python run. Questions: do the partial deep events at 82 and 219 yr persist; does the full-to-full interval converge. If both resolutions agree with each other but not with 250 m, 250 m was underresolved. If all three agree, the box is the cause of the mismatch with the BP1 reference.

6. **Box size.** At the resolution found adequate in item 5, run a 60 x 100 km box and compare the catalogue. This tests the Dirichlet-loading hypothesis directly. If the partial events disappear and the recurrence approaches 78 yr, the remaining differences from the BP1 reference are the half-space truncation and the first-order state treatment.

7. **Iteration scaling at high resolution.** Record iterations per step (interseismic and coseismic) at 125 m, 60 m, and 30 m. The linear scaling with n observed at 80 and 200 nodes is the main cost driver of the method; confirming or breaking it at large n decides whether a multigrid-type acceleration is needed.

8. **Tolerance at high resolution.** Repeat the 1 Pa vs 0.1 Pa comparison at 60 m from a common post-event state. The 250 m result may not carry over when stress gradients along the fault are steeper.

9. **Time integration order.** The state is first order in time. A predictor-corrector (Heun) on slip and state would double the PT solves per step. Test at 125 m whether it changes the catalogue; if not, keep first order.

Items 5 and 6 are the ones that determine whether the method reproduces BP1. Items 2, 3, and 7 determine what it costs.

## 5. Open assessment

The elliptic PT solve was never the limiting factor in these runs: it converged reliably with warm starts, and the tolerance question that looked serious at 1 km turned out to be an unstable-initial-condition effect. The costs that remain are the ones shared with any volume method for rate-and-state faults (resolution set by h*, time steps set by Dc / V) multiplied by the PT method's linear iteration scaling with resolution. Whether that product is acceptable is what the GPU tests in section 4 will show. The structural advantage of the approach, switching to physical inertia during rupture within the same code, has not been exercised yet.
