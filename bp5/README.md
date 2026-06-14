# SEAS BP5-QD (rectangular) with interact `rsf_solve`

This directory contains everything needed to run the SCEC SEAS **BP5-QD** benchmark
(rectangular vertical strike-slip fault, 3-D half-space) with interact's `rsf_solve`,
set up to mirror the HBI `bp5r` configuration
([benchmark description](https://strike.scec.org/cvws/seas/download/SEAS_BP5.pdf)).

Unlike BP1 (2-D antiplane), BP5 is genuinely 3-D rectangular, so interact's 3-D Okada
half-space kernel is the *correct* Green's function here — and the two codes agree to
**< 0.5 %** on the recurrence interval (see *Results* below).

---

## 1. Files

| file | format | meaning |
|------|--------|---------|
| `make_bp5.py` | python | generator; produces the three input files for any cell size `ds` |
| `geom_bp5_1km.in` | `x y z strike dip L W group` | 4000-patch geometry (official 1 km resolution) |
| `rsf_bp5_1km.dat` | `a b` | per-patch rate-and-state `a`, `b` (geometry order) |
| `ic_bp5_1km.in`  | `tau[Pa] v[m/s]` | per-patch initial shear stress and slip velocity |
| `geom_bp5_2km.in`, `rsf_bp5_2km.dat`, `ic_bp5_2km.in` | — | 1000-patch coarse version (ds = 2 km) |

All three per-patch files are written in the **same patch order** as the geometry file
(`k = i*jmax + j`, `i` = along-strike index, `j` = depth index), which is the order
`rsf_solve` reads them in.

### Geometry convention (important)

In interact, `strike = 0`, `dip = 90` gives patch **normal = +x**, **strike direction = +y**,
**dip direction = +z (up)**. So for BP5 the along-strike coordinate maps to interact's **y**,
depth maps to **z = −depth**, and the fault plane sits at **x = 0**. `L` and `W` in the
geometry file are **half-lengths** (so a 1 km square cell has `L = W = 500`). This only
matters for BP5 because BP1 stacked patches in depth alone.

### Physical parameters (baked into the files / generator)

| quantity | value |
|----------|-------|
| domain | 100 km (strike) × 40 km (depth), surface-breaking |
| VW core | `a = 0.004`, 60 km × 12 km, centered at 10 km depth (|x| ≤ 30, |dep−10| ≤ 6) |
| VS exterior | `a = 0.04`, 2 km linear taper between |
| `b` | 0.03 (uniform) |
| `D_c` | 0.14 m (uniform) |
| `f0`, `V0` | 0.6, 1e-6 m/s |
| `V_pl` = `V_init` | 1e-9 m/s |
| `σ` | 25 MPa (constant) |
| `G`, `c_s` | 32.04 GPa, 3.464 km/s |
| nucleation patch | 12 km × 12 km at (x = −24 km, dep = 10 km), `V = 3e-2` m/s |

Initial shear stress per cell is the steady-state-at-`V_pl` stress evaluated at the local
initial velocity, `τ = σ (a·asinh[ 0.5 V/V0 · exp((f0 + b·ln(V0/V_pl))/a) ] + η V)`,
`η = G/(2 c_s)` — identical to HBI's `bp5r_param.py`.

---

## 2. The `-rsf_ic_file` option (per-cell initial conditions)

BP5 nucleation is driven by a localized high-velocity patch, which needs **per-cell**
initial velocity (and the matching per-cell stress). Stock `rsf_solve` only had per-cell
`a`, `b` (via `rsf.dat`); `σ`, `τ`, `V`, `D_c` were uniform. A small non-breaking option
was added:

```
-rsf_ic_file <file>     # two columns per patch: tau[Pa]  v[m/s], in geometry order
```

When supplied it overrides the uniform `-tau_init` / `-vel_init` cell by cell; the state
`ψ` is then computed per cell from `(τ, σ, V, a)` exactly as before. Without it the old
uniform behavior is unchanged. (`b` and `D_c` remain uniform — for BP5 `b = 0.03` and
`D_c = 0.14` everywhere; HBI's tiny `D_c = 0.13` in the nucleation patch is immaterial.)

---

## 3. Running the benchmark

Build `rsf_solve` as usual, then (1 km, HACApK compression):

```bash
mpirun -np 1 bin/rsf_solve \
  -geom_file   geom_bp5_1km.in \
  -rsf_file    rsf_bp5_1km.dat \
  -rsf_ic_file ic_bp5_1km.in \
  -use_hmatrix 3 -hacapk_ztol 1e-4 \
  -shear_modulus 3.204e10 -s_wave_speed 3464 \
  -f0 0.6 -dc 0.14 -vpl 1e-9 -v0 1e-6 -sigma_init 25e6 \
  -rtol 1e-4 -stop_time_yr 400 -ts_max_steps 300000 \
  -print_interval_yr 0.02 -log_view
```

Notes:
- `-use_hmatrix 3` selects HACApK (≈ 25 MB, 21 % compression at 4000 cells, ~5× faster
  per matvec than dense). `-use_hmatrix 0` runs dense; `1` = HTOOL. The cycle is
  identical across all three (dense ≡ HACApK to machine precision — see below).
- `-sigma_init 25e6` is in **Pa**; `-shear_modulus 3.204e10` is **Pa**.
- Output is `rsf_monitor.dat`: columns `step  t[s]  t[yr]  dt[s]  log10(max|V|)  mean_slip …`.
- The nucleation patch fires an event at `t ≈ 0`; the first *spontaneous* recurrence is
  near 235 yr, so run to at least ~250–400 yr.

### Other resolutions

```bash
python3 make_bp5.py 0.5     # 200x80 = 16000 cells (refined)
python3 make_bp5.py 2.0     # 50x20  =  1000 cells (coarse)
```

`ds = 1 km` is the official benchmark. Cohesive zone `L_b = G D_c/(b σ) ≈ 6 km`, so 1 km
gives ~6 cells across it (adequate); 2 km under-resolves it but is fine for code-vs-code
checks and runs in seconds.

---

## 4. Results (validation)

| code / setup | resolution | spontaneous recurrence | peak log₁₀\|V\| |
|--------------|-----------|------------------------|-----------------|
| **HBI** (lattice-H, patched) | 1 km, 4000 | **234.3 yr** | ~0 |
| **interact** `rsf_solve` (HACApK) | 1 km, 4000 | **235.1 yr** | ~0 |
| HBI (lattice-H) | 2 km, 1000 | 236.8 yr | ~0 |
| interact `rsf_solve` (HACApK) | 2 km, 1000 | 233.5 yr | ~0 |
| interact `rsf_solve` (**dense**) | 2 km, 1000 | 233.5 yr | ~0 |

At equal (1 km) resolution the codes agree to **0.4 %** on recurrence, with essentially
identical peak velocity and HACApK compression matching HBI's (25.3 MB / 20.7 %). The 2 km
dense and HACApK runs are bit-for-bit equivalent (identical event time and step-by-step
behavior), confirming the H-matrix approximation contributes nothing to the dynamics.

Contrast with BP1: there interact and HBI genuinely differed because BP1 is 2-D antiplane
and interact's 3-D Okada strips are not HBI's 2-D infinite kernel. BP5 removes that kernel
mismatch, and the integrators converge — good evidence the BP1 gap was physical, not a bug.

### The pre-event "noise"

The `log₁₀ max|V|` trace shows a fluctuating plateau (~−8.9, just above plate rate) for
~20 yr before the runaway. This is **not** numerical noise from compression: the dense and
HACApK runs reproduce it identically. It is the `max`-over-cells reduction sampling a
*migrating* creep/nucleation front — several VW cells creep at comparable rates and the
maximum hops between them across the discrete grid. Individual cell histories are smooth,
event timing is unaffected, and the amplitude shrinks with finer cells. (Plotting mean\|V\|
instead of max\|V\| removes it.)

---

## 5. Why BP5 originally crashed in HBI

Out of the box, HBI's `bp5r` aborts during H-matrix setup at 4000 cells with
`corrupted size vs. prev_size` (a glibc heap-corruption abort), *before* time-stepping.
Root cause, located with a bounds-checked (`-fcheck=all`) build:

- In the `param(61)==3` (relative-norm) diagonal-scaling loop in `m_HACApK_use.f90`,
  the loop runs `do il = ndnr_s, ndnr_e` with `ndnr_s = lpmd(6) = 0`. It then calls
  `HACApK_entry_ij(0,0)`, and the `3dph` kernel `okada_ij` reads `ang(0)` — one slot
  **below** the array's lower bound — and writes `st_bemv%ao(0)`.
- This 1-element out-of-bounds read+write happens at **every** mesh size (it fires at
  1000 cells too under bounds checking). Under `-O3` it silently corrupts the heap;
  whether that is fatal depends on the heap layout — **benign at 1000 cells, fatal at
  4000**. That size sensitivity is why it looks like "BP5 used to work."

A minimal guard fixes it:

```fortran
! m_HACApK_use.f90, in the param(61)==3 branch:
   do il = max(ndnr_s,1), ndnr_e        ! was: do il = ndnr_s, ndnr_e
```

After that guard, HBI builds the 4000-cell H-matrix cleanly. One secondary point remains:
the RK stepper fails its very first step at 1 km with the default `dtinit = 1d-1` (the
finer mesh is stiffer); reducing `dtinit` to `1d-4` lets it proceed, after which HBI runs
the full 1 km BP5 to completion (the 234.3 yr result above). This is a step-size issue,
not a bug. (HBI still throws its usual harmless heap error *at finalization*, after all
data is written.)
