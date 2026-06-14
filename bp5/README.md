# SEAS BP5-QD (rectangular) with interact `rsf_solve`

This directory contains everything needed to run the SCEC SEAS **BP5-QD** benchmark
(rectangular vertical strike-slip fault, 3-D half-space) with interact's `rsf_solve`,
set up to mirror the HBI `bp5r` configuration and the published benchmark
(Jiang et al., 2022, *JGR Solid Earth*, doi:10.1029/2021JB023519;
[description](https://strike.scec.org/cvws/seas/download/SEAS_BP5.pdf)).

Unlike BP1 (2-D antiplane), BP5 is genuinely 3-D rectangular, so interact's 3-D Okada
half-space kernel is the *correct* Green's function — and interact and HBI agree on the
earthquake recurrence to **~0.01 yr** at both 1 km and 2 km resolution (see *Results*).

---

## 1. Files

| file | format | meaning |
|------|--------|---------|
| `make_bp5.py` | python | generator; emits the four input files for any cell size `ds` |
| `geom_bp5_1km.in` | `x y z strike dip L W group` | 4000-patch geometry (official 1 km resolution) |
| `rsf_bp5_1km.dat` | `a b` | per-patch rate-and-state `a`, `b` (geometry order) |
| `ic_bp5_1km.in`  | `tau[Pa] v[m/s]` | per-patch initial shear stress and slip velocity |
| `dc_bp5_1km.in`  | `dc[m]` | per-patch `D_RS` (0.14 m bulk, 0.13 m in the nucleation patch) |
| `*_2km.*` | — | coarse version (ds = 2 km, 1000 patches) |

All per-patch files are written in the **same patch order** as the geometry file
(`k = i*jmax + j`, `i` = along-strike index, `j` = depth index).

### Geometry convention

`strike = 0`, `dip = 90` ⇒ patch **normal = +x**, **strike direction = +y**, **dip = +z (up)**.
So BP5's along-strike coordinate maps to interact's **y**, depth to **z = −depth**, and the
fault plane sits at **x = 0**. `L`,`W` in the geometry file are **half-lengths** (a 1 km square
cell ⇒ `L = W = 500`).

### Physical parameters (baked into the generator)

| quantity | value |
|----------|-------|
| domain | 100 km (strike) × 40 km (depth), surface-breaking |
| VW core | `a = 0.004`, 60 km × 12 km, centered at 10 km depth |
| VS exterior | `a = 0.04`, 2 km linear taper |
| `b` | 0.03 (uniform) |
| `D_RS` | 0.14 m bulk, **0.13 m in the nucleation patch** |
| `f0`, `V0` | 0.6, 1e-6 m/s |
| `V_pl` = `V_init` | 1e-9 m/s |
| `sigma` | 25 MPa (constant) |
| `G`, `c_s` | 32.04 GPa, 3.464 km/s |
| nucleation patch | 12 km × 12 km at (x = −24 km, dep = 10 km), `V = 3e-2` m/s, `D_RS = 0.13` |

`tau0` per cell is the steady-state-at-`V_pl` stress evaluated at the local initial velocity
(Eq 16 of Jiang et al.). Note `tau0`, `V`, and the regularized state
`psi = a·ln(2V*/V·sinh(tau/(sigma a)))` are all **independent of `D_RS`** — `D_RS` enters only
the aging-law *rate* (`dpsi/dt ∝ 1/D_RS`), which is why the reduced-`D_RS` patch is supplied as
a separate file rather than changing the initial conditions.

---

## 2. Per-cell overrides in `rsf_solve` (two added options)

Stock `rsf_solve` read only per-cell `a`, `b` (via `rsf.dat`); `sigma`, `tau`, `V`, `D_c` were
uniform. Two small, non-breaking options were added so BP5 can be run faithfully:

```
-rsf_ic_file <file>    # per patch: tau[Pa]  v[m/s]   (initial conditions; nucleation patch)
-rsf_dc_file <file>    # per patch: D_c[m]            (overrides uniform -dc; nucleation D_RS)
```

When `-rsf_ic_file` is given it overrides the uniform `-tau_init`/`-vel_init` cell by cell; the
state `psi` is then computed per cell from `(tau, sigma, V, a)`. When `-rsf_dc_file` is given,
the per-cell `D_c` replaces the scalar `-dc` everywhere it is used (the aging law and the
cohesive-zone check). Omit either and the old uniform behavior is recovered exactly. `b` and
`sigma` remain uniform (BP5 uses `b = 0.03`, `sigma = 25 MPa` everywhere).

---

## 3. Running the benchmark

```bash
mpirun -np 1 bin/rsf_solve \
  -geom_file   geom_bp5_1km.in \
  -rsf_file    rsf_bp5_1km.dat \
  -rsf_ic_file ic_bp5_1km.in \
  -rsf_dc_file dc_bp5_1km.in \
  -use_hmatrix 3 -hacapk_ztol 1e-4 \
  -shear_modulus 3.204e10 -s_wave_speed 3464 \
  -f0 0.6 -dc 0.14 -vpl 1e-9 -v0 1e-6 -sigma_init 25e6 \
  -rtol 1e-4 -stop_time_yr 1800 -ts_max_steps 2000000 \
  -print_interval_yr 0.02 -log_view
```

Notes:
- `-use_hmatrix 3` = HACApK (≈25 MB, 21 % compression at 4000 cells). `0` = dense, `1` = HTOOL;
  the cycle is identical across all three (dense ≡ HACApK to machine precision).
- `-dc 0.14` sets the bulk value; `-rsf_dc_file` overrides it in the nucleation patch.
- `-sigma_init` and `-shear_modulus` are in **Pa**.
- Output `rsf_monitor.dat`: `step  t[s]  t[yr]  dt[s]  log10(max|V|)  mean_slip …`.
- The nucleation patch fires an event at `t ≈ 0`; the first spontaneous recurrence is near
  235 yr, and the benchmark runs to `t_f = 1800 yr` (~8 events). The `D_RS = 0.13` patch pins
  re-nucleation to the same spot, so the sequence stays periodic.
- The 0.02-yr monitor cadence resolves the long-term cycle but **not** the ~30 s coseismic
  phase; for per-event rupture duration / stress drop, monitor far more finely across one event.

### Other resolutions

```bash
python3 make_bp5.py 0.5     # 200x80 = 16000 cells (refined)
python3 make_bp5.py 2.0     #  50x20 =  1000 cells (coarse)
```

`ds = 1 km` is the official benchmark (cohesive zone `L_b = G·D_RS/(b·sigma) ≈ 6 km`, ~6 cells
across it). 2 km under-resolves it but is fine for fast code-vs-code checks.

---

## 4. Results (validation)

Spontaneous recurrence of the first natural event (the nucleation-driven event is at `t ≈ 0`):

| code / setup | resolution | recurrence | note |
|--------------|-----------|------------|------|
| **HBI** (lattice-H, patched) | 1 km | **234.28 yr** | reference |
| **interact** `rsf_solve` (HACApK) | 1 km | **234.29 yr** | matches HBI to 0.01 yr |
| **HBI** | 2 km | 236.81 yr | reference |
| **interact** | 2 km | 236.81 yr | exact match |

The published benchmark (Jiang et al., Fig 9a) reports inter-event times scattering around
**~235 yr** over the 1800-yr window, with the boundary-element codes (BICyclE, ESAM, FDRA, HBI,
TriBIE, Unicycle) clustered tightly and the two volume-discretized codes (EQsimu, GARNET)
sitting apart near ~250 yr. interact and HBI both land in the BEM cluster, as expected (same 3-D
Okada half-space kernel, same backslip loading).

**The `D_RS = 0.13` patch matters even for the first recurrence.** Without it (uniform 0.14 m),
interact gives 235.1 yr (1 km) and 233.5 yr (2 km); adding it brings both into ~0.01 yr agreement
with HBI. Robust coseismic quantities reported by the benchmark — rupture duration ~30 s, stress
drop ~5 MPa (surface) / ~10 MPa (VW core), peak slip rate ~1 m/s — are reproduced here too.

### The pre-event "noise" is a grid effect, not numerics

The `log10 max|V|` trace shows a fluctuating plateau (~−8.9, just above plate rate) for ~20 yr
before each event. We isolated its cause:

- **Not the H-matrix.** Dense and HACApK runs are bit-for-bit identical (same event time, same
  step-by-step trace).
- **Not the time integrator.** Tightening `rtol` from 1e-4 to 1e-6 (100×) leaves the detrended
  plateau amplitude unchanged (RMS ≈ 0.11, 0.18 reversals/yr at 2 km) — finer sampling, same
  fluctuation.
- **It is the spatial grid.** Refining 2 km → 1 km drops the plateau RMS 3.6× (0.109 → 0.030)
  and peak-to-peak 2.8×, while reversals/yr *rises* (0.18 → 0.61).

This is the signature of the `max`-over-cells reduction sampling a creep/nucleation front that
migrates across discrete cells: each cell-crossing is one hop, so finer cells give more, smaller
hops. Individual cell histories are smooth and event timing is unaffected. This is consistent
with Jiang et al., who find the analogous pre-event transients vanish with finer **grid** (their
BP4 at 500 m). The lever is grid spacing, not solver tolerance; plotting mean|V| instead of
max|V| also removes it.

---

## 5. Why BP5 originally crashed in HBI

Out of the box, HBI's `bp5r` aborts during H-matrix setup at 4000 cells with
`corrupted size vs. prev_size` (a glibc heap-corruption abort), *before* time-stepping. Root
cause, located with a bounds-checked (`-fcheck=all`) build:

- In the `param(61)==3` diagonal-scaling loop in `m_HACApK_use.f90`, the loop runs
  `do il = ndnr_s, ndnr_e` with `ndnr_s = lpmd(6) = 0` (the lattice-H path sets a 0-based
  index range). It then calls `HACApK_entry_ij(0,0)`, and the `3dph` kernel `okada_ij` reads
  `ang(0)` — one below the array's lower bound — and writes `st_bemv%ao(0)`.
- This 1-element out-of-bounds read+write fires at every mesh size; under `-O3` it silently
  corrupts the heap — benign at 1000 cells, fatal at 4000.

Minimal guard (`hbi_bp5_fix.patch`):

```fortran
   do il = max(ndnr_s,1), ndnr_e        ! was: do il = ndnr_s, ndnr_e
```

After that, HBI builds the 4000-cell H-matrix cleanly. One secondary point: the RK stepper
fails its first step at 1 km with the default `dtinit = 1d-1` (finer mesh is stiffer); reducing
`dtinit` to `1d-4` lets it proceed, after which HBI runs the full 1 km BP5 to completion (the
234.3 yr result above). The published Ozawa et al. (2021) HBI ran 1 km BP5 fine, so this is a
post-publication regression — worth flagging upstream.

**Note:** interact's own HACApK copy (`HACApK/v.1.0.0/C_interface`) has the identical loop and
takes the same `param(61)==3` branch, but its standard (non-lattice) generate path gives
`lpmd(6) = 1`, so it never hits index 0 (verified with a bounds-checked build). The same one-line
`max(ndnr_s,1)` guard is applied there defensively (`interact_m_HACApK_use.f90`).
