# SCEC SEAS benchmark problems with `rsf_solve`

Community code-verification benchmarks from the SCEC initiative on
Sequences of Earthquakes and Aseismic Slip (SEAS), set up for
interact's `rsf_solve`.  Problem descriptions and the community
comparison platform: https://strike.scec.org/cvws/seas/ .

Each case directory is self-contained (generator, driver, notes);
the scripts at this level are shared helpers used mainly by the BP5
workflow.  All drivers take positional arguments with defaults and
no environment-variable indirection.

## Case directories

| directory | problem | geometry | status |
|----|----|----|----|
| `bp5/` | BP5-QD | 3-D rectangular vertical strike-slip, half space | reference case; matches HBI recurrence to ~0.01 yr at 1 and 2 km |
| `bp4/` | BP4-QD | 3-D rectangular, whole space (BP5's whole-space sibling) | generator + driver, no notes file yet |
| `bp3/` | BP3-QD | 2-D in-plane (plane strain), dipping fault, half space | dip 90/60/30, thrust and normal; see README_bp3.md |
| `bp_thrust/` | not a SEAS benchmark | curved (listric) surface-breaking thrust, Ozawa et al. (2023) in spirit | exercises `-rsf_slip_mode 1` and `-calc_sigma_dot` |

**BP1-QD is not here.**  The 2-D ANTIPLANE benchmark is part of the
antiplane cycle suite, `../test_twod/anti_plane_cycles/`, because
that directory is also the general multi-fault antiplane experiment
machinery and shares its generator (`gen_cycles2d.py profile bp1`).
Community anchor there: recurrence 78.2 yr at fine resolution
against Erickson et al. (SRL 2020).

## Quick start

    cd bp5   && less README.md          # the fullest documented case
    cd bp3   && ./run_bp3 0.025 60 1e-9 bp3_d60_thrust 1500 1 \
                  "-calc_sigma_dot -limit_sigma 1 -min_sigma 1e6 -max_sigma 2e8"
             && python3 plot_bp3.py bp3_d60_thrust
    cd bp4   && ./run_bp4.sh
    cd bp_thrust && ./test_thrust.sh

`bp3/run_bp3` and the BP5 workflow are idempotent: re-running the
same output directory resumes from the last checkpoint, and a
completed run exits as a no-op.

## Shared helper scripts (this directory)

| script | what it does |
|----|----|
| `create_run` | assemble a BP5-style run directory (geometry, friction, initial conditions, driver options) |
| `run_many` | run a list of BP5 variants (evolution laws, IMEX vs RK) in sequence |
| `event_stats` | print the Mw >= mw_min event sequence and its recurrence statistics for a run |
| `check_events` | sanity-check an event catalog |
| `combine_groups` | merge per-group monitor files `rsf_monitor.gGGG.dat` into one time series |
| `compare_evol_laws` | compare BP5 results across state-evolution laws |
| `find_directory` | locate the run directory matching a parameter set |

## Benchmark-specific `rsf_solve` features

These options exist largely because of the benchmarks; see
`../rsf_solve.md` for the full reference.

- `-rsf_ic_file` / `-rsf_dc_file` / `-rsf_sigma_file`: per-cell
  initial conditions, `D_RS`, and normal stress (all benchmarks
  prescribe cell-dependent values).
- `-rsf_slip_mode 1` and `-calc_sigma_dot`: dip-slip source
  direction and normal-stress evolution, needed by BP3 and the
  thrust tests (BP3's effective normal stress includes the
  slip-induced change).
- `-rsf_stations FILE`, `-rsf_station_dt_yr`, `-rsf_station_vdense`:
  SEAS-style per-station time series `fltst_NAME.dat` with the
  benchmark's field order and sign conventions, written densely
  through events and sparsely between them.  `bp3/gen_bp3.py`
  writes the station file for the 12 official BP3 stations.
- `-rsf_checkpoint`, `-rsf_checkpoint_wall`, `-rsf_restart`: long
  benchmark runs are chained; the wallclock cadence is usually what
  you want.

## Notes on running these

- **Resolution matters more for some problems than others.**  BP5
  and BP1 are converged at their suggested cell sizes; BP3 is not
  converged until the suggested 25 m, and intermediate resolutions
  produce spurious period-2 event patterns (documented in
  `bp3/README_bp3.md`).
- **Integrator.**  The explicit RK path (`-ts_rk_type 3bs`) is
  usually fastest and most robust.  The exception is BP3 with
  normal-stress coupling, where trial stages can drive sigma
  through zero at the near-trace cells during large events; `-imex`
  integrates through those stretches (or raise
  `-domain_check_max_reject`).
- **Comparing against the community.**  The station files are
  directly comparable to the submissions on the SEAS platform,
  column for column.  Our slip sign follows the segment tangent, so
  check the sense convention first if curves come out mirrored.

## Viscoelastic extensions

The benchmarks are elastic.  Viscoelastic cycle work built on the
same machinery lives outside this directory:

- `../test_twod/anti_plane/` antiplane viscoelastic cycles,
  validation ledger and loading-condition map, Kato (2002) and
  Miyake & Noda (2019) replications.
- `../test_relax/inplane_ve_proto/` plane-strain viscoelastic
  kernels with gravity (the BP3 geometry over a Maxwell
  half-space), their analytic and literature gates, and the
  generator that writes `-ve_prony_file` kernels for
  `rsf_solve -ve_mode 3`.
- `../test_twod/ve_cycle_demos/` runnable demos reproducing
  published viscoelastic cycle systematics, including a BP3-like
  dipping thrust over a relaxing substrate.
