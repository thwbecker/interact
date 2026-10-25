# Hand-off: fully dynamic extension of rsf_solve

Written 2026-10-25 so that the work can be resumed from a clean machine.
The current state, with numbers, is in `fd_plan.md` section 0; user-facing
documentation in `../rsf_solve_fd.md`; design notes of the C code in
`fd_c_design.md`; the sandbox build recipe in `c_port/SANDBOX.md`.

## 1. Where things are

- C code: `src/rsf_solve/fd_*.c` (2-D marcher, kernels, tables, FFT, band
  solves), `src/rsf_solve/fd3d.c` (3-D marcher: tables, mirror, opening
  layer, full surface, long steps), `src/rsf_solve/rsf_fd.c` (hand-over
  with the QD integrator, setup from interact's geometry, options), headers
  `src/includes/fd*.h`. All FD code is inactive without `-fd_mode 1`.
- Tests: `src/testing/rsf_solve/run_tests.sh` (unit level, about 5 min:
  tables, FFT, segments, TPV205 and Harris-Day against the prototype, 3-D
  kernels, mirror, opening layer, full surface, surface stability,
  quasi-2-D dipping thrust), `test_fd_sequence.sh` (2-D cycle),
  `test_fd3d_gate.sh` (3-D statics against Is), `test_fd3d_dip.sh`
  (dipping thrust in rsf_solve), and the QD bit-identity harness
  `fully_dynamic/c_port/regress_rsf.sh` (needs a reference binary built
  from the master, `REFBIN=...`). Reference series in
  `src/testing/rsf_solve/ref/`.
- Python prototype: `fully_dynamic/*.py` (2-D: `shells.py`, `stepper.py`,
  `inplane*.py`, `antiplane.py`, `densesys.py` for arbitrary element sets
  with a virtual surface; 3-D: `tada_quad.py`, `fs3d_*.py`). Every C
  component was first built and checked here; the prototype is the
  reference for roundoff comparisons.
- Benchmarks and scripts: `seas_bp_tests/bp1`, `bp4`, `bp5` (generators and
  run scripts), `fully_dynamic/c_port/run_*.sh`, `fully_dynamic/ref/`
  (SEAS reference series, untracked).

## 2. How to build and check

    export PETSC_DIR=... PETSC_ARCH=...
    make bin/rsf_solve rsf_testing
    src/testing/rsf_solve/run_tests.sh
    src/testing/rsf_solve/test_fd_sequence.sh
    src/testing/rsf_solve/test_fd3d_gate.sh
    src/testing/rsf_solve/test_fd3d_dip.sh
    REFBIN=/path/to/master/rsf_solve fully_dynamic/c_port/regress_rsf.sh x 2

All of these passed on 2026-10-25. On the sandbox the build needs the
configuration described in `c_port/SANDBOX.md`; neither it nor the
makefile change is part of the code.

## 3. Typical runs

- 2-D BP1-FD: `fully_dynamic/c_port/run_bp1_fd.sh 50 440 6000`
- 2-D stepover sequences: `fully_dynamic/c_port/run_step_fd.sh`,
  `run_many_stepover`, `summarize_rsf_step.py`
- 3-D BP5-FD first event (mirror): see `seas_bp_tests/bp5`; add
  `-fd_mode 1 -fd_back_eps 0.05 -fd_stations <marcher indices>`; with
  `-fd_layer 2` the opening layer; with `-fd_surface 1` the full surface.
- 3-D dipping fault: a planar geometry of equal squares in a half space
  with `-rsf_slip_mode 1` (dip slip) or 0 (strike slip), `-fd_mode 1`,
  `-calc_sigma_dot 1` for normal-stress coupling; the surface is built
  automatically (`-fd_surface_ext`, `-fd_surface_m`, `-fd_surface_grade`,
  `-fd_surface_extend`). Station files `fd_phase_<n>_st<i>.dat`.

## 4. Conventions that bite

- The 3-D marcher's element order is `i*nz + j` (i along grid axis 1,
  j rows by depth); `-fd_stations` takes these indices. In surface mode axis
  1 is the strike; otherwise axis 1 is the slip direction.
- Normal stress: interact's `In` is compression positive; the tables'
  normal traction `tn` is tension positive; the marcher uses
  `sigma = sigma0 - tn` with `tn_static = -In delta`.
- Dip slip: in the surface code positive slip moves the hanging wall up
  dip (thrust); interact's sign is mapped by `fsign` in `rsf_fd.c`.
- Surface rows: widths must keep collocation points off other elements'
  edge lines unless exactly on them (powers of two of the finest width do).
- Time scheme: `-fd_nl 1` (default) is second order; the first-order
  scheme is stable only near cS dt/dx = 0.5 in mode II. The 3-D surface
  needs cS dt/hs away from 0.85 to 1.05 (default 2).
- Patches: every change was delivered as `git apply --index` patches
  against the GitHub master; the last cumulative one is
  `fd_fully_dynamic_2026-10-24.patch` (in the master as of 2026-10-25),
  followed by `fd_fs3d_b7.patch`.

## 5. What to do next (see fd_plan.md 0.4)

1. Fault-resolution study of a surface-breaking thrust at the trace,
   quasi-2-D against the 2-D prototype; then TPV10 when the SCEC data are
   available.
2. Newton-Krylov long-step solver in 3-D, then the coarse history.
3. Stepover ensembles with the second-order scheme.

## 6. Known weak points

- 3-D sequences are slow (fixed fine steps; the surface multiplies the
  cost by 2 to 12).
- Surface-breaking dipping faults carry a trace error that shrinks only
  with fault resolution.
- The 2-D C code has no non-parallel segments and no in-plane free
  surface (the prototype has both).
- TPV101 and TPV10 station comparisons are still owed.
