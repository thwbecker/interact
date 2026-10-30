# Hand-off: fully dynamic extension of rsf_solve

Updated 2026-10-29 (H7, D1, H8) so that the work can be resumed in a new
session or on another machine. The record of every step, with numbers, is
`fd_plan.md` (the sections G1 to H10 near its end describe the most recent
work); user-facing
documentation is in `../rsf_solve_fd.md`; design notes of the C code in
`fd_c_design.md`; the sandbox build recipe in `c_port/SANDBOX.md`.

## 1. Where things are

- C code: `src/rsf_solve/fd_*.c` (2-D marcher, kernels, tables, FFT, band
  solves), `src/rsf_solve/fd3d.c` (3-D: tables, one fault, images, virtual
  surface; several planes (G1); several faults with one virtual surface
  (G3b, `fd3d_msurf_*`); per-fault surface strips (H6, `fd3d_strips_*`,
  `fd3d_elem_resp`, compressed dense couplings `fd3d_dtab`)),
  `src/rsf_solve/rsf_fd.c` (hand-over with the QD integrator, setup from
  interact's geometry, options), headers `src/includes/fd*.h`. All FD code
  is inactive without `-fd_mode 1`.
- Tests: `src/testing/rsf_solve/run_tests.sh` (unit level, about 10 min,
  including the 3-D checks `fd3d_test multi`, `msurf`, `h6frames`,
  `h6strips`, `h6step`), `test_fd_sequence.sh`, `test_fd3d_gate.sh`,
  `test_fd3d_dip.sh`, `test_fd3d_pair.sh`, `test_fd3d_multi.sh`,
  `test_fd3d_msurf.sh`, `test_fd3d_mpi.sh` (rsf_solve on 1, 2 and 4 ranks,
  station files compared byte for byte), and the QD bit-identity harness
  `fully_dynamic/c_port/regress_rsf.sh` (needs a reference binary built from
  the master, `REFBIN=...`).
- Study drivers in `fd3d_test`: `h6strike <angle> <steps> [m2] [m1]`
  (resolution study of strips of different strikes), `tpv14 <dx> [14|15]
  [t_end] [gap_m]` (SCEC TPV14/15 with free-surface images).
- Python: prototypes in `fully_dynamic/*.py`; the gate studies of the H
  series in `fully_dynamic/h_gates/` (domain partitioning, compression,
  surface layouts).

## 2. How to build and check

    export PETSC_DIR=... PETSC_ARCH=...
    make bin/rsf_solve rsf_testing
    src/testing/rsf_solve/run_tests.sh
    src/testing/rsf_solve/test_fd_sequence.sh     (and the other test_fd3d_*.sh)
    REFBIN=/path/to/master/rsf_solve fully_dynamic/c_port/regress_rsf.sh x 2

All of these passed on 2026-10-29 (H7). Results at round-off level can change with
the compiler, the flags or a reorganisation of the code (fused multiply-add);
`-ffp-contract=off` makes builds reproduce each other exactly.

## 3. Half space in 3-D: which method is used

- One plane: images (vertical strike slip) or the virtual surface (dipping).
- Several vertical strike-slip planes: images (G3a).
- Several dipping planes sharing one strike: one virtual surface (G3b).
- Several planes of different strikes that need a surface: a strip per
  fault (H6), chosen automatically; `-fd_strips 1` forces it also for a
  shared strike. Options: `-fd_surface_band <m>` (half width of a strip,
  default the fault's depth extent), `-fd_surface_m` (cells dx/m).
- Any geometry: `-fd_hs_static 1` (whole-space dynamics, exact half-space
  statics; adequate for buried faults, not for slip rates near a trace).

## 4. Conventions that bite

- Element order in the 3-D marchers: `i*nz + j` (i along strike, j down
  dip); `-fd_stations` takes these indices.
- Normal stress: interact's `In` is compression positive; the tables' `tn`
  is tension positive; `sigma = sigma0 - tn`.
- Frames of the several-plane code: e1 strike, e2 up dip, e3 = e1 x e2; the
  marcher's slip is that of the +e3 side along e1 or e2 (in TPV14's frames
  right-lateral is negative).
- Surface cells: cS dt/hs must stay away from 0.85 to 1.05, so cells of
  dx/2 at the default time step are unstable; use dx/3, dx/4 or finer.
  Traces and fault ends must lie on cell boundaries.
- Strips of different strikes: cells are kept for the nearest trace, and
  where turned cells overlap the one closer to the boundary is masked.
  Faults closer than a cell or so (crossing, branching) are not covered by
  the strip checks; `h6strike` refuses such set-ups.
- C style of the code: declarations at the top of each routine, one
  statement per line, loops on separate lines.
- Dense lag tables are partitioned (D1, default on): changes to how they
  are built must keep the per-entry verification against the full
  increments; FD3D_DP=0 gives the unpartitioned reference, and
  src/testing/rsf_solve/fd3d_dumpdiff.py compares exact dumps.
- MPI: every 3-D path splits its tables over the ranks; the dense couplings
  by receiver rows (H7). A receiver row is always computed whole on its
  owner and the results are gathered, so results do not depend on the
  number of ranks; any new code that touches a dense table must keep that
  rule (own rows only, then gather) or the 1/2/4-rank checks will fail.
  `FD3D_DUMP=<file>` makes fd3d_test write exact dumps for these checks.
- Patches: every change is delivered as a `git apply --index` patch against
  the current GitHub master.

## 5. What to do next (fd_plan.md, H7, D1, H8)

Done: MPI for the dense couplings and the strip solve (H7); domain
partitioning of the dense lag tables (D1); band block-Jacobi preconditioner
for large strip sets (H8).
1. Faults that meet in a half space (H9): corners behave like separated
   faults; branches from the middle of a trace converge poorly in shear at
   the junction cell. Next: the gap rule for the strip path, a TPV14-like
   dipping branch in rsf_solve.
1b. Free surface for deep faults (H10): surface strips need bands of the
   order of the fault depth (graded rows); the mirror image alone makes the
   rupture reach the top 3 km about 0.5 s late at 1 km (TPV14 main fault),
   so it is not a reference for near-surface comparisons.
1c. H11, measurements on a multi-core machine:
   fully_dynamic/run_h11_measurements.sh (regression, memory of the
   partitioned tables at 1000 to 250 m, MPI speed-up at 500 m, graded strips
   at 750 and 500 m, the band preconditioned strip solve); its summary.txt
   is what the next session needs.
2. Memory and speed-up measurements on a multi-core machine: D1 at 500 and
   250 m, H7 scaling, H8 above 6000 strip unknowns (the sandbox has one core).
3. Benchmarks: TPV14/15 at 250 m, TPV22/23 (needs forced nucleation),
   TPV10 and TPV36/37 at 250 m, TPV102.

## 6. Known weak points

- Memory and build time of the dense couplings: TPV14 at 500 m needs about
  10 GB; a small strip pair at dx/4 did not build within the sandbox limits.
- The strip LU is dense and assembled on every rank (na^2 doubles twice
  during assembly); beyond test sizes the GMRES path needs a better
  preconditioner.
- TPV14 is not yet reproduced (the branch ruptures fully or not at all,
  depending on the junction gap); the resolution so far is 5 to 10 times
  coarser than the benchmark's.
- Faults that meet (junctions) have no special treatment of the shared edge.
- 3-D sequences are slow (fixed fine steps).
