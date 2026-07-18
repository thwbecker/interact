# BigWham backend: build + shim scaffolding (set up like HMMVP)

Base commit: `0e555b9`.

This is the staging step for the BigWham full-space H-matrix backend
(`use_hmatrix 5`): the repo subdirectory, a compile script, the C shim, and the
opt-in build wiring. It deliberately does NOT yet add the `use_hmatrix==5`
dispatch or the geometry / slip-projection mapping (that is the next step,
below). Everything here is off by default, so a normal interact build is
unchanged.

## Where each piece goes

- `BigWham/compile_for_interact.sh`  -> place inside the BigWham checkout at
  `interact/BigWham/`. You clone the library yourself:

  ```
  cd interact
  git clone https://github.com/GeoEnergyLab-EPFL/BigWham.git
  cd BigWham && ./compile_for_interact.sh      # -> build/libBigWham.a
  ```

  The script runs the CMake build (Python/Julia/Mathematica interfaces off,
  OpenMP on, `BigWhamStatic` target). Dependencies: CMake, a C++17 compiler,
  OpenMP, and a BLAS providing `cblas.h` + `lapacke.h` (Ubuntu/Debian:
  `apt-get install libopenblas-dev liblapacke-dev`; TACC/MKL:
  `IL_MATH_VENDOR=MKL_sequential ./compile_for_interact.sh`).

- `src/la_and_geo/bigwham_shim.cc`  -> the `extern "C"` wrapper around
  `BigWhamIO`, living next to `hmmvp_c_shim.cpp`. Unlike the hmmvp shim (which
  compresses interact's own Okada kernel via the `ckernel_func` callback),
  BigWham owns its kernel, so the shim takes a mesh `(coor, conn)` + elastic
  constants and exposes `cbigwham_create / _build / _mvp / _get_diagonal /
  _get_info / _delete`. The operator is BigWham's native full-space 3N x 3N
  rectangular DDM operator (3 slip DoF per element).

- `config/makefile.petsc`, `makefile`, `src/includes/petsc_prototypes.h`
  -> the build wiring, mirroring the HMMVP variables and rules. Apply with
  `git apply bigwham_build_wiring.patch` from the repo root (or use the whole
  files provided here).

## Enabling the build

BigWham is GPL v3, OpenMP-only, and not part of interact, so it is OFF by
default. To turn it on, uncomment the `BIGWHAM_*` lines in
`config/makefile.petsc` (after cloning + building BigWham as above). If
`cblas.h` is not on the default include path, add its directory to
`BIGWHAM_INC` (e.g. `-I/usr/include/x86_64-linux-gnu` on Debian/Ubuntu).

## Verified in a clean container

- The shim compiles against `libBigWham.a` and, driven through its C API on a
  6x6 `3DR0-H` mesh (108 DoF), reproduces the direct `BigWhamIO` results
  exactly (compression 1.75, ||A*ones|| 1.76871, diag[0] 0.420148).
- With BigWham OFF (default), `compress_interaction_matrix` builds and links
  unchanged.
- With BigWham ON (variables set), the shim builds through the new makefile
  rule and links into the binary (the `cbigwham_*` symbols are present).

These are specific to the tested configuration (this compiler, OpenBLAS,
BigWham v0.2.0); other versions or BLAS vendors may behave differently.

## Licensing

BigWham is GPL v3. Linking it into interact makes the combined, distributed
work subject to GPL v3. The backend is kept compile-time optional
(`USE_BIGWHAM`, off by default) so a non-GPL interact build simply omits it.

## Next step (not in this delivery)

Wire `use_hmatrix==5` into `compress_interaction_matrix.c` (then
`petsc_interact` / `rsf_solve`):

1. Build `(coor, conn)` for BigWham `3DR0` rectangles from interact's
   patch geometry (centre + strike/dip + half-lengths -> 4 corner nodes),
   and convert elastic constants (`E = 2 G (1 + nu)`).
2. Resolve the convention mapping between interact's single slip-mode /
   stress-mode (N) and BigWham's 3-component per-element slip / traction
   (3N), using `GetRotationMatrix` / `GetElementNormals` to align frames.
3. Wrap `cbigwham_mvp` in a PETSc `MATSHELL` and use `KSP` for the solve,
   with `cbigwham_get_diagonal` for a Jacobi preconditioner.
4. Validate against the native full-space Okada operator using
   `full_vs_half_test.sh` with `ENABLE_BIGWHAM=1`: in FULL mode the
   `relerr` column for bigwham is the apples-to-apples error vs the
   `-full_space 1` dense reference. BigWham is full-space only and must
   never be compared to the half-space operator.
