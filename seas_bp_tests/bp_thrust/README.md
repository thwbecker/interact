# Dip-slip thrust test for the normal-stress path

This directory sets up a surface-breaking, listric (curved) thrust for
`rsf_solve`, used to exercise and validate the dip-slip mode (`-rsf_slip_mode 1`)
and the normal-stress evolution (`-calc_sigma_dot`). It follows Ozawa et al.
(2023) in spirit rather than reproducing their run exactly: a thrust about 50 km
along strike and 20 km down the curved dip path, with the dip tapering from 30
degrees at the surface to 10 degrees at depth.

## Files

- `make_thrust.py` generates the inputs. All parameters are hardcoded at the top
  of the script. It writes:
  - `geom_thrust.in` (x y z strike dip L W group; meters and degrees; L and W are
    half-lengths),
  - `rsf_thrust.dat` (a, b per patch),
  - `ic_thrust.in` (tau[Pa], v[m/s] per patch, for `-rsf_ic_file`),
  - `sigma_thrust.in` (a depth-dependent sigma0[Pa] per patch, for the optional
    `-rsf_sigma_file` demo).
- `test_thrust.sh` builds the geometry and runs three checks, described below.
  All parameters are hardcoded at the top of the script.

## What the test shows

Running `bash test_thrust.sh` (with `PETSC_DIR` and `PETSC_ARCH` set, and after
`bin/rsf_solve` has been rebuilt with the Task 2 changes) does three things:

- A. With `-calc_sigma_dot` the `In` (normal-stress) matrix is built and the
  normal stress moves away from its uniform initial value once the seeded patch
  slips. On the coarse 1 km mesh the spread is around 1 MPa, with both increases
  and decreases, which is the behavior expected for dip slip on a bending,
  surface-breaking fault.
- B. Without the flag the `In` matrix is not built and the normal stress stays
  fixed at its initial value for the whole run, confirming the path is inert
  when off.
- C. With `-rsf_sigma_file` a depth-dependent initial normal stress is read and
  honored at the start of the run.

## Caveats

The default resolution (1 km cells) with an Ozawa-like `D_c` only marginally
resolves the cohesive zone, so this setup validates the dip-slip and
normal-stress machinery rather than providing a converged benchmark. For
quantitative work, reduce `ds` in `make_thrust.py` and match the parameter set
you want to compare against. Results quoted here are specific to this
configuration; a different resolution or kernel may shift them.
