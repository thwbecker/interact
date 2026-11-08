# Compile notes (internal)

Notes on building `rsf_solve` for the scripts in this directory, as
done so far. Only the Ubuntu builds below have been tested; other
machines (e.g. the TACC makefiles in `config/`) have not been tried
with these scripts.

## Requirements

- PETSc 3.19 or newer (`PetscCall`, `PETSC_SUCCESS`); 3.15 is too old.
  For PETSc older than 3.21, `src/rsf_solve/rsf_fd.c` falls back from
  `TSSetPostEventStep` to `TSSetPostEventIntervalStep`. That fallback
  only affects the 3-D fully dynamic path, which these scripts do not
  use, and it is untested there. The 3.21 cutoff has not been checked
  against the PETSc changelog.
- MPI, gfortran, BLAS and LAPACK.
- PETSc's own H-matrix packages (HTOOL, H2OPUS) are not needed. HACApK
  and hmmvp, which come with interact, are optional and can be built
  as below.
- Python 3 with numpy, matplotlib, scipy; h5py for the planned HDF5
  export.
- GNU time (`/usr/bin/time`), which `run_model` uses.

## Ubuntu 24.04 with the distribution PETSc (sandbox, root)

    apt-get update
    apt-get install -y --no-install-recommends petsc-dev gfortran \
        libopenmpi-dev openmpi-bin liblapack-dev libblas-dev time
    pip install h5py --break-system-packages

This gives PETSc 3.19.6 in
`/usr/lib/petscdir/petsc3.19/x86_64-linux-gnu-real` and OpenMPI 4.1.6.

Optional H-matrix libraries, from the top of the interact tree:

    (cd hmat/HACApK/v.1.0.0/C_interface && mkdir -p objects && make)
    (cd hmat/hmmvp && mkdir -p lib bin && make mode=mpi libhmmvp)

HACApK's test `test_hacpak_c` should run (H against dense solve
difference about 2.5e-5 at ztol 1e-6). It does not call
`MPI_Finalize`, so under mpirun it reports an abnormal exit; run it
directly. The hmmvp changes for newer MPI described in
`hmat/hmmvp/README_hmmvp_mpi_build.md` are already in the repository.

`rsf_solve`, without PETSc's H-matrix packages:

    P=/usr/lib/petscdir/petsc3.19/x86_64-linux-gnu-real
    make bin/rsf_solve PETSC_DIR=$P PETSC_ARCH= PETSC_HTOOL_USED= \
        SCARGS="-Wno-discarded-qualifiers" SFARGS="-ffixed-line-length-132" CXXFLAGS=

The overrides are needed with the Debian/Ubuntu PETSc package. Its
pkg-config file sets `cflags_extra`, `fflags_extra` and
`cxxflags_extra` to the distribution build flags, including
`-flto=auto -O2 -g` and hardening options. With LTO the link fails
with a TLS mismatch on the OpenMP threadprivate common block `/C0/` in
`src/green/dc3dmod.F`. Overriding the three variables removes those
flags; `OPTIM_FLAGS` from `config/makefile.gcc_petsc` still applies.
Note that in `config/makefile.gcc_petsc` the pkg-config value of
`SFARGS` replaces the file's own `SFARGS` in any PETSc build.

If HACApK or hmmvp were not built, also pass empty values for their
variables, e.g. `HACAPK_DEFINES= HACAPK_LIBS= HMMVP_DEFINES=
HMMVP_INC= HMMVP_LIBS= HMMVP_OBJS=` (not tested in this form).

Check:

    bin/rsf_solve -h 2>&1 | grep -E "rsf_series|sigma_dot_interact"

## Ubuntu 22.04 with PETSc from source (first session)

The distribution PETSc 3.15 is too old. PETSc 3.23.4 was cloned from
GitHub and configured with `--with-cxx=0 --with-fc=0
--with-debugging=0` against the system BLAS/LAPACK and an MPI. In
`config/makefile.gcc_petsc`, `PETSC_HTOOL_USED` was emptied and the
HACApK and hmmvp variables commented out. In a sandbox without root,
gfortran and OpenMPI were extracted from downloaded packages with
`dpkg-deb -x` and used through wrapper scripts; PETSc's `make all`
had to be split (`make -j2 libs`, twice) to fit the time limit per
call.

## Running

`run_model` sources an optional `config.sh` next to it, for example

    rsf_solve=/path/to/interact/bin/rsf_solve
    mpirun="mpirun --allow-run-as-root --oversubscribe"

(the mpirun options were needed in the sandbox, which runs as root on
one core). In the sandbox, long runs started in the background have to
be detached with `setsid`; a plain `nohup ... &` job was killed when
the shell call that started it ended. Check free disk space before
large runs (`df -h`); the per-cell series can be hundreds of MB per
member (README.md, Series output).

## Reference checks

- `./run_model s1 3000 1 --rot 0`: 9 events, first onset 69.587 yr,
  recurrence 327.04 yr from the fourth event on, Mw 7.394. Two builds
  (PETSc 3.23.4 and 3.19.6) agreed to these digits, with RHS
  evaluation counts that differed by 0.3 percent.
- `examples/make_examples 3000 4 rundir`: single gives 14 events, Mw
  7.2, recurrence 218 yr.
- With all new options at their defaults, the series and catalog
  files are byte-identical to those written before the options were
  added.
