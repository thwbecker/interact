# Sandbox build used for the C port (hand-off notes)

How the development environment was set up (2026-10-17) so that every patch
could be compiled, unit tested and run through the bit-identity harness
before delivery. One core, 3 GB, Ubuntu 24.04, root.

## PETSc 3.25.5

    apt-get install gfortran liblapack-dev libblas-dev libopenmpi-dev openmpi-bin make bison flex \
        python3-numpy python3-scipy python3-sympy
    pip download --no-deps -d /tmp/x petsc==3.25.5      # the sdist is the plain PETSc tarball
    tar xzf /tmp/x/petsc-3.25.5.tar.gz; mv petsc-3.25.5 /opt/petsc/src
    cd /opt/petsc/src
    python3 ./configure --prefix=/opt/petsc/build --with-debugging=0 --with-fc=0 --with-cxx=0 \
        --with-mpi-dir=/usr --with-blas-lapack-lib="-llapack -lblas" COPTFLAGS=-O2 \
        --with-shared-libraries=1 --download-hwloc=0
    make PETSC_DIR=/opt/petsc/src PETSC_ARCH=arch-linux-c-opt all      # about 20 min on one core
    make PETSC_DIR=/opt/petsc/src PETSC_ARCH=arch-linux-c-opt install

No Fortran bindings (interact's Fortran parts use gfortran directly), no
external packages, no Htool/H2OPUS.

## interact

    export PETSC_DIR=/opt/petsc/build PETSC_ARCH=
    cp config/makefile.gcc_petsc config/makefile.gcc_petsc_sandbox

and in the copy: `HACAPK_DEFINES`, `HACAPK_LIBS`, `HMMVP_DEFINES`,
`HMMVP_INC`, `HMMVP_LIBS`, `HMMVP_OBJS` emptied, `PETSC_HTOOL_USED` emptied
(`petsc_matrix_operations.c` otherwise needs `MatHtoolKernelFn`, which this
PETSc does not provide), `-march=x86-64-v3` replaced by `-march=native`;
then the `include config/makefile.gcc_petsc` line of the makefile pointed at
the copy. Neither the copy nor the makefile line is part of any patch.

    make bin/rsf_solve          # about 4 min
    make rsf_testing            # bin/test_fd_fft, test_fd_segments, test_fd_tables, fd_dynamic
    src/testing/rsf_solve/run_tests.sh
    src/testing/rsf_solve/test_fd_sequence.sh

## Harness

A clean master is built in a second tree and its binary kept as the
reference, then

    REFBIN=/path/to/rsf_solve.ref fully_dynamic/c_port/regress_rsf.sh x 2

(the worktree route of the harness also works once `hmat/` is linked in,
which the script does). Identity has passed on every patch; the timing
column is noise at the run lengths that fit the sandbox (0.5 to 3 s).

## Test inputs used

- `/tmp/bp1c`: BP1 geometry, friction and initial-condition files written by
  the snippet in `rsf_solve_fd.md`; QD and FD runs to 520 yr.
- `/tmp/cyc`: the prototype's single-fault case (the sequence regression
  script writes the same files).
- `/tmp/step`: two parallel segments (`gen_step.py` in the script of
  2026-10-18, now `src/testing/rsf_solve/gen_step.py`).

## Prototype

`fully_dynamic/*.py` runs with python3, numpy, scipy, sympy, matplotlib;
`fully_dynamic/ref/` (untracked) holds the reference series. The prototype
dumps used by the C tests are produced by `dump_tables.py` and by
`run_tpv205.py`/`fd_results` (`ref_tpv205_dx100.txt`).

## Reset of 2026-10-24

The sandbox lost all working trees, PETSc and run directories between two
turns; only the delivered patches survived. Rebuilding with the recipe
above took about 35 minutes (PETSc 25 min on one core, interact 4 min);
the reference binary for the harness is a second tree at the master
commit with the same configuration. Work not yet in a patch (the graded
surface rows) had to be re-implemented, so: patch early, patch often.
