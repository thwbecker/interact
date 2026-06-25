# One-step scaling test, 40000 patches

This note documents the strong-scaling test in `test_matrix_inverse_scaling/` and
summarizes the 40000 patch results collected on the node `theo5`. The numbers below
are specific to that machine and to this build and PETSc configuration; a different
node, compiler, or PETSc options file may give different absolute times and
efficiencies.

## What is tested

Each run is a single one-step `interact` problem on a planar fault discretized into
40000 boundary elements, solved across increasing MPI core counts (1, 2, 4, 8, 16, 32)
on a single node. Two geometries and two problem types are compared:

- Geometry: rectangular (quad) elements, and triangular elements (two triangles per
  quad). Both are run at the same total element count of 40000.
- Problem type: a forward case and an inverse case.
  - Forward (`ntest=1`): slip is prescribed on every patch, so the run evaluates the
    full patch-to-patch interaction (a Green's function evaluation) to resolve the
    fault-local fields. There is no linear solve. This is effectively a dense
    matrix-times-vector cost that grows with the number of patches and parallelizes
    directly.
  - Inverse (`ntest=2`): a stress boundary condition is prescribed, so the run
    assembles the dense interaction matrix and solves the linear system for slip.

## How the runs are produced

The driver is `run_new_benchmark`. For each combination it does the following.

Geometry is built with `makefault`. Quads use `makefault -n 200 -m 200`. Triangles use
`makefault -n 200 -m 100 | patchquad2patchtri 2`, which yields the same 40000 element
total. A one-step `bc.in` selects the case: `-1 -1 -1 0 1` for the forward (prescribed
slip) case and `-1 -1 -1 10 -2e3` for the inverse (stress boundary condition) case.

`interact` is then run under `mpirun --bind-to core -np $ncore` with three options:

- `-fpetsc` forces the PETSc solver path at every core count. Without it `interact`
  uses a serial direct solver for serial runs, which would make the single-core
  baseline a different solver from the parallel runs and distort the scaling.
- `-npsfse` skips the post-solve resolved-stress evaluation. That step is a serial
  Green's function evaluation that is not needed when only the slip solution is
  wanted, and leaving it in previously dominated the inverse runtime and hid the
  scaling of assembly and solve.
- `-log_view` emits the PETSc performance table.

`OMP_NUM_THREADS` is set to 1 so the sweep measures pure MPI scaling. The PETSc options
file used here selects an `fgmres` Krylov method with a `jacobi` preconditioner.

## How timing is recorded

`stats.theo5.log` carries, per run, a total wall time and a per-phase breakdown.

- `total_wall_s` is the whole process time from `/usr/bin/time`.
- `asm_wall_s` and `solve_wall_s` come from `interact`'s own wall-clock markers
  ("parallel assembly done" and "parallel solve completed") in the PETSc solve path.
- `ksp_solve_s`, `mat_mult_s`, `pc_apply_s`, `pc_setup_s` are the corresponding PETSc
  `-log_view` event times, which break the solve into the Krylov solve, the dense
  matrix-vector products, and the preconditioner.

## What the results show

On 32 cores, relative to the single-core baseline on theo5:

| case            | T(1) (s) | T(32) (s) | speedup | efficiency |
|-----------------|---------:|----------:|--------:|-----------:|
| forward, quad   |  1628.8  |    59.9   |  27.2x  |    85%     |
| forward, tri    |  7771.6  |   256.3   |  30.3x  |    95%     |
| inverse, quad   |  1795.6  |    71.6   |  25.1x  |    78%     |
| inverse, tri    |  7962.5  |   273.4   |  29.1x  |    91%     |
| assembly, quad  |  1727.9  |    60.7   |  28.5x  |    89%     |
| assembly, tri   |  7878.7  |   260.0   |  30.3x  |    95%     |

Three points stand out for this configuration.

The cost is assembly-dominated and the assembly scales well. The matrix build is about
96 percent of the inverse runtime at one core, and it scales close to ideal, reaching
roughly 89 percent efficiency for quads and 95 percent for triangles on 32 cores. The
total inherits this, so forward and inverse track each other closely.

Triangles scale slightly better than quads. The triangular kernel carries more work per
element, so the ratio of compute to communication and fixed overhead is higher, and the
efficiency is correspondingly higher (about 91 to 95 percent versus 78 to 89 percent on
32 cores here).

The linear solve is a small but weaker-scaling piece. The solve is almost entirely the
dense matrix-vector product inside `fgmres`: `mat_mult_s` and `ksp_solve_s` are nearly
equal, and the preconditioner is negligible. That matvec scales to about 6.5x on 32
cores (near 20 percent efficiency). Because assembly scales and the solve does not, the
solve share grows with core count, from about 4 percent of the quad inverse at one core
to about 15 percent at 32 cores, although it remains a minor fraction of the total at
this problem size.

In short, for 40000 patches on theo5 the assembly parallelizes well and sets the
overall good scaling of both the forward and the inverse problem, while the remaining
and expected limiter is the dense matrix-vector product in the iterative solve. If the
solve fraction is to be reduced further at larger problem sizes, the matvec is the part
to target, for instance through the hierarchical-matrix backends.
