# One-step scaling test, 40000 patches

This note documents the strong-scaling test in `test_matrix_inverse_scaling/` and
summarizes the 40000 patch results collected on the node `walter`. The numbers below
are specific to that machine and to this build and PETSc configuration; a different
node, compiler, or PETSc options file may give different absolute times and
efficiencies. An earlier version of this note reported the same test on `theo5` up to
32 cores; the walter sweep below extends it to 64 cores and supersedes those numbers,
though the qualitative picture is unchanged.

## What is tested

Each run is a single one-step `interact` problem on a planar fault discretized into
40000 boundary elements, solved across increasing MPI core counts (1, 2, 4, 8, 16, 32,
64) on a single node. Two geometries and two problem types are compared:

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

`stats.walter.log` carries, per run, a total wall time and a per-phase breakdown.

- `total_wall_s` is the whole process time from `/usr/bin/time`.
- `asm_wall_s` and `solve_wall_s` come from `interact`'s own wall-clock markers
  ("parallel assembly done" and "parallel solve completed") in the PETSc solve path.
- `ksp_solve_s`, `mat_mult_s`, `pc_apply_s`, `pc_setup_s` are the corresponding PETSc
  `-log_view` event times, which break the solve into the Krylov solve, the dense
  matrix-vector products, and the preconditioner.

## What the results show

On 64 cores, relative to the single-core baseline on walter:

| case            | T(1) (s) | T(64) (s) | speedup | efficiency |
|-----------------|---------:|----------:|--------:|-----------:|
| forward, quad   |  1632.9  |    30.5   |  53.6x  |    84%     |
| forward, tri    |  7804.4  |   130.0   |  60.1x  |    94%     |
| inverse, quad   |  1795.8  |    36.9   |  48.6x  |    76%     |
| inverse, tri    |  8012.2  |   138.4   |  57.9x  |    90%     |
| assembly, quad  |  1725.1  |    30.7   |  56.2x  |    88%     |
| assembly, tri   |  7924.5  |   130.5   |  60.7x  |    95%     |

Three points stand out for this configuration, and they hold steady from 32 to 64
cores.

The cost is assembly-dominated and the assembly scales well. The matrix build is about
96 percent of the inverse runtime at one core, and it scales close to ideal, reaching
roughly 88 percent efficiency for quads and 95 percent for triangles on 64 cores (86 and
95 percent on 32). The total inherits this, so forward and inverse track each other
closely.

Triangles scale slightly better than quads. The triangular kernel carries more work per
element, so the ratio of compute to communication and fixed overhead is higher, and the
efficiency is correspondingly higher (about 90 to 95 percent versus 76 to 88 percent on
64 cores here).

The linear solve is a small but weaker-scaling piece, and on 64 cores its effect on the
quad total is now visible. The solve is almost entirely the dense matrix-vector product
inside `fgmres`: `mat_mult_s` and `ksp_solve_s` are nearly equal, and the preconditioner
is negligible. That matvec scales to about 13x on 64 cores (near 20 percent efficiency)
for both geometries. Because assembly scales and the solve does not, the solve share
grows with core count. For quads it rises from about 4 percent of the inverse at one core
to about 15 percent at 64, which is what pulls the quad inverse total efficiency (76
percent) below the quad assembly efficiency (88 percent). For triangles the assembly is
heavy enough per element that the solve stays under 5 percent even at 64 cores, so the
tri total stays close to the tri assembly.

In short, for 40000 patches on walter the assembly parallelizes well to 64 cores and
sets the overall good scaling of both the forward and the inverse problem, while the
remaining and expected limiter is the dense matrix-vector product in the iterative solve,
which begins to cap the quad inverse total around 64 cores. If the solve fraction is to
be reduced further at larger problem sizes or higher core counts, the matvec is the part
to target, for instance through the hierarchical-matrix backends.
