# Solving the UCERF one-step interaction problem: tests and recommendations

Status as of 2026-09-06, interact master with the near-field
preconditioner (`-near_pc_rfac`) and the single assembly path in
compress_interaction_matrix. Problem considered throughout: strike slip
on all patches from a unit strike stress drop (bc code 10, value 1),
i.e. one N x N operator, strike slip -> strike shear stress. Numbers
are from one machine (theo6, 16 MPI ranks) and one sandbox (2 ranks);
they should transfer in ratios, not in absolute seconds.

## 1. Tools and what they solve

- `interact`: dense operator, assembled row by row on the fly
  (`par_assemble_a_matrix`, N^2/np Okada evaluations), PETSc KSP solve,
  optional post-slip stress evaluation (another N^2/np pass). General
  boundary conditions (any slip/stress mix, Coulomb correction).
- `compress_interaction_matrix`: H-matrix (or dense) operator through
  `calc_petsc_Isn_matrices`, solve of A x = 1 only, slip written to
  flt.dat, no stress evaluation. Solver options need the `htool_`
  prefix (all backends), e.g. `-htool_ksp_type gmres`. Un-prefixed
  `-ksp_*` options only affect the optional `-nsolve` timing solves.
- `-near_pc_rfac f` (both programs): sparse matrix of the operator
  entries with source-receiver distance below f mean patch lengths
  (sqrt of the mean patch area), used as PETSc preconditioning matrix.
  `-near_pc_radius r` gives the radius directly. Combine with
  `-pc_type asm -sub_pc_type lu` (interact) or `-htool_pc_type asm
  -htool_sub_pc_type lu` (compress). Do not derive the radius from
  column 6 of the geometry file: it is negative for triangles.

Test drivers: `ucerf/load_test/run_model` (interact, modes 1 to 7) and
`ucerf/solve_test/run_test` (compress, modes 1 to 4), both writing
`log.<mode>.dat`, `flt.<mode>.dat` and one line per mode to
`summary.dat`.

## 2. Measurements

### 2.1 Dense operator, interact, N = 100000 (UCERF subset), 16 ranks, rtol 1e-6

Assembly 3760 to 3780 s in every run. Solve:

    mode  KSP / PC                                   its   solve s  slip diff to mode 1
    1     fgmres(2000), Jacobi (previous setting)    490    733     -
    2     fgmres(2000), none                         501    749     8e-7
    3     gmres(30), near-field asm / sub lu          20     33     2e-6
    4     gmres(30), near-field asm / sub ilu         86    132     2e-6
    5     bcgs,      near-field asm / sub lu          14     34     2e-6
    6     gmres(30), near-field bjacobi / sub lu      73    113     2e-6
    7     as 3, with post-slip stress evaluation      20     33     (stress pass 3723 s)

Near-field matrix at rfac 4: 4.5 M nonzeros, 45 per row, 0.045 percent
of dense; per-rank sparse LU 0.44 s, one application 4 ms. One dense
matvec 0.74 s. The whole 33 s solve is 41 matvecs.

Reading: Jacobi does nothing (the diagonal is nearly uniform). The
near-field ASM/LU preconditioner cuts the iterations by 25x and the
solve time by 22x. ASM overlap matters (block Jacobi is 4x worse), ILU
is 4x worse than LU and the LU is cheap, so there is no reason to use
ILU. bcgs and gmres(30) are equivalent in time. After this, the solve
is 0.4 percent of a run with stress evaluation; the two N^2 Okada
passes are everything.

### 2.2 H matrix, compress_interaction_matrix, same geometry, 16 ranks

hmmvp, block-relative (BREM) tolerance 1e-3: assembly 112 s
(33x faster than dense), compression ratio 62, 1.2 GB instead of
80 GB.

    mode  KSP / PC                       its     solve CPU s
    1     fgmres(2000), none             501     16
    2     gmres(30), none              20000     544 (not converged, 37 percent off)

The unpreconditioned iteration count equals the dense one, i.e. the
operator is well approximated; one iteration costs about 30 ms
instead of 740 ms. Modes 3 and 4 (near-field ASM) were not yet run on
this geometry; on the 4000 patch case below they match the dense
counts (21 iterations).

### 2.3 Sandbox case, N = 4000, three faults with a 4 degree intersection, 2 ranks

Chosen because it reproduces the UCERF solver behaviour: gmres(30)
stagnates, fgmres(2000) needs about 400 iterations. Serial LAPACK LU
available as reference.

Operator accuracy (H solve to rtol 1e-8, rel L2 error of slip vs LU):

    hmmvp global norm (inorm 3) 1e-4    2.0e-2   ratio 6.8   (old run_test setting)
    hmmvp global norm (inorm 3) 1e-5    9.2e-4   ratio 4.9
    hmmvp block norm  (inorm 1) 1e-3    2.0e-4   ratio 4.2
    hmmvp block norm  (inorm 1) 1e-4    2.0e-5   ratio 3.1
    HACApK ztol 1e-4                    1.8e-7   ratio 1.45

The global-norm setting is the only one that gave a visibly smoother
slip field than LU, and it is the least accurate. The exact
(LU) solution itself has sign-reversed patches at the shallow
intersection; the smoothing seen with global-norm 1e-4 is operator
error damping the near-singular modes, not a better answer.

Dense solve accuracy vs LU: rtol 1e-4 gives 8e-3, 1e-6 gives 5e-5,
1e-8 gives 2e-7 (all fgmres/jacobi).

Preconditioner sweep (dense, rtol 1e-6): none 375 its; near-field
asm/lu at 2, 3, 4, 6 patch lengths: 94, 31, 21, 9 its; the same 21 at
4 patch lengths for hmmvp (BREM 1e-3) and HACApK (1e-4) operators.

### 2.4 Scaling to the full set (265k patches)

From the 100k measurements and N^2 scaling: dense assembly about
7.3 h on 16 ranks or 1.8 h on 64 (consistent with the 6835 s measured
earlier on 64), the same again for the stress pass; the dense matrix
is 560 GB and has to fit in distributed memory. H-matrix assembly
about 13 min on 16 ranks (a few minutes on 64, in line with the
150 to 400 s in ucerf/hmat_storage.dat); solve with the near-field
preconditioner a few seconds to a minute.

## 3. Best practices

- Always use the near-field preconditioner with ASM and sub-block LU.
  rfac 4 is a good default; rfac 6 halves the iterations again at a
  few times the (small) memory. Never bjacobi (no overlap) or ILU.
- gmres(30) is sufficient with the preconditioner; the restart-2000
  Krylov memory is not needed. bcgs is an equivalent alternative.
- Judge convergence in the unpreconditioned (true residual) norm:
  `-ksp_norm_type unpreconditioned` (or `-htool_ksp_norm_type ...`),
  and use rtol 1e-6 or tighter before interpreting fine-scale slip.
  rtol 1e-4 leaves percent-level errors in the slip.
- Always pass `-ksp_converged_reason` / `-htool_ksp_converged_reason`
  and check the reason string in the log; a non-converged solve still
  writes flt.dat.
- H matrix: hmmvp with `-hmmvp_inorm 1 -hmmvp_tol 1e-3` (block
  relative), or HACApK with `-hacapk_ztol 1e-4`. Do not use the
  global-norm hmmvp tolerance (inorm 3) for inverse problems.
- compress_interaction_matrix solver options need the `htool_`
  prefix. Set `-nsolve 0` unless the random-RHS timing is wanted.
- The dense route (interact) is needed for anything other than the
  uniform b = 1 strike problem (mixed BCs, Coulomb, dip). There, run
  with `-npsfse` unless stresses are needed; the stress pass costs as
  much as the assembly.
- Check the first lines of the solve in the log: "near-field radius
  ... times mean patch length" and "... nonzeros, ... per row" have
  to be present, otherwise the preconditioner is not active (the code
  now aborts on a missing or non-positive radius).

## 4. Recommendation for production solves of the full UCERF set

Primary: H-matrix solve with compress_interaction_matrix

    compress_interaction_matrix -geom_file geom.in -skip_dense -nrandom 0 -nsolve 0 -test_forward 0 \
      -use_hmatrix 4 -hmmvp_inorm 1 -hmmvp_tol 1e-3 \
      -near_pc_rfac 4 -htool_ksp_type gmres -htool_pc_type asm -htool_sub_pc_type lu \
      -htool_ksp_norm_type unpreconditioned -htool_ksp_rtol 1e-6 -htool_ksp_max_it 2000 \
      -htool_ksp_converged_reason

Expected at 265k on 64 ranks: minutes of assembly, tens of iterations,
seconds of solve, operator error about 1e-4 in slip.

Backup: dense solve with interact, only if the H solve does not
converge or the problem needs general boundary conditions

    interact -npsfse -geom geom.in -near_pc_rfac 4 -ksp_type gmres -pc_type asm -sub_pc_type lu \
      -ksp_norm_type unpreconditioned -ksp_rtol 1e-6 -ksp_max_it 2000 -ksp_converged_reason

Expected at 265k on 64 ranks: about 2 h assembly (plus 2 h stress
pass without -npsfse), 560 GB distributed, a minute of solve.

Independent check (recommended once per geometry): repeat the H solve
with `-use_hmatrix 3 -hacapk_ztol 1e-4` and compare flt.dat; the two
operators are approximated differently, so agreement to about 1e-4 in
the relative L2 slip norm indicates the operator error is under
control. Disagreement points to the tolerance, not the solver.

`ucerf/run_production` implements this sequence (H solve, optional
HACApK cross-check, dense fallback) with the convergence checks.

## 5. Open items

- H matrix inside interact's general one-step solve (mixed BCs,
  Coulomb) is not implemented; the dense route remains for those.
- The post-slip stress evaluation re-evaluates all Green's functions;
  for single-component BCs the solved component could be obtained
  from the assembled operator by one MatMult.
- The near-field pair search in calc_petsc_Isn_matrices is an
  N^2/np distance test (about a minute at 265k on 64 ranks). A kd-tree
  would remove it if it becomes noticeable.
- ASM iteration counts were measured up to 100k (20 iterations); the
  265k count is expected to be similar but is not yet measured.
