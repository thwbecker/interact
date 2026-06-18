# `scaling_tests`: single-fault H-matrix backend benchmarks

Performance benchmarks for `rsf_solve`'s H-matrix backends (HTOOL, HACApK,
hmmvp, and dense as reference) on the single SEAS BP5 fault. The goal is to pick
backend and tolerance settings and to understand how cost and compression scale
with MPI ranks and with problem size. The two-fault accuracy sweep (which turned
out to be confounded by chaotic recurrence) lives separately in `hmatrix_test/`,
and is not part of this directory.

See also `../rsf_solve.md` for the backend performance summary and recommended
defaults, and `../bp5/README.md` for the benchmark problem setup and validation.

## Scripts

| script | what it varies | what it measures |
|--------|----------------|------------------|
| `rsf_solve_scaling_test.sh` | MPI rank count, at a fixed matched ~1e-6 tolerance band | assembly, matvec, total, steps, memory vs ranks (the parallel-scaling view) |
| `rsf_solve_compression_test.sh` | each backend's tolerance (`SWEEP=tol`) or admissibility eta (`SWEEP=eta`), at a fixed rank count, optionally across a list of resolutions | compression ratio, memory, assembly, matvec, total, steps (the compression-vs-speed and compression-vs-N view) |
| `hmat_scaling_test.sh` | MPI rank count, at the matched band | the same scaling done on the forward operator via `compress_interaction_matrix` (block apply), with forward error |
| `bigwham_omp_vs_mpi_scaling.sh` | thread/rank count for BigWham vs the MPI backends | full-space matvec scaling; documents BigWham's OpenMP matvec anti-scaling |
| `full_vs_half_test.sh` | full-space vs half-space kernel | consistency check between kernels |

All four `rsf_solve` backends run under MPI on `procs` ranks (hmmvp included; it
was moved from OpenMP to MPI). Tolerances are each package's own knob and are not
comparable across packages at equal nominal value; the matched ~1e-6 band is
HTOOL `epsilon` 3e-5, HACApK `ztol` 1e-1, hmmvp `tol` 1e-7.

## Notes

| file | content |
|------|---------|
| `rsf_solve_compression.md` | compression-vs-speed and compression-vs-N findings, with the matched-band reading and the relation to `hmat_scaling_test.sh` |
| `compress_interaction_matrix.md` | the tolerance-to-forward-error mapping (at N=14400) that the matched band rests on |
| `bigwham_status_note.md` | measured status note on the BigWham full-space backend (correct, assembles well, matvec does not scale as shipped; paused) |

## Results and figures

- `compression_vs_N.png`, `compression_speed_tradeoff.png`: from the compression
  sweep (companions to `rsf_solve_compression.md`).
- `compress_interaction_matrix_scaling.png`, `bp5_0p5km_scaling_clean_hosts.png`:
  forward-operator and multi-host scaling figures.
- `rsf_solve_compression_tol_*.csv`: tolerance-sweep results per resolution.
- `rsf_solve_scaling_*.csv`, `hmat_scaling.*`: rank-scaling results per host.

The `.walter` / `.theo2` suffixes denote the host the numbers were measured on;
timing is host and build dependent, so read cross-host numbers as indicative.
