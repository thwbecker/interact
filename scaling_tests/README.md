# `scaling_tests`: H-matrix backend benchmarks and scaling

Performance and scaling benchmarks for `rsf_solve` and `compress_interaction_matrix`,
covering the H-matrix backends HTOOL, HACApK, and hmmvp (with the dense operator
as a self-checking reference, and BigWham as an experimental full-space option).
The goal is to choose backend and tolerance settings, to see how cost and
compression scale with problem size and rank count, and to keep the comparisons
at matched accuracy rather than matched tolerance.

All numbers produced here are specific to the geometry, accuracy band, machine,
and build that produced them; a different size, tolerance, host, or library
version may shift them. Where a backend looks weak in a particular regime that is
noted as configuration-specific, not as a property of the library.

## What gets measured, and along which axis

The scripts differ along three axes, which is the quickest way to pick one:

- single BP5 fault versus two interacting faults,
- the full earthquake-cycle solve versus the forward operator (assembly and
  matvec) alone,
- and what is swept: MPI rank count, a per-backend accuracy tolerance, the
  admissibility constant eta, or the cluster construction.

## Scripts in this directory

### `rsf_solve_scaling_test.sh [RES] [stop_yr]`
Parallel-scaling benchmark for the full single-fault BP5 cycle. Fixes the
backend tolerances at the matched ~1e-6 band and sweeps MPI rank count times
backend, running `rsf_solve` with `-log_view` and parsing assembly, matvec, total
time, and step count. Use it to see how the cycle scales with cores at fixed
accuracy. It is also the script that answers whether a given backend's assembly
or matvec actually scale with ranks (hold RES fixed and read the per-rank
columns). Example: `./rsf_solve_scaling_test.sh 1km 1000`. Writes
`rsf_solve_scaling_<RES>.<host>.csv`.

### `rsf_solve_compression_test.sh [RES] [np] [stop_yr]`
Compression-versus-speed sweep for the single-fault cycle, the companion to the
scaling script: it fixes the rank count and sweeps a compression knob.
`SWEEP=tol` (default) varies each backend's accuracy tolerance; `SWEEP=eta` holds
tolerance at the matched band and varies the admissibility constant. `RES` can be
a space-separated list (for example `"2km 1km 0.5km"`) to read the
compression-versus-N behaviour across sizes. Use it to choose a backend and
tolerance and to map memory and compression against accuracy. Example:
`./rsf_solve_compression_test.sh 0.25km 16 50`. Writes
`rsf_solve_compression_<sweep>_<RES>.<host>.csv` and feeds `rsf_solve_compression.md`.

### `hmat_scaling_test.sh`
Forward-operator scaling: benchmarks the `compress_interaction_matrix` backends
against core count, with no ODE solve. Each run rebuilds the MPI-distributed
dense reference, reports the relative error of `b = A x`, the H-matrix assembly
time, and dense and H-matrix matvec timings, so every line is self-checking. Use
it to isolate assembly and matvec scaling from the integrator; this is the script
to probe the HTOOL assembly cost. Parameters (the core list, fault size, the
per-backend options, and which backends to loop) are set as plain assignments at
the top. Writes `hmat_scaling.dat`.

#### Isolating one backend's assembly scaling

To dig into a single backend's assembly cost versus rank count (for example to
check whether HTOOL assembly actually parallelizes), `compress_interaction_matrix`
prints `H matrix assembly took` separately from the dense build and the matvec,
so the assembly time is clean even on the normal path. A minimal sweep, HTOOL
only, at the 0.5 km size:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    # ~16000 patches, 0.5 km cells; -z -20 keeps all patches below the z=0 free
    # surface (a dip-90 fault with half-width 20 centered at z=0 would straddle
    # it, which is invalid for the half-space Okada kernel)
    makefault -n 200 -m 80 -l 50 -w 20 -z -20 -strike 0 -dip 90 > /tmp/g16k.in
    for P in 1 2 4 8 16 24 48; do
      asm=$(mpirun --bind-to core --map-by core -np $P \
            compress_interaction_matrix -geom_file /tmp/g16k.in \
            -use_hmatrix 1 -mat_htool_epsilon 3e-5 -mat_htool_eta 10 -nrandom 0 2>&1 \
            | awk '/H matrix assembly took/{print $(NF-1)}')
      echo "np=$P  htool_assembly_s=$asm"
    done

Each run also rebuilds the MPI-distributed dense reference, which costs time but
is parallel and does not contaminate the assembly number; set `-nrandom 0` to
skip the matvec timing, or raise it to keep the `b = A x` sanity check. If the
assembly time stays roughly flat across ranks the assembly is effectively serial
at scale despite the even row split; if it falls about as 1/P it is parallel and
simply expensive at this N. Rerunning with `-use_hmatrix 3` (HACApK) and
`-use_hmatrix 4` (hmmvp) gives controls, and adding `-log_view` to a large-rank
run lets the per-event max/min time ratio reveal load imbalance (one rank holding
most blocks) rather than a uniformly slow assembly. The eta value may matter at
large N, so it is worth sweeping `-mat_htool_eta` (for example 100, 10, 2) as a
second axis.

### `bigwham_omp_vs_mpi_scaling.sh [n] [m] [procs] [nrandom]`
Full-space assembly and matvec scaling for BigWham, comparing OpenMP and MPI
parallelism. Experimental and currently set aside; see `bigwham_status_note.md`
for the standing assessment (strong full-space assembler, weaker matvec as
shipped, licensing to settle before any default-on use). Use only when revisiting
BigWham. Example: `./bigwham_omp_vs_mpi_scaling.sh 60 40 "1 2 4 8 16 32 48" 100`.

### `full_vs_half_test.sh`
Correctness check rather than a performance test: compares the interact Okada
operator in half-space (the standard default) against the native full-space
operator built with `-full_space 1`, reporting the dense forward norms and the
backend-versus-native error. Use it to validate the full-space path, not to time
backends.

## Companion notes in this directory

- `rsf_solve_compression.md`: the single-fault compression-versus-speed results
  (2 km through 0.25 km), the findings and backend-choice guidance, and the
  investigation of the HTOOL assembly cost at large N.
- `compress_interaction_matrix.md`: the forward-operator accuracy and compression
  sweep (around N = 14400) that fixes the matched ~1e-6 accuracy band the cycle
  scripts rely on.
- `bigwham_status_note.md`: BigWham status and the reasons it is paused.

The `.png` figures and `.csv` / `.dat` files are saved outputs of the scripts
above; the CSVs are tagged by host because timing is machine-specific while
compression and accuracy are not.

## Related tests in other directories

### `hmatrix_test/` (two interacting faults, and clustering)
- `sweep_rsf_hmatrix.sh [ds_km] [stop_yr] [procs]`: two-fault accuracy-versus-speed
  cycle sweep. Builds a coupled two-fault BP5-style problem (`make_two_fault.py`),
  sweeps each backend's tolerance, and scores accuracy against the dense solve
  (`rsf_accuracy.py`) alongside speed. Use it to check whether compression
  corrupts a coupled trajectory, not just how much memory it saves.
- `cluster_compression_test.sh [SCALE] [np]` with `cluster_compare.py`: tests
  whether splitting patches by fault before clustering compresses better than a
  single geometric tree, and runs a quick backend compressibility survey on a
  makefault panel. Uses the `-dump_matrix` / `-dump_coords` flags of
  `compress_interaction_matrix`. See `hmatrix_test/cluster_compression.md`.

### `hbi_tests/` (head-to-head with HBI)
- `hbi_bp5_scaling_test.sh [RES] [tmax_yr] [eps_h_label]`: MPI-scaling benchmark
  for HBI (sozawa94/hbi, HACApK-only) on the same BP5 problem, so interact and HBI
  can be compared rank for rank. See `hbi_tests/rsf_solve_vs_hbi_scaling.md`.

## Which to run for what

- Parallel scaling of the cycle at fixed accuracy: `rsf_solve_scaling_test.sh`.
- Memory and compression versus accuracy or versus N, single fault:
  `rsf_solve_compression_test.sh`.
- Assembly and matvec scaling of a backend, isolated from the ODE solve:
  `hmat_scaling_test.sh`.
- Does compression stay accurate when faults couple: `sweep_rsf_hmatrix.sh`.
- Fault-aware versus joint clustering, or a quick compressibility check:
  `cluster_compression_test.sh`.
- Half-space versus full-space physics: `full_vs_half_test.sh`.
- Comparison against HBI: `hbi_bp5_scaling_test.sh`.

A typical workflow is to pick a backend and tolerance with
`rsf_solve_compression_test.sh`, confirm it parallelizes with
`rsf_solve_scaling_test.sh` (or `hmat_scaling_test.sh` for the operator alone),
and validate accuracy once on a coupled problem with `sweep_rsf_hmatrix.sh`.

## Conventions

- Thread hygiene applies to all timing runs: `OMP_NUM_THREADS=1` and the BLAS
  thread pins (`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`) keep MPI ranks from
  oversubscribing and stop BLAS from nesting threads inside block-level
  parallelism; hmmvp raises its own pool internally from `-hmmvp_nthreads`.
- Resolution tags (`2km`, `1km`, `0.5km`, `0.25km`) set the BP5 cell size and
  therefore N (roughly 1000, 4000, 16000, 64000 patches).
- Cross-backend statements are made at matched accuracy, not matched tolerance;
  the matched ~1e-6 band is documented in the two `.md` notes above.
- Parameters are positional or set as plain assignments at the top of each
  script; the runs do not read tunable parameters from the environment, so they
  are deterministic.
