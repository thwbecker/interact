# MPI interface for hmmvp in interact (compress_interaction_matrix)

This documents switching interact's hmmvp backend (`-use_hmatrix 4`) from the
previous OpenMP/in-memory interface to hmmvp's MPI interface, for distributed
assembly and matvec, focused on `compress_interaction_matrix`. Scope note: the
correctness checks below are from a single-core container (MPI oversubscribed),
so they validate the interface, not scaling; the scaling numbers need a real
multi-core run on your side.

## What hmmvp's MPI mode does (and how it differs from OpenMP)

The two hmmvp paths are structurally different:

- **OpenMP / in-memory (previous):** `Compressor::CompressInMemory` builds the
  whole H-matrix in memory on a single rank; `Hmat::Mvp` applies it on full
  vectors. Under `mpirun -np P` every rank redundantly built the whole matrix.
- **MPI (new):** `Compressor::CompressToFile` assembles the H-matrix with the
  blocks distributed across ranks and writes one file (collective);
  `hmmvp::MpiHmat` loads that file (each rank holds its share). In
  `MpiHmat::Mvp` the full `x` need be valid only on the root and the full `y`
  is valid only on the root afterwards, with the block products distributed and
  reduced. So the matvec is a root-gather / distributed-compute / root-scatter
  pattern (similar in spirit to HACApK's gather-and-reduce).

## interact-side changes

- `hmmvp_c_shim.cpp`: added `chmmvp_compress_to_file` (distributed assembly to a
  file) and, under `#ifdef UTIL_MPI`, `chmmvp_mpi_load` / `chmmvp_mpi_mvp` /
  `chmmvp_mpi_get_info` / `chmmvp_mpi_delete` wrapping `hmmvp::MpiHmat<double>`.
  The in-memory functions are retained for the non-MPI build. The file also
  defines `OMPI_SKIP_MPICXX`/`MPICH_SKIP_MPICXX` before any MPI include (see
  build notes).
- `petsc_interact.c`: `MatMult_hmmvp` now has an MPI path (gather `x` to the root
  with a `VecScatterToZero`, call the collective `MpiHmat` matvec into a
  full-length buffer, scatter `y` back) under `USE_HMMVP_MPI`, with the old
  in-memory path kept under `#else`. `calc_petsc_Isn_matrices` case 4 compresses
  to a temporary file, loads the `MpiHmat`, and removes the file after load.
- `compress_interaction_matrix.c`: the inline (default, non-`-make_matrix_externally`)
  case 4 was given the same MPI path; this is the code the test program actually
  runs by default.
- `petsc_prototypes.h`: prototypes for the new shim functions.
- `makefile.petsc`: `HMMVP_LIBS` links `-lhmmvp_mpi`; `HMMVP_DEFINES` adds
  `-DUSE_HMMVP_MPI`. `makefile`: the shim is compiled with `-DUTIL_MPI`.
- `hmat_scaling_test.sh`: the hmmvp loop now runs under `mpirun -np P` like HTOOL
  and HACApK (previously OpenMP via `-hmmvp_nthreads`).

## Building hmmvp in MPI mode

See `hmmvp_mpi_patches/` for two small compatibility edits to hmmvp needed on
modern OpenMPI/MPICH (removed MPI-1 error-handler calls; suppress the removed
MPI C++ bindings). Then `make clean && make mode=mpi libhmmvp` produces
`lib/libhmmvp_mpi.a`. These edits are MPI-version compatibility fixes; on an
older MPI stack hmmvp may build in MPI mode unchanged.

## Correctness (single-core container, MPI oversubscribed)

`compress_interaction_matrix` on a 1600-patch makefault geometry, comparing the
H-matrix product against the dense reference, `-hmmvp_tol 1e-6`:

| np | relative error |b - b_h|/|b| |
|---:|---|
| 1 | 1.92e-6 |
| 2 | 1.91e-6 |
| 4 | 1.90e-6 |

The error tracks the requested tolerance and is consistent across rank counts
(the tiny variation is reduction-order in the distributed matvec). Compression
ratio and stored-scalar count are identical across np, as expected for the same
H-matrix assembled in a distributed fashion. A 3600-patch case at `-hmmvp_tol
1e-5` gave ~9e-5, again consistent with tolerance.

## What is NOT yet measured

No scaling numbers: this container has one physical core, so MPI runs are
oversubscribed (correctness only). The `--bind-to core --map-by core` options in
`hmat_scaling_test.sh` will (correctly) fail to oversubscribe here but are right
for a real node. On your 24-core box, `hmat_scaling_test.sh` should now produce
hmmvp lines under MPI alongside HTOOL/HACApK. Two things worth watching, stated
as expectations rather than results:

- the root-gather/scatter in the matvec funnels the vectors through rank 0 each
  product, which may limit matvec scaling at high rank counts (as it can for
  HACApK), independent of how well the distributed block products scale;
- `CompressToFile` does real disk I/O for the H-matrix; on a fast problem the
  file write/read may be a visible part of assembly, and on a shared filesystem
  it is worth placing the temp file on node-local storage.

These are observations to verify with the benchmark, not conclusions.
