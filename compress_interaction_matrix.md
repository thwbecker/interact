# Testing H-matrix backends with `compress_interaction_matrix`

`compress_interaction_matrix` builds the (strike-stress, by default)
fault interaction operator twice — once dense, once with a hierarchical
matrix backend constructed from the *same* Green's function kernel —
and compares the two for the forward product `b = A x` (always) and
optionally the inverse solve `x = A\b`. It is the fast way to choose
H-matrix settings for a given geometry before committing to long
`rsf_solve` runs.

All results below: two vertical fault groups offset by one x unit
(`makefault -n .. -m ..`, see below), half-space Okada rectangles,
serial, June 2026. Errors are relative two-norm differences against the
dense operator.

## Build requirements

| backend | `-use_hmatrix` | needs |
|---|---|---|
| dense | (always built) | PETSc |
| HTOOL | 1 (default) | PETSc configured `--download-htool`, compile with `-DUSE_PETSC_HMAT` (`PETSC_HTOOL_USED` in `makefile.petsc`) |
| H2OPUS | 2 | PETSc configured `--download-h2opus --download-thrust` (CPU), same define |
| HACApK | 3 | `cd HACApK/v.1.0.0/C_interface; make; ar rcs libhacapk.a m_*.o HACApK_*.o`, then `HACAPK_DEFINES = -DUSE_HACAPK` and `HACAPK_LIBS` in `makefile.petsc`. **Build the library with `-O3`** (`OPT` in its Makefile): the unoptimized build has a ~13x slower matvec. **The library MUST be compiled with the same MPI implementation as PETSc** (otherwise `HACApK_init` fails with `MPI_Comm_size failed`): use the wrappers PETSc reports, `pkg-config --variable=fcompiler PETSc`, or `$PETSC_DIR/$PETSC_ARCH/bin/mpifort` for `--download-mpich` builds. |

PETSc example configure for everything at once:

    ./configure --with-debugging=0 --with-cc=mpicc --with-cxx=mpicxx \
      --with-fc=0 --download-htool --download-h2opus --download-thrust \
      COPTFLAGS=-O2 CXXOPTFLAGS=-O2

## Test geometry

    makefault -n 20 -m 10           >  geom.in   # fault group 0
    makefault -n 20 -m 10 -x 1 -grp 1 >> geom.in # group 1, offset
    # -n 20 -m 10  -> N=400; 40/20 -> 1600; 80/40 -> 6400; 100/50 -> 10000

## Command lines

Forward accuracy (error is always printed; `-mat_view ::ascii_info`
adds the htool compression report):

    compress_interaction_matrix -geom_file geom.in -use_hmatrix 1 \
        -mat_htool_epsilon 3e-5 -mat_htool_eta 10        # HTOOL
    compress_interaction_matrix -geom_file geom.in -use_hmatrix 2 \
        -mat_h2opus_maxrank 256                          # H2OPUS
    compress_interaction_matrix -geom_file geom.in -use_hmatrix 3 \
        -hacapk_ztol 1e-5                                # HACApK

Matvec timing (adds n random-vector multiplies, timed for dense and H):

    compress_interaction_matrix -geom_file geom.in -use_hmatrix 3 \
        -hacapk_ztol 1e-5 -nrandom 300

Inverse solve test (dense LU vs unpreconditioned KSP on the H operator):

    compress_interaction_matrix -geom_file geom.in -use_hmatrix 1 \
        -mat_htool_epsilon 3e-5 -mat_htool_eta 10 -test_forward false

Assembly times for the dense and the H matrix are printed on stderr
("dense assembly took", "H matrix assembly took").

## Accuracy: error vs tolerance setting

Forward `b = Ax`, error tracks the tolerance for HTOOL and HACApK:

| setting | N=400 | N=1600 | N=6400 | N=10000 |
|---|---|---|---|---|
| HTOOL eps 1e-3 (ACA) | 5.8e-4* | 5.8e-4 | — | — |
| HTOOL eps 1e-4 (ACA, eta 100) | — | 7.0e-5 | 2.1e-4 | — |
| HTOOL eps 3e-5 (ACA, eta 10) | — | 3.2e-7 | 5.8e-7 | 8.8e-6 |
| HTOOL eps 1e-8 | 1e-15 (= dense) | — | — | — |
| H2OPUS (maxrank 256) | 2.3e-3 | 2.3e-3 | 2.6e-3 | — |
| HACApK ztol 1e-4 | 2.0e-6 | 4.7e-6 | — | 1.8e-5 |
| HACApK ztol 1e-5 | — | 2.5e-7 | 3.5e-7 | 6.2e-7 |

(*N=400 value at eta=100; the eps 1e-3 row is representative.)

Inverse `x = A\b` (N=1600, vs dense LU): HTOOL 5.1e-6 (recommended
config), HACApK works through the same KSP/MATSHELL path, H2OPUS
8.1e-4 (floor, see below).

## Performance (300 matvecs; assembly of one H matrix)

| N=10000 | err | H assembly | matvec total | vs dense matvec |
|---|---|---|---|---|
| dense | — | ~85 s | 25.6 s | 1x |
| HTOOL (ACA eps 3e-5, eta 10) | 8.8e-6 | ~19 s | 6.6 s | 3.9x |
| **HACApK ztol 1e-5** | **6.2e-7** | **11.7 s** | **3.8 s** | **6.7x** |
| HACApK ztol 1e-4 | 1.8e-5 | 9.0 s | 2.9 s | 8.8x |

Speedups grow with N (HTOOL matvec: 1.7x at N=1600, 3.5-8x at 6400);
dense matvec scales N^2, the H backends ~N log N.

## Recommended settings

- **HACApK (`-use_hmatrix 3`) is the best fit for this operator**:
  index-based kernel interface, robust default ACA (no reliability
  tuning needed; `param(61)=1` normalization is set by the interface),
  fastest assembly and matvec at matched error. Use `-hacapk_ztol`
  ~3-10x below the error target (1e-5 gives <1e-6 everywhere tested).
- **HTOOL (`-use_hmatrix 1`)**: use `sympartialACA` (default) with
  **eta = 10 (default)** and eps ~3x below the target. Caution: large
  eta (e.g. 100) compresses better but partial-pivoting ACA can
  mis-converge on big marginally-separated blocks (the kernel's
  quadrant sign structure produces near-zero pivot rows): at N=10000,
  eps 3e-5 gave 3.1e-4 with eta=100 vs 8.8e-6 with eta=10. The
  reliable eta=10 costs more assembly (70 s vs 10 s at N=6400). The
  fullACA/SVD compressors are reliable and compress best but have
  O(N^2) assembly — only for operators reused over >~1e5 matvecs.
- **H2OPUS (`-use_hmatrix 2`)**: construction is via sampling of the
  assembled dense operator (`MatCreateH2OpusFromMat`), because the
  kernel interface (Chebyshev interpolation at arbitrary points) is
  incompatible with patch-pair Green's functions. **This h2opus only
  implements sampling construction for SYMMETRIC matrices** (a warning
  is printed): the result approximates (K+K^T)/2, and the measured
  operator asymmetry |K-K^T|/|K| ~ 4-8e-4 is an irreducible error
  floor (~2.5e-3 in the matvec). Needs `-mat_h2opus_maxrank 256` for
  N >= 1600 (the default 64 truncates and *wastes* memory). Fine for
  symmetric problems; not for this operator.
- Context for the error target: in quasi-dynamic `rsf_solve` cycles,
  operator errors ~1e-4 shift event onsets by O(100 s) and can
  restructure two-fault rupture sequences; ~1e-6 reproduces dense
  onsets to ~5e-4 s. The H2OPUS floor is therefore not
  production-viable for earthquake cycles, while HTOOL at eps 1e-6
  and HACApK at ztol 1e-5..1e-6 are.
