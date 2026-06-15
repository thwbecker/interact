# HACApK (v1.0.0) — H-matrix backend for `interact`

This subdirectory is **not** original HACApK development. It contains a **copy of
HACApK v1.0.0** together with a thin **C interface** (added for `interact`, T. W. Becker)
that lets `interact` call HACApK as one of its hierarchical-matrix (H-matrix) backends for
the boundary-element fault-interaction (Okada) kernels. All credit for HACApK itself
belongs to its original authors; only the C wrapper here is local work.

## Provenance — what this is and where it came from

HACApK is the hierarchical-matrix library developed within the **ppOpen-HPC** project
(Ida et al.). The sources here follow the **v1.0.0** release. HACApK is also the H-matrix
engine used by **HBI** (So Ozawa and co-workers, <https://github.com/sozawa94/hbi>); the
copy bundled here is consistent with the HACApK that HBI ships — in particular
`HACApK_lib.f90`
(<https://github.com/sozawa94/hbi/blob/master/HACApK_lib.f90>) — which `interact` also
uses as a reference point for the rate-and-state solver comparisons and for the
thread-safe Okada validation.

This is a *copy*, pinned to v1.0.0; a newer or differently configured HACApK release may
behave differently from what is included here. The original attribution and citations are
preserved verbatim in the `COPYRIGHT` file alongside this README.

## Contents

- `v.1.0.0/C_interface/` — the **C interface** (the part added for `interact`):
  `HACApK_c_interface.f90` and `HACApK_c_interface.h` (the Fortran↔C bridge),
  `test_hacpak_c.c` (a C driver/test), the HACApK Fortran modules
  (`m_HACApK_base.f90`, `m_HACApK_solve.f90`, `m_HACApK_calc_entry_ij.f90`,
  `HACApK_lib.f90`, …), and a `Makefile`.
- `v.1.0.0/fortran_interface/` — the original HACApK Fortran interface/test, kept for
  reference.
- top level — the HACApK v1.0.0 Fortran sources, the upstream manual
  (`manual_Hacapk_1.0.0.pdf`), and the original `COPYRIGHT`.

The C interface exposes routines such as `cinit_hacapk_struct`,
`cset_hacapk_struct_coord`, … to build an N×N H-matrix from a user-supplied kernel
callback (`ckernel_func`) and apply it; see `HACApK_c_interface.h` for the full set.

## Building (within `interact`)

The C-interface library is built from `v.1.0.0/C_interface/`, e.g.

```
cd v.1.0.0/C_interface
make F90=mpifort F77=mpifort CC=mpicc LD=mpifort libhacapk.a
```

`interact` then links `libhacapk.a` and selects this backend at run time with
`-use_hmatrix 3` (HACApK), with the ACA tolerance set via `-hacapk_ztol`.

## Attribution and references

Please cite the original HACApK / ppOpen-HPC work when using this backend. References
(as recorded in `COPYRIGHT`):

- Standard H-matrices:
  <https://www.jstage.jst.go.jp/article/ipsjjip/22/4/22_642/_article/-char/en>
- Lattice H-matrices:
  <https://ieeexplore.ieee.org/abstract/document/8425193>
- So Ozawa, Akihiro Ida, Tetsuya Hoshino, Ryosuke Ando (2023), "Large-scale earthquake
  sequence simulations of 3D geometrically complex faults using the boundary element
  method accelerated by lattice H-matrices," *Geophysical Journal International*, 232 (3),
  1471–1481, <https://doi.org/10.1093/gji/ggac386>.
- So Ozawa, Yuyun Yang, Eric M. Dunham (2024), "Fault Valve Instability: A mechanism for
  slow slip events," *Journal of Geophysical Research: Solid Earth*, 129,
  <https://doi.org/10.1029/2024JB029165>.

HBI (So Ozawa et al., <https://github.com/sozawa94/hbi>) is the boundary-integral
earthquake-cycle code used here as a reference point for the rate-and-state solver
comparisons and the thread-safe Okada validation.

## License / terms

HACApK is third-party software included here as a copy; its use is governed by the terms
of the upstream ppOpen-HPC HACApK project — see `manual_Hacapk_1.0.0.pdf` and the original
`COPYRIGHT` file in this directory. The C interface added for `interact` is part of
`interact` and is distributed under `interact`'s terms. If you redistribute or build on
this, retain the upstream HACApK attribution.
