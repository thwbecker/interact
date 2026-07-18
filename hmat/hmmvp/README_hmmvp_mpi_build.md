# Building hmmvp (github.com/ambrad/hmmvp) in MPI mode on a modern MPI

hmmvp predates two MPI changes that affect recent OpenMPI (>= 5) / MPICH (>= 4).
Two small, local edits let `mode=mpi` compile and link; neither is needed on
older MPI stacks, and both are MPI-version compatibility fixes rather than
problems with hmmvp itself.

## 1. `src/Mpi.cpp` — removed MPI-1 error-handler calls

`MPI_Errhandler_create` / `MPI_Errhandler_set` were deprecated in MPI-2 and
removed in MPI-3. Replace with the drop-in equivalents (identical arguments):

```
-  MPI_Errhandler_create(MPI_Handler_Crash, &eh);
-  MPI_Errhandler_set(MPI_COMM_WORLD, eh);
+  MPI_Comm_create_errhandler(MPI_Handler_Crash, &eh);
+  MPI_Comm_set_errhandler(MPI_COMM_WORLD, eh);
```

(in `util::mpi::Init`, inside the `#ifdef UTIL_MPI` block).

## 2. `make.inc` — suppress the removed MPI C++ bindings

OpenMPI >= 5 / MPICH >= 4 dropped the deprecated `MPI::` C++ bindings, but
`mpi.h` still pulls in `mpicxx.h` unless told not to. hmmvp uses only the C MPI
API, so define the standard skip guards in the C++ flags:

```
-CPPFLAGS = -std=c++14
+CPPFLAGS = -std=c++14 -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX
```

## 3. Build the MPI library

```
cd hmmvp
# make.inc should define MPICPP = mpic++   (plus fortranint/BLAS/LAPACK as usual)
make clean
make mode=mpi libhmmvp      # -> lib/libhmmvp_mpi.a
```

Only the library is needed by interact; the `hmmvpbuild`/`mvp` driver targets
pull in extra flags and are not required here.
