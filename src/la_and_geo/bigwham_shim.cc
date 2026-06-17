/*

  bigwham_shim.cc

  extern "C" wrapper around BigWham's BigWhamIO class so the C / Fortran +
  PETSc interact code can build and apply a BigWham hierarchical matrix the
  same way it uses the hmmvp and HACApK interfaces (a MATSHELL whose MatMult
  calls cbigwham_mvp, with PETSc KSP for the solve).

  IMPORTANT difference from the hmmvp shim: the hmmvp shim compresses
  interact's OWN Okada kernel via the ckernel_func callback. BigWham instead
  supplies its own elastic kernel (for example "3DR0-H": a 3D rectangular
  displacement-discontinuity element, collocation at the centre, 3 slip DoF
  per element), so this shim does NOT take a kernel callback. It is handed a
  mesh (node coordinates + element connectivity) plus elastic constants,
  builds the hierarchical matrix, and exposes matvec / diagonal / info /
  destroy.

  Scope and conventions (read before wiring use_hmatrix==5):
    - BigWham is FULL-SPACE only and OpenMP only (no MPI). It must therefore
      be compared against the native full-space Okada operator
      (compress / rsf_solve run with -full_space 1), never the half-space one.
    - The operator built here is BigWham's native 3N x 3N traction-from-slip
      operator, with the 3 DoF per element being slip components in BigWham's
      per-element local frame (N = number of elements). Mapping interact's
      rectangular patches to (coor, conn), and projecting interact's single
      slip-mode / stress-mode N-vectors onto the 3-component BigWham 3N
      vectors (and back), is the job of the caller (the use_hmatrix==5
      wiring), not of this shim.
    - cbigwham_mvp uses BigWhamIO::MatVec, which takes and returns vectors in
      natural (user) DoF ordering and permutes internally. For a hot matvec
      loop one can later switch to the permuted / in-place variants.

  build: see the BIGWHAM_* variables in config/makefile.petsc. Needs a C++17
  compiler, OpenMP, BigWham/src and BigWham/il on the include path with the IL
  math-vendor defines, and links against BigWham/build/libBigWham.a. Build the
  static library first with BigWham/compile_for_interact.sh.

  This wrapper makes no claim about BigWham's accuracy or performance; those
  are properties of the library and the chosen accuracy parameters
  (max_leaf, eta, eps_aca) and may differ across versions and configurations.

*/

#include <vector>
#include <string>
#include <cstring>
#include <cstddef>

#include "io/bigwham_io.h"   /* from BigWham/src (set by BIGWHAM_INC) */

extern "C" {
  void *cbigwham_create(const double *coor, int n_coor,
                        const int *conn, int n_conn,
                        const char *kernel, double E, double nu,
                        int n_threads);
  void  cbigwham_build(void *h, int max_leaf, double eta, double eps_aca);
  void  cbigwham_mvp(void *h, const double *x, double *y);
  void  cbigwham_get_diagonal(void *h, double *diag);
  void  cbigwham_get_info(void *h, int *m, int *n, double *compression);
  void  cbigwham_delete(void *h);
}

/*
  create a BigWham operator from a mesh.

    coor   : flat array of node coordinates, length n_coor (= 3 * n_nodes)
    conn   : flat array of element connectivity, length n_conn
             (= vertices_per_element * n_elements; 4 for "3DR0")
    kernel : BigWham kernel string, e.g. "3DR0-H"
    E, nu  : Young's modulus and Poisson ratio (note interact carries shear
             modulus G; convert with E = 2 G (1 + nu) in the caller)
    n_threads : OpenMP threads for BigWham

  returns an opaque handle (BigWhamIO *); the hierarchical matrix is not built
  until cbigwham_build is called.
*/
void *cbigwham_create(const double *coor, int n_coor,
                      const int *conn, int n_conn,
                      const char *kernel, double E, double nu,
                      int n_threads)
{
  std::vector<double> vcoor(coor, coor + n_coor);
  std::vector<int>    vconn(conn, conn + n_conn);
  std::vector<double> props = {E, nu};
  BigWhamIO *io = new BigWhamIO(vcoor, vconn, std::string(kernel), props, n_threads);
  return (void *)io;
}

/* build the hierarchical matrix (block cluster tree + ACA) */
void cbigwham_build(void *h, int max_leaf, double eta, double eps_aca)
{
  BigWhamIO *io = (BigWhamIO *)h;
  io->BuildHierarchicalMatrix(max_leaf, eta, eps_aca);
}

/* y = A x, both length n = MatrixSize(1) = 3 * n_elements, natural ordering */
void cbigwham_mvp(void *h, const double *x, double *y)
{
  BigWhamIO *io = (BigWhamIO *)h;
  int n = io->MatrixSize(1);
  std::vector<double> vx(x, x + n);
  std::vector<double> vy = io->MatVec(vx);
  std::memcpy(y, vy.data(), sizeof(double) * static_cast<std::size_t>(n));
}

/* diagonal of A, length n = MatrixSize(1) (useful for a Jacobi preconditioner) */
void cbigwham_get_diagonal(void *h, double *diag)
{
  BigWhamIO *io = (BigWhamIO *)h;
  std::vector<double> d;
  io->GetDiagonal(d);
  std::memcpy(diag, d.data(), sizeof(double) * d.size());
}

/* m, n (rows, cols) and the compression ratio (dense_storage / hmat_storage) */
void cbigwham_get_info(void *h, int *m, int *n, double *compression)
{
  BigWhamIO *io = (BigWhamIO *)h;
  *m = io->MatrixSize(0);
  *n = io->MatrixSize(1);
  *compression = io->GetCompressionRatio();
}

void cbigwham_delete(void *h)
{
  if (h) delete (BigWhamIO *)h;
}
