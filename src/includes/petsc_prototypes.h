#include <petscksp.h>
#include "petscts.h"

/* mode descriptions */
#define IHMAT_TYPE_DENSE 0
#define IHMAT_TYPE_HTOOLS 1
#define IHMAT_TYPE_H2OPUS 2
#define IHMAT_TYPE_HACAPK 3
#define IHMAT_TYPE_HMMVP 4
#define IHMAT_TYPE_BIGWHAM 5

PetscErrorCode rsf_ODE_RHSFunction(TS, PetscReal, Vec, Vec, void*);
PetscErrorCode rsf_TS_Monitor(TS, PetscInt, PetscReal, Vec, void*);
PetscErrorCode rsf_domain_check(TS, PetscReal, Vec, PetscBool*);

PetscErrorCode print_petsc_matrix(Mat , PetscInt , PetscInt ,char *);
/* kernel function */


PetscErrorCode GenKEntries_petsc(PetscInt , PetscInt , PetscInt ,const PetscInt *,
				 const PetscInt *, PetscScalar *, void *);

void report_hmat_storage(struct med *, const char *,PetscInt , PetscInt , long );
PetscReal vel_from_rsf(PetscReal, PetscReal, PetscReal, PetscReal,PetscReal,
		       PetscReal *, PetscReal *, PetscReal *,struct med *);
const char *hmat_backend_name(int );

PetscErrorCode interact_petsc_initialize(int *, char ***);
PetscErrorCode calc_petsc_Isn_matrices(struct med *, struct flt *,PetscInt ,PetscReal, int, Mat *,hacapk_shell_ctx *);
PetscErrorCode set_hmat_defaults_and_options(struct med *,int);


#ifdef USE_PETSC_HMAT
PetscScalar GenKEntries_h2opus(PetscInt, PetscReal [], PetscReal [], void *);

#if PETSC_VERSION_LT(3,22,0)
/* 
   compatibility with older PETSc: the Mat*Kernel callback types were
   function-pointer typedefs before the 3.22 rename to function-type
   *Fn typedefs; define the function-type names so declarations like
   "MatHtoolKernelFn *kernel" work for both
*/
typedef PetscErrorCode MatHtoolKernelFn(PetscInt, PetscInt, PetscInt, const PetscInt *, const PetscInt *, PetscScalar *, void *);
typedef PetscScalar MatH2OpusKernelFn(PetscInt, PetscReal[], PetscReal[], void *);
#endif


void setup_kdtree(int ,int ,PetscReal *,struct flt *,struct med *);

#endif


#ifdef USE_HACAPK
/* 
   HACApK support via the C interface in HACApK/v.1.0.0/C_interface:
   index-based kernel (a natural fit for patch-pair Green's functions,
   unlike interpolation-based approaches), wrapped as a PETSc MATSHELL
   so all comparison and timing machinery works unchanged.
   build the library there first (make in that directory, or ar the
   objects into libhacapk.a) and set HACAPK_DEFINES/HACAPK_LIBS, see
   makefile.petsc
*/
extern void *cinit_hacapk_struct(int, void *);
extern void cdeallocate_hacapk_struct(void *);
extern void cset_hacapk_struct_coord(void *, double *, double *, double *);
extern void cset_hacapk_eta(void *, double);
extern void cset_hacapk_inorm(void *, int);
extern void cmake_hacapk_struct_hmat(void *, double);
/* return this rank's number of stored scalars in the assembled HACApk
   H-matrix (dense leaves ndl*ndt, low-rank leaves kt*(ndl+ndt)); the caller
   sums across ranks for the global count. */
extern long cget_hacapk_nnz(void *);
extern void chacapk_mult_Ax_H(void *, double *, double *);
PetscErrorCode MatMult_HACApK(Mat , Vec , Vec );
#endif

#ifdef USE_HMMVP
/* 
   hmmvp (A.M. Bradley) support via hmmvp_c_shim.cpp: index-based
   block kernel (same ckernel_func as HACApK), in-memory compression
   with a WHOLE-MATRIX relative Frobenius tolerance -hmmvp_tol
   (||B-A||_F <= tol ||B||_F, i.e. tol bounds exactly the error this
   tool measures), OpenMP-threaded construction and matvec. wrapped
   as a MATSHELL like HACApK; the Mvp needs global vectors, so the
   same scatter-to-all context is used (at np>1 every rank holds the
   full H matrix and computes the full product - correct but
   redundant; intended for serial/threaded use)
*/
extern void *chmmvp_compress_in_memory(int, double *, double *, double *,
				       double, double, int, void *);
extern void chmmvp_mvp(void *, double *, double *);
extern void chmmvp_get_info(void *, int *, int *, long *);
extern void chmmvp_delete(void *);
/* MPI interface (distributed assembly to file + MpiHmat matvec). These
   MUST be declared here: without the prototypes C assumes chmmvp_mpi_load
   returns int and truncates the 64-bit MpiHmat* handle to 32 bits, which
   then segfaults on first use. */
extern int   chmmvp_compress_to_file(int, double *, double *, double *,
				     double, double, void *, const char *);
extern void *chmmvp_mpi_load(const char *, int);
extern void  chmmvp_mpi_mvp(void *, double *, double *);
extern void  chmmvp_mpi_get_info(void *, int *, int *, long *);
extern void  chmmvp_mpi_delete(void *);


PetscErrorCode MatMult_hmmvp(Mat , Vec , Vec );
#endif

#ifdef USE_BIGWHAM
/* BigWham full-space H-matrix shim (src/la_and_geo/bigwham_shim.cc). BigWham
   owns its kernel, so create takes a mesh (coor,conn) + elastic constants, not
   a kernel callback. matvec/diagonal are length 3*n_elements (natural order). */
extern void *cbigwham_create(const double *, int, const int *, int,
			     const char *, double, double, int);
extern void  cbigwham_build(void *, int, double, double);
extern void  cbigwham_mvp(void *, const double *, double *);
extern void  cbigwham_get_diagonal(void *, double *);
extern void  cbigwham_get_info(void *, int *, int *, double *);
extern void  cbigwham_delete(void *);
/* interact-side BigWham driver (petsc_interact.c) */
PetscErrorCode MatMult_bigwham(Mat, Vec, Vec);
PetscErrorCode setup_bigwham_matshell(struct med *, struct flt *, PetscReal, int,
				      Mat *, hacapk_shell_ctx **);
#endif


