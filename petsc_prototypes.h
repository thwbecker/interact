#include <petscksp.h>
#include "petscts.h"

PetscErrorCode rsf_ODE_RHSFunction(TS, PetscReal, Vec, Vec, void*);
PetscErrorCode rsf_TS_Monitor(TS, PetscInt, PetscReal, Vec, void*);
PetscErrorCode rsf_domain_check(TS, PetscReal, Vec, PetscBool*);

/* kernel function */


PetscErrorCode GenKEntries_htools(PetscInt , PetscInt , PetscInt ,const PetscInt *,
				  const PetscInt *, PetscScalar *, void *);
PetscScalar GenKEntries_h2opus(PetscInt, PetscReal [], PetscReal [], void *);

PetscReal vel_from_rsf(PetscReal, PetscReal, PetscReal, PetscReal,PetscReal,
		       PetscReal *, PetscReal *, PetscReal *,struct med *);


PetscErrorCode calc_petsc_Isn_matrices(struct med *, struct flt *,PetscInt ,PetscReal, int, Mat *,hacapk_shell_ctx *);

#ifdef USE_PETSC_HMAT
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



PetscErrorCode set_htools_defaults_and_options(struct med *);
PetscErrorCode set_h2opus_defaults_and_options(struct med *);

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
extern void cmake_hacapk_struct_hmat(void *, double);
extern void chacapk_mult_Ax_H(void *, double *, double *);
PetscErrorCode MatMult_HACApK(Mat , Vec , Vec );
PetscErrorCode set_hacapk_defaults_and_options(struct med *);
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
PetscErrorCode set_hmmvp_defaults_and_options(struct med *);

PetscErrorCode MatMult_hmmvp(Mat , Vec , Vec );
#endif


