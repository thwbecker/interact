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


PetscErrorCode calc_petsc_Isn_matrices(struct med *, struct flt *,PetscInt ,PetscReal, int, Mat *);

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
#endif
