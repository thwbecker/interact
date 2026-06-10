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
