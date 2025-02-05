#include <petscksp.h>
#include "petscts.h"

PetscErrorCode rsf_ODE_RHSFunction(TS, PetscReal, Vec, Vec, void*);

/* kernel function */


PetscErrorCode GenKEntries(PetscInt , PetscInt , PetscInt ,const PetscInt *,
			   const PetscInt *, PetscScalar *, void *);

PetscReal vel_from_rsf(PetscReal, PetscReal, PetscReal, PetscReal,PetscReal,
		       PetscReal *, PetscReal *, PetscReal *, PetscReal *);


PetscErrorCode calc_petsc_Isn_matrices(struct med *, struct flt *,PetscBool ,PetscReal, int, Mat *);
