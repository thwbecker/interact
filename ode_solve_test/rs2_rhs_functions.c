#include "interact.h"
/*

  solving ordinary differential equations originally based on
  https://www.cs.usask.ca/~spiteri/CMPT851/notes/PETSc.pdf

  this examples solves a 3-D rate state friction example, see Becker
  (2000)

  this holds RHS versions

*/
#ifdef USE_PETSC


#include "petscts.h"
#include "ode_headers.h"

/* 
   two state variable RHS ODE - x = log(v/v0), y = (tau-tau0)/a, z = b2 log(v0*theta2/dc2), steady state = {0,0,0}
*/
PetscErrorCode RHSFunction3D(TS ts,PetscReal time,Vec X,Vec F,void *ptr)
{
  PetscScalar *f;
  const PetscScalar *x;
  struct AppCtx *par;
  PetscScalar expx;
  PetscFunctionBeginUser;
  par = (struct AppCtx *)ptr;
  PetscCall(VecGetArrayRead(X,&x));PetscCall(VecGetArray(F,&f));
  expx = PetscExpReal(x[0]);
  if(!isfinite(expx)){
    /* a stage state went out of range; emit a non-finite RHS so the adaptive
       controller rejects this step and retries with a smaller dt, instead of
       aborting the whole run. myDomainError() normally catches this at the
       accepted-state level before it happens. */
    f[0] = f[1] = f[2] = PETSC_INFINITY;
    PetscCall(VecRestoreArrayRead(X,&x));PetscCall(VecRestoreArray(F,&f));
    PetscFunctionReturn(PETSC_SUCCESS);
  } /* keep this order, we are using f[1] and f[2] for f[0]! */
  f[1] =  (1.0 - expx) * par->k;
  f[2] = -expx * par->r * (par->b2 * x[0] + x[2]);
  f[0] = expx * ((par->b1 - 1.0) * x[0] + x[1] - x[2]) + f[1] -  f[2];
  PetscCall(VecRestoreArrayRead(X,&x));PetscCall(VecRestoreArray(F,&f));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* state + tangent RHS; tangent columns v1=x[3..5], v2=x[6..8], v3=x[9..11],
   dv/dt = J(x) v with the analytical Jacobian given in the header */
PetscErrorCode RHSFunction3DTangent(TS ts,PetscReal time,Vec X,Vec F,void *ptr)
{
  PetscScalar *f;
  const PetscScalar *x;
  struct AppCtx *par = (struct AppCtx *)ptr;
  PetscReal E,f0,f1,f2,J[3][3];
  PetscInt iv,i,j,off;
  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(X,&x));PetscCall(VecGetArray(F,&f));
  E = PetscExpReal(PetscRealPart(x[0]));
  if(!isfinite(E)){
    for(i=0;i<BASIN_NFULL;i++)
      f[i] = PETSC_INFINITY;
    PetscCall(VecRestoreArrayRead(X,&x));PetscCall(VecRestoreArray(F,&f));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  f1 =  (1.0 - E) * par->k;
  f2 = -E * par->r * (par->b2 * PetscRealPart(x[0]) + PetscRealPart(x[2]));
  f0 =  E * ((par->b1 - 1.0) * PetscRealPart(x[0]) + PetscRealPart(x[1]) - PetscRealPart(x[2])) + f1 - f2;
  f[0] = f0; f[1] = f1; f[2] = f2;
  /* Jacobian */
  J[0][0] = f0 - f1 + E * (par->b1 - 1.0 - par->k + par->r * par->b2);
  J[0][1] = E;
  J[0][2] = E * (par->r - 1.0);
  J[1][0] = -par->k * E;
  J[1][1] = 0.0;
  J[1][2] = 0.0;
  J[2][0] = f2 - par->r * par->b2 * E;
  J[2][1] = 0.0;
  J[2][2] = -par->r * E;
  for(iv=0;iv < 3;iv++){	/* three tangent vectors */
    off = BASIN_NSTATE + iv*3;
    for(i=0;i < 3;i++){
      f[off+i] = 0.0;
      for(j=0;j < 3;j++)
	f[off+i] += J[i][j] * x[off+j];
    }
  }
  PetscCall(VecRestoreArrayRead(X,&x));PetscCall(VecRestoreArray(F,&f));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RHSFunction4D(TS ts,PetscReal time,Vec X,Vec F,void *ptr) /* 4D
									    version
									    with
									    exponential
									    solve */
{
  PetscScalar *f;
  const PetscScalar *x;
  struct AppCtx *par;
  PetscFunctionBeginUser;
  par = (struct AppCtx *)ptr;
  PetscCall(VecGetArrayRead(X,&x));PetscCall(VecGetArray(F,&f));
  if(!isfinite(PetscRealPart(x[3]))){
    /* exp-state out of range; force a step rejection (see RHSFunction3D) */
    f[0] = f[1] = f[2] = f[3] = PETSC_INFINITY;
    PetscCall(VecRestoreArrayRead(X,&x));PetscCall(VecRestoreArray(F,&f));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  /* 
     
  */  
  f[1] =  (1.0 - x[3]) * par->k;
  f[2] = -x[3] * par->r * (par->b2 * x[0] + x[2]);
  f[0] =  x[3] * ((par->b1 - 1.0) * x[0] + x[1] - x[2]) + f[1] -  f[2];
  f[3] = f[0] * x[3]; // x[3] = exp(x), \dot{x[3]} = x[3] \dot{x} 
  PetscCall(VecRestoreArrayRead(X,&x));PetscCall(VecRestoreArray(F,&f));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif
