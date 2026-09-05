/*

  solving ordinary differential equations originally based on
  https://www.cs.usask.ca/~spiteri/CMPT851/notes/PETSc.pdf

  this examples solves a 3-D rate state friction example, see Becker
  (2000)

  this holds RHS versions

*/
#ifdef USE_PETSC


#include "ode_headers.h"

PetscErrorCode  init_stiff_par(struct AppCtx *par, PetscReal knd_default)
{
  PetscReal kcr1,kcr2;
  par->knd = knd_default;  
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-knd",&par->knd,NULL));


  par->b1 = 1.0;
  par->b2 = 0.84;
  par->r  = 0.048;
  
  kcr1 = par->b1 - 1.0;
  kcr2 = (kcr1 + par->r*(2.*par->b1 + (par->b2 - 1.)*(2. + par->r)) +
	  sqrt(4.*par->r*par->r*(kcr1 + par->b2) +
	       pow(kcr1 + par->r*par->r*(par->b2 - 1.),2)))/(2. + 2.*par->r);
  par->k = par->knd * kcr2;
  PetscFunctionReturn(PETSC_SUCCESS);
}


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

/* plain state RHS evaluation for crossing interpolation */
void eval_rhs(const struct AppCtx *par,const PetscReal x[3],PetscReal f[3])
{
  PetscReal E = PetscExpReal(x[0]);
  f[1] =  (1.0 - E) * par->k;
  f[2] = -E * par->r * (par->b2 * x[0] + x[2]);
  f[0] =  E * ((par->b1 - 1.0) * x[0] + x[1] - x[2]) + f[1] - f[2];
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


/*
   orthonormal basis {d1, d2} of the most unstable invariant plane of
   the Jacobian at the steady sliding fixed point {0,0,0}.

   the eigenvalues of J(0) are found from the characteristic cubic
   (closed form, cf. Numerical Recipes sec. 5.6). for a complex
   conjugate pair (the generic Hopf-unstable case for knd < 1, cf.
   Rice & Ruina 1983; Gu et al. 1984), the plane is spanned by the
   real and imaginary parts of the corresponding eigenvector; if all
   eigenvalues are real, the eigenvectors of the two largest are used.
   eigenvectors are obtained as complex cross products of two rows of
   (J - lambda I). returns Re and Im of the leading eigenvalue.
   signs are fixed such that the largest-magnitude component of each
   basis vector is positive.
*/
#include <complex.h>
PetscErrorCode unstable_plane(const struct AppCtx *par,PetscReal d1[3],PetscReal d2[3],
			      PetscReal *re_lam,PetscReal *im_lam)
{
  PetscReal J[3][3],a,b,c,Q,R,theta,A,B,sq,nrm,dot,amax;
  double complex lam[3],rows[3][3],v[3],w[3],tmp;
  PetscInt i,j,imax,i1,i2;
  PetscFunctionBeginUser;
  /* Jacobian at the origin (E = 1, f = 0) */
  J[0][0] = par->b1 - 1.0 - par->k + par->r * par->b2;
  J[0][1] = 1.0;
  J[0][2] = par->r - 1.0;
  J[1][0] = -par->k;  J[1][1] = 0.0; J[1][2] = 0.0;
  J[2][0] = -par->r * par->b2; J[2][1] = 0.0; J[2][2] = -par->r;
  /* characteristic cubic lam^3 + a lam^2 + b lam + c = 0 */
  a = -(J[0][0] + J[1][1] + J[2][2]);
  b =  (J[0][0]*J[1][1] - J[0][1]*J[1][0]) +
       (J[0][0]*J[2][2] - J[0][2]*J[2][0]) +
       (J[1][1]*J[2][2] - J[1][2]*J[2][1]);
  c = -(J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
      - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
      + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]));
  Q = (a*a - 3.0*b)/9.0;
  R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
  if(R*R < Q*Q*Q){		/* three real roots */
    theta = acos(R/PetscSqrtReal(Q*Q*Q));
    for(i=0;i < 3;i++)
      lam[i] = -2.0*PetscSqrtReal(Q)*cos((theta + 2.0*M_PI*(PetscReal)i)/3.0) - a/3.0;
  }else{			/* one real root, complex pair */
    A = -copysign(pow(fabs(R) + PetscSqrtReal(R*R - Q*Q*Q),1.0/3.0),R);
    B = (A != 0.0) ? Q/A : 0.0;
    sq = 0.5*sqrt(3.0)*(A - B);
    lam[0] = (A + B) - a/3.0;			  /* real */
    lam[1] = -0.5*(A + B) - a/3.0 + sq*I;
    lam[2] = -0.5*(A + B) - a/3.0 - sq*I;
  }
  /* sort by descending real part (i1: leading, i2: second) */
  i1 = 0;
  for(i=1;i < 3;i++)
    if(creal(lam[i]) > creal(lam[i1]))
      i1 = i;
  i2 = (i1 == 0) ? 1 : 0;
  for(i=0;i < 3;i++)
    if((i != i1) && (creal(lam[i]) > creal(lam[i2])))
      i2 = i;
  *re_lam = creal(lam[i1]);
  *im_lam = fabs(cimag(lam[i1]));
  /* eigenvector of lam[i1] via complex cross product of two rows of (J - lam I),
     picking the pair with the largest norm for robustness */
  for(i=0;i < 3;i++)
    for(j=0;j < 3;j++)
      rows[i][j] = J[i][j] - ((i == j) ? lam[i1] : 0.0);
  {
    double complex cr[3][3];
    PetscReal cn,cnmax = -1.0;
    PetscInt ic,jc,ibest = 0,jbest = 1;
    for(ic=0;ic < 3;ic++)
      for(jc=ic+1;jc < 3;jc++){
	cr[0][0] = rows[ic][1]*rows[jc][2] - rows[ic][2]*rows[jc][1];
	cr[0][1] = rows[ic][2]*rows[jc][0] - rows[ic][0]*rows[jc][2];
	cr[0][2] = rows[ic][0]*rows[jc][1] - rows[ic][1]*rows[jc][0];
	cn = 0.0;
	for(i=0;i < 3;i++)
	  cn += creal(cr[0][i]*conj(cr[0][i]));
	if(cn > cnmax){
	  cnmax = cn; ibest = ic; jbest = jc;
	  for(i=0;i < 3;i++)
	    v[i] = cr[0][i];
	}
      }
    (void)ibest;(void)jbest;
  }
  if(cimag(lam[i1]) != 0.0){
    /* complex pair: plane from Re(v), Im(v) */
    for(i=0;i < 3;i++){
      d1[i] = creal(v[i]);
      d2[i] = cimag(v[i]);
    }
  }else{
    /* all real: second direction from the eigenvector of lam[i2] */
    for(i=0;i < 3;i++){
      d1[i] = creal(v[i]);
      for(j=0;j < 3;j++)
	rows[i][j] = J[i][j] - ((i == j) ? lam[i2] : 0.0);
    }
    w[0] = rows[0][1]*rows[1][2] - rows[0][2]*rows[1][1];
    w[1] = rows[0][2]*rows[1][0] - rows[0][0]*rows[1][2];
    w[2] = rows[0][0]*rows[1][1] - rows[0][1]*rows[1][0];
    tmp = w[0]*conj(w[0]) + w[1]*conj(w[1]) + w[2]*conj(w[2]);
    if(creal(tmp) < 1e-20){	/* degenerate pair of rows, use other rows */
      w[0] = rows[1][1]*rows[2][2] - rows[1][2]*rows[2][1];
      w[1] = rows[1][2]*rows[2][0] - rows[1][0]*rows[2][2];
      w[2] = rows[1][0]*rows[2][1] - rows[1][1]*rows[2][0];
    }
    for(i=0;i < 3;i++)
      d2[i] = creal(w[i]);
  }
  /* orthonormalize {d1, d2} by Gram-Schmidt */
  nrm = 0.0;
  for(i=0;i < 3;i++)
    nrm += d1[i]*d1[i];
  nrm = PetscSqrtReal(nrm);
  for(i=0;i < 3;i++)
    d1[i] /= nrm;
  dot = 0.0;
  for(i=0;i < 3;i++)
    dot += d2[i]*d1[i];
  for(i=0;i < 3;i++)
    d2[i] -= dot*d1[i];
  nrm = 0.0;
  for(i=0;i < 3;i++)
    nrm += d2[i]*d2[i];
  nrm = PetscSqrtReal(nrm);
  for(i=0;i < 3;i++)
    d2[i] /= nrm;
  /* sign convention: largest-magnitude component positive */
  imax = 0; amax = fabs(d1[0]);
  for(i=1;i < 3;i++)
    if(fabs(d1[i]) > amax){
      amax = fabs(d1[i]); imax = i;
    }
  if(d1[imax] < 0.0)
    for(i=0;i < 3;i++)
      d1[i] = -d1[i];
  imax = 0; amax = fabs(d2[0]);
  for(i=1;i < 3;i++)
    if(fabs(d2[i]) > amax){
      amax = fabs(d2[i]); imax = i;
    }
  if(d2[imax] < 0.0)
    for(i=0;i < 3;i++)
      d2[i] = -d2[i];
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif
