#include "interact.h"
#include "properties.h"
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#include "rsf.h"




/* 
   compute the sliding velocity from tau, sigma, and state

   v = 2 v0 exp(-state/a) sinh(tau/(sigma a))

   and a bunch of helper factors to be reused
   mu = tau/sigma
   scaled_tau = mu/a
   exp_fac = exp(-state/a) 
 */
PetscReal vel_from_rsf(PetscReal tau, PetscReal sigma, PetscReal state, PetscReal a,PetscReal v0,
		       PetscReal *mu, PetscReal *scaled_tau, PetscReal *exp_fac, 
		       struct med *medium)
{
  PetscReal vel;
  *mu = tau/sigma;
  *scaled_tau = (*mu)/a;
  *exp_fac = PetscExpReal(-state/a);
  vel =  2.*v0 * (*exp_fac) * PetscSinhReal(*scaled_tau);
  return vel;
}



/* 
   compute time-derivatives for the state vector
   x = (psi1,tau1,sigma1,s1,psi2,tau2,sigma2,s2,....) of size medium->rsf->dim*n

   this implements the same ODE as HBI's derivs()/deriv() for the
   aging law without fluid pressure or viscous flow, see comments at
   top of file

   the local portion of X/F on each rank corresponds to patches
   [medium->rs, medium->re) of the interaction matrix row layout, by
   construction of x in the driver above
*/
PetscErrorCode rsf_ODE_RHSFunction(TS ts,PetscReal time,Vec X,Vec F,void *ptr)
{
  PetscScalar *f,*vel,cosh_fac,mu,scaled_tau,exp_fac,pre_fac;
  PetscScalar dvdtau,dvdsigma,dvdstate,a,b,sdot;
  const PetscScalar *x,*tau_dot,*sigma_dot,*velr;
  struct med *medium;struct flt *fault;
  PetscInt i,j,k,ln;
  struct interact_ctx *par;
  struct rsf_vars *rsf;
  
  PetscFunctionBeginUser;
  par = (struct interact_ctx *)ptr;
  medium = par->medium;
  fault = par->fault;
  rsf = medium->rsf;
  
  /* consistency check of the layout alignment (cheap) */
  PetscCall(VecGetLocalSize(X, &ln));
  if(ln != rsf->dim*medium->rn){
    fprintf(stderr,"rsf_ODE_RHSFunction: layout mismatch, local %i vs %i x %i\n",
	    (int)ln,rsf->dim,(int)medium->rn);
    exit(-1);
  }
  PetscCall(VecGetArrayRead(X,&x));
  /* 
     compute the slip rate vector
     v = 2 v0 exp(-psi/a) sinh(tau/(sigma a)) 
  */
  /* 
     NOTE (possible optimization, intentionally NOT applied here):
     this loop and the derivative loop below share per-cell factors.
     vel_from_rsf returns mu = tau/sigma, scaled_tau = mu/a and
     exp_fac = exp(-psi/a); the derivative loop then recomputes all three
     and additionally calls cosh(scaled_tau).  One could instead form
     exp(scaled_tau) ONCE here, get both sinh and cosh from it
     (sinh=(e-1/e)/2, cosh=(e+1/e)/2), and stash the single combination the
     derivative loop needs, cosh_fac = cosh(scaled_tau)*exp(-psi/a), in a
     scratch array sized to medium->rn for reuse below.  On a BP5 2 km dense
     serial test that made the RHS *arithmetic* ~1.6-1.7x faster, machine-
     precision-identical (first-event time unchanged to ~1e-9 yr).  The
     end-to-end effect is small here because the matvec dominates the step;
     it would matter more where the matvec is cheap (H-matrix, high core
     count).  It is left explicit on purpose: the gain is modest and a
     precomputed/stashed form is easy to get subtly wrong if the friction
     formulation is later changed, so clarity is preferred for now.
  */
  PetscCall(VecGetArray(rsf->vel,&vel));
  /* i global patch, j local x offset, k local patch */
  for (i = medium->rs, j=0, k=0; i < medium->re; i++, j+=rsf->dim, k++) {
    /* local x layout: state = x[j], tau = x[j+1], sigma = x[j+2], slip = x[j+3] */
    vel[k] = vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,rsf->v0,
			  &mu, &scaled_tau, &exp_fac, medium);
  }
  PetscCall(VecRestoreArray(rsf->vel,&vel));
  /* 
     stressing rates from all slipping faults, this works for dense
     and H matrix; background backslip loading (sinc) is added below
  */
  PetscCall(MatMult(medium->Is, rsf->vel, rsf->tau_dot));
  if(rsf->calc_sigma_dot)
    PetscCall(MatMult(medium->In, rsf->vel, rsf->sigma_dot));
  /* 
     compute derivatives 
  */
  /* unpack */
  PetscCall(VecGetArrayRead(rsf->vel,&velr));
  PetscCall(VecGetArrayRead(rsf->tau_dot,&tau_dot));
  if(rsf->calc_sigma_dot)
    PetscCall(VecGetArrayRead(rsf->sigma_dot,&sigma_dot));
  /*  */
  PetscCall(VecGetArray(F,&f));
  for (i = medium->rs, j=0, k=0; i < medium->re; i++, j+=rsf->dim, k++) {
    /* a = fault[].mu_s, b = fault[].mu_d as read by read_rsf */
    a = fault[i].mu_s;b = fault[i].mu_d;
    /* 
       d psi/dt = b/dc (v0 exp((f0-psi)/b) - |v|)
       
       aging law; HBI uses |v|, important for the sinh regularization
       which allows v < 0
    */
    if(b != 0.0)		/* b == 0 guard as in HBI */
      f[j] = (b/((rsf->dc_vec)?(rsf->dc_vec[i]):(rsf->dc))) *
	(rsf->v0 * PetscExpReal((rsf->f0 - x[j])/b) - fabs(velr[k]));
    else
      f[j] = 0.0;
    /* 
       d sigma/dt, compression positive (In was scaled by -1)
    */
    if(rsf->calc_sigma_dot)
      sdot = sigma_dot[k] + fault[i].sinc[1];
    else{
      sdot = 0.;
    }
    if(rsf->limit_sigma){	/* HBI limitsigma behavior */
      if((x[j+2] < rsf->min_sigma)||(x[j+2] > rsf->max_sigma))
	sdot = 0.0;
    }
    f[j+2] = sdot;
    /* 
       d tau/dt from the quasi-dynamic strength = stress condition,
       eta = G/(2cs)
    */
    mu = x[j+1]/x[j+2];				/* tau/sigma */
    scaled_tau = mu/a;				/* tau/(sigma a) */
    exp_fac = PetscExpReal(-x[j]/a);		/* exp(-psi/a) */
    pre_fac =  2.0*rsf->v0/(a * x[j+2]);	/* 2v0/(a sigma) */
    cosh_fac  = PetscCoshReal(scaled_tau) * exp_fac;
    dvdtau   =   pre_fac * cosh_fac;
    dvdsigma =  -pre_fac * cosh_fac * mu;
    dvdstate = -velr[k]/a;
    f[j+1]  = (tau_dot[k] + fault[i].sinc[0] 
	       - rsf->shear_mod_over_2cs_si * (dvdsigma * f[j+2] + dvdstate * f[j]));
    f[j+1] /= (1.0 + rsf->shear_mod_over_2cs_si * dvdtau);
    /* d slip/dt */
    f[j+3] = velr[k];
  }
  /*  */
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(F,&f));
  PetscCall(VecRestoreArrayRead(rsf->vel,&velr));
  PetscCall(VecRestoreArrayRead(rsf->tau_dot,&tau_dot));
  if(rsf->calc_sigma_dot)
    PetscCall(VecRestoreArrayRead(rsf->sigma_dot,&sigma_dot));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* 
   domain check: a stage/solution state is acceptable only if all
   entries are finite, normal stress is positive, and the friction
   arguments will not overflow exp/sinh; collective so that all ranks
   agree
*/
PetscErrorCode rsf_domain_check(TS ts,PetscReal time,Vec X,PetscBool *accept)
{
  const PetscScalar *x;
  struct med *medium;struct flt *fault;
  PetscInt i,j,lok,gok;
  PetscReal a,b;
  const PetscReal arg_max = 600.; /* exp/sinh overflow guard */
  struct rsf_vars *rsf;

  PetscFunctionBeginUser;
  medium = rsf_par_static->medium;
  rsf = medium->rsf;
  fault = rsf_par_static->fault;

  lok = 1;
  PetscCall(VecGetArrayRead(X,&x));
  for (i = medium->rs, j=0; (i < medium->re) && lok; i++, j+=rsf->dim) {
    a = fault[i].mu_s;b = fault[i].mu_d;
    if((!PetscIsNormalReal(x[j+2])) || (x[j+2] <= 0.0)){lok = 0;break;} /* sigma */
    if(PetscIsInfOrNanReal(x[j]) || PetscIsInfOrNanReal(x[j+1]) || PetscIsInfOrNanReal(x[j+3])){lok = 0;break;}
    if(fabs(x[j+1]/(x[j+2]*a)) > arg_max){lok = 0;break;} /* sinh(tau/(sigma a)) */
    if(fabs(x[j]/a) > arg_max){lok = 0;break;}		  /* exp(-psi/a) */
    if((b != 0.0) && (fabs((rsf->f0 - x[j])/b) > arg_max)){lok = 0;break;} /* aging law exp */
  }
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCallMPI(MPI_Allreduce(&lok,&gok,1,MPIU_INT,MPI_MIN,PETSC_COMM_WORLD));
  *accept = (gok)?(PETSC_TRUE):(PETSC_FALSE);
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
