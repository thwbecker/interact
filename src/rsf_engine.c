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
  PetscScalar dvdtau,dvdsigma,dvdstate,a,b,sdot,vabs,psi_ss,epsi,omega,gate,lomega;
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
       state evolution, in the psi variable used throughout:

       aging law (state_law 0, default):
         d psi/dt = b/dc (v0 exp((f0-psi)/b) - |v|)

       slip law (state_law 1):
         d psi/dt = -(|v|/dc) (psi - psi_ss),  psi_ss = f0 - b ln(|v|/v0)

       PRZ law (state_law 2; Perrin, Rice and Zheng, 1995), in the normalization
       with Omega = |v| theta/dc and Omega_ss = 1:
         d theta/dt = 1/2 ( 1 - (|v| theta/dc)^2 )
       Substituting theta = (dc/v0) exp((psi-f0)/b) gives
         d psi/dt = b/(2 dc) ( v0 exp((f0-psi)/b) - v^2/v0 exp((psi-f0)/b) )

       All three laws share the steady state psi_ss = f0 - b ln(|v|/v0), i.e. the
       usual f_ss = f0 + (a-b) ln(|v|/v0), and all three have the same
       linearization about it (d theta/dt ~ -(Omega-1)), so they agree on linear
       stability and differ only away from steady state. Note that theta carries a
       constant-factor arbitrariness: PRZ is sometimes written
       d theta/dt = 1 - (|v| theta/(2 dc))^2, which is the same law under
       theta -> theta/2 but puts theta_ss at 2 dc/|v|, shifting psi_ss by b ln 2
       and making it NOT directly comparable to aging and slip at equal dc. The
       Omega_ss = 1 form above is used here so that all three laws are comparable.

       Sato-type law (state_law 3) and the Kato and Tullis composite law
       (state_law 4): both are the slip law plus the aging law's healing term
       gated so that healing only operates near stationary contact,
         d theta/dt = gate - Omega ln Omega,   Omega = |v| theta/dc
       with gate = exp(-Omega/sato_beta) for the Sato-type law and
       gate = exp(-|v|/kt_vc) for Kato and Tullis. In psi,
         d psi/dt = b/dc ( v0 exp((f0-psi)/b) gate - |v| ln Omega )
       whose second term is exactly the slip law. sato_beta and kt_vc are held
       fixed (see rsf_init.c). Steady state: the Sato gate is negligible at
       Omega ~ 1 so Omega_ss = 1 as for the other laws, whereas the Kato and
       Tullis gate does not vanish for |v| << kt_vc, which shifts psi_ss up by
       b W(1) = 0.567 b there (the composite law's known steady-state offset);
       see rsf_solve.md.

       |v| is used (not v) for the sinh regularization, which allows v < 0. In
       the slip law |v| is floored at vmin_state so ln(|v|/v0) stays finite as
       v -> 0 (where the slip law's d psi/dt -> 0 in any case).
    */
    if(b != 0.0){		/* b == 0 guard as in HBI */
      if(rsf->state_law == RSF_SLIP_LAW){
	vabs = fabs(velr[k]);
	if(vabs < rsf->vmin_state)
	  vabs = rsf->vmin_state;
	psi_ss = rsf->f0 - b * PetscLogReal(vabs/rsf->v0);
	f[j] = -(vabs/((rsf->dc_vec)?(rsf->dc_vec[i]):(rsf->dc))) * (x[j] - psi_ss);
      }else if(rsf->state_law == RSF_PRZ_LAW){
	vabs = fabs(velr[k]);
	epsi  = PetscExpReal((x[j] - rsf->f0)/b); /* exp((psi-f0)/b) */
	f[j] = (b/(2.0*((rsf->dc_vec)?(rsf->dc_vec[i]):(rsf->dc)))) *
	  (rsf->v0/epsi - (vabs*vabs/rsf->v0)*epsi);
      }else if((rsf->state_law == RSF_SATO_LAW) || (rsf->state_law == RSF_KT_LAW)){
	/*
	  the two gated (composite) laws. Both are the slip law plus the aging
	  law's healing term multiplied by a gate that shuts the healing off away
	  from stationary contact; they differ only in what the gate keys on:

	    Sato-type (state_law 3): gate = exp(-Omega/sato_beta)   (keys on Omega)
	    Kato and Tullis (state_law 4, composite law):
	                             gate = exp(-|v|/kt_vc)         (keys on |v|)

	  In theta form, d theta/dt = gate - Omega ln Omega, with Omega = |v| theta/dc.
	*/
	vabs = fabs(velr[k]);
	if(vabs < rsf->vmin_state)
	  vabs = rsf->vmin_state;
	epsi  = PetscExpReal((x[j] - rsf->f0)/b);       /* exp((psi-f0)/b)   */
	omega = (vabs/rsf->v0) * epsi;			/* Omega = |v| th/dc */
	gate  = (rsf->state_law == RSF_SATO_LAW)?
	  (PetscExpReal(-omega/rsf->sato_beta)):
	  (PetscExpReal(-vabs/rsf->kt_vc));
	lomega = PetscLogReal(vabs/rsf->v0) + (x[j] - rsf->f0)/b;  /* ln Omega */
	f[j] = (b/((rsf->dc_vec)?(rsf->dc_vec[i]):(rsf->dc))) *
	  ((rsf->v0/epsi)*gate - vabs*lomega);
      }else{			/* aging law */
	f[j] = (b/((rsf->dc_vec)?(rsf->dc_vec[i]):(rsf->dc))) *
	  (rsf->v0 * PetscExpReal((rsf->f0 - x[j])/b) - fabs(velr[k]));
      }
    }else
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
    /* the exp(+-(f0-psi)/b) overflow test applies to the aging law and to PRZ,
       both of which exponentiate (psi-f0)/b; the slip law has no such term */
    if((rsf->state_law != RSF_SLIP_LAW) && (b != 0.0) &&
       (fabs((rsf->f0 - x[j])/b) > arg_max)){lok = 0;break;} /* aging/PRZ exp */
  }
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCallMPI(MPI_Allreduce(&lok,&gok,1,MPIU_INT,MPI_MIN,PETSC_COMM_WORLD));
  *accept = (gok)?(PETSC_TRUE):(PETSC_FALSE);
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
