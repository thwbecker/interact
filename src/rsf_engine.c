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
   shared state evolution rate S = d psi/dt for global patch i, given
   psi and the SIGNED velocity vin, with the same expressions, guards,
   and floors as previously inlined in rsf_ODE_RHSFunction (see the
   long comment there for the laws and references).  If dSdpsi/dSdvabs
   are non-NULL, also return the partial derivatives of S with respect
   to psi (at fixed |v|) and with respect to |v| (at fixed psi); the
   floored-velocity branches return dSdvabs = 0 there, consistent with
   the floor.
*/
PetscReal rsf_state_rate(PetscInt i, PetscReal psi, PetscReal vin,
			 struct interact_ctx *par,
			 PetscReal *dSdpsi, PetscReal *dSdvabs)
{
  PetscReal S,vabs,psi_ss,epsi,omega,gate,lomega,b,D;
  PetscBool floored = PETSC_FALSE;
  struct rsf_vars *rsf;
  rsf = par->medium->rsf;
  b = par->fault[i].mu_d;
  if(b == 0.0){			/* b == 0 guard as in HBI */
    if(dSdpsi)*dSdpsi = 0.0;
    if(dSdvabs)*dSdvabs = 0.0;
    return 0.0;
  }
  D = (rsf->dc_vec)?(rsf->dc_vec[i]):(rsf->dc);
  if(rsf->state_law == RSF_SLIP_LAW){
    vabs = fabs(vin);
    if(vabs < rsf->vmin_state){
      vabs = rsf->vmin_state;floored = PETSC_TRUE;
    }
    psi_ss = rsf->f0 - b * PetscLogReal(vabs/rsf->v0);
    S = -(vabs/D) * (psi - psi_ss);
    if(dSdpsi)
      *dSdpsi = -vabs/D;
    if(dSdvabs)
      *dSdvabs = (floored)?(0.0):(-(psi - psi_ss)/D - b/D);
  }else if(rsf->state_law == RSF_PRZ_LAW){
    vabs = fabs(vin);
    epsi  = PetscExpReal((psi - rsf->f0)/b); /* exp((psi-f0)/b) */
    S = (b/(2.0*D)) * (rsf->v0/epsi - (vabs*vabs/rsf->v0)*epsi);
    if(dSdpsi)
      *dSdpsi = -(1.0/(2.0*D)) * (rsf->v0/epsi + (vabs*vabs/rsf->v0)*epsi);
    if(dSdvabs)
      *dSdvabs = -(b/D) * (vabs/rsf->v0) * epsi;
  }else if((rsf->state_law == RSF_SATO_LAW) || (rsf->state_law == RSF_KT_LAW)){
    vabs = fabs(vin);
    if(vabs < rsf->vmin_state){
      vabs = rsf->vmin_state;floored = PETSC_TRUE;
    }
    epsi  = PetscExpReal((psi - rsf->f0)/b);       /* exp((psi-f0)/b)   */
    omega = (vabs/rsf->v0) * epsi;		   /* Omega = |v| th/dc */
    gate  = (rsf->state_law == RSF_SATO_LAW)?
      (PetscExpReal(-omega/rsf->sato_beta)):
      (PetscExpReal(-vabs/rsf->kt_vc));
    lomega = PetscLogReal(vabs/rsf->v0) + (psi - rsf->f0)/b;  /* ln Omega */
    S = (b/D) * ((rsf->v0/epsi)*gate - vabs*lomega);
    if(dSdpsi){
      /* d/dpsi[(v0/epsi) gate] = (v0/epsi)(dgate/dpsi - gate/b);
	 Sato: dgate/dpsi = -gate omega/(b sato_beta); KT: 0 */
      *dSdpsi = (b/D) * ((rsf->v0/epsi) *
			 (((rsf->state_law == RSF_SATO_LAW)?
			   (-gate*omega/(b*rsf->sato_beta)):(0.0)) - gate/b)
			 - vabs/b);
    }
    if(dSdvabs){
      if(floored)
	*dSdvabs = 0.0;
      else
	*dSdvabs = (b/D) * ((rsf->v0/epsi) *
			    ((rsf->state_law == RSF_SATO_LAW)?
			     (-gate*epsi/(rsf->v0*rsf->sato_beta)):
			     (-gate/rsf->kt_vc))
			    - lomega - 1.0);
    }
  }else{			/* aging law */
    epsi = PetscExpReal((rsf->f0 - psi)/b); /* exp((f0-psi)/b), NOT the epsi above */
    S = (b/D) * (rsf->v0 * epsi - fabs(vin));
    if(dSdpsi)
      *dSdpsi = -(rsf->v0/D) * epsi;
    if(dSdvabs)
      *dSdvabs = -b/D;
  }
  return S;
}
/*
   fill rsf->vel from the state vector X and apply the interaction
   matvec(s) into rsf->tau_dot (and rsf->sigma_dot if enabled); shared
   between the monolithic RHS and the IMEX explicit RHS.  Includes the
   layout consistency check formerly at the top of rsf_ODE_RHSFunction.
*/
PetscErrorCode rsf_compute_vel_and_stressing(Vec X, struct interact_ctx *par)
{
  PetscScalar *vel,mu,scaled_tau,exp_fac;
  const PetscScalar *x;
  struct med *medium;struct flt *fault;
  PetscInt i,j,k,ln;
  struct rsf_vars *rsf;
  PetscFunctionBeginUser;
  medium = par->medium;fault = par->fault;rsf = medium->rsf;
  PetscCall(VecGetLocalSize(X, &ln));
  if(ln != rsf->dim*medium->rn){
    fprintf(stderr,"rsf_compute_vel_and_stressing: layout mismatch, local %i vs %i x %i\n",
	    (int)ln,rsf->dim,(int)medium->rn);
    exit(-1);
  }
  PetscCall(VecGetArrayRead(X,&x));
  PetscCall(VecGetArray(rsf->vel,&vel));
  /* i global patch, j local x offset, k local patch */
  for (i = medium->rs, j=0, k=0; i < medium->re; i++, j+=rsf->dim, k++) {
    /* local x layout: state = x[j], tau = x[j+1], sigma = x[j+2], slip = x[j+3] */
    vel[k] = vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,rsf->v0,
			  &mu, &scaled_tau, &exp_fac, medium);
  }
  PetscCall(VecRestoreArray(rsf->vel,&vel));
  PetscCall(VecRestoreArrayRead(X,&x));
  /*
     stressing rates from all slipping faults, dense or H matrix;
     background backslip loading (sinc) is added by the caller
  */
  PetscCall(MatMult(medium->Is, rsf->vel, rsf->tau_dot));
  if(rsf->calc_sigma_dot)
    PetscCall(MatMult(medium->In, rsf->vel, rsf->sigma_dot));
  PetscFunctionReturn(PETSC_SUCCESS);
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
  PetscScalar *f,cosh_fac,mu,scaled_tau,exp_fac,pre_fac;
  PetscScalar dvdtau,dvdsigma,dvdstate,a,sdot;
  const PetscScalar *x,*tau_dot,*sigma_dot,*velr;
  struct med *medium;struct flt *fault;
  PetscInt i,j,k;
  struct interact_ctx *par;
  struct rsf_vars *rsf;
  
  PetscFunctionBeginUser;
  par = (struct interact_ctx *)ptr;
  medium = par->medium;
  fault = par->fault;
  rsf = medium->rsf;
  
  /* layout consistency is checked inside rsf_compute_vel_and_stressing */
  PetscCall(rsf_compute_vel_and_stressing(X,par));
  PetscCall(VecGetArrayRead(X,&x));
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
    a = fault[i].mu_s;	/* b = fault[].mu_d is used inside rsf_state_rate */
    /* 
       state evolution, in the psi variable used throughout
       (implemented in rsf_state_rate above, shared with the IMEX path):

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
    /* state evolution rate; laws and derivatives implemented once in
       rsf_state_rate (shared with the IMEX path) */
    f[j] = rsf_state_rate(i,x[j],velr[k],par,NULL,NULL);
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
