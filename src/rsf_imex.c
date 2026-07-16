#include "interact.h"
#include "properties.h"
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#include "rsf.h"
/*

   IMEX (additive Runge-Kutta) variant of the rate-and-state ODE
   system solved in rsf_solve, selected with -imex.

   Motivation: the state evolution term is local (per cell) but can be
   very stiff away from steady state, most severely for the PRZ law
   (state_law 2), whose theta rate grows like Omega^2 above steady
   state.  At rupture fronts the state lags steady state by many
   e-folds regardless of law, and for PRZ the stable EXPLICIT step
   then collapses like 2 dc/(|v| Omega) (observed ~3e-8 s in BP5
   tests at 1-2 km).  The nonlocal stress transfer (the matvec), by
   contrast, is not stiff on those scales.  This is the classic IMEX
   split: integrate the local stiff terms implicitly, keep the matvec
   explicit, so the step size is set by the elastic coupling rather
   than by the state law.

   Split (per cell, X = (psi,tau,sigma,slip), v = 2 v0 exp(-psi/a)
   sinh(tau/(sigma a)), eta = G/(2 cs), denom = 1 + eta dv/dtau):

     EXPLICIT G(t,X)                        IMPLICIT F(t,X,Xdot) = Xdot - Fimpl(X)
     G_psi   = 0                            Fimpl_psi   = S(psi,|v|)
     G_sigma = sigma_dot_mv + sinc[1]       Fimpl_sigma = 0
     G_tau   = (tau_dot_mv + sinc[0]        Fimpl_tau   = -eta (dv/dpsi) S / denom
                - eta (dv/dsigma) G_sigma)
               / denom
     G_slip  = |v| -> v                     Fimpl_slip  = 0

   G + Fimpl reproduces rsf_ODE_RHSFunction exactly (the tau equation
   there is (tau_dot_mv + sinc[0] - eta(dv/dsigma f_sigma + dv/dpsi
   f_psi))/denom with f_psi = S), so -imex changes the integrator, not
   the ODE.  The implicit part is purely local: the IFunction does no
   matvec and the IJacobian is block diagonal, one 4x4 block per cell
   (rows psi and tau couple to columns psi,tau,sigma; sigma and slip
   rows are diagonal).  Block-Jacobi ILU(0) is then an EXACT solve for
   the stage systems, so the implicit stages cost per-cell work plus
   one reduction, no additional matvecs.

   The old explicit path (rsf_ODE_RHSFunction + TSRK) is untouched and
   remains the default; -imex switches the TS to TSARKIMEX with the
   functions below.  The state law itself lives in rsf_state_rate(),
   shared by both paths, so a new evolution law is added in one place.

   Jacobian note: the per-cell block Jacobian is analytic and exact,
   including the derivatives of the radiation-damping prefactor in the
   tau row (a frozen-prefactor quasi-Newton variant was tried first and
   caused frequent SNES line-search failures during coseismic phases,
   where the state rate S is large and the dropped S dQ/dX terms
   dominate).
*/

/*
   IMEX explicit right hand side G(t,X): the nonlocal stressing rates
   plus the non-stiff local terms; see the split in the header comment.
   Contains everything of rsf_ODE_RHSFunction EXCEPT the state rate S
   and the -eta (dv/dpsi) S/denom radiation-damping correction, which
   live in the IFunction.
*/
PetscErrorCode rsf_IMEX_RHSFunction(TS ts,PetscReal time,Vec X,Vec G,void *ptr)
{
  PetscScalar *g,cosh_fac,mu,scaled_tau,exp_fac,pre_fac;
  PetscScalar dvdtau,dvdsigma,a,sdot;
  const PetscScalar *x,*tau_dot,*sigma_dot,*velr;
  struct med *medium;struct flt *fault;
  PetscInt i,j,k;
  struct interact_ctx *par;
  struct rsf_vars *rsf;
  PetscFunctionBeginUser;
  par = (struct interact_ctx *)ptr;
  medium = par->medium;fault = par->fault;rsf = medium->rsf;
  PetscCall(rsf_compute_vel_and_stressing(X,par));
  PetscCall(VecGetArrayRead(X,&x));
  PetscCall(VecGetArrayRead(rsf->vel,&velr));
  PetscCall(VecGetArrayRead(rsf->tau_dot,&tau_dot));
  if(rsf->calc_sigma_dot)
    PetscCall(VecGetArrayRead(rsf->sigma_dot,&sigma_dot));
  PetscCall(VecGetArray(G,&g));
  for (i = medium->rs, j=0, k=0; i < medium->re; i++, j+=rsf->dim, k++) {
    a = fault[i].mu_s;
    /* state evolution is fully implicit */
    g[j] = 0.0;
    /* d sigma/dt, compression positive (In was scaled by -1) */
    if(rsf->calc_sigma_dot)
      sdot = sigma_dot[k] + fault[i].sinc[1];
    else
      sdot = 0.;
    if(rsf->limit_sigma){	/* HBI limitsigma behavior */
      if((x[j+2] < rsf->min_sigma)||(x[j+2] > rsf->max_sigma))
	sdot = 0.0;
    }
    g[j+2] = sdot;
    /* explicit part of d tau/dt (no dv/dpsi S term here) */
    mu = x[j+1]/x[j+2];				/* tau/sigma */
    scaled_tau = mu/a;				/* tau/(sigma a) */
    exp_fac = PetscExpReal(-x[j]/a);		/* exp(-psi/a) */
    pre_fac =  2.0*rsf->v0/(a * x[j+2]);	/* 2v0/(a sigma) */
    cosh_fac  = PetscCoshReal(scaled_tau) * exp_fac;
    dvdtau   =   pre_fac * cosh_fac;
    dvdsigma =  -pre_fac * cosh_fac * mu;
    g[j+1]  = (tau_dot[k] + fault[i].sinc[0]
	       - rsf->shear_mod_over_2cs_si * (dvdsigma * sdot));
    g[j+1] /= (1.0 + rsf->shear_mod_over_2cs_si * dvdtau);
    /* d slip/dt */
    g[j+3] = velr[k];
  }
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(G,&g));
  PetscCall(VecRestoreArrayRead(rsf->vel,&velr));
  PetscCall(VecRestoreArrayRead(rsf->tau_dot,&tau_dot));
  if(rsf->calc_sigma_dot)
    PetscCall(VecRestoreArrayRead(rsf->sigma_dot,&sigma_dot));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*
   IMEX implicit function F(t,X,Xdot) = Xdot - Fimpl(X); purely local,
   NO matvec: the sliding velocity is recomputed per cell from X.
*/
PetscErrorCode rsf_IMEX_IFunction(TS ts,PetscReal time,Vec X,Vec Xdot,Vec F,void *ptr)
{
  PetscScalar *f,cosh_fac,mu,scaled_tau,exp_fac,pre_fac;
  PetscScalar dvdtau,dvdstate,a,v,S,denom;
  const PetscScalar *x,*xd;
  struct med *medium;struct flt *fault;
  PetscInt i,j;
  struct interact_ctx *par;
  struct rsf_vars *rsf;
  PetscFunctionBeginUser;
  par = (struct interact_ctx *)ptr;
  medium = par->medium;fault = par->fault;rsf = medium->rsf;
  PetscCall(VecGetArrayRead(X,&x));
  PetscCall(VecGetArrayRead(Xdot,&xd));
  PetscCall(VecGetArray(F,&f));
  for (i = medium->rs, j=0; i < medium->re; i++, j+=rsf->dim) {
    a = fault[i].mu_s;
    v = vel_from_rsf(x[j+1],x[j+2],x[j],a,rsf->v0,
		     &mu,&scaled_tau,&exp_fac,medium);
    S = rsf_state_rate(i,x[j],v,par,NULL,NULL);
    pre_fac  = 2.0*rsf->v0/(a * x[j+2]);
    cosh_fac = PetscCoshReal(scaled_tau) * exp_fac;
    dvdtau   = pre_fac * cosh_fac;
    dvdstate = -v/a;
    denom = 1.0 + rsf->shear_mod_over_2cs_si * dvdtau;
    f[j]   = xd[j]   - S;
    f[j+1] = xd[j+1] + (rsf->shear_mod_over_2cs_si * dvdstate * S)/denom;
    f[j+2] = xd[j+2];
    f[j+3] = xd[j+3];
  }
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArrayRead(Xdot,&xd));
  PetscCall(VecRestoreArray(F,&f));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*
   IMEX implicit Jacobian J = shift dF/dXdot + dF/dX = shift I -
   dFimpl/dX; block diagonal, one block per cell: rows psi and tau
   have columns (psi,tau,sigma), rows sigma and slip are diagonal.
   The psi row is exact; the tau row freezes the prefactor Q (see
   header comment).
*/
PetscErrorCode rsf_IMEX_IJacobian(TS ts,PetscReal time,Vec X,Vec Xdot,PetscReal shift,
				  Mat A,Mat B,void *ptr)
{
  PetscScalar cosh_fac,mu,scaled_tau,exp_fac,pre_fac;
  PetscScalar dvdtau,dvdsigma,dvdstate,a,v,vabs,sgnv,S,denom,Q,E;
  PetscScalar dSdpsi,dSdvabs,Spsi,Stau,Ssig,vals[3];
  PetscScalar ddvdtau_psi,ddvdtau_tau,ddvdtau_sig,Qpsi,Qtau,Qsig;
  const PetscScalar *x;
  PetscInt i,j,gj,row,cols[3];
  struct med *medium;struct flt *fault;
  struct interact_ctx *par;
  struct rsf_vars *rsf;
  PetscFunctionBeginUser;
  par = (struct interact_ctx *)ptr;
  medium = par->medium;fault = par->fault;rsf = medium->rsf;
  PetscCall(VecGetArrayRead(X,&x));
  for (i = medium->rs, j=0; i < medium->re; i++, j+=rsf->dim) {
    a = fault[i].mu_s;
    v = vel_from_rsf(x[j+1],x[j+2],x[j],a,rsf->v0,
		     &mu,&scaled_tau,&exp_fac,medium);
    S = rsf_state_rate(i,x[j],v,par,&dSdpsi,&dSdvabs);
    vabs = fabs(v);sgnv = (v < 0.0)?(-1.0):(1.0);
    pre_fac  = 2.0*rsf->v0/(a * x[j+2]);
    cosh_fac = PetscCoshReal(scaled_tau) * exp_fac;
    dvdtau   =  pre_fac * cosh_fac;
    dvdsigma = -pre_fac * cosh_fac * mu;
    dvdstate = -v/a;
    E = rsf->shear_mod_over_2cs_si;
    denom = 1.0 + E * dvdtau;
    Q = E * dvdstate / denom;	/* Fimpl_tau = -Q S, i.e. F_tau = xdot_tau + Q S */
    /* total derivatives of S(psi,|v(psi,tau,sigma)|):
       d|v|/dpsi = -|v|/a, d|v|/dtau = sgn(v) dv/dtau, d|v|/dsigma = sgn(v) dv/dsigma */
    Spsi = dSdpsi + dSdvabs * (-vabs/a);
    Stau = dSdvabs * sgnv * dvdtau;
    Ssig = dSdvabs * sgnv * dvdsigma;
    /* derivatives of Q = E dvdstate/denom, with
       d(dvdstate)/dX = -(1/a) dv/dX and
       d(dvdtau)/dpsi = -dvdtau/a, d(dvdtau)/dtau = v/(a^2 sigma^2),
       d(dvdtau)/dsigma = -dvdtau/sigma - v tau/(a^2 sigma^3) */
    ddvdtau_psi = -dvdtau/a;
    ddvdtau_tau =  v/(a*a * x[j+2]*x[j+2]);
    ddvdtau_sig = -dvdtau/x[j+2] - v*x[j+1]/(a*a * x[j+2]*x[j+2]*x[j+2]);
    Qpsi = (E/(denom*denom)) * ((-dvdstate/a) /* = -(1/a) dv/dpsi = v/a^2 */ * denom
				- dvdstate * E * ddvdtau_psi);
    Qtau = (E/(denom*denom)) * ((-dvdtau/a)   * denom - dvdstate * E * ddvdtau_tau);
    Qsig = (E/(denom*denom)) * ((-dvdsigma/a) * denom - dvdstate * E * ddvdtau_sig);
    gj = i * rsf->dim;		/* global row offset of this cell (rank-contiguous layout) */
    cols[0] = gj;cols[1] = gj+1;cols[2] = gj+2;
    /* psi row: shift - dS/dX */
    row = gj;
    vals[0] = shift - Spsi;vals[1] = -Stau;vals[2] = -Ssig;
    PetscCall(MatSetValues(B,1,&row,3,cols,vals,INSERT_VALUES));
    /* tau row: F_tau = xdot_tau + Q S, so J = shift + Q dS/dX + S dQ/dX (exact) */
    row = gj+1;
    vals[0] = Q*Spsi + S*Qpsi;
    vals[1] = shift + Q*Stau + S*Qtau;
    vals[2] = Q*Ssig + S*Qsig;
    PetscCall(MatSetValues(B,1,&row,3,cols,vals,INSERT_VALUES));
    /* sigma and slip rows: diagonal shift */
    row = gj+2;
    PetscCall(MatSetValue(B,row,row,shift,INSERT_VALUES));
    row = gj+3;
    PetscCall(MatSetValue(B,row,row,shift,INSERT_VALUES));
  }
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY));
  if(A != B){
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*
   switch the given TS to the IMEX (TSARKIMEX) formulation: set the
   explicit RHS, the implicit function and its block-diagonal
   Jacobian, and solver defaults that make the implicit stage solve a
   direct per-cell operation (block Jacobi + ILU(0) is an exact
   factorization for this sparsity pattern).  Everything else
   (tolerances, adaptivity, domain check, events, monitor) is attached
   by the driver exactly as for the explicit path.  The created
   Jacobian matrix is returned for the driver to destroy.
*/
PetscErrorCode rsf_IMEX_setup(TS ts, struct interact_ctx *par, Mat *Jimex)
{
  Mat J;
  PetscInt *d_nnz,ln,k,c;
  struct med *medium;
  struct rsf_vars *rsf;
  PetscFunctionBeginUser;
  medium = par->medium;rsf = medium->rsf;
  ln = rsf->dim * medium->rn;
  /* block-diagonal Jacobian: per cell rows (psi,tau) have 3 entries,
     (sigma,slip) are diagonal; no off-rank coupling */
  PetscCall(PetscMalloc1(ln,&d_nnz));
  for(k=0;k < medium->rn;k++){
    for(c=0;c < rsf->dim;c++)
      d_nnz[k*rsf->dim+c] = (c < 2)?(3):(1);
  }
  PetscCall(MatCreateAIJ(PETSC_COMM_WORLD,ln,ln,PETSC_DETERMINE,PETSC_DETERMINE,
			 0,d_nnz,0,NULL,&J));
  PetscCall(PetscFree(d_nnz));
  PetscCall(MatSetOption(J,MAT_NO_OFF_PROC_ENTRIES,PETSC_TRUE));
  /*  */
  PetscCall(TSSetType(ts,TSARKIMEX));
  /* tableau default: ARKIMEX L2 (L-stable, 2nd order, stage order 2).
     The nominally more accurate ARKIMEX3 suffers stiff stage-order
     reduction on the rate-state stage problems: its embedded error
     estimator underestimates the state-variable error through fast
     transitions, so the controller accepts steps whose true error is
     far above the requested tolerance.  On the single-patch slider
     benchmark at rtol 1e-6 this cost a factor of about 50 in
     event-time accuracy relative to L2 (about 8e-2 versus 2e-3 yr RMS,
     uniformly across the evolution laws), with ARKIMEX4 worse still
     and ars443 numerically quenching the stick-slip cycle outright, so
     tableau choice here affects correctness, not just efficiency.
     L2 needs more steps (2nd order), which is the acceptable
     direction of the trade-off for a default; -ts_arkimex_type
     overrides for experimentation.  See slider/README.md. */
  PetscCall(TSARKIMEXSetType(ts,TSARKIMEXL2));
  PetscCall(TSSetRHSFunction(ts,NULL,rsf_IMEX_RHSFunction,par));
  PetscCall(TSSetIFunction(ts,NULL,rsf_IMEX_IFunction,par));
  PetscCall(TSSetIJacobian(ts,J,J,rsf_IMEX_IJacobian,par));
  /* the IMEX stage solvers get their own options prefix so that
     UNQUALIFIED global solver options (command line or the auto-read
     petsc_settings.yaml, which is tuned for the dense solves of the
     interact main program) cannot reach them; without this, a yaml
     containing e.g. ksp_type fgmres / pc_type jacobi / ksp_rtol 1e-8
     silently replaces the exact per-cell block solve with hundreds of
     poorly preconditioned Krylov iterations per stage, each with a
     global reduction, and degrades the Newton directions enough to
     cause step rejections (observed: an aging BP5 2 km run at 8.1e6
     linear iterations, 594 nonlinear failures, 146 s on 32 cores,
     versus seconds when configured as intended).  Deliberate overrides
     remain available as -imex_snes_*, -imex_ksp_*, -imex_pc_* */
  {
    SNES snes;
    PetscCall(TSGetSNES(ts,&snes));
    PetscCall(SNESSetOptionsPrefix(snes,"imex_"));
  }
  /* the block-diagonal stage solver (preonly + block Jacobi/ILU) is installed
     by rsf_IMEX_set_stage_solver, which the driver calls again after
     TSSetFromOptions so the options pass cannot silently revert it to an
     iterative Krylov solve over the full elastic operator */
  PetscCall(rsf_IMEX_set_stage_solver(ts));
  HEADNODE{
    fprintf(stderr,"rsf_IMEX_setup: IMEX (ARKIMEX, tableau l2) time integration: implicit local state terms, explicit stress transfer\n");
  }
  *Jimex = J;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Install the stage linear solver for the IMEX implicit systems. Those systems
  are block diagonal (one 4x4 block per cell, no inter-cell coupling), so the
  exact solve is preonly + block Jacobi with a per-block ILU(0), and it needs no
  global matvec. This is set both in rsf_IMEX_setup and again from the driver
  after TSSetFromOptions, because that call resets the SNES KSP/PC to PETSc
  defaults. To keep a deliberate override working, only force the types when
  the user has NOT set -imex_ksp_type / -imex_pc_type explicitly (the stage
  solvers live under the "imex_" options prefix, see rsf_IMEX_setup).
*/
PetscErrorCode rsf_IMEX_set_stage_solver(TS ts)
{
  SNES snes;KSP ksp;PC pc;PetscBool kset,pcset;
  PetscFunctionBeginUser;
  PetscCall(TSGetSNES(ts,&snes));
  PetscCall(SNESGetKSP(snes,&ksp));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-imex_ksp_type",&kset));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-imex_pc_type",&pcset));
  if(!kset)
    PetscCall(KSPSetType(ksp,KSPPREONLY));
  if(!pcset)
    PetscCall(PCSetType(pc,PCBJACOBI));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
