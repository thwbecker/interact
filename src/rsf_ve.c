#include "interact.h"
#include "properties.h"
/*

  visco-elastic hereditary stressing for rsf_solve: Prony
  representation of the step-response interaction kernel

    K(t) = C_inf + sum_{p=1..np} C_p exp(-t/tau_p),  K(0) = K^elastic

  carried by per-patch memory states h_p with

    dh_p/dt = -h_p/tau_p + [C_p (v - vpl)],
    tau_dot = [K^el (v - vpl)] - sum_p h_p/tau_p

  (the relaxed operator C_inf never appears explicitly: it is
  implicitly K^el - sum_p C_p, so the elastic response stays exactly
  the assembled Is matrix).  This is the generator-agnostic
  rsf_ve_implementation_plan.md machinery: the dynamics below consume
  only rates and amplitude operators, so 3-D kernels (e.g. from the
  plate-over-Maxwell propagator, via a future -ve_prony_file) plug in
  without touching the time stepping.

  h treatment: the memory-light "step" mode (plan step 5): h lives
  OUTSIDE the TS state.  Between accepted steps the exact-exponential
  update

    h_p <- h_p e^(-dt/tau_p) + [C_p (vbar - vpl)] tau_p (1 - e^(-dt/tau_p))

  is applied in a TS monitor, with vbar = (slip increment)/dt the
  EXACT step-average velocity (slip is part of the TS state).  Inside
  a step the RHS sees the decay-exact, forcing-lagged states
  h_p(t) = h_p(t0) e^(-(t-t0)/tau_p), applied as a per-pole VecAXPY
  in rsf_compute_vel_and_stressing.  The lag error is second order in
  the step and negligible wherever VE matters (tau_p >> coseismic
  time scales).

  kernel generators (option -ve_mode):

    1  uniform Maxwell medium.  Implemented for 2-D antiplane
       geometries, where the kernel depends on the moduli only
       through G so the representation is a SINGLE exponential with
       tau_1 = t_M and C_1 = K^el, exactly (no sampling, no fit).
       (General 3-D geometries have the three-pole homogeneous
       representation of prony_kernel.c; not wired here yet.)

    2  elastic plate over Maxwell half space, 2-D antiplane: the
       validated layered image kernel (both image families, see
       ve_rybicki_antiplane_uz) sampled in the Laplace domain through
       the actual elements and fitted per pair on a geometric rate
       ladder [b/45, 1.5 b], b = g1/((g1+g2) t_M), samples
       [b/100, 5 b], with two held-out samples verified; the fitted
       constant is discarded (relaxed operator implicit).  This is
       the operator-level version of the ve_sp_cycle testbed.

  options (all parsed here):
    -ve_mode <1|2>         generator (required to switch VE on)
    -ve_tmaxwell_yr <t>    Maxwell time of the relaxing medium [yr]
    -ve_plate_h <H>        elastic plate thickness [m] (mode 2)
    -ve_g2fac <f>          substrate/plate rigidity ratio (mode 2, default 1)
    -ve_np <n>             exponential terms (mode 2, default 6, max 6)
    -ve_h_stage <0|1>      1 (default) stage-consistent sink forcing;
                           0 frozen decay-only sink (see rsf_ve_apply_sink)
    -ve_h_init <0|1>       initial memory states: 1 (default) backslip
                           steady state (no secular VE spin-up),
                           0 virgin medium (zero memory)

  the sampling image depth and the held-out abort threshold are
  compile-time constants (VE_NIMG_DEF, VE_FIT_TOL_DEF in
  properties.h): good values are geometry independent, see there

  restart: the h states, slip_prev, and t0 travel in the checkpoint
  (version 1.1); a VE run cannot restart from an elastic checkpoint
  and vice versa.

  part of interact

*/
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#include "rsf.h"
#include "prony_kernel.h"

/*
   parse options and, if VE is requested, build rates, amplitude
   operators, and state vectors.  Called from rsf_solve_run AFTER the
   Is matrix and the backslip products exist.  negvpl_in holds -vpl
   per patch (the slip-rate vector used for the sinc products).
   scale converts interact-internal stress to SI, as for Is.
*/
PetscErrorCode rsf_ve_setup(struct interact_ctx *par, Vec negvpl_in,
			    PetscReal scale)
{
  struct med *medium = par->medium;
  struct rsf_vars *rsf = medium->rsf;
  struct flt *fault = par->fault;
  struct prony_spec spec;
  PetscReal tm_yr=0.0,plate_h=0.0,g2fac=1.0;
  const PetscReal fit_tol = VE_FIT_TOL_DEF;
  PetscReal gamma0,rate_b,worst_res;
  PetscInt np=6,mode=0;
  const PetscInt n_img = VE_NIMG_DEF;
  PetscInt i,j,k,n,it,p,jn;
  PetscBool flg;
  PetscFunctionBeginUser;
  rsf->ve_np = 0;rsf->ve_t0 = 0.0;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ve_mode",&mode,&flg));
  if((!flg) || (mode == 0))
    PetscFunctionReturn(PETSC_SUCCESS); /* VE off: nothing changes */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-ve_tmaxwell_yr",&tm_yr,&flg));
  if((!flg)||(tm_yr <= 0.0))
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
	    "rsf_ve_setup: -ve_mode %ld requires -ve_tmaxwell_yr > 0",(long)mode);
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ve_np",&np,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-ve_plate_h",&plate_h,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-ve_g2fac",&g2fac,NULL));
  /* both implemented generators require pure 2-D antiplane geometry */
  for(i=0;i < medium->nrflt;i++)
    if(fault[i].type != TWO_DIM_ANTIPLANE)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,
	      "rsf_ve_setup: -ve_mode %ld: patch %ld type %i, only 2-D antiplane elements are wired up (3-D kernels: future -ve_prony_file)",
	      (long)mode,(long)i,(int)fault[i].type);
  if(rsf->calc_sigma_dot)
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,
	    "rsf_ve_setup: VE normal-stress coupling not implemented (antiplane has none); drop -calc_sigma_dot");
  rsf->ve_mode = mode;
  rsf->ve_hstage = 1;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ve_h_stage",&rsf->ve_hstage,NULL));
  if(mode == 1){
    /* uniform Maxwell, antiplane: exactly one pole, C_1 = K^el */
    rsf->ve_np = 1;
    rsf->ve_tau[0] = tm_yr*SEC_PER_YEAR;
    PetscCall(MatDuplicate(medium->Is,MAT_COPY_VALUES,&rsf->ve_C[0]));
    HEADNODE
      fprintf(stderr,"rsf_ve_setup: mode 1 (uniform Maxwell, antiplane): single EXACT pole, tau = %g yr, C_1 = Is\n",
	      tm_yr);
  }else if(mode == 2){
    /* layered antiplane: sampled image kernel, per-pair Prony fit */
    struct flt img[1];
    COMP_PRECISION slip[3],u[3],sm[3][3],gam[VE_MAX_NS],wfac;
    COMP_PRECISION *srow,smp[VE_MAX_NS],amp[VE_MAX_NP+1],res,scl,elastic_dev,val;
    PetscScalar *crow;
    PetscInt *colidx;
    if((np < 2)||(np > VE_MAX_NP))
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"rsf_ve_setup: -ve_np %ld out of [2, %d]",(long)np,VE_MAX_NP);
    if(plate_h <= 0.0)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"rsf_ve_setup: -ve_mode 2 requires -ve_plate_h [m] > 0");
    for(i=0;i < medium->nrflt;i++)
      if(-(fault[i].x[INT_Y] - fault[i].l) > plate_h + 1e-3)
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
		"rsf_ve_setup: patch %ld reaches depth %g m below the plate (H = %g m)",
		(long)i,(double)(-(fault[i].x[INT_Y]-fault[i].l)),(double)plate_h);
    ve_layered_gamma_pars(medium->elastic.shear,g2fac*medium->elastic.shear,
			  tm_yr*SEC_PER_YEAR,&gamma0,&rate_b);
    /* spec: ladder [b/45, 1.5 b], samples [b/100, 5 b], as validated
       in the ve_sp_cycle testbed */
    spec.g0 = medium->elastic.shear;spec.nu0 = medium->elastic.poisson;
    spec.t_M = tm_yr*SEC_PER_YEAR;
    spec.bulk0 = bulk_mod_from_G_nu(spec.g0,spec.nu0);
    spec.np = np;
    for(i=0;i < np;i++)
      spec.tau[i] = 1.0/(rate_b*(1.0/45.0)*
			 pow(67.5,(COMP_PRECISION)i/((COMP_PRECISION)(np-1))));
    spec.has_const = TRUE;	/* fitted, then discarded (implicit) */
    spec.nterm = np + 1;
    spec.ns = spec.nterm + VE_NHELD;
    for(i=0;i < spec.nterm;i++)
      spec.sk[i] = rate_b*(1.0/100.0)*
	pow(500.0,(COMP_PRECISION)i/((COMP_PRECISION)spec.nterm - 1.0));
    spec.sk[spec.nterm]   = rate_b/30.0;
    spec.sk[spec.nterm+1] = rate_b*0.6;
    ve_solve_weights(&spec);
    for(k=0;k < spec.ns;k++)
      gam[k] = ve_layered_gamma_of_s(gamma0,rate_b,spec.sk[k]);
    rsf->ve_np = np;
    for(p=0;p < np;p++){
      rsf->ve_tau[p] = spec.tau[p];
      PetscCall(MatCreate(PETSC_COMM_WORLD,&rsf->ve_C[p]));
      PetscCall(MatSetSizes(rsf->ve_C[p],medium->rn,PETSC_DECIDE,
			    medium->nrflt,medium->nrflt));
      PetscCall(MatSetType(rsf->ve_C[p],MATDENSE));
      PetscCall(MatSetUp(rsf->ve_C[p]));
    }
    /* sample rows for the owned receivers: srow[k*n + j] holds
       sample k of the kernel row; the image families are walked once
       per (i,j) and accumulated into all samples */
    srow = (COMP_PRECISION *)malloc((size_t)spec.ns*medium->nrflt*sizeof(COMP_PRECISION));
    colidx = (PetscInt *)malloc((size_t)medium->nrflt*sizeof(PetscInt));
    crow = (PetscScalar *)malloc((size_t)medium->nrflt*sizeof(PetscScalar));
    if((!srow)||(!colidx)||(!crow))MEMERROR("rsf_ve_setup");
    for(j=0;j < medium->nrflt;j++)
      colidx[j] = j;
    slip[STRIKE] = 1.0;slip[DIP] = slip[NORMAL] = 0.0;
    worst_res = elastic_dev = 0.0;
    HEADNODE
      fprintf(stderr,"rsf_ve_setup: mode 2 (plate over Maxwell, antiplane): H %g m g2/g1 %g t_M %g yr, np %ld, n_img %ld, sampling...\n",
	      (double)plate_h,(double)g2fac,(double)tm_yr,(long)np,(long)n_img);
    for(i=medium->rs;i < medium->re;i++){
      int iret;
      memset(srow,0,(size_t)spec.ns*medium->nrflt*sizeof(COMP_PRECISION));
      for(j=0;j < medium->nrflt;j++){
	for(n=0;n <= n_img;n++){
	  /* image family n: plain shift and (n >= 1) the surface
	     mirror, each with its own auto free-surface image via
	     the half-plane kernel; see ve_rybicki_antiplane_uz */
	  for(jn=0;jn < ((n == 0)?(1):(2));jn++){
	    if(n == 0){
	      eval_green_at_receiver(fault,i,j,slip,u,sm,&iret,
				     GC_STRESS_ONLY,TRUE,FALSE,
				     medium->elastic);
	    }else{
	      *img = fault[j];
	      img->x[INT_Y] = ((jn)?(-fault[j].x[INT_Y]):(fault[j].x[INT_Y]))
		- 2.0*(COMP_PRECISION)n*plate_h;
	      eval_green(fault[i].x,img,slip,u,sm,&iret,
			 GC_STRESS_ONLY,FALSE,FALSE,medium->elastic);
	    }
	    val = resolve_stress_on_fault(sm,(fault+i),STRIKE);
	    if(n == 0){
	      for(k=0;k < spec.ns;k++)
		srow[k*medium->nrflt+j] += val;
	    }else{
	      for(k=0;k < spec.ns;k++){
		wfac = pow(gam[k],(COMP_PRECISION)n);
		srow[k*medium->nrflt+j] += wfac*val;
	      }
	    }
	  }
	  /* all samples converged geometrically? largest |gam| decides */
	  if(pow(fabs(gam[0]) > fabs(gam[spec.ns-1])?fabs(gam[0]):fabs(gam[spec.ns-1]),
		 (COMP_PRECISION)(n+1)) < 1e-16)
	    break;
	}
      }
      /* per-pair fit and fill of the np amplitude rows */
      for(p=0;p < np;p++){
	for(j=0;j < medium->nrflt;j++){
	  for(k=0;k < spec.ns;k++)
	    smp[k] = srow[k*medium->nrflt+j];
	  for(it=0;it < spec.nterm;it++){
	    amp[it] = 0.0;
	    for(k=0;k < spec.nterm;k++)
	      amp[it] += spec.W[it][k]*smp[k];
	  }
	  if(p == 0){		/* held-out and elastic checks, once per pair */
	    scl = 0.0;
	    for(k=0;k < spec.ns;k++)
	      if(fabs(smp[k]) > scl)scl = fabs(smp[k]);
	    if(scl > 0.0){
	      for(k=spec.nterm;k < spec.ns;k++){
		res = 0.0;
		for(it=0;it < spec.nterm;it++)
		  res += amp[it]*ve_basis(&spec,it,spec.sk[k]);
		res = fabs(res - smp[k])/scl;
		if(res > worst_res)worst_res = res;
	      }
	    }
	  }
	  crow[j] = (PetscScalar)(amp[p]*scale);
	}
	PetscCall(MatSetValues(rsf->ve_C[p],1,&i,medium->nrflt,colidx,crow,INSERT_VALUES));
      }
    }
    free(srow);free(colidx);free(crow);
    for(p=0;p < np;p++){
      PetscCall(MatAssemblyBegin(rsf->ve_C[p],MAT_FINAL_ASSEMBLY));
      PetscCall(MatAssemblyEnd(rsf->ve_C[p],MAT_FINAL_ASSEMBLY));
    }
    /* elastic-consistency: || (sum_p C_p) x ||: the exponential sum
       may deviate from Is by the implicit relaxed operator C_inf,
       which for a fault NOT reaching through the whole plate is
       genuinely nonzero; report the norm ratio for information */
    {
      Vec xr,y1,y2;
      PetscReal n1,n2;
      PetscCall(MatCreateVecs(medium->Is,&xr,&y1));
      PetscCall(VecDuplicate(y1,&y2));
      PetscCall(VecSet(xr,1.0));
      PetscCall(MatMult(medium->Is,xr,y1));
      PetscCall(VecSet(y2,0.0));
      for(p=0;p < np;p++){
	PetscCall(MatMultAdd(rsf->ve_C[p],xr,y2,y2));
      }
      PetscCall(VecNorm(y1,NORM_2,&n1));
      PetscCall(VecAXPY(y2,-1.0,y1)); /* -(implicit C_inf) x */
      PetscCall(VecNorm(y2,NORM_2,&n2));
      HEADNODE
	fprintf(stderr,"rsf_ve_setup: worst held-out residual %.2e (tol %.1e); ||C_inf 1||/||K^el 1|| = %.3e (relaxed operator, uniform probe)\n",
		(double)worst_res,(double)fit_tol,(double)(n2/((n1 > 0.0)?n1:1.0)));
      PetscCall(VecDestroy(&xr));PetscCall(VecDestroy(&y1));PetscCall(VecDestroy(&y2));
    }
    if(worst_res > fit_tol)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,
	      "rsf_ve_setup: held-out residual %.2e above -ve_fit_tol %.2e",
	      (double)worst_res,(double)fit_tol);
  }else{
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"rsf_ve_setup: -ve_mode %ld undefined (1 or 2)",(long)mode);
  }
  /* state vectors, Is row layout */
  for(p=0;p < rsf->ve_np;p++){
    PetscCall(MatCreateVecs(medium->Is,NULL,&rsf->ve_h[p]));
    PetscCall(VecSet(rsf->ve_h[p],0.0));
  }
  PetscCall(MatCreateVecs(medium->Is,&rsf->ve_vrel,&rsf->ve_work));
  PetscCall(VecDuplicate(rsf->ve_vrel,&rsf->ve_slip_prev));
  PetscCall(VecSet(rsf->ve_slip_prev,0.0));
  PetscCall(VecDuplicate(rsf->ve_vrel,&rsf->ve_negvpl));
  PetscCall(VecCopy(negvpl_in,rsf->ve_negvpl));
  /*
     initial memory states (-ve_h_init, default 1): 1 = the backslip
     steady state h_p = C_p (-vpl) tau_p, i.e. the medium has been
     loaded by steady plate motion since forever, which removes the
     secular VE spin-up (time scale 3x the slowest ladder pole, about
     3 x 90 t_M for the layered default) from cycle studies; 0 = a
     virgin (unrelaxed) medium with zero memory, appropriate for
     relaxation-from-scratch tests such as the locked-fault
     hereditary-loading check.  NOTE for -ve_mode 1 (uniform Maxwell):
     the spun state has ZERO net fault loading (the memory sink
     exactly cancels the backslip products, the C_inf = 0
     saturation), so steady-state uniform-Maxwell backslip cycles do
     not exist; the events seen with -ve_h_init 0 are the (possibly
     very long lived, ~t_M) transient of the virgin start
  */
  {
    PetscInt hinit = 1;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-ve_h_init",&hinit,NULL));
    if(hinit == 1){
      for(p=0;p < rsf->ve_np;p++){
	PetscCall(MatMult(rsf->ve_C[p],rsf->ve_negvpl,rsf->ve_h[p]));
	PetscCall(VecScale(rsf->ve_h[p],rsf->ve_tau[p]));
      }
      HEADNODE
	fprintf(stderr,"rsf_ve_setup: memory states initialized at the backslip steady state (-ve_h_init 1)\n");
    }else{
      HEADNODE
	fprintf(stderr,"rsf_ve_setup: memory states initialized to zero (virgin medium, -ve_h_init 0)\n");
    }
  }
  HEADNODE{
    fprintf(stderr,"rsf_ve_setup: VE ON, np %ld, rates ladder [yr]:",(long)rsf->ve_np);
    for(p=0;p < rsf->ve_np;p++)
      fprintf(stderr," %.3g",(double)(rsf->ve_tau[p]/SEC_PER_YEAR));
    fprintf(stderr,"; h: exact-exponential, %s sink forcing\n",
	    (rsf->ve_hstage)?("stage-consistent"):("frozen (lagged)"));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*
   apply the memory sink to tau_dot at stage time `time` (no-op when VE
   is off); called from rsf_compute_vel_and_stressing with the stage
   state X.

   Default (-ve_h_stage 1), stage-consistent forcing: the sink is the
   exact-exponential h propagated from the last accepted step USING THE
   STAGE STATE'S OWN SLIP INCREMENT as the within-step forcing,

     tau_dot -= sum_p [ h_p(t0) e^(-dt/tau_p) / tau_p
                        + C_p qbar (1 - e^(-dt/tau_p)) ],
     qbar = (slip(t) - slip(t0))/dt - vpl,   dt = time - t0,

   which is exact for velocity constant within the step and, because
   qbar depends on the stage state, is seen by the TS error controller
   (the frozen variant is not, and its lagged forcing accumulates a
   first-order bias during rapidly varying afterslip that can shift
   recurrence by tens of percent in relaxation-critical problems, e.g.
   the uniform-Maxwell patch of Miyake and Noda 2019; the bias was
   found by cross-checking against an exact-Erlang-chain integrator).
   Cost: np matvecs per RHS evaluation.

   -ve_h_stage 0 restores the frozen decay-only sink
   tau_dot -= sum_p h_p(t0) e^(-dt/tau_p)/tau_p (np matvecs per
   ACCEPTED STEP only).
*/
PetscErrorCode rsf_ve_apply_sink(struct interact_ctx *par, Vec X, PetscReal time)
{
  struct med *medium = par->medium;
  struct rsf_vars *rsf = medium->rsf;
  const PetscScalar *x;
  const PetscScalar *sp;
  PetscScalar *vr;
  PetscReal dt,efac;
  PetscInt p,i,j,k;
  PetscFunctionBeginUser;
  if(rsf->ve_np <= 0)
    PetscFunctionReturn(PETSC_SUCCESS);
  dt = time - rsf->ve_t0;
  if(dt < 0.0)dt = 0.0;
  if(rsf->ve_hstage && (dt > 1e-3)){
    /* stage-average relative velocity from the stage state's slip */
    PetscCall(VecGetArrayRead(X,&x));
    PetscCall(VecGetArrayRead(rsf->ve_slip_prev,&sp));
    PetscCall(VecGetArray(rsf->ve_vrel,&vr));
    for(i=medium->rs,j=0,k=0;i < medium->re;i++,j+=rsf->dim,k++)
      vr[k] = (x[j+3] - sp[k])/dt;
    PetscCall(VecRestoreArray(rsf->ve_vrel,&vr));
    PetscCall(VecRestoreArrayRead(rsf->ve_slip_prev,&sp));
    PetscCall(VecRestoreArrayRead(X,&x));
    PetscCall(VecAXPY(rsf->ve_vrel,1.0,rsf->ve_negvpl));
    for(p=0;p < rsf->ve_np;p++){
      efac = PetscExpReal(-dt/rsf->ve_tau[p]);
      PetscCall(VecAXPY(rsf->tau_dot,-efac/rsf->ve_tau[p],rsf->ve_h[p]));
      PetscCall(MatMult(rsf->ve_C[p],rsf->ve_vrel,rsf->ve_work));
      PetscCall(VecAXPY(rsf->tau_dot,-(1.0-efac),rsf->ve_work));
    }
  }else{
    for(p=0;p < rsf->ve_np;p++){
      efac = PetscExpReal(-dt/rsf->ve_tau[p]);
      PetscCall(VecAXPY(rsf->tau_dot,-efac/rsf->ve_tau[p],rsf->ve_h[p]));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*
   accepted-step h update (TS monitor): exact-exponential recurrence
   with the step-average relative velocity from the slip increment
*/
PetscErrorCode rsf_ve_monitor(TS ts, PetscInt step, PetscReal time, Vec X, void *ptr)
{
  struct interact_ctx *par = (struct interact_ctx *)ptr;
  struct med *medium = par->medium;
  struct rsf_vars *rsf = medium->rsf;
  const PetscScalar *x;
  PetscScalar *vr,*sp;
  PetscReal dt,efac;
  PetscInt i,j,k,p;
  PetscFunctionBeginUser;
  if(rsf->ve_np <= 0)
    PetscFunctionReturn(PETSC_SUCCESS);
  dt = time - rsf->ve_t0;
  if(dt <= 0.0)
    PetscFunctionReturn(PETSC_SUCCESS);
  /* vbar - vpl = (slip - slip_prev)/dt + negvpl */
  PetscCall(VecGetArrayRead(X,&x));
  PetscCall(VecGetArray(rsf->ve_vrel,&vr));
  PetscCall(VecGetArray(rsf->ve_slip_prev,&sp));
  for(i=medium->rs,j=0,k=0;i < medium->re;i++,j+=rsf->dim,k++){
    vr[k] = (x[j+3] - sp[k])/dt;
    sp[k] = x[j+3];
  }
  PetscCall(VecRestoreArray(rsf->ve_slip_prev,&sp));
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(rsf->ve_vrel,&vr));
  PetscCall(VecAXPY(rsf->ve_vrel,1.0,rsf->ve_negvpl));
  for(p=0;p < rsf->ve_np;p++){
    efac = PetscExpReal(-dt/rsf->ve_tau[p]);
    PetscCall(MatMult(rsf->ve_C[p],rsf->ve_vrel,rsf->ve_work));
    PetscCall(VecScale(rsf->ve_h[p],efac));
    PetscCall(VecAXPY(rsf->ve_h[p],rsf->ve_tau[p]*(1.0-efac),rsf->ve_work));
  }
  rsf->ve_t0 = time;
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
