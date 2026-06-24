#include "interact.h"
#include "properties.h"
/*
  
  quasi-dynamic rate-and-state earthquake sequence solver using PETSc
  TS time steppers, set up to solve the same ODE system as HBI (Ozawa
  et al., 2023, https://github.com/sozawa94/hbi, main_LH.f90) in the
  aging-law, no-fluid, backslip-loading case:

  state vector per patch:  (psi, tau, sigma, slip), where

  v         = 2 v0 exp(-psi/a) sinh(tau/(sigma a))     (regularized RSF)
  dpsi/dt   = (b/dc) (v0 exp((f0-psi)/b) - |v|)        (aging law, as in HBI)
  dsigma/dt = [K_n (v - vpl)]_i                        (compression positive)
  dtau/dt   = ([K_s (v - vpl)]_i 
               - eta (dv/dpsi dpsi/dt + dv/dsigma dsigma/dt))
              / (1 + eta dv/dtau)                       (quasi-dynamic,
                                                        eta = G/(2 c_s))
  dslip/dt  = v

  unlike HBI, slip is part of the ODE state (HBI integrates slip
  outside the RK step with a trapezoidal update); slip is effectively
  excluded from the step-size control via a loose absolute tolerance

  time integration: adaptive embedded Runge-Kutta (default
  Dormand-Prince 5(4), HBI uses Cash-Karp 5(4)) with a pure-relative
  infinity error norm and rtol = eps_r as in HBI's rkqs, step growth
  limited to 2x, shrink to 0.5x per rejection

  units are SI: Pa, m, s; geometry input (geom.in) must hence be in meters

*/
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#include "rsf.h"
/* single definition of the domain-check context pointer declared in rsf.h
   (the callback registered with TSSetFunctionDomainError has no user pointer) */
struct interact_ctx *rsf_par_static = NULL;
#endif
/*
   thin driver: gather settings, then run the solver; the two phases
   live in rsf_get_settings and rsf_solve_run below
*/
int main(int argc,char **argv)
{
#ifdef USE_PETSC
  struct interact_ctx par[1];	/* fault geometry and medium */
  struct med *medium;
  struct rsf_solve_settings set[1];
  FILE *tst;
  char *home_dir = getenv("HOME");
  char par_file[STRLEN];
  /* only read the YAML defaults file if it exists */
  snprintf(par_file,STRLEN,"%s/progs/src/interact/petsc_settings.yaml",(home_dir)?(home_dir):("."));
  tst = fopen(par_file,"r");
  if(tst){
    HEADNODE
      fprintf(stderr,"%s: found and using Petsc options in %s\n",argv[0],par_file);
    fclose(tst);
  }else{
    par_file[0]='\0';
  }
  PetscInitialize(&argc,&argv,(par_file[0])?(par_file):(NULL),NULL); /* initialize from setting file and or command line */
  /* holds the parameets */
  par->medium=(struct med *)calloc(1,sizeof(struct med)); /* init as all zeros */
  medium = par->medium;
  
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  /* 
     all defaults and command-line/-YAML overrides 
  */
  PetscCall(rsf_get_settings(argc,argv,par,set));
  /* 
     geometry, interaction matrices, initial state, time integration, output 
  */
  PetscCall(rsf_solve_run(argc,argv,par,set));
  PetscCall(PetscFinalize());
#endif
  exit(0);
}
#ifdef USE_PETSC

/*
   the solver proper: read geometry, build the interaction matrices,
   set the initial state, configure the TS integrator and the output,
   integrate to the stop time, and clean up.  All configuration arrives
   through struct rsf_solve_settings; the body below is the original driver
   logic with the settings copied back into locals
*/
PetscErrorCode rsf_solve_run(int argc,char **argv,struct interact_ctx *par,
			     struct rsf_solve_settings *set)
{
  TS ts;			/* timestepping context */
  TSAdapt adapt;
  PetscReal *values = NULL;
  const PetscReal sec_per_year = 365.25*24.*60.*60.;
  VecScatter ctx;
  Vec xout,x,vatol,islip_rate_vec,stress_rate;
  struct med *medium = par->medium;
  struct rsf_vars *rsf = medium->rsf;
  struct flt *fault;
  PetscInt m,i,n,j,use_hmatrix = medium->use_hmatrix;
  PetscRandom prand;
  PetscReal state,patch_l,stable_l,rand_fac,dummy[3];
  PetscInt event_direction[1] = {0};		/* detect both crossing directions */
  PetscBool event_terminate[1] = {PETSC_FALSE};
  PetscLogDouble tb0,tb1;
  struct rsf_out_ctx uc[1];
  FILE *iin;
  PetscInt ii;
  PetscReal *ic_tau=NULL,*ic_vel=NULL;
  PetscBool warned = PETSC_FALSE,flg;
  /* per-group field output (rank 0); built and owned here, freed in finalize */
  struct rsf_group_grid *rsf_groups=NULL;
  int rsf_ngroup=0;
  double *rsf_vbuf=NULL;
  /* settings copied back into locals so the solver body reads unchanged */
  PetscReal shear_modulus_si = set->shear_modulus_si, s_wave_speed_si = set->s_wave_speed_si;
  PetscReal sigma_init = set->sigma_init, tau_init = set->tau_init, vel_init = set->vel_init;
  PetscReal rand_amp = set->rand_amp, rtol = set->rtol, atol_slip = set->atol_slip;
  PetscReal dt_init = set->dt_init, dt_max = set->dt_max;
  PetscReal dt_monitor = set->dt_monitor, rdx_monitor = set->rdx_monitor;
  PetscReal adx_monitor = set->adx_monitor, monitor_tmin = set->monitor_tmin;
  PetscReal vel_event = set->vel_event, vel_event_hyst = set->vel_event_hyst, event_tmin = set->event_tmin;
  PetscBool track_events = set->track_events, have_ic = set->have_ic, have_dc = set->have_dc;
  
  char geom_file[STRLEN],rsf_file[STRLEN],rsf_ic_file[STRLEN],rsf_dc_file[STRLEN];
  PetscFunctionBeginUser;
  strncpy(geom_file,  set->geom_file,  STRLEN);
  strncpy(rsf_file,   set->rsf_file,   STRLEN);
  strncpy(rsf_ic_file,set->rsf_ic_file,STRLEN);
  strncpy(rsf_dc_file,set->rsf_dc_file,STRLEN);
  /* get the geometry */
  read_geometry(geom_file,&medium,&par->fault,TRUE,FALSE,FALSE,FALSE);
  fault = par->fault;
  n = medium->nrflt;
  HEADNODE{
    fprintf(stderr,"%s: read geometry from %s, %i patches, on %i cores\n",
	    argv[0],geom_file,n,medium->comm_size);
    fprintf(stderr,"%s: G %.6e Pa cs %g m/s eta %.6e f0 %g dc %g m v0 %.3e vpl %.3e m/s\n",
	    argv[0],shear_modulus_si,s_wave_speed_si,rsf->shear_mod_over_2cs_si,
	    rsf->f0,rsf->dc,rsf->v0,rsf->vpl);
    fprintf(stderr,"%s: sigma0 %.6e tau0 %.6e Pa vinit %.3e m/s rand_amp %g\n",
	    argv[0],sigma_init,tau_init,vel_init,rand_amp);
    fprintf(stderr,"%s: rtol %.1e dt_init %g s dt_max %g s stop %g yr\n",
	    argv[0],rtol,dt_init,dt_max,medium->stop_time/sec_per_year);
  }
  /* 
     now, read in a,b variations (stored in fault[].mu_s, fault[].mu_d)
  */
  read_rsf(rsf_file,medium,fault);
  /* optional per-cell initial tau,vel (e.g. BP5 nucleation patch) */
  if(have_ic){
    ic_tau=(PetscReal *)malloc((size_t)n*sizeof(PetscReal));
    ic_vel=(PetscReal *)malloc((size_t)n*sizeof(PetscReal));
    if((!ic_tau)||(!ic_vel)){
      fprintf(stderr,"%s: per-cell IC alloc failed\n",argv[0]);
      exit(-1);
    }
    iin = myopen(rsf_ic_file,"r");
    for(ii=0;ii < n;ii++)
      if(fscanf(iin,"%lf %lf",(ic_tau+ii),(ic_vel+ii))!=2){
	fprintf(stderr,"%s: error reading tau vel for patch %i from %s\n",
		argv[0],ii,rsf_ic_file);
	exit(-1);
      }
    fclose(iin);
    HEADNODE
      fprintf(stderr,"%s: read per-cell initial tau,vel from %s\n",argv[0],rsf_ic_file);
  }
  /* optional per-cell D_c (e.g. reduced D_RS in the BP5 nucleation patch) */
  if(have_dc){
    rsf->dc_vec=(PetscReal *)malloc((size_t)n*sizeof(PetscReal));
    if(!rsf->dc_vec){
      fprintf(stderr,"%s: per-cell dc alloc failed\n",argv[0]);
      exit(-1);
    }
    iin=myopen(rsf_dc_file,"r");
    for(ii=0;ii < n;ii++)
      if(fscanf(iin,"%lf",(rsf->dc_vec+ii))!=1){
	fprintf(stderr,"%s: error reading dc for patch %i from %s\n",
		argv[0],ii,rsf_dc_file);
	exit(-1);
      }
    fclose(iin);
    HEADNODE
      fprintf(stderr,"%s: read per-cell D_c from %s\n",argv[0],rsf_dc_file);
  }
  /* 
     create and calculate interaction matrices, scaled from interact's
     internal shear modulus to SI

     NOTE on conventions: interact uses the physics (extension
     positive) convention for stress.  the RSF formulation, like HBI,
     needs sigma compression positive, hence the NEGATIVE scale factor
     for the normal stress interaction matrix.  the shear (strike)
     convention is consistent as is: positive slip reduces the shear
     traction on the slipping patch (negative diagonal), so backslip
     -vpl produces positive loading
  */
  /* instrument the interaction-matrix build so its wall time can be
     compared directly against the MatAssemblyEnd event in -log_view.
     the two should agree; if the build wall is small but -log_view
     still charges a large MatAssemblyEnd, the cost is deferred. */
  
  PetscCall(PetscBarrier(NULL));
  PetscCall(PetscTime(&tb0));
  calc_petsc_Isn_matrices(medium,fault,use_hmatrix,
			  shear_modulus_si/SHEAR_MODULUS,0,&medium->Is,medium->Is_hctx); /* shear stress */
  PetscCall(PetscBarrier(NULL));
  PetscCall(PetscTime(&tb1));
  HEADNODE
    fprintf(stderr,"rsf_solve: Is (shear) interaction matrix build wall = %12.4f s\n",(double)(tb1-tb0));
  if(rsf->calc_sigma_dot){
    PetscCall(PetscBarrier(NULL));
    PetscCall(PetscTime(&tb0));
    calc_petsc_Isn_matrices(medium,fault,use_hmatrix,
			    -shear_modulus_si/SHEAR_MODULUS,1,&medium->In,medium->In_hctx); /* normal stress, compression positive */
    PetscCall(PetscBarrier(NULL));
    PetscCall(PetscTime(&tb1));
    HEADNODE
      fprintf(stderr,"rsf_solve: In (normal) interaction matrix build wall = %12.4f s\n",(double)(tb1-tb0));
  }

  if(use_hmatrix)		/* this should only print simple info,
				   not full matrix for H matrix */
    PetscCall(MatView(medium->Is,PETSC_VIEWER_STDOUT_WORLD));
  /* 
     preallocate the RHS work vectors with the matrix row layout
  */
  PetscCall(MatCreateVecs(medium->Is,&rsf->vel,&rsf->tau_dot));
  if(rsf->calc_sigma_dot)
    PetscCall(VecDuplicate(rsf->tau_dot,&rsf->sigma_dot));
  /*
     compute backslip stressing rates sinc[0] (shear) and sinc[1]
     (normal, compression positive) from slip rate -vpl, made
     available on all nodes
  */
  PetscCall(MatCreateVecs(medium->Is,&islip_rate_vec,&stress_rate));
  PetscCall(VecSet(islip_rate_vec,-rsf->vpl));
  for(i=0;i < ((rsf->calc_sigma_dot)?(2):(1));i++){
    if(i==0)
      PetscCall(MatMult(medium->Is,islip_rate_vec,stress_rate)); /* shear */
    else
      PetscCall(MatMult(medium->In,islip_rate_vec,stress_rate)); /* normal */
    /* assign to the faults for every node */
    PetscCall(VecScatterCreateToAll(stress_rate,&ctx,&xout));
    PetscCall(VecScatterBegin(ctx,stress_rate,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(  ctx,stress_rate,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGetArray(xout,&values));
    for(j=0;j < medium->nrflt;j++) /* store the backslip loading rates */
      fault[j].sinc[i] = values[j];
    PetscCall(VecRestoreArray(xout,&values));
    PetscCall(VecScatterDestroy(&ctx));
    PetscCall(VecDestroy(&xout));
  }
  PetscCall(VecDestroy(&stress_rate));
  PetscCall(VecDestroy(&islip_rate_vec));

  /* 
     solution vector: local size has to be rsf->dim times the local
     matrix row range so that local patch i of the matrix layout owns
     entries [i*rsf->dim, (i+1)*rsf->dim) of x on the same rank
  */
  m=rsf->dim*n;
  PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
  PetscCall(VecSetSizes(x,rsf->dim*medium->rn,m));
  PetscCall(VecSetFromOptions(x));
  /* 
     initialize x vector: uniform tau,sigma; state from inverting the
     velocity relation at vel_init, as in HBI initcond():
     psi = a log(2 v0/v sinh(tau/(sigma a)))
  */
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&prand));
  PetscCall(PetscRandomSetType(prand,PETSCRAND48));
  if(rand_amp > 0)
    PetscCall(PetscRandomSetInterval(prand,1.0-rand_amp,1.0+rand_amp)); 
  for (i = medium->rs,j=i*rsf->dim; i < medium->re; i++,j+=rsf->dim){
    PetscCall(PetscRandomGetValue(prand,&rand_fac));
    fault[i].s[NORMAL] = sigma_init;	/* compression positive */
    if(have_ic){
      fault[i].s[STRIKE] = ic_tau[i];	/* per-cell initial shear stress [Pa] */
      fault[i].u[0]      = ic_vel[i];	/* per-cell initial slip velocity [m/s] */
    }else{
      fault[i].s[STRIKE] = tau_init;
      fault[i].u[0] = vel_init;
    }
    /* state consistent with v = vel_init */
    state = fault[i].mu_s*log(2.0*rsf->v0/fault[i].u[0]*sinh(fault[i].s[STRIKE]/fault[i].s[NORMAL]/fault[i].mu_s));
    if(rand_amp > 0)
      state *= rand_fac;
    PetscCall(VecSetValue(x, (PetscInt)(j+0),state, INSERT_VALUES));
    PetscCall(VecSetValue(x, (PetscInt)(j+1),fault[i].s[STRIKE], INSERT_VALUES)); /* shear stress */
    PetscCall(VecSetValue(x, (PetscInt)(j+2),fault[i].s[NORMAL], INSERT_VALUES)); /* normal stress */
    PetscCall(VecSetValue(x, (PetscInt)(j+3),0.0, INSERT_VALUES)); /* total slip */
    /* check discretization against the quasi-static cohesive zone
       size Lb = G dc/(b sigma), using the smallest patch dimension */
    patch_l = 2.0*((fault[i].l < fault[i].w)?(fault[i].l):(fault[i].w));
    stable_l = shear_modulus_si*(rsf->dc_vec?rsf->dc_vec[i]:rsf->dc)/fault[i].mu_d/fault[i].s[NORMAL];
    if((patch_l > stable_l/3.)&&(!warned)){
      fprintf(stderr,"%s: WARNING: patch %05i: size %.3e m vs. Lb = G dc/(b sigma) %.3e m, coarser than Lb/3\n",
	      argv[0],i,patch_l,stable_l);
      fprintf(stderr,"%s: a %g b %g tau %.4e sigma %.4e psi %.4e v %.4e (suppressing further warnings)\n",
	      argv[0],fault[i].mu_s,fault[i].mu_d,fault[i].s[STRIKE],fault[i].s[NORMAL],state,
	      vel_from_rsf(fault[i].s[STRIKE], fault[i].s[NORMAL], state,
			   fault[i].mu_s,rsf->v0,dummy,(dummy+1),(dummy+2),medium));
      warned = PETSC_TRUE;
    }
  }
  PetscCall(VecAssemblyBegin(x));PetscCall(VecAssemblyEnd(x));
  PetscCall(PetscRandomDestroy(&prand));
  /* 
     absolute tolerance vector: essentially zero for psi,tau,sigma
     (pure relative control as in HBI's |yerr/y| norm), loose for slip
     to exclude it from step size control
  */
  PetscCall(VecDuplicate(x,&vatol));
  for (i = medium->rs,j=i*rsf->dim; i < medium->re; i++,j+=rsf->dim){
    PetscCall(VecSetValue(vatol,(PetscInt)(j+0),1e-12, INSERT_VALUES)); /* psi */
    PetscCall(VecSetValue(vatol,(PetscInt)(j+1),1e-3,  INSERT_VALUES)); /* tau [Pa] */
    PetscCall(VecSetValue(vatol,(PetscInt)(j+2),1e-3,  INSERT_VALUES)); /* sigma [Pa] */
    PetscCall(VecSetValue(vatol,(PetscInt)(j+3),atol_slip,INSERT_VALUES)); /* slip [m] */
  }
  PetscCall(VecAssemblyBegin(vatol));PetscCall(VecAssemblyEnd(vatol));
  /*
    Create timestepper context
  */
  PetscCall(TSCreate(PETSC_COMM_WORLD,&ts));
  PetscCall(TSSetProblemType(ts,TS_NONLINEAR));
  PetscCall(TSSetSolution(ts,x));
  PetscCall(TSSetRHSFunction(ts, NULL, rsf_ODE_RHSFunction,par));
  /*
    adaptive embedded Runge-Kutta; HBI uses Cash-Karp RK5(4), which
    PETSc does not provide, Dormand-Prince 5(4) is the equivalent
    default here (override with -ts_rk_type)
  */
  PetscCall(TSSetType(ts,TSRK));
  PetscCall(TSRKSetType(ts,TSRK5DP));
  /*
     NOTE on the integrator order: after the step-size controller (see the
     longer note further down), the RK order is the other large production
     lever.  5DP is kept as the default because its 5(4) order matches HBI's
     Cash-Karp for the cross-code comparison.  A lower-order embedded pair
     does fewer stage matvecs per accepted step -- 3BS (Bogacki-Shampine
     RK3(2)) costs 3 vs 6 for 5DP -- and, although it takes more (and more
     frequently rejected) steps, in our tests it still nets FEWER matvecs.
     Judged by the physically meaningful metric, the event RECURRENCE
     interval (which is far better converged than the absolute event phase:
     the phase drifts by ~rtol run-to-run, the interval does not),
       -ts_rk_type 3bs
     reproduced the 5DP recurrence to within ~0.005-0.02 yr (<~0.01% of the
     ~230 yr interval) at every resolution we tried -- BP5 2 km dense, 1 km
     dense AND HACApK, 0.5 km HACApK -- while reducing matvecs by very roughly
     24% (2 km), 35% (1 km) and 41% (0.5 km).  The saving appears to grow with
     resolution because finer meshes take more steps, so the cheaper-per-step
     method compounds.  Guidance (all specific to these single-host serial
     tests and one event sequence, so confirm per case): keep 5DP for anything
     compared against HBI; consider -ts_rk_type 3bs for production where a
     recurrence error of order 0.02 yr is acceptable (it stacks with
     -ts_adapt_type dsp below).  Do NOT loosen rtol below ~1e-4: at 1e-3 the
     rejection rate climbed enough that the run was both less accurate and not
     faster, and over long runs could make very slow progress.  In our runs
     2a (RK2(1)) needed far more matvecs (tiny steps) and 5bs did not behave
     well on the stiff coseismic phase; higher order (8vr) only repaid its
     per-step cost for very tight ABSOLUTE event times, not long multi-cycle
     runs.  A different problem, mesh, or machine could of course rank these
     differently.
  */
  /*
    Set the initial time and the initial timestep given above.
  */
  PetscCall(TSSetTime(ts,medium->time));	   /* initial time, set above */
  PetscCall(TSSetTimeStep(ts,dt_init));		   /* initial timestep */
  PetscCall(TSSetMaxSteps(ts,(PetscInt)1000000000)); /* do not stop on step count */
  /* HBI's rkqs retries shrinking steps without limit; PETSc defaults
     to a max of 10 consecutive rejections which is far too few for
     the interseismic -> coseismic dt collapse */
  PetscCall(TSSetMaxStepRejections(ts,-1));
  PetscCall(TSSetMaxSNESFailures(ts,-1)); /* domain errors count as solver failures */
  /* output times should not constrain the natural step size selection,
     allow rough final time as in ode_solve_test.c */
  PetscCall(TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER));
  /* HBI rkqs uses a pure relative infinity norm error measure */
  PetscCall(TSSetTolerances(ts,0.0,vatol,rtol,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-ts_adapt_wnormtype",&flg));
  if(!flg)
    PetscCall(PetscOptionsSetValue(NULL,"-ts_adapt_wnormtype","infinity"));
  /* HBI rkqs: step growth limited to 2x, shrink to >= 0.5x */
  PetscCall(TSGetAdapt(ts,&adapt));
  PetscCall(TSAdaptSetClip(adapt,0.5,2.0));
  /* 
     NOTE on the step-size controller: the default PETSc RK controller
     ("basic") is used here because it mirrors HBI's rkqs (the clip and the
     infinity norm above), which matters for the cross-code comparison.  For
     PRODUCTION speed, the digital-signal-processing controller is often a
     bigger lever than anything in the RHS, because rejected steps each waste
     a full set of stage matvecs: on a BP5 2 km dense serial test,
       -ts_adapt_type dsp
     cut rejected steps by ~63% (188 -> 70) for the SAME accepted-step count,
     i.e. ~12% fewer matvecs and ~9% less time in TSStep, at equal tolerance.
     It does change the accepted step sequence, so event/recurrence times
     shift by ~rtol (e.g. first event 236.81 -> 236.85 yr at rtol=1e-4) and it
     no longer matches HBI's rkqs.  Recommendation: keep the default for
     anything compared against HBI; pass -ts_adapt_type dsp for production
     runs where reproducibility to << rtol is not required.  Numbers are
     specific to that test (resolution, tolerance, machine); confirm per case,
     since a different problem could trade rejections for more accepted steps.
  */
  /* cap the step at the monitor interval, cf. ode_solve_test.c, so
     the change-triggered monitor cannot be stepped over */
  PetscCall(TSAdaptSetStepLimits(adapt,0.0,(dt_max < dt_monitor)?(dt_max):(dt_monitor)));
  /* 
     reject steps whose stages leave the physical domain (overflow,
     non-positive normal stress); the equivalent of HBI rkqs' NaN/Inf
     guard, without this an overshooting stage can produce a NaN error
     estimate which the controller would silently accept
  */
  rsf_par_static = par;		/* need to assign this for the bounds
				   check */
  /*  */
  PetscCall(TSSetFunctionDomainError(ts,rsf_domain_check));
  PetscCall(TSSetFromOptions(ts));
  /*
     compact slip-rate field output: build the per-group grids and write
     each group's static geometry once (rank 0 only).  n fault groups
     give n separate grids, written to rsf_geom.gGGG.dat.  The group
     structures and scratch are handed to the monitor context and freed
     in the finalize routine
  */
  if(set->field_enable && (medium->comm_rank == 0)){
    int ig;
    char gfile[STRLEN];
    rsf_ngroup = rsf_build_groups(medium,fault,&rsf_groups);
    rsf_vbuf = (double *)malloc((size_t)n*sizeof(double));
    if((rsf_ngroup < 1)||(!rsf_groups)||(!rsf_vbuf)){
      fprintf(stderr,"%s: field-output setup failed, disabling field output\n",argv[0]);
      if(rsf_vbuf){free(rsf_vbuf);rsf_vbuf=NULL;}
      rsf_free_groups(rsf_groups,rsf_ngroup);rsf_groups=NULL;rsf_ngroup=0;
      set->field_enable = PETSC_FALSE;
    }else{
      for(ig=0;ig < rsf_ngroup;ig++){
	snprintf(gfile,STRLEN,"rsf_geom.g%03d.dat",rsf_groups[ig].id);
	rsf_write_group_geometry(rsf_groups+ig,fault,sigma_init,gfile);
      }
      fprintf(stderr,"%s: field output on, %i fault group(s), one frame every %i accepted steps\n",
	      argv[0],rsf_ngroup,set->field_step_interval);
    }
  }
  /* 
     set up the change/time triggered state monitor, the periodic
     field output, and the velocity threshold event tracker,
     cf. hmatrix_test/ode_solve_test.c
  */
  PetscCall(rsf_init_monitor_and_event(uc,par,dt_monitor,adx_monitor,rdx_monitor,
				       monitor_tmin,event_tmin,vel_event,vel_event_hyst,
				       track_events,medium->time,x,vel_init,
				       set->field_enable,set->field_step_interval,set->field_tmin,
				       rsf_groups,rsf_ngroup,rsf_vbuf));
  PetscCall(TSMonitorSet(ts,rsf_TS_Monitor,(void *)uc,NULL));
  if(track_events){
    PetscCall(TSSetEventHandler(ts,1,event_direction,event_terminate,
				rsf_event_function,rsf_post_event,(void *)uc));
    /* locate the crossing to within 0.01 in log10(v) */
    PetscCall(TSSetEventTolerances(ts,1e-2,NULL));
  }
  /* 
     single TSSolve to the stop time

     NOTE: do NOT chop the integration into output intervals with
     repeated TSSolve calls and TS_EXACTFINALTIME_INTERPOLATE - each
     restart resumes from the interpolated rather than the true
     trajectory state, which introduces a systematic, tolerance-
     independent drift (verified against a reference integration);
     all output is instead handled in the monitor at accepted steps
     and by the event handler
  */
  PetscCall(TSSetMaxTime(ts,medium->stop_time));
  PetscCall(TSSolve(ts, x));
  PetscCall(TSGetSolveTime(ts,&(medium->time)));
  PetscCall(rsf_finalize_monitor_and_event(uc));
  
  /*View information about the time-stepping method and the solution at the end time.*/
  PetscCall(TSView(ts, PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(MatDestroy(&medium->Is));
  
  PetscCall(VecDestroy(&rsf->vel));
  PetscCall(VecDestroy(&rsf->tau_dot));
  if(rsf->calc_sigma_dot){
    PetscCall(MatDestroy(&medium->In));
    PetscCall(VecDestroy(&rsf->sigma_dot));
  }
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&vatol));

  /*  */
  PetscCall(TSDestroy(&ts)); 
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* 
   event function, cf. myEventFunction in ode_solve_test.c: locate
   crossings of the maximum slip velocity through the vel_event
   threshold (in log space); evaluated collectively
*/
PetscErrorCode rsf_event_function(TS ts,PetscReal t,Vec X,PetscScalar *fvalue,void *ctx)
{
  const PetscScalar *x;
  struct rsf_out_ctx *uc;
  struct med *medium;struct flt *fault;struct rsf_vars *rsf;
  PetscInt i,j;
  PetscReal v,lvmax,gvmax,d1,d2,d3;
  PetscFunctionBeginUser;
  uc = (struct rsf_out_ctx *)ctx;
  medium = uc->par->medium;
  rsf = medium->rsf;
  fault = uc->par->fault;
  lvmax = 0.0;
  PetscCall(VecGetArrayRead(X,&x));
  for (i = medium->rs, j=0; i < medium->re; i++, j+=rsf->dim) {
    v = fabs(vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,rsf->v0,&d1,&d2,&d3,medium));
    if(v > lvmax)lvmax = v;
  }
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCallMPI(MPI_Allreduce(&lvmax,&gvmax,1,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD));
  /* hysteresis: once slipping, arrest only below vel_event * hyst,
     which debounces phantom crossings from dense-output oscillation
     through the near-singular coseismic acceleration */
  fvalue[0] = log10((gvmax > 1e-300)?(gvmax):(1e-300)) -
    log10(uc->vel_event * ((uc->slipping)?(uc->vel_event_hyst):(1.0)));
  PetscFunctionReturn(PETSC_SUCCESS);
}


#endif
