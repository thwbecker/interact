#include "interact.h"
#include "properties.h"
/*
  
  quasi-dynamic rate-and-state earthquake sequence solver using PETSc
  TS time steppers, set up to solve the same ODE system as HBI
  (Ozawa et al., 2023, https://github.com/sozawa94/hbi, main_LH.f90)
  in the aging-law, no-fluid, backslip-loading case:

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

  time integration: adaptive embedded Runge-Kutta (default Dormand-Prince
  5(4), HBI uses Cash-Karp 5(4)) with a pure-relative infinity error
  norm and rtol = eps_r as in HBI's rkqs, step growth limited to 2x,
  shrink to 0.5x per rejection

  units are SI: Pa, m, s; geometry input (geom.in) must hence be in meters

*/
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#endif

#define INT_RSF_DIM 4		/* psi,tau,sigma,u(slip) per patch */

#ifdef USE_PETSC
/* 
   parameters for optional normal stress limiter, cf. HBI limitsigma
   (file scope since they are only used by the driver and RHS here) 
*/
static PetscBool rsf_limit_sigma = PETSC_FALSE;
static PetscReal rsf_min_sigma = 1e6, rsf_max_sigma = 300e6; /* Pa */
/* monitor output */
static FILE *rsf_monitor_out = NULL;
static PetscReal rsf_monitor_last_time = -1;
/* context access for the domain check */
static struct interact_ctx *rsf_par_static = NULL;
#endif

int main(int argc,char **argv)
{
#ifdef USE_PETSC

  TS ts; /* timestepping context */
  TSAdapt adapt;
  PetscReal *values = NULL;
  /* material parameters */
  PetscReal shear_modulus_si = 32.04e9, s_wave_speed_si = 3.464e3;
  const PetscReal sec_per_year = 365.25*24.*60.*60.;
  /*  */
  VecScatter ctx;
  Vec xout,x,vatol,islip_rate_vec,stress_rate;
  FILE *fout1=NULL,*fout2;
  struct med *medium;struct flt *fault;
  PetscInt m,i,n,j,field_out,ierr,use_hmatrix=0;
  PetscRandom prand;
  PetscReal state,slip,sum[3],patch_l,stable_l,rand_fac,dummy[3],vmean,vstd,vmin,vmax,smean;
  PetscReal sigma_init,tau_init,vel_init,rtol,atol_slip,dt_init,dt_max,rand_amp,tmp;
  struct interact_ctx par[1]; /* user-defined work context */
  char geom_file[STRLEN]="geom.in",rsf_file[STRLEN]="rsf.dat";
  PetscBool read_value,warned = PETSC_FALSE,flg;
  char *home_dir = getenv("HOME");char par_file[STRLEN],vel_file[STRLEN];
  snprintf(par_file,STRLEN,"%s/progs/src/interact/petsc_settings.yaml",(home_dir)?(home_dir):("."));
  /* set up structure */
  par->medium=(struct med *)calloc(1,sizeof(struct med));
  medium = par->medium;
  /* 
     default frictional and loading parameters, can be overridden via
     options below; defaults follow the SEAS BP1 benchmark / HBI
     conventions where applicable
  */
  medium->f0  = 0.6; 		/* f_0 reference friction */
  medium->dc = 0.008;		/* D_c [m] */
  medium->vpl = 1e-9;		/* plate motion [m/s] */
  medium->v0 = 1e-6;		/* reference speed [m/s] (HBI: vref) */
  /* 
     timestepping 
  */
  medium->time = medium->slip_line_time = 0.;
  medium->stop_time = 3000*sec_per_year;	      /* stop time */
  /* output */
  medium->print_interval = 0.1*sec_per_year;	      /* for average property output */
  medium->slip_line_dt = sec_per_year;		      /* for output of velocity grid */
  
  /* start up MPI/PETSc, only read the YAML defaults file if it exists */
  PetscFunctionBegin;
  {
    FILE *tst = fopen(par_file,"r");
    if(tst){fclose(tst);}else{par_file[0]='\0';}
  }
  PetscInitialize(&argc,&argv,(par_file[0])?(par_file):(NULL),NULL);
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));

  /* options for this code: 0 dense 1 htools 2 h2opus */
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-use_hmatrix", &use_hmatrix,&read_value));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-rsf_file", rsf_file, STRLEN,&read_value));
  /* physical parameters */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-shear_modulus",&shear_modulus_si,NULL)); /* G [Pa] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-s_wave_speed",&s_wave_speed_si,NULL));   /* c_s [m/s] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-f0",&medium->f0,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-dc",&medium->dc,NULL));     /* [m] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-vpl",&medium->vpl,NULL));   /* [m/s] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-v0",&medium->v0,NULL));     /* reference velocity [m/s] */
  /* G/(2 c_s) radiation damping factor */
  medium->shear_mod_over_2cs_si = shear_modulus_si /(2.0*s_wave_speed_si);
  /* initial conditions */
  sigma_init = 50e6;		/* [Pa] */
  tau_init = -1;		/* [Pa], if < 0, will use f0 * sigma_init */
  vel_init = medium->vpl;	/* [m/s] */
  rand_amp = 0;			/* random multiplier amplitude for initial state, e.g. 0.05 */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-sigma_init",&sigma_init,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-tau_init",&tau_init,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-vel_init",&vel_init,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-rand_amp",&rand_amp,NULL));
  if(tau_init < 0)
    tau_init = medium->f0 * sigma_init;
  /* normal stress limiter as in HBI's limitsigma (default off here,
     irrelevant for planar faults) */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-limit_sigma",&rsf_limit_sigma,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-min_sigma",&rsf_min_sigma,NULL)); /* [Pa] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-max_sigma",&rsf_max_sigma,NULL)); /* [Pa] */
  /* time stepping controls */
  rtol = 1e-4;			/* HBI eps_r default */
  atol_slip = 1e-3;		/* [m] absolute tolerance for the slip entries,
				   effectively excluding slip from error control as in HBI */
  dt_init = 1.0;		/* [s], HBI dtinit default */
  dt_max = 1e10;		/* [s], HBI dtmax default */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-rtol",&rtol,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-atol_slip",&atol_slip,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-dt_init",&dt_init,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-dt_max",&dt_max,NULL));
  tmp = medium->stop_time/sec_per_year;
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-stop_time_yr",&tmp,NULL));
  medium->stop_time = tmp * sec_per_year;
  tmp = medium->print_interval/sec_per_year;
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-print_interval_yr",&tmp,NULL));
  medium->print_interval = tmp * sec_per_year;
  tmp = medium->slip_line_dt/sec_per_year;
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-slip_line_dt_yr",&tmp,NULL));
  medium->slip_line_dt = tmp * sec_per_year;

  /* get the geometry */
  read_geometry(geom_file,&medium,&par->fault,TRUE,FALSE,FALSE,FALSE);
  fault = par->fault;
  n = medium->nrflt;
  HEADNODE{
    fprintf(stderr,"%s: read geometry from %s, %i patches, on %i cores\n",
	    argv[0],geom_file,n,medium->comm_size);
    fprintf(stderr,"%s: G %.6e Pa cs %g m/s eta %.6e f0 %g dc %g m v0 %.3e vpl %.3e m/s\n",
	    argv[0],shear_modulus_si,s_wave_speed_si,medium->shear_mod_over_2cs_si,
	    medium->f0,medium->dc,medium->v0,medium->vpl);
    fprintf(stderr,"%s: sigma0 %.6e tau0 %.6e Pa vinit %.3e m/s rand_amp %g\n",
	    argv[0],sigma_init,tau_init,vel_init,rand_amp);
    fprintf(stderr,"%s: rtol %.1e dt_init %g s dt_max %g s stop %g yr\n",
	    argv[0],rtol,dt_init,dt_max,medium->stop_time/sec_per_year);
  }
  /* 
     now, read in a,b variations (stored in fault[].mu_s, fault[].mu_d)
  */
  read_rsf(rsf_file,medium,fault);
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
  calc_petsc_Isn_matrices(medium,fault,use_hmatrix, shear_modulus_si/SHEAR_MODULUS,0,&medium->Is); /* shear stress */
  calc_petsc_Isn_matrices(medium,fault,use_hmatrix,-shear_modulus_si/SHEAR_MODULUS,1,&medium->In); /* normal stress, compression positive */
  if(use_hmatrix)
    PetscCall(MatView(medium->Is,PETSC_VIEWER_STDOUT_WORLD));
  /* 
     preallocate the RHS work vectors with the matrix row layout
  */
  PetscCall(MatCreateVecs(medium->Is,&medium->rsf_vel,&medium->rsf_tau_dot));
  PetscCall(VecDuplicate(medium->rsf_tau_dot,&medium->rsf_sigma_dot));
  /* 
     compute backslip stressing rates sinc[0] (shear) and sinc[1]
     (normal, compression positive) from slip rate -vpl, made
     available on all nodes
  */
  PetscCall(MatCreateVecs(medium->Is,&islip_rate_vec,&stress_rate));
  PetscCall(VecSet(islip_rate_vec,-medium->vpl));
  for(i=0;i < 2;i++){
    if(i==0)
      PetscCall(MatMult(medium->Is,islip_rate_vec,stress_rate)); /* this is shear */
    else
      PetscCall(MatMult(medium->In,islip_rate_vec,stress_rate)); /* this is normal */
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
     solution vector: local size has to be INT_RSF_DIM times the local
     matrix row range so that local patch i of the matrix layout owns
     entries [i*INT_RSF_DIM, (i+1)*INT_RSF_DIM) of x on the same rank
  */
  m=INT_RSF_DIM*n;
  PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
  PetscCall(VecSetSizes(x,INT_RSF_DIM*medium->rn,m));
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
  for (i = medium->rs,j=i*INT_RSF_DIM; i < medium->re; i++,j+=INT_RSF_DIM){
    PetscCall(PetscRandomGetValue(prand,&rand_fac));
    fault[i].s[NORMAL] = sigma_init;	/* compression positive */
    fault[i].s[STRIKE] = tau_init;
    fault[i].u[0] = vel_init;
    /* state consistent with v = vel_init */
    state = fault[i].mu_s*log(2.0*medium->v0/fault[i].u[0]*sinh(fault[i].s[STRIKE]/fault[i].s[NORMAL]/fault[i].mu_s));
    if(rand_amp > 0)
      state *= rand_fac;
    PetscCall(VecSetValue(x, (PetscInt)(j+0),state, INSERT_VALUES));
    PetscCall(VecSetValue(x, (PetscInt)(j+1),fault[i].s[STRIKE], INSERT_VALUES)); /* shear stress */
    PetscCall(VecSetValue(x, (PetscInt)(j+2),fault[i].s[NORMAL], INSERT_VALUES)); /* normal stress */
    PetscCall(VecSetValue(x, (PetscInt)(j+3),0.0, INSERT_VALUES)); /* total slip */
    /* check discretization against the quasi-static cohesive zone
       size Lb = G dc/(b sigma), using the smallest patch dimension */
    patch_l = 2.0*((fault[i].l < fault[i].w)?(fault[i].l):(fault[i].w));
    stable_l = shear_modulus_si*medium->dc/fault[i].mu_d/fault[i].s[NORMAL];
    if((patch_l > stable_l/3.)&&(!warned)){
      fprintf(stderr,"%s: WARNING: patch %05i: size %.3e m vs. Lb = G dc/(b sigma) %.3e m, coarser than Lb/3\n",
	      argv[0],i,patch_l,stable_l);
      fprintf(stderr,"%s: a %g b %g tau %.4e sigma %.4e psi %.4e v %.4e (suppressing further warnings)\n",
	      argv[0],fault[i].mu_s,fault[i].mu_d,fault[i].s[STRIKE],fault[i].s[NORMAL],state,
	      vel_from_rsf(fault[i].s[STRIKE], fault[i].s[NORMAL], state,
			   fault[i].mu_s,medium->v0,dummy,(dummy+1),(dummy+2),medium));
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
  for (i = medium->rs,j=i*INT_RSF_DIM; i < medium->re; i++,j+=INT_RSF_DIM){
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
  /* output times should not constrain the natural step size selection */
  PetscCall(TSSetExactFinalTime(ts,TS_EXACTFINALTIME_INTERPOLATE));
  /* HBI rkqs uses a pure relative infinity norm error measure */
  PetscCall(TSSetTolerances(ts,0.0,vatol,rtol,NULL));
  PetscCall(PetscOptionsHasName(NULL,NULL,"-ts_adapt_wnormtype",&flg));
  if(!flg)
    PetscCall(PetscOptionsSetValue(NULL,"-ts_adapt_wnormtype","infinity"));
  /* HBI rkqs: step growth limited to 2x, shrink to >= 0.5x */
  PetscCall(TSGetAdapt(ts,&adapt));
  PetscCall(TSAdaptSetClip(adapt,0.5,2.0));
  PetscCall(TSAdaptSetStepLimits(adapt,0.0,dt_max));
  /* 
     reject steps whose stages leave the physical domain (overflow,
     non-positive normal stress); the equivalent of HBI rkqs' NaN/Inf
     guard, without this an overshooting stage can produce a NaN error
     estimate which the controller would silently accept
  */
  rsf_par_static = par;
  PetscCall(TSSetFunctionDomainError(ts,rsf_domain_check));
  PetscCall(TSSetFromOptions(ts));
  /* per-step monitor, comparable to HBI's monitor.dat */
  HEADNODE{
    rsf_monitor_out = fopen("rsf_monitor.dat","w");
    fprintf(rsf_monitor_out,"# step time[s] time[yr] dt[s] log10(max|v|[m/s]) mean_slip[m] mean_mu max_sigma[Pa] min_sigma[Pa]\n");
  }
  PetscCall(TSMonitorSet(ts,rsf_TS_Monitor,par,NULL));

  field_out = 0;
  HEADNODE{
    fout1 = fopen("rsf_stats.dat","w");
    fprintf(fout1,"# time[yr] mean_vel std_vel min_vel max_vel mean_slip\n");
  }
  /* for sending to head node */
  PetscCall(VecScatterCreateToZero(x,&ctx,&xout));		/* solution vector */
  while(medium->time < medium->stop_time){	/* advance until next output */
    PetscCall(TSSetMaxTime(ts,(medium->time + medium->print_interval))); /* advance solution to next output time */
    PetscCall(TSSolve(ts, x));
    PetscCall(TSGetSolveTime(ts,&(medium->time))); /* get the current time */
    /* 
       distribute to zero node - collective on all ranks!
    */
    PetscCall(VecScatterBegin(ctx,x,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx,x,xout,INSERT_VALUES,SCATTER_FORWARD));
    HEADNODE{
      PetscCall(VecGetArray(xout,&values));
      /* all patches loop */
      for(sum[0]=sum[1]=sum[2]=0.,vmin=1e20,vmax=-1e20,
	    i=0,j=0;i < n;i++,j += INT_RSF_DIM){
	state = values[j];
	fault[i].s[STRIKE] = values[j+1]; /* shear stress */
	fault[i].s[NORMAL] = values[j+2]; /* normal stress */
	slip = values[j+3];
	/* velocity */
	fault[i].u[STRIKE] = vel_from_rsf(fault[i].s[STRIKE],fault[i].s[NORMAL],
					  state,fault[i].mu_s,medium->v0,dummy,(dummy+1),(dummy+2),medium);
	/* compute average slip rate */
	if(fault[i].u[STRIKE] < vmin)
	  vmin = fault[i].u[STRIKE];
	if(fault[i].u[STRIKE] > vmax)
	  vmax = fault[i].u[STRIKE];
	sum[0] += fault[i].u[STRIKE];
	sum[1] += fault[i].u[STRIKE] * fault[i].u[STRIKE];
	sum[2] += slip;				     /* mean slip */
      }
      /* output */
      vstd = sqrt(((COMP_PRECISION)n * sum[1] - sum[0] * sum[0]) / ((COMP_PRECISION)(n*(n-1))));
      vmean = sum[0]/(COMP_PRECISION)n;
      smean = sum[2]/(COMP_PRECISION)n;
      /* print some stats
	 
	 time[yr] mean_vel std_vel min_vel max_vel mean_slip 
      */
      fprintf(fout1,"%.8f %e %e %e %e %e\n",medium->time/sec_per_year,vmean,vstd,vmin,vmax,smean);
      /* some more output */
      if(medium->time - medium->slip_line_time > medium->slip_line_dt){
	if(!field_out){
	  ierr = system("mkdir -p tmp_rsf");
	  if(ierr){
	    fprintf(stderr,"%s: error making output directory\n",argv[0]);
	    exit(-1);
	  }
	}
	/* 
	   print slip rate field to file 
	*/
	snprintf(vel_file,STRLEN,"tmp_rsf/vel-%012.5e-%06i-gmt",medium->time/sec_per_year,field_out);
	fout2 = fopen(vel_file,"w");
	for(i=0;i < n;i++)
	  print_patch_geometry_and_bc(i,fault,PSXYZ_STRIKE_DISP_OUT_MODE,medium->time,FALSE,fout2,FALSE,dummy);
	fclose(fout2);
	field_out++;
	medium->slip_line_time = medium->time;
      }
      PetscCall(VecRestoreArray(xout,&values));
    }
  }
  HEADNODE{
    fclose(fout1);
    if(rsf_monitor_out)
      fclose(rsf_monitor_out);
  }
  PetscCall(VecScatterDestroy(&ctx));
  
  /*View information about the time-stepping method and the solution at the end time.*/
  PetscCall(TSView(ts, PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(MatDestroy(&medium->Is));
  PetscCall(MatDestroy(&medium->In));
  PetscCall(VecDestroy(&medium->rsf_vel));
  PetscCall(VecDestroy(&medium->rsf_tau_dot));
  PetscCall(VecDestroy(&medium->rsf_sigma_dot));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&vatol));
  PetscCall(VecDestroy(&xout));

  /*  */
  PetscCall(TSDestroy(&ts)); 
  PetscCall(PetscFinalize());

#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);
}

#ifdef USE_PETSC
/* 
   compute time-derivatives for the state vector
   x = (psi1,tau1,sigma1,s1,psi2,tau2,sigma2,s2,....) of size INT_RSF_DIM*n

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
  PetscFunctionBeginUser;
  par = (struct interact_ctx *)ptr;medium = par->medium;fault = par->fault;
  /* consistency check of the layout alignment (cheap) */
  PetscCall(VecGetLocalSize(X, &ln));
  if(ln != INT_RSF_DIM*medium->rn){
    fprintf(stderr,"rsf_ODE_RHSFunction: layout mismatch, local %i vs %i x %i\n",
	    (int)ln,INT_RSF_DIM,(int)medium->rn);
    exit(-1);
  }
  PetscCall(VecGetArrayRead(X,&x));
  /* 
     compute the slip rate vector
     v = 2 v0 exp(-psi/a) sinh(tau/(sigma a)) 
  */
  PetscCall(VecGetArray(medium->rsf_vel,&vel));
  /* i global patch, j local x offset, k local patch */
  for (i = medium->rs, j=0, k=0; i < medium->re; i++, j+=INT_RSF_DIM, k++) {
    /* local x layout: state = x[j], tau = x[j+1], sigma = x[j+2], slip = x[j+3] */
    vel[k] = vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,medium->v0,
			  &mu, &scaled_tau, &exp_fac, medium);
  }
  PetscCall(VecRestoreArray(medium->rsf_vel,&vel));
  /* 
     stressing rates from all slipping faults, this works for dense
     and H matrix; background backslip loading (sinc) is added below
  */
  PetscCall(MatMult(medium->Is, medium->rsf_vel, medium->rsf_tau_dot));
  PetscCall(MatMult(medium->In, medium->rsf_vel, medium->rsf_sigma_dot));
  /* 
     compute derivatives 
  */
  PetscCall(VecGetArrayRead(medium->rsf_vel,&velr));
  PetscCall(VecGetArrayRead(medium->rsf_tau_dot,&tau_dot));
  PetscCall(VecGetArrayRead(medium->rsf_sigma_dot,&sigma_dot));
  PetscCall(VecGetArray(F,&f));
  for (i = medium->rs, j=0, k=0; i < medium->re; i++, j+=INT_RSF_DIM, k++) {
    /* a = fault[].mu_s, b = fault[].mu_d as read by read_rsf */
    a = fault[i].mu_s;b = fault[i].mu_d;
    /* 
       d psi/dt = b/dc (v0 exp((f0-psi)/b) - |v|)
       
       aging law; HBI uses |v|, important for the sinh regularization
       which allows v < 0
    */
    if(b != 0.0)		/* b == 0 guard as in HBI */
      f[j] = (b/medium->dc) * (medium->v0 * PetscExpReal((medium->f0 - x[j])/b) - fabs(velr[k]));
    else
      f[j] = 0.0;
    /* 
       d sigma/dt, compression positive (In was scaled by -1)
    */
    sdot = sigma_dot[k] + fault[i].sinc[1];
    if(rsf_limit_sigma){	/* HBI limitsigma behavior */
      if((x[j+2] < rsf_min_sigma)||(x[j+2] > rsf_max_sigma))
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
    pre_fac =  2.0*medium->v0/(a * x[j+2]);	/* 2v0/(a sigma) */
    cosh_fac  = PetscCoshReal(scaled_tau) * exp_fac;
    dvdtau   =   pre_fac * cosh_fac;
    dvdsigma =  -pre_fac * cosh_fac * mu;
    dvdstate = -velr[k]/a;
    f[j+1]  = (tau_dot[k] + fault[i].sinc[0] 
	       - medium->shear_mod_over_2cs_si * (dvdsigma * f[j+2] + dvdstate * f[j]));
    f[j+1] /= (1.0 + medium->shear_mod_over_2cs_si * dvdtau);
    /* d slip/dt */
    f[j+3] = velr[k];
  }
  /*  */
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(F,&f));
  PetscCall(VecRestoreArrayRead(medium->rsf_vel,&velr));
  PetscCall(VecRestoreArrayRead(medium->rsf_tau_dot,&tau_dot));
  PetscCall(VecRestoreArrayRead(medium->rsf_sigma_dot,&sigma_dot));
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
  PetscFunctionBeginUser;
  medium = rsf_par_static->medium;fault = rsf_par_static->fault;
  lok = 1;
  PetscCall(VecGetArrayRead(X,&x));
  for (i = medium->rs, j=0; (i < medium->re) && lok; i++, j+=INT_RSF_DIM) {
    a = fault[i].mu_s;b = fault[i].mu_d;
    if((!PetscIsNormalReal(x[j+2])) || (x[j+2] <= 0.0)){lok = 0;break;} /* sigma */
    if(PetscIsInfOrNanReal(x[j]) || PetscIsInfOrNanReal(x[j+1]) || PetscIsInfOrNanReal(x[j+3])){lok = 0;break;}
    if(fabs(x[j+1]/(x[j+2]*a)) > arg_max){lok = 0;break;} /* sinh(tau/(sigma a)) */
    if(fabs(x[j]/a) > arg_max){lok = 0;break;}		  /* exp(-psi/a) */
    if((b != 0.0) && (fabs((medium->f0 - x[j])/b) > arg_max)){lok = 0;break;} /* aging law exp */
  }
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCallMPI(MPI_Allreduce(&lok,&gok,1,MPIU_INT,MPI_MIN,PETSC_COMM_WORLD));
  *accept = (gok)?(PETSC_TRUE):(PETSC_FALSE);
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* 
   per accepted time step monitor, output comparable to HBI's
   monitor.dat: max slip velocity, mean slip, mean friction, normal
   stress extrema
*/
PetscErrorCode rsf_TS_Monitor(TS ts,PetscInt step,PetscReal time,Vec X,void *ptr)
{
  const PetscScalar *x;
  struct med *medium;struct flt *fault;
  struct interact_ctx *par;
  PetscInt i,j;
  PetscReal v,lsum[2],gsum[2],lminmax[3],gminmax[3],dt,d1,d2,d3;
  const PetscReal sec_per_year = 365.25*24.*60.*60.;
  PetscFunctionBeginUser;
  par = (struct interact_ctx *)ptr;medium = par->medium;fault = par->fault;
  if(time == rsf_monitor_last_time) /* skip repeat calls at interpolated output times */
    PetscFunctionReturn(PETSC_SUCCESS);
  rsf_monitor_last_time = time;
  PetscCall(TSGetTimeStep(ts,&dt));
  PetscCall(VecGetArrayRead(X,&x));
  lsum[0]=lsum[1]=0.0;
  lminmax[0] = -1e30;		/* max |v| */
  lminmax[1] = -1e30;		/* max sigma */
  lminmax[2] = -1e30;		/* -min sigma */
  for (i = medium->rs, j=0; i < medium->re; i++, j+=INT_RSF_DIM) {
    v = vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,medium->v0,
		     &d1,&d2,&d3,medium);
    v = fabs(v);
    if(v > lminmax[0])lminmax[0] = v;
    if( x[j+2] > lminmax[1])lminmax[1] =  x[j+2];
    if(-x[j+2] > lminmax[2])lminmax[2] = -x[j+2];
    lsum[0] += x[j+3];		/* slip */
    lsum[1] += x[j+1]/x[j+2];	/* mu */
  }
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCallMPI(MPI_Reduce(lsum,gsum,2,MPIU_REAL,MPI_SUM,0,PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(lminmax,gminmax,3,MPIU_REAL,MPI_MAX,0,PETSC_COMM_WORLD));
  if((medium->comm_rank == 0) && rsf_monitor_out){
    fprintf(rsf_monitor_out,"%9i %20.8e %17.10f %15.8e %12.7f %15.8e %10.6f %15.8e %15.8e\n",
	    (int)step,time,time/sec_per_year,dt,
	    log10(gminmax[0]),gsum[0]/(PetscReal)medium->nrflt,
	    gsum[1]/(PetscReal)medium->nrflt,
	    gminmax[1],-gminmax[2]);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

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

#endif
