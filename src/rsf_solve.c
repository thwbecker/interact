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


#define INT_RSF_DIM 4		/* psi,tau,sigma,u(slip) per patch */


/* 
   parameters for optional normal stress limiter, cf. HBI limitsigma
   (file scope since they are only used by the driver and RHS here) 
*/
static PetscBool rsf_limit_sigma = PETSC_FALSE;
static PetscReal rsf_min_sigma = 1e6, rsf_max_sigma = 300e6; /* Pa */
/* optional per-cell D_c [m] (geometry order); NULL => use uniform
   medium->dc.  Used for the BP5 nucleation patch (reduced D_RS) and
   read in the RHS, so it lives at file scope alongside the statics above */
static PetscReal *rsf_dc_vec = NULL;
/* 
   context for the monitor and event functions, following the layout
   in hmatrix_test/ode_solve_test.c; the RHS keeps using the plain
   interact_ctx which is shared with the matrix assembly routines
*/
struct rsf_out_ctx{
  struct interact_ctx *par;	/* fault geometry and medium */
  /* change/time triggered state monitor */
  PetscReal old_time,dt_monitor,adx_monitor,rdx_monitor,monitor_tmin;
  Vec Xold;
  FILE *fout_monitor;
  /* periodic full-field output (stats line, velocity snapshots) */
  PetscReal next_print_time;
  int field_out;
  VecScatter gather;
  Vec gathered;
  FILE *fout_stats;
  /* slip velocity threshold crossing events */
  PetscBool track_events,slipping;
  PetscReal vel_event,vel_event_hyst,event_tmin;
  int nevent;
  FILE *fout_event;
};
static PetscErrorCode rsf_init_monitor_and_event(struct rsf_out_ctx *,struct interact_ctx *,
						 PetscReal,PetscReal,PetscReal,PetscReal,
						 PetscReal,PetscReal,PetscReal,PetscBool,
						 PetscReal,Vec,PetscReal);
static PetscErrorCode rsf_finalize_monitor_and_event(struct rsf_out_ctx *);
static PetscErrorCode rsf_event_function(TS,PetscReal,Vec,PetscScalar *,void *);
static PetscErrorCode rsf_post_event(TS,PetscInt,PetscInt[],PetscReal,Vec,PetscBool,void *);
/* context access for the domain check, whose callback has no user pointer */
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
  PetscBool use_full_space=PETSC_FALSE;
  /*  */
  VecScatter ctx;
  PetscInt tvmode = 0;		/* triangle evaluattion mode */
  Vec xout,x,vatol,islip_rate_vec,stress_rate;
  struct med *medium;struct flt *fault;
  PetscInt m,i,n,j,use_hmatrix=IHMAT_TYPE_DENSE;
  PetscRandom prand;
  PetscReal state,patch_l,stable_l,rand_fac,dummy[3];
  PetscReal sigma_init,tau_init,vel_init,rtol,atol_slip,dt_init,dt_max,rand_amp,tmp;
  PetscReal dt_monitor,rdx_monitor,adx_monitor,monitor_tmin,vel_event,vel_event_hyst,event_tmin;
  PetscBool track_events;
  PetscInt event_direction[1] = {0}; /* detect both crossing directions */
  PetscBool event_terminate[1] = {PETSC_FALSE};
  struct rsf_out_ctx uc[1];
  struct interact_ctx par[1]; /* user-defined work context */
  char geom_file[STRLEN]="geom.in",rsf_file[STRLEN]="rsf.dat",rsf_ic_file[STRLEN]="",rsf_dc_file[STRLEN]="";
  FILE *iin;
  PetscInt ii;
  PetscReal *ic_tau=NULL,*ic_vel=NULL;
  PetscBool have_ic=PETSC_FALSE,have_dc=PETSC_FALSE; /* I/O */
  PetscBool read_value,warned = PETSC_FALSE,flg;
  char *home_dir = getenv("HOME");char par_file[STRLEN];
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
  /*  */
  medium->calc_sigma_dot = FALSE; /* compute normal stress change? */
  /* initial conditions, can overwrite below */
  sigma_init = 50e6;		/* [Pa] */
  tau_init = -1;		/* [Pa], if < 0, will use f0 * sigma_init */
  vel_init = medium->vpl;	/* [m/s] */
  rand_amp = 0;			/* random multiplier amplitude for initial state, e.g. 0.05 */
  
  /* time stepping controls, can override below */
  rtol = 1e-4;			/* HBI eps_r default */
  atol_slip = 1e-3;		/* [m] absolute tolerance for the slip entries,
				   effectively excluding slip from error control as in HBI */
  dt_init = 1.0;		/* [s], HBI dtinit default */
  dt_max = 1e10;		/* [s], HBI dtmax default */
  /* 
     timestepping 
  */
  dt_monitor = 5.0;		/* [yr], also caps the step size */
  rdx_monitor = 1e-4;		/* relative state change trigger */
  adx_monitor = 0.0;		/* absolute state change trigger, <= 0: off */
  monitor_tmin = 0.0;		/* [yr] suppress monitor output before */
  vel_event = 1e-3;		/* [m/s] event velocity threshold, cf. HBI velth */
  vel_event_hyst = 0.5;		/* arrest at vel_event * hyst, debounces spurious
				   crossings from interpolant oscillation during
				   the rapid coseismic acceleration */
  event_tmin = 0.0;		/* [yr] suppress event output before */

  medium->time = medium->slip_line_time = 0.;
  medium->stop_time = 3000*sec_per_year;	      /* stop time */
  /* output */
  medium->print_interval = 0.1*sec_per_year;	      /* for average property output */
  /*
     velocity-field snapshots to tmp_rsf/ (written in rsf_TS_Monitor): a
     full-field dump for slip-evolution visualization, NOT required for the
     BP5 benchmark quantities, so it defaults OFF here.  The earlier default
     (1 yr) writes one frame per model-year; over a multi-century BP5 run
     that is ~10^3 frames, each roughly the size of the on-fault field (~1 MB
     at 0.5 km / N=16000), and the monitor does not purge tmp_rsf between
     runs -- so repeated runs ACCUMULATE frames and can fill the working
     disk.  A full disk then makes the snapshot writes (and the
     rsf_monitor.dat writes) fail, which is easy to misread as a solver
     stall.  Opt in when you want the frames, e.g. -slip_line_dt_yr 1 for the
     fine SEAS-style interseismic cadence or a coarser value (say 10-20 yr)
     for a cycle-scale view, and clean tmp_rsf between runs either way.
  */
  medium->slip_line_dt = 1.0e9*sec_per_year;	      /* OFF; opt in via -slip_line_dt_yr */
  
  /* start up MPI/PETSc, only read the YAML defaults file if it exists */
  PetscFunctionBegin;
  {
    FILE *tst = fopen(par_file,"r");
    if(tst){fclose(tst);}else{par_file[0]='\0';}
  }
  PetscInitialize(&argc,&argv,(par_file[0])?(par_file):(NULL),NULL);
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));

  /* options for this code: 0 dense 1 htools 2 h2opus 3 hacapk 4 hmmvp 5 BIGWHAM, see petsc_prototypes.h */
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-use_hmatrix", &use_hmatrix,&read_value));
  /* HTOOL (use_hmatrix==IHMAT_TYPE_HTOOLS ) compressor default: prefer sympartialACA over
     PETSc's built-in SVD default.  On the SEAS BP5 4000-cell test (1 core),
     the SVD compressor spent ~48 s assembling the H-matrix vs ~6.5 s for
     HACApK and ~8 s for sympartialACA, for no measurable accuracy gain
     (recurrence and the full max|V| trace agree with the dense solution to
     <0.01 yr and <0.005 in log10|V| at epsilon<=1e-4 -- the error floor is
     set by the ODE rtol, not the H-matrix tolerance).  Override on the
     command line with -mat_htool_compressor {SVD,fullACA,partialACA}. See
     rsf_solve.md for the full HTOOL/HACApK/dense performance comparison. */
  medium->use_hmatrix = use_hmatrix;
  switch(medium->use_hmatrix){
  case IHMAT_TYPE_DENSE:
    break;
  case IHMAT_TYPE_HTOOLS: /* this now defaults to sympartialACA, not SVD */
#ifdef USE_PETSC_HMAT		/* htools and H2opus  */
    set_htools_defaults_and_options(medium);
#else
    fprintf(stderr,"%s: h matrix type %i not compiled in, see makefile.petsc \n",argv[0]);
    exit(-1);
#endif
    break;
  case IHMAT_TYPE_H2OPUS: /* this did not used to get called */
#ifdef USE_PETSC_HMAT		/* htools and H2opus  */
    set_h2opus_defaults_and_options(medium);
#else
    fprintf(stderr,"%s: h matrix type %i not compiled in, see makefile.petsc \n",argv[0]);
    exit(-1);
#endif
    break;
  case IHMAT_TYPE_HACAPK: /* this did not used to get called */
#ifdef USE_HACAPK
    set_hacapk_defaults_and_options(medium);
#else
    fprintf(stderr,"%s: h matrix type %i not compiled in, see makefile.petsc \n",argv[0]);
    exit(-1);
#endif
   break;
  case IHMAT_TYPE_HMMVP: /* this did not used to get called */
#ifdef USE_HMMVP
    set_hmmvp_defaults_and_options(medium);
#else
    fprintf(stderr,"%s: h matrix type %i not compiled in, see makefile.petsc \n",argv[0]);
    exit(-1);
#endif
    break;
  case IHMAT_TYPE_BIGWHAM: /* this did not used to get called */
#ifdef USE_BIGWHAM
    set_bigwham_defaults_and_options(medium);
#else
    fprintf(stderr,"%s: h matrix type %i not compiled in, see makefile.petsc \n",argv[0]);
    exit(-1);
#endif
    break;
  default:
    fprintf(stderr,"%s: h matrix type %i undefined\n",argv[0],medium->use_hmatrix);
    exit(-1);
    break;
  }
  
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-full_space", &use_full_space,&read_value)); /* use read in or default */
  medium->full_space = (my_boolean)use_full_space;
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
 
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-sigma_init",&sigma_init,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-tau_init",&tau_init,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-vel_init",&vel_init,NULL));
  /* optional per-cell initial conditions (tau[Pa] vel[m/s] per patch,
     one row per patch in geometry order) for e.g. the SEAS BP5
     nucleation patch; when given, overrides uniform tau_init/vel_init */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-rsf_ic_file",rsf_ic_file,STRLEN,&have_ic));
  /* optional per-cell D_c file (one D_c[m] per patch, geometry order);
     overrides the uniform -dc per cell, e.g. the reduced D_RS in the
     SEAS BP5 nucleation patch */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-rsf_dc_file",rsf_dc_file,STRLEN,&have_dc));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-rand_amp",&rand_amp,NULL));
  if(tau_init < 0)
    tau_init = medium->f0 * sigma_init;
  /* normal stress limiter as in HBI's limitsigma (default off here,
     irrelevant for planar faults) */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-limit_sigma",&rsf_limit_sigma,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-min_sigma",&rsf_min_sigma,NULL)); /* [Pa] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-max_sigma",&rsf_max_sigma,NULL)); /* [Pa] */
  
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
  /* 
     monitor and event controls, cf. hmatrix_test/ode_solve_test.c:
     the state monitor logs when the relative or absolute state change
     since the last logged state exceeds the limits, or at least every
     dt_monitor; events are slip velocity threshold crossings located
     by the TS event handler
  */
 
  track_events = PETSC_TRUE;
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-dt_monitor_yr",&dt_monitor,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-rdx_monitor",&rdx_monitor,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-adx_monitor",&adx_monitor,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-monitor_tmin_yr",&monitor_tmin,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-vel_event",&vel_event,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-vel_event_hyst",&vel_event_hyst,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-event_tmin_yr",&event_tmin,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-track_events",&track_events,NULL));
  /* triangular patch evaluation scheme, cf. -tv of interact */
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-tv",&tvmode,NULL));
  medium->tri_eval_mode = (MODE_TYPE)tvmode;

  dt_monitor *= sec_per_year;monitor_tmin *= sec_per_year;event_tmin *= sec_per_year;

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
    rsf_dc_vec=(PetscReal *)malloc((size_t)n*sizeof(PetscReal));
    if(!rsf_dc_vec){
      fprintf(stderr,"%s: per-cell dc alloc failed\n",argv[0]);
      exit(-1);
    }
    iin=myopen(rsf_dc_file,"r");
    for(ii=0;ii < n;ii++)
      if(fscanf(iin,"%lf",(rsf_dc_vec+ii))!=1){
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
  calc_petsc_Isn_matrices(medium,fault,use_hmatrix,
			  shear_modulus_si/SHEAR_MODULUS,0,&medium->Is,medium->Is_hctx); /* shear stress */
  if(medium->calc_sigma_dot)
    calc_petsc_Isn_matrices(medium,fault,use_hmatrix,
			    -shear_modulus_si/SHEAR_MODULUS,1,&medium->In,medium->In_hctx); /* normal
											       stress,
											       compression
											       positive */
  if(use_hmatrix)
    PetscCall(MatView(medium->Is,PETSC_VIEWER_STDOUT_WORLD));
  /* 
     preallocate the RHS work vectors with the matrix row layout
  */
  PetscCall(MatCreateVecs(medium->Is,&medium->rsf_vel,&medium->rsf_tau_dot));
  if(medium->calc_sigma_dot)
    PetscCall(VecDuplicate(medium->rsf_tau_dot,&medium->rsf_sigma_dot));
  /* 
     compute backslip stressing rates sinc[0] (shear) and sinc[1]
     (normal, compression positive) from slip rate -vpl, made
     available on all nodes
  */
  PetscCall(MatCreateVecs(medium->Is,&islip_rate_vec,&stress_rate));
  PetscCall(VecSet(islip_rate_vec,-medium->vpl));
  for(i=0;i < ((medium->calc_sigma_dot)?(2):(1));i++){
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
    if(have_ic){
      fault[i].s[STRIKE] = ic_tau[i];	/* per-cell initial shear stress [Pa] */
      fault[i].u[0]      = ic_vel[i];	/* per-cell initial slip velocity [m/s] */
    }else{
      fault[i].s[STRIKE] = tau_init;
      fault[i].u[0] = vel_init;
    }
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
    stable_l = shear_modulus_si*(rsf_dc_vec?rsf_dc_vec[i]:medium->dc)/fault[i].mu_d/fault[i].s[NORMAL];
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
  rsf_par_static = par;
  PetscCall(TSSetFunctionDomainError(ts,rsf_domain_check));
  PetscCall(TSSetFromOptions(ts));
  /* 
     set up the change/time triggered state monitor, the periodic
     field output, and the velocity threshold event tracker,
     cf. hmatrix_test/ode_solve_test.c
  */
  PetscCall(rsf_init_monitor_and_event(uc,par,dt_monitor,adx_monitor,rdx_monitor,
				       monitor_tmin,event_tmin,vel_event,vel_event_hyst,
				       track_events,medium->time,x,vel_init));
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
  
  PetscCall(VecDestroy(&medium->rsf_vel));
  PetscCall(VecDestroy(&medium->rsf_tau_dot));
  if(medium->calc_sigma_dot){
    PetscCall(MatDestroy(&medium->In));
    PetscCall(VecDestroy(&medium->rsf_sigma_dot));
  }
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&vatol));

  /*  */
  PetscCall(TSDestroy(&ts)); 
  PetscCall(PetscFinalize());
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
  if(medium->calc_sigma_dot)
    PetscCall(MatMult(medium->In, medium->rsf_vel, medium->rsf_sigma_dot));
  /* 
     compute derivatives 
  */
  PetscCall(VecGetArrayRead(medium->rsf_vel,&velr));
  PetscCall(VecGetArrayRead(medium->rsf_tau_dot,&tau_dot));
  if(medium->calc_sigma_dot)
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
      f[j] = (b/(rsf_dc_vec?rsf_dc_vec[i]:medium->dc)) * (medium->v0 * PetscExpReal((medium->f0 - x[j])/b) - fabs(velr[k]));
    else
      f[j] = 0.0;
    /* 
       d sigma/dt, compression positive (In was scaled by -1)
    */
    if(medium->calc_sigma_dot)
      sdot = sigma_dot[k] + fault[i].sinc[1];
    else{
      sdot = 0.;
    }
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
  if(medium->calc_sigma_dot)
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
   monitoring function, cf. myMonitor in hmatrix_test/ode_solve_test.c:
   logs the global state summary (max slip velocity, mean slip, mean
   friction, normal stress extrema) when the relative or absolute
   change of the state vector since the last logged state exceeds the
   limits, or when more than dt_monitor has elapsed; also handles the
   periodic full-field output (stats line and velocity snapshots);
   unlike ode_solve_test.c the reference state is only updated on
   output, so the criteria decimate relative to accepted steps
*/
PetscErrorCode rsf_TS_Monitor(TS ts,PetscInt step,PetscReal time,Vec X,void *ptr)
{
  const PetscScalar *x;
  struct med *medium;struct flt *fault;
  struct rsf_out_ctx *uc;
  Vec DX;
  PetscInt i,j;
  PetscBool bail;
  PetscReal v,lsum[2],gsum[2],lminmax[3],gminmax[3],dt,d1,d2,d3,dx_norm,x_norm;
  const PetscReal sec_per_year = 365.25*24.*60.*60.;
  PetscFunctionBeginUser;
  uc = (struct rsf_out_ctx *)ptr;
  medium = uc->par->medium;fault = uc->par->fault;
  if(step < 0)	      /* negative indicates an interpolated solution */
    PetscFunctionReturn(PETSC_SUCCESS);
  if(time >= uc->monitor_tmin){
    /* change since the last logged state */
    PetscCall(VecNorm(X,NORM_2,&x_norm));
    PetscCall(VecDuplicate(X,&DX));
    PetscCall(VecWAXPY(DX,-1.0,X,uc->Xold)); /* dx = x_old - x */
    PetscCall(VecNorm(DX,NORM_2,&dx_norm));
    PetscCall(VecDestroy(&DX));
    bail = PETSC_FALSE;
    if((x_norm > 1e-15) && (dx_norm/x_norm > uc->rdx_monitor))
      bail = PETSC_TRUE;
    if((uc->adx_monitor > 0) && (dx_norm > uc->adx_monitor))
      bail = PETSC_TRUE;
    if(fabs(time - uc->old_time) > uc->dt_monitor)
      bail = PETSC_TRUE;
    if(bail){
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
      if((medium->comm_rank == 0) && uc->fout_monitor){
	fprintf(uc->fout_monitor,"%9i %20.8e %17.10f %15.8e %12.7f %15.8e %10.6f %15.8e %15.8e\n",
		(int)step,time,time/sec_per_year,dt,
		log10(gminmax[0]),gsum[0]/(PetscReal)medium->nrflt,
		gsum[1]/(PetscReal)medium->nrflt,
		gminmax[1],-gminmax[2]);
      }
      /* store last logged state */
      uc->old_time = time;
      PetscCall(VecCopy(X,uc->Xold));
    }
  }
  /* 
     periodic full-field output (stats line and, less frequently,
     velocity snapshots); the decision is identical on all ranks
     and the scatter is collective
  */
  if(time >= uc->next_print_time){
    while(uc->next_print_time <= time) /* may cross several intervals in one step */
      uc->next_print_time += medium->print_interval;
    PetscCall(VecScatterBegin(uc->gather,X,uc->gathered,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(uc->gather,X,uc->gathered,INSERT_VALUES,SCATTER_FORWARD));
    if(medium->comm_rank == 0){
      PetscScalar *values;
      PetscReal sum[3],vmin,vmax,vmean,vstd,smean,du1,du2,du3;
      FILE *fout2;
      char vel_file[STRLEN];
      int n = medium->nrflt,ierr2;
      PetscCall(VecGetArray(uc->gathered,&values));
      for(sum[0]=sum[1]=sum[2]=0.,vmin=1e20,vmax=-1e20,
	    i=0,j=0;i < n;i++,j += INT_RSF_DIM){
	fault[i].s[STRIKE] = values[j+1]; /* shear stress */
	fault[i].s[NORMAL] = values[j+2]; /* normal stress */
	fault[i].u[STRIKE] = vel_from_rsf(values[j+1],values[j+2],values[j],fault[i].mu_s,
					  medium->v0,&du1,&du2,&du3,medium);
	if(fault[i].u[STRIKE] < vmin)vmin = fault[i].u[STRIKE];
	if(fault[i].u[STRIKE] > vmax)vmax = fault[i].u[STRIKE];
	sum[0] += fault[i].u[STRIKE];
	sum[1] += fault[i].u[STRIKE] * fault[i].u[STRIKE];
	sum[2] += values[j+3];				     /* slip */
      }
      vstd = (n > 1)?(sqrt(((PetscReal)n * sum[1] - sum[0]*sum[0])/((PetscReal)n*(PetscReal)(n-1)))):(0.0);
      vmean = sum[0]/(PetscReal)n;
      smean = sum[2]/(PetscReal)n;
      if(uc->fout_stats)
	fprintf(uc->fout_stats,"%.8f %e %e %e %e %e\n",time/sec_per_year,vmean,vstd,vmin,vmax,smean);
      if(time - medium->slip_line_time > medium->slip_line_dt){
	if(!uc->field_out){
	  ierr2 = system("mkdir -p tmp_rsf");
	  if(ierr2)
	    fprintf(stderr,"rsf_TS_Monitor: WARNING: error making tmp_rsf output directory\n");
	}
	snprintf(vel_file,STRLEN,"tmp_rsf/vel-%012.5e-%06i-gmt",time/sec_per_year,uc->field_out);
	fout2 = fopen(vel_file,"w");
	if(fout2){
	  for(i=0;i < n;i++)
	    print_patch_geometry_and_bc(i,fault,PSXYZ_STRIKE_DISP_OUT_MODE,time,FALSE,fout2,FALSE,&du1);
	  fclose(fout2);
	}
	uc->field_out++;
	medium->slip_line_time = time;
      }
      PetscCall(VecRestoreArray(uc->gathered,&values));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* 
   event function, cf. myEventFunction in ode_solve_test.c: locate
   crossings of the maximum slip velocity through the vel_event
   threshold (in log space); evaluated collectively
*/
static PetscErrorCode rsf_event_function(TS ts,PetscReal t,Vec X,PetscScalar *fvalue,void *ctx)
{
  const PetscScalar *x;
  struct rsf_out_ctx *uc;
  struct med *medium;struct flt *fault;
  PetscInt i,j;
  PetscReal v,lvmax,gvmax,d1,d2,d3;
  PetscFunctionBeginUser;
  uc = (struct rsf_out_ctx *)ctx;
  medium = uc->par->medium;fault = uc->par->fault;
  lvmax = 0.0;
  PetscCall(VecGetArrayRead(X,&x));
  for (i = medium->rs, j=0; i < medium->re; i++, j+=INT_RSF_DIM) {
    v = fabs(vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,medium->v0,&d1,&d2,&d3,medium));
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
/* 
   log the event, cf. myPostEventFunction in ode_solve_test.c: write
   time, onset/arrest flag, and the global state summary
*/
static PetscErrorCode rsf_post_event(TS ts,PetscInt nevents,PetscInt event_list[],
				     PetscReal t,Vec X,PetscBool forwardsolve,void *ctx)
{
  const PetscScalar *x;
  struct rsf_out_ctx *uc;
  struct med *medium;struct flt *fault;
  PetscInt i,j;
  PetscReal v,lsum[2],gsum[2],lvmax,gvmax,d1,d2,d3;
  const PetscReal sec_per_year = 365.25*24.*60.*60.;
  PetscFunctionBeginUser;
  uc = (struct rsf_out_ctx *)ctx;
  medium = uc->par->medium;fault = uc->par->fault;
  uc->slipping = (uc->slipping)?(PETSC_FALSE):(PETSC_TRUE); /* threshold crossed */
  if(t >= uc->event_tmin){
    lsum[0]=lsum[1]=0.0;lvmax = 0.0;
    PetscCall(VecGetArrayRead(X,&x));
    for (i = medium->rs, j=0; i < medium->re; i++, j+=INT_RSF_DIM) {
      v = fabs(vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,medium->v0,&d1,&d2,&d3,medium));
      if(v > lvmax)lvmax = v;
      lsum[0] += x[j+3];	/* slip */
      lsum[1] += x[j+1]/x[j+2];	/* mu */
    }
    PetscCall(VecRestoreArrayRead(X,&x));
    PetscCallMPI(MPI_Allreduce(&lvmax,&gvmax,1,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Reduce(lsum,gsum,2,MPIU_REAL,MPI_SUM,0,PETSC_COMM_WORLD));
    uc->nevent++;
    if((medium->comm_rank == 0) && uc->fout_event){
      fprintf(uc->fout_event,"%20.8e %17.10f %2i %12.7f %15.8e %10.6f\n",
	      t,t/sec_per_year,(uc->slipping)?(1):(-1),
	      log10((gvmax > 1e-300)?(gvmax):(1e-300)),
	      gsum[0]/(PetscReal)medium->nrflt,gsum[1]/(PetscReal)medium->nrflt);
      fflush(uc->fout_event);
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* 
   set up the monitor and event tracking environment,
   cf. init_monitor_and_event in ode_solve_test.c
*/
static PetscErrorCode rsf_init_monitor_and_event(struct rsf_out_ctx *uc,struct interact_ctx *par,
						 PetscReal dt_monitor,PetscReal adx_monitor,
						 PetscReal rdx_monitor,PetscReal monitor_tmin,
						 PetscReal event_tmin,PetscReal vel_event,
						 PetscReal vel_event_hyst,PetscBool track_events,
						 PetscReal t_init,Vec X0,PetscReal vel_init)
{
  struct med *medium;
  PetscFunctionBeginUser;
  medium = par->medium;
  uc->par = par;
  uc->dt_monitor = dt_monitor;
  uc->adx_monitor = adx_monitor;
  uc->rdx_monitor = rdx_monitor;
  uc->monitor_tmin = monitor_tmin;
  uc->event_tmin = event_tmin;
  uc->vel_event = vel_event;
  uc->vel_event_hyst = vel_event_hyst;
  uc->track_events = track_events;
  uc->slipping = (vel_init > vel_event)?(PETSC_TRUE):(PETSC_FALSE);
  uc->nevent = 0;
  uc->field_out = 0;
  uc->next_print_time = t_init;	/* fields at start, then every print_interval */
  /* force the first monitor call to log */
  uc->old_time = t_init - 2.0*dt_monitor;
  uc->fout_monitor = uc->fout_stats = uc->fout_event = NULL;
  if(medium->comm_rank == 0){
    uc->fout_monitor = fopen("rsf_monitor.dat","w");
    fprintf(uc->fout_monitor,"# step time[s] time[yr] dt[s] log10(max|v|[m/s]) mean_slip[m] mean_mu max_sigma[Pa] min_sigma[Pa]\n");
    uc->fout_stats = fopen("rsf_stats.dat","w");
    fprintf(uc->fout_stats,"# time[yr] mean_vel std_vel min_vel max_vel mean_slip\n");
    if(uc->track_events){
      uc->fout_event = fopen("rsf_events.dat","w");
      fprintf(uc->fout_event,"# time[s] time[yr] onset(1)/arrest(-1) log10(max|v|[m/s]) mean_slip[m] mean_mu, |v| threshold %.3e m/s\n",
	      uc->vel_event);
    }
  }
  PetscCall(VecDuplicate(X0,&uc->Xold));
  PetscCall(VecCopy(X0,uc->Xold));
  /* gather context for the full-field output */
  PetscCall(VecScatterCreateToZero(X0,&uc->gather,&uc->gathered));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* close output files and free work space, cf. ode_solve_test.c */
static PetscErrorCode rsf_finalize_monitor_and_event(struct rsf_out_ctx *uc)
{
  struct med *medium;
  PetscFunctionBeginUser;
  medium = uc->par->medium;
  if(medium->comm_rank == 0){
    if(uc->fout_monitor)fclose(uc->fout_monitor);
    if(uc->fout_stats)fclose(uc->fout_stats);
    if(uc->fout_event){
      fclose(uc->fout_event);
      fprintf(stderr,"rsf_finalize_monitor_and_event: tracked %i events (|v| through %.3e m/s)\n",
	      uc->nevent,uc->vel_event);
    }
  }
  PetscCall(VecDestroy(&uc->Xold));
  PetscCall(VecScatterDestroy(&uc->gather));
  PetscCall(VecDestroy(&uc->gathered));
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
