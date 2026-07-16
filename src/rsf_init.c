#include "interact.h"
#include "properties.h"
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#include "rsf.h"

void init_medium_rsf(struct med *medium)
{
  struct rsf_vars *rsf;
  /* the rsf parameter block hangs off medium and is allocated here,
     before any field is set; medium itself was calloc-ed in main */
  medium->rsf = (struct rsf_vars *)calloc(1,sizeof(struct rsf_vars));
  if(!medium->rsf){
    fprintf(stderr,"init_medium_rsf: rsf_vars allocation failed\n");
    exit(-1);
  }
  rsf = medium->rsf;
  /* 
     default frictional and loading parameters, can be overridden via
     options below; defaults follow the SEAS BP1 benchmark / HBI
     conventions where applicable
  */
  rsf->f0  = 0.6; 		/* f_0 reference friction */
  rsf->dc = 0.008;		/* D_c [m] */
  rsf->vpl = 1e-9;		/* plate motion [m/s] */
  rsf->v0 = 1e-6;		/* reference speed [m/s] (HBI: vref) */
  /*  */
  rsf->calc_sigma_dot = FALSE; /* compute normal stress change? */
  /* optional per-cell D_c [m] (geometry order); NULL => use uniform
   rsf->dc.  Used for the BP5 nucleation patch (reduced D_RS) and
   read in the RHS, so it lives at file scope alongside the statics above */
  rsf->dc_vec = NULL;
  /* optional per-cell initial normal stress sigma0 [Pa] (geometry order);
     NULL => uniform sigma_init.  Read in the driver alongside dc_vec */
  rsf->sigma_vec = NULL;
  /* slip direction: 0 = strike (default, unchanged), 1 = dip (thrust / normal
     fault).  Parsed below as -rsf_slip_mode */
  rsf->slip_mode = STRIKE;
  /* 
   parameters for optional normal stress limiter, cf. HBI limitsigma
   (file scope since they are only used by the driver and RHS here) 
  */
  rsf->limit_sigma = PETSC_FALSE;
  rsf->state_law = RSF_AGING_LAW;  /* 1 aging (default), 2 slip law 3 PRZ 4: Sato 5: Tullis */

  rsf->vmin_state = 1e-16;         /* |v| floor in the slip law's ln(|v|/v0) */

  /* fixed constants of the two gated laws, at the values of the MATLAB reference
     implementation (Omega threshold 0.01, and Vc = v0/100). Not exposed as
     options for now; edit here if they need to change */
  rsf->sato_beta = 1e-2;
  rsf->kt_vc     = -1.0;		/* set from v0 below, after -v0 has been read */
  rsf->min_sigma = 1e6;
  rsf->max_sigma = 300e6; /* Pa */

  rsf->dim = 4;/* psi,tau,sigma,u(slip) per patch */
}
/*
  print a grouped summary of the rsf_solve command-line options and a few
  relevant PETSc options, then the caller exits.  triggered by -h.  kept here
  next to rsf_get_settings so the two stay in step: when an option is added to
  the parser above, add its line here too.  defaults shown match the values set
  in init_medium_rsf and rsf_get_settings.
*/
void rsf_print_help(const char *prog)
{
  fprintf(stderr,"\n");
  fprintf(stderr,"rsf_solve: quasi-dynamic rate-and-state earthquake-cycle solver (PETSc TS)\n");
  fprintf(stderr,"usage: %s [options]   (options may also be placed in petsc_settings.yaml)\n",prog);
  fprintf(stderr,"\n");

  fprintf(stderr,"geometry and input files\n");
  fprintf(stderr,"  -geom_file <file>       fault geometry, one patch per row (default geom.in)\n");
  fprintf(stderr,"  -rsf_file <file>        per-cell a b friction parameters (default rsf.dat)\n");
  fprintf(stderr,"  -rsf_ic_file <file>     per-cell initial tau[Pa] vel[m/s]; overrides uniform IC\n");
  fprintf(stderr,"  -rsf_dc_file <file>     per-cell D_c[m]; overrides uniform -dc\n");
  fprintf(stderr,"  -rsf_sigma_file <file>  per-cell initial sigma0[Pa]; overrides uniform -sigma_init\n");
  fprintf(stderr,"  -full_space <bool>      whole-space Green functions (default 0 = half space)\n");
  fprintf(stderr,"  -tv <int>               triangular patch evaluation mode (default 0)\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"elastic and rate-and-state parameters\n");
  fprintf(stderr,"  -shear_modulus <Pa>     shear modulus G (default 32.04e9)\n");
  fprintf(stderr,"  -s_wave_speed <m/s>     shear wave speed c_s, sets radiation damping (default 3464)\n");
  fprintf(stderr,"  -f0 <val>               reference friction f0 (default 0.6)\n");
  fprintf(stderr,"  -dc <m>                 characteristic slip distance D_c (default 0.008)\n");
  fprintf(stderr,"  -v0 <m/s>               reference velocity (default 1e-6)\n");
  fprintf(stderr,"  -vpl <m/s>              plate loading rate (default 1e-9)\n");
  fprintf(stderr,"  -rsf_slip_mode <0|1>    slip direction: 0 strike (default), 1 dip (thrust/normal)\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"initial conditions\n");
  fprintf(stderr,"  -sigma_init <Pa>        uniform initial normal stress (default 50e6)\n");
  fprintf(stderr,"  -tau_init <Pa>          uniform initial shear stress (default f0*sigma_init + eta*vel_init)\n");
  fprintf(stderr,"  -vel_init <m/s>         uniform initial slip rate (default vpl)\n");
  fprintf(stderr,"  -rand_amp <val>         random initial-state multiplier amplitude (default 0)\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"normal-stress evolution and limiter (for dip slip / nonplanar faults)\n");
  fprintf(stderr,"  -calc_sigma_dot <bool>  evolve normal stress via the In matrix (default 0 = off)\n");
  fprintf(stderr,"  -limit_sigma <bool>     clamp sigma to [min,max], as HBI limitsigma (default 0)\n");
  fprintf(stderr,"  -min_sigma <Pa>         limiter floor (default 1e6)\n");
  fprintf(stderr,"  -max_sigma <Pa>         limiter ceiling (default 300e6)\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"state evolution\n");
  fprintf(stderr,"  -state_law n            state evolution law: %i aging (default), %i slip, %i PRZ,\n",
	 RSF_AGING_LAW,RSF_SLIP_LAW,RSF_PRZ_LAW);
  fprintf(stderr,"                          %i Sato-type, %i Kato and Tullis composite\n",RSF_SATO_LAW,RSF_KT_LAW);
  fprintf(stderr,"                          aging: d psi/dt = b/dc (v0 exp((f0-psi)/b) - |v|)\n");
  fprintf(stderr,"                          slip:  d psi/dt = -(|v|/dc) (psi - psi_ss)\n");
  fprintf(stderr,"                          with psi_ss = f0 - b ln(|v|/v0); both share this\n");
  fprintf(stderr,"                          steady state, i.e. f_ss = f0 + (a-b) ln(|v|/v0)\n");
  fprintf(stderr,"                          PRZ (Perrin, Rice and Zheng 1995), theta form\n");
  fprintf(stderr,"                          d theta/dt = 1/2 (1 - (|v| theta/dc)^2):\n");
  fprintf(stderr,"                          d psi/dt = b/(2 dc) (v0 exp((f0-psi)/b)\n");
  fprintf(stderr,"                                              - v^2/v0 exp((psi-f0)/b))\n");
  fprintf(stderr,"                          aging, slip and PRZ share psi_ss = f0 - b ln(|v|/v0)\n");
  fprintf(stderr,"                          and the same linearization: comparable at equal dc\n");
  fprintf(stderr,"                          Sato and Kato-Tullis are the slip law plus a\n");
  fprintf(stderr,"                          gated aging healing term, d theta/dt = gate\n");
  fprintf(stderr,"                          - Omega ln Omega, Omega = |v| theta/dc, with\n");
  fprintf(stderr,"                          gate = exp(-Omega/beta), beta = 1e-2      (Sato)\n");
  fprintf(stderr,"                          gate = exp(-|v|/Vc),     Vc   = v0/100    (Kato-Tullis)\n");
  fprintf(stderr,"                          beta and Vc are fixed in rsf_init.c for now. NOTE the\n");
  fprintf(stderr,"                          Kato-Tullis gate does not vanish at |v| << Vc, which\n");
  fprintf(stderr,"                          raises psi_ss there by b*W(1) = 0.567 b (its known\n");
  fprintf(stderr,"                          steady-state offset); Sato has no such shift\n");
  fprintf(stderr,"  -vmin_state <m/s>       |v| floor in the slip law's ln(|v|/v0) (default 1e-16)\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"time stepping\n");
  fprintf(stderr,"  -rtol <val>             ODE relative tolerance (default 1e-4)\n");
  fprintf(stderr,"  -atol_slip <m>          absolute tolerance for slip entries (default 1e-3)\n");
  fprintf(stderr,"  -dt_init <s>            initial step size (default 1)\n");
  fprintf(stderr,"  -dt_max <s>             maximum step size (default 1e10)\n");
  fprintf(stderr,"  -imex                   IMEX (ARKIMEX) integration: local state terms implicit,\n");
  fprintf(stderr,"                          stress transfer explicit; helps stiff laws (e.g. PRZ)\n");
  fprintf(stderr,"                          during coseismic phases (default off = explicit RK).\n");
  fprintf(stderr,"                          The stage solvers use the options prefix imex_ (e.g.\n");
  fprintf(stderr,"                          -imex_ksp_type); unprefixed ksp/pc/snes options,\n");
  fprintf(stderr,"                          including those from petsc_settings.yaml, do not\n");
  fprintf(stderr,"                          affect them.\n");
  fprintf(stderr,"  -stop_time_yr <yr>      integration stop time (default 3000)\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"monitor and event detection\n");
  fprintf(stderr,"  -print_interval_yr <yr> averaged-property output cadence (default 0.1)\n");
  fprintf(stderr,"  -dt_monitor_yr <yr>     monitor cadence, also caps the step (default 5)\n");
  fprintf(stderr,"  -rdx_monitor <val>      relative state-change monitor trigger (default 1e-4)\n");
  fprintf(stderr,"  -adx_monitor <val>      absolute state-change trigger, <=0 off (default 0)\n");
  fprintf(stderr,"  -monitor_tmin_yr <yr>   suppress monitor output before this time (default 0)\n");
  fprintf(stderr,"  -track_events <bool>    locate slip-rate threshold crossings (default 1 = on)\n");
  fprintf(stderr,"  -vel_event <m/s>        event onset threshold (default 1e-3)\n");
  fprintf(stderr,"  -vel_event_hyst <val>   arrest at vel_event*hyst, debounces (default 0.5)\n");
  fprintf(stderr,"  -event_tmin_yr <yr>     suppress event output before this time (default 0)\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"optional outputs (all default off)\n");
  fprintf(stderr,"  -rsf_catalog            write rsf_catalog.dat (per-event slip, drop, M0, Mw)\n");
  fprintf(stderr,"  -rsf_rupture_time       write rsf_rupture_time.dat (first event front times)\n");
  fprintf(stderr,"  -slip_budget            write rsf_slip_budget.dat (slip vs plate-rate reference)\n");
  fprintf(stderr,"  -rupture_vth <m/s>      rupture-front threshold for the two above (default vel_event)\n");
  fprintf(stderr,"  -field_step_interval <n> slip-rate field frame every n accepted steps (0 = off)\n");
  fprintf(stderr,"  -field_tmin_yr <yr>     suppress field frames before this time (default 0)\n");
  fprintf(stderr,"  -slip_line_dt_yr <yr>   tmp_rsf slip-line snapshot cadence (default off)\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"interaction-matrix backend\n");
  fprintf(stderr,"  -use_hmatrix <0..5>     0 dense (default) 1 HTOOL 2 H2OPUS 3 HACApK 4 hmmvp 5 BigWham\n");
  fprintf(stderr,"  -hacapk_ztol <val>      HACApK compression tolerance\n");
  fprintf(stderr,"  -hacapk_eta <val>       HACApK admissibility eta\n");
  fprintf(stderr,"  -hacapk_inorm <int>     HACApK tolerance norm (1 = block-local)\n");
  fprintf(stderr,"  -hmmvp_tol <val>        hmmvp compression tolerance\n");
  fprintf(stderr,"  -hmmvp_eta <val>        hmmvp admissibility eta\n");
  fprintf(stderr,"  -hmmvp_inorm <int>      hmmvp tolerance norm\n");
  fprintf(stderr,"  -hmmvp_nthreads <int>   hmmvp compression threads\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"relevant PETSc options (a full list is printed by -help)\n");
  fprintf(stderr,"  -ts_rk_type <3bs|5dp|...>   Runge-Kutta variant (rsf uses an explicit RK)\n");
  fprintf(stderr,"  -ts_adapt_type <basic|none> step-size adaptation\n");
  fprintf(stderr,"  -ts_max_steps <n>           cap the number of accepted steps\n");
  fprintf(stderr,"  -ts_monitor                 print time and step at every accepted step\n");
  fprintf(stderr,"  -log_view                   PETSc performance and timing summary at the end\n");
  fprintf(stderr,"  -mat_htool_compressor <SVD|fullACA|partialACA>  HTOOL compressor choice\n");
  fprintf(stderr,"  -options_file <file>        read options from a file (petsc_settings.yaml auto-read)\n");
  fprintf(stderr,"  -help                       PETSc's full registered-option dump\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"output files\n");
  fprintf(stderr,"  rsf_monitor.dat         time series on the adaptive monitor cadence (-dt_monitor,\n");
  fprintf(stderr,"                          -adx_monitor, -rdx_monitor): step, time[s], time[yr], dt[s],\n");
  fprintf(stderr,"                          log10(max|v|), mean_slip, mean_mu, max_sigma, min_sigma.\n");
  fprintf(stderr,"                          Flushed as it is written\n");
  fprintf(stderr,"  rsf_vel.times           one row per field frame (-field_step_interval): frame, step,\n");
  fprintf(stderr,"                          time[yr], time[s], log10(max|v|), mean|v|, std|v|, min|v|,\n");
  fprintf(stderr,"                          mean_slip. The |v| statistics are slip SPEEDS (v is signed).\n");
  fprintf(stderr,"                          This is the index for the tmp_rsf/rsf_vel.gGGG.NNNNNN.bin\n");
  fprintf(stderr,"                          frames and replaces the former rsf_stats.dat\n");
  fprintf(stderr,"  rsf_events.dat          one row per slip-rate threshold crossing (-vel_event)\n");
  fprintf(stderr,"  rsf_catalog.dat         one row per completed event (with -rsf_catalog)\n");
  fprintf(stderr,"  rsf_geom.gGGG.dat       per-group patch geometry for the field frames\n");
  fprintf(stderr,"\n");

  fprintf(stderr,"examples\n");
  fprintf(stderr,"  strike-slip BP5 at 2 km, dense, with an event catalog:\n");
  fprintf(stderr,"    %s -geom_file geom_bp5_2km.in -rsf_file rsf_bp5_2km.dat \\\n",prog);
  fprintf(stderr,"       -rsf_ic_file ic_bp5_2km.in -rsf_dc_file dc_bp5_2km.in \\\n");
  fprintf(stderr,"       -sigma_init 25e6 -dc 0.14 -stop_time_yr 240 -rsf_catalog\n");
  fprintf(stderr,"  dip-slip thrust with evolving normal stress, hmmvp backend:\n");
  fprintf(stderr,"    %s -geom_file geom_thrust.in -rsf_file rsf_thrust.dat -rsf_ic_file ic_thrust.in \\\n",prog);
  fprintf(stderr,"       -sigma_init 58e6 -dc 0.02 -rsf_slip_mode 1 -calc_sigma_dot -limit_sigma \\\n");
  fprintf(stderr,"       -use_hmatrix 4 -hmmvp_tol 1e-4 -rsf_catalog\n");
  fprintf(stderr,"\n");
}

/*
   gather all run settings: defaults, then PetscOptions overrides (which
   also pick up the optional petsc_settings.yaml read in main).  Fills
   struct rsf_solve_settings and the medium fields that have a home there, so
   the solver routine sees a fully prepared configuration
*/
PetscErrorCode rsf_get_settings(int argc,char **argv,struct interact_ctx *par,
				struct rsf_solve_settings *set)
{
  struct med *medium = par->medium;
  const PetscReal sec_per_year = 365.25*24.*60.*60.;
  /* material parameters */
  PetscReal shear_modulus_si = 32.04e9, s_wave_speed_si = 3.464e3;
  PetscBool use_full_space=PETSC_FALSE;
  PetscInt tvmode = 0;		/* triangle evaluation mode */
  PetscInt use_hmatrix=IHMAT_TYPE_DENSE;
  PetscReal sigma_init,tau_init,vel_init,rtol,atol_slip,dt_init,dt_max,rand_amp,tmp;
  PetscReal dt_monitor,rdx_monitor,adx_monitor,monitor_tmin,vel_event,vel_event_hyst,event_tmin;
  PetscBool track_events;
  char geom_file[STRLEN]="geom.in",rsf_file[STRLEN]="rsf.dat",rsf_ic_file[STRLEN]="",
    rsf_dc_file[STRLEN]="",rsf_sigma_file[STRLEN]="";
  PetscBool have_ic=PETSC_FALSE,have_dc=PETSC_FALSE,have_sigma=PETSC_FALSE;
  PetscBool read_value;
  struct rsf_vars *rsf;
  PetscInt field_step_interval=0;PetscReal field_tmin_yr=0.0;PetscBool fset=PETSC_FALSE;
  /* output defaults */
  PetscBool cat_enable=PETSC_FALSE,rup_enable=PETSC_FALSE,budget_enable=PETSC_FALSE;
  PetscReal rupture_vth=-1.0;
  PetscFunctionBeginUser;
  /* -h prints the grouped option summary and exits (before any setup).  PETSc's
     own -help still works and dumps the full registered-option list. */
  {
    PetscBool help=PETSC_FALSE;
    PetscCall(PetscOptionsHasName(NULL,NULL,"-h",&help));
    if(help){
      int rank=0;
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      if(rank==0)
	rsf_print_help(argv[0]);
      fflush(stdout);
      PetscCall(PetscFinalize());
      exit(0);
    }
  }
  /* 
     default frictional and loading parameters, can be overridden via
     options below; defaults follow the SEAS BP1 benchmark / HBI
     conventions where applicable
  */
  init_medium_rsf(medium);
  rsf = medium->rsf;
  
  /* initial conditions, can overwrite below */
  sigma_init = 50e6;		/* [Pa] */
  tau_init = -1;		/* [Pa], if < 0, will use f0 * sigma_init */
  vel_init = rsf->vpl;	        /* [m/s] */
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
  /* 
     field output settings 
  */
  field_step_interval=0;
  field_tmin_yr=0.0;
  
  /* 
     H matrix

     options for this code: 
     0 dense 1 htools 2 h2opus 
     3 hacapk 4 hmmvp 5 BIGWHAM, 
     see petsc_prototypes.h 
  */
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-use_hmatrix", &use_hmatrix,&read_value));
  /* HTOOL (use_hmatrix==IHMAT_TYPE_HTOOLS ) 
     compressor default: prefer sympartialACA over
     PETSc's built-in SVD default.  On the SEAS BP5 4000-cell test (1 core),
     the SVD compressor spent ~48 s assembling the H-matrix vs ~6.5 s for
     HACApK and ~8 s for sympartialACA, for no measurable accuracy gain
     (recurrence and the full max|V| trace agree with the dense solution to
     <0.01 yr and <0.005 in log10|V| at epsilon<=1e-4 -- the error floor is
     set by the ODE rtol, not the H-matrix tolerance).  Override on the
     command line with -mat_htool_compressor {SVD,fullACA,partialACA}. See
     rsf_solve.md for the full HTOOL/HACApK/dense performance comparison. */
  medium->use_hmatrix = use_hmatrix;
  if(medium->use_hmatrix)
    set_hmat_defaults_and_options(medium,medium->use_hmatrix);

  
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-full_space", &use_full_space,&read_value)); /* use read in or default */
  medium->full_space = (my_boolean)use_full_space;
  PetscCall(PetscOptionsGetString(NULL, NULL, "-rsf_file", rsf_file, STRLEN,&read_value));
  /* physical parameters */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-shear_modulus",&shear_modulus_si,NULL)); /* G [Pa] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-s_wave_speed",&s_wave_speed_si,NULL));   /* c_s [m/s] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-f0",&rsf->f0,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-dc",&rsf->dc,NULL));     /* [m] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-vpl",&rsf->vpl,NULL));   /* [m/s] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-v0",&rsf->v0,NULL));     /* reference velocity [m/s] */
  /* G/(2 c_s) radiation damping factor */
  rsf->shear_mod_over_2cs_si = shear_modulus_si /(2.0*s_wave_speed_si);
 
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
  /* optional per-cell initial normal stress file (one sigma0[Pa] per patch,
     geometry order); overrides uniform -sigma_init per cell.  Natural companion
     to -calc_sigma_dot for a dipping or nonplanar fault where sigma0 varies */
  PetscCall(PetscOptionsGetString(NULL,NULL,"-rsf_sigma_file",rsf_sigma_file,STRLEN,&have_sigma));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-rand_amp",&rand_amp,NULL));
  if(tau_init < 0){
    /* default uniform prestress from the quasi-dynamic strength balance,
       tau = sigma f + eta v, with f = f0 at v = v0. The radiation-damping term
       eta*vel_init is added here for consistency; it is tiny at plate rate
       (eta*vpl ~ 5e-3 Pa) but matters if vel_init is set to a seismic rate */
    tau_init = rsf->f0 * sigma_init + rsf->shear_mod_over_2cs_si * vel_init;
  }
  /* evolve normal stress (dsigma/dt from the In interaction matrix)?  Off by
     default; needed for nonplanar or dipping faults where slip changes the
     normal stress.  Previously only settable in code */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-calc_sigma_dot",&rsf->calc_sigma_dot,NULL));
  /* slip direction for the rate-and-state solve: 0 strike (default), 1 dip.
     Dip slip is what a thrust or normal fault needs and is the case that
     exercises the normal-stress path meaningfully */
  {
    PetscInt sm = rsf->slip_mode;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-rsf_slip_mode",&sm,NULL));
    if((sm != STRIKE) && (sm != DIP)){
      fprintf(stderr,"rsf_get_settings: -rsf_slip_mode must be %i (strike) or %i (dip); got %i\n",
	      STRIKE,DIP,(int)sm);
      exit(-1);
    }
    rsf->slip_mode = (int)sm;
  }
  /* normal stress limiter as in HBI's limitsigma (default off here,
     irrelevant for planar faults) */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-limit_sigma",&rsf->limit_sigma,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-min_sigma",&rsf->min_sigma,NULL)); /* [Pa] */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-max_sigma",&rsf->max_sigma,NULL)); /* [Pa] */
  /* state evolution law: 1,2,3,4,5, 1 default ageing */
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-state_law",&rsf->state_law,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-vmin_state",&rsf->vmin_state,NULL)); /* [m/s] */
  if((rsf->state_law < RSF_AGING_LAW) || (rsf->state_law > RSF_KT_LAW)){
    if(medium->comm_rank == 0)
      fprintf(stderr,"rsf_get_settings: -state_law has to be %i (aging), %i (slip), %i (PRZ), %i (Sato) or %i (Kato-Tullis), not %i\n",
	      RSF_AGING_LAW,RSF_SLIP_LAW,RSF_PRZ_LAW,RSF_SATO_LAW,RSF_KT_LAW,(int)rsf->state_law);
    exit(-1);
  }
  /* Kato and Tullis Vc scales with the reference velocity, so it has to be set
     after -v0 has been parsed */
  if(rsf->kt_vc < 0)
    rsf->kt_vc = 1e-2 * rsf->v0;
  
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
  set->use_imex = PETSC_FALSE;
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-imex",&set->use_imex,NULL));
  /* other outputs */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-rsf_catalog",&cat_enable,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-rsf_rupture_time",&rup_enable,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-slip_budget",&budget_enable,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-rupture_vth",&rupture_vth,NULL));

  /* triangular patch evaluation scheme, cf. -tv of interact */
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-tv",&tvmode,NULL));
  medium->tri_eval_mode = (MODE_TYPE)tvmode;

  dt_monitor *= sec_per_year;monitor_tmin *= sec_per_year;event_tmin *= sec_per_year;
  /*
     optional compact slip-rate field output for later GMT
     visualization; off unless -field_step_interval > 0 is given.  A
     frame is written every field_step_interval accepted steps, so the
     cadence follows the solver's own step density (fine through
     nucleation and rupture, coarse through the interseismic) rather
     than model time.  -field_tmin_yr optionally suppresses early frames
  */
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-field_step_interval",&field_step_interval,&fset));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-field_tmin_yr",&field_tmin_yr,NULL));
  set->field_step_interval = (fset && (field_step_interval > 0))?(field_step_interval):(0);
  set->field_enable = (set->field_step_interval > 0)?(PETSC_TRUE):(PETSC_FALSE);
  set->field_tmin = field_tmin_yr * sec_per_year;
  /* cat output */
  set->cat_enable = cat_enable;
  set->rup_enable = rup_enable;
  set->budget_enable = budget_enable;
  set->rupture_vth = rupture_vth;
  /* 
     hand the gathered settings to the solver routine 
  */
  set->shear_modulus_si = shear_modulus_si;
  set->s_wave_speed_si  = s_wave_speed_si;
  set->sigma_init = sigma_init;
  set->tau_init   = tau_init;
  set->vel_init   = vel_init;
  set->rand_amp   = rand_amp;
  set->rtol       = rtol;
  set->atol_slip  = atol_slip;
  set->dt_init    = dt_init;
  set->dt_max     = dt_max;
  set->dt_monitor   = dt_monitor;
  set->rdx_monitor  = rdx_monitor;
  set->adx_monitor  = adx_monitor;
  set->monitor_tmin = monitor_tmin;
  set->vel_event      = vel_event;
  set->vel_event_hyst = vel_event_hyst;
  set->event_tmin     = event_tmin;
  set->track_events   = track_events;
  set->have_ic = have_ic;
  set->have_dc = have_dc;
  set->have_sigma = have_sigma;
  strncpy(set->geom_file,  geom_file,  STRLEN);
  strncpy(set->rsf_file,   rsf_file,   STRLEN);
  strncpy(set->rsf_ic_file,rsf_ic_file,STRLEN);
  strncpy(set->rsf_dc_file,rsf_dc_file,STRLEN);
  strncpy(set->rsf_sigma_file,rsf_sigma_file,STRLEN);
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif
