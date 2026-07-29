#ifndef __READ_RSF_HEADER__
/* 
   rate-and-state state evolution laws, for rsf->state_law (-state_law) 
*/
#define RSF_AGING_LAW 1
#define RSF_SLIP_LAW  2
#define RSF_PRZ_LAW   3
#define RSF_SATO_LAW  4
#define RSF_KT_LAW    5


/* 

   default values

*/
#define RSF_SHEAR_MODULUS_DEF 3.204e10
#define RSF_S_WAVE_SPEED_DEF 3464.0
#define RSF_DEFAULT_FRICTION 0.6
#define RSF_DEFAULT_DC 0.008
#define RSF_DEFAULT_VPL 1e-9
#define RSF_DEFAULT_V0 1e-6
#define RSF_DEFAULT_NORMAL_STRESS  5.0e7

/* 

   defaultt absolute error tolerance for ODE integration

*/
#define RSF_AERROR_PSI 1e-12 /* psi */
#define RSF_AERROR_TAU 1e-3 /* tau [Pa] */
#define RSF_AERROR_SIGMA 1e-3 /* sigma [Pa] */

/* 
   context for the monitor and event functions, following the layout
   in hmatrix_test/ode_solve_test.c; the RHS keeps using the plain
   interact_ctx which is shared with the matrix assembly routines
*/
/*
   one regular grid per fault group (the last column of the geometry
   input).  An n-group simulation produces n separate grids; each
   group's coordinates are measured in that group's own strike/dip
   frame, so faults with different orientations still grid cleanly
*/
struct rsf_group_grid{
  int id;			/* fault group id from the geometry */
  int np;			/* number of patches in this group */
  int *idx;			/* global patch indices, geometry order */
  double *xs,*ys;		/* on-fault coordinates [m], group frame, min-shifted */
  float *buf;			/* 3*np scratch for the xyz frame writes */
  double xmin,xmax,ymin,ymax,dx,dy;
  int nx,ny,regular;
};

struct rsf_out_ctx{
  struct interact_ctx *par;	/* fault geometry and medium */
  /* change/time triggered state monitor */
  PetscReal old_time,dt_monitor,adx_monitor,rdx_monitor,monitor_tmin;
  Vec Xold;
  FILE *fout_monitor;
  /*
     optional per fault group monitor (-rsf_monitor_by_group): the same
     columns as fout_monitor, on the same cadence, but reduced over the
     patches of one group of the geometry input instead of over the whole
     fault.  mgrp_n == 0 means the feature is off and none of the arrays
     below are allocated.  mgrp_map is rank local, length medium->rn, and
     holds the group slot of each owned patch; mgrp_np holds the GLOBAL
     patch count per group, so the per group means need no extra reduction
  */
  int mgrp_n;			/* number of distinct groups, 0: off */
  int *mgrp_id;			/* group ids, first seen order, length mgrp_n */
  int *mgrp_np;			/* global patches per group, length mgrp_n */
  int *mgrp_map;		/* owned patch -> group slot, length medium->rn */
  PetscReal *mgrp_lsum,*mgrp_gsum; /* 2*mgrp_n: sum slip, sum mu */
  PetscReal *mgrp_lmx,*mgrp_gmx;   /* 3*mgrp_n: max|v|, max sigma, -min sigma */
  FILE **mgrp_fout;		/* rank 0: one file per group, length mgrp_n */
  /* periodic full-field output (stats line, velocity snapshots) */
  PetscReal next_print_time;
  int field_out;
  VecScatter gather;
  Vec gathered;
  /* slip velocity threshold crossing events */
  PetscBool track_events,slipping;
  PetscReal vel_event,vel_event_hyst,event_tmin;
  int nevent;
  FILE *fout_event;
  /*
     compact per-fault slip-rate field output for later GMT
     visualization (see rsf_build_groups / rsf_TS_Monitor).  Geometry
     is written once per group to rsf_geom.gGGG.dat; the per-frame
     log10|v| field goes to tmp_rsf/rsf_vel.gGGG.NNNNNN.bin as
     xyz2grd-ready float triples, and rsf_vel.times indexes frames by
     step and time.  The cadence is a fixed number of accepted steps
     rather than model time, so frames densify automatically through
     nucleation and rupture (many small steps) and thin out through
     the long interseismic (few large steps)
  */
  PetscBool field_enable;
  PetscInt field_step_interval;	/* write a frame every this many accepted steps */
  PetscReal field_tmin;		/* [s] suppress field output before this model time */
  int field_frame;
  FILE *fout_field_times;
  struct rsf_group_grid *groups;	/* rank 0: one per fault group */
  int ngroup;
  double *vbuf;			/* rank 0: per-patch |v| scratch, length nrflt */
  /*
    Task 1: SEAS-style event catalog, rupture-time field, and slip-budget
    diagnostic.  All three are opt-in and default off; when off none of the
    arrays below are allocated and the monitor/event hooks are no-ops.  The
    per-cell arrays are rank-local, length medium->rn (owned patches), indexed
    by local k with global patch i = medium->rs + k.
  */
  /*
     checkpoint/restart: the full evolving state is the TS solution
     vector (psi,tau,sigma,slip per patch), so a checkpoint is that
     vector plus time/step/dt and validation metadata, written as a
     PETSc binary file (rank-count portable through VecLoad). Written
     every ckpt_every accepted steps (0: off) and once more after a
     regular TSSolve return; the previous checkpoint is kept as
     <file>.prev. Restart validates nrflt/dim and warns on law or
     slip-mode changes; output files are opened in append mode on
     restart with a marker comment. Note the in-progress event tracker
     is not checkpointed: an event spanning the restart is split, and
     -ts_max_steps counts absolute step numbers, so raise it when
     chaining runs.
  */
  PetscInt ckpt_every;
  char ckpt_file[300];
  PetscBool restarted;
  PetscInt ckpt_dim,ckpt_slip_mode,ckpt_law; /* stashed for the metadata */
  PetscBool cat_enable;		/* -rsf_catalog: write rsf_catalog.dat */
  PetscBool rup_enable;		/* -rsf_rupture_time: write rsf_rupture_time.dat (event 1) */
  PetscBool budget_enable;	/* -slip_budget: write rsf_slip_budget.dat */
  PetscReal rupture_vth;	/* [m/s] rupture-front threshold (default = vel_event) */
  PetscReal shear_modulus;	/* [Pa] G, for the seismic-moment sum */
  PetscReal vpl;		/* [m/s] plate rate, for the slip budget */
  PetscReal total_area;		/* [m^2] summed patch area, for the slip budget */
  FILE *fout_catalog;
  FILE *fout_budget;
  PetscReal *snap_tau0;		/* [Pa] onset shear-stress snapshot, length rn */
  PetscReal *snap_slip0;	/* [m]  onset slip snapshot, length rn */
  int       *cell_ruptured;	/* 1 if cell exceeded rupture_vth this event, length rn */
  PetscReal *rup_time;		/* [s] first-crossing time after event-1 onset, or -1; length rn */
  PetscReal peakv_local;	/* running max |v| this event on this rank [m/s] */
  PetscReal onset_time;		/* [s] onset time of the current event */
  PetscBool ev_open;		/* an onset snapshot is captured; event is open */
  int ncat;			/* completed catalog events (arrests written) */
  PetscBool rup_armed;		/* event-1 rupture-time capture is active */
  PetscBool rup_done;		/* event 1 finished; rup_time is frozen */
};
/*
   all run settings gathered up front by rsf_get_settings and consumed
   by rsf_solve_run, so that parameter handling and the solver proper
   live in separate routines; quantities that already have a natural
   home in struct med (f0, dc, vpl, v0, stop_time, print_interval,
   slip_line_dt, use_hmatrix, full_space, tri_eval_mode) are stored
   there and not duplicated here
*/

struct rsf_solve_settings{
  PetscReal shear_modulus_si,s_wave_speed_si;	    /* G [Pa], c_s [m/s] */
  PetscReal sigma_init,tau_init,vel_init,rand_amp;  /* initial conditions */
  PetscReal rtol,atol_slip,dt_init,dt_max;	    /* time stepping */
  PetscReal dt_monitor,rdx_monitor,adx_monitor,monitor_tmin; /* monitor */
  PetscReal vel_event,vel_event_hyst,event_tmin;    /* event tracking */
  PetscBool track_events;
  PetscInt field_step_interval;			    /* compact field output cadence [steps] */
  PetscReal field_tmin;				    /* [s] field output time floor */
  PetscBool field_enable;
  PetscBool monitor_by_group;	/* -rsf_monitor_by_group: per group monitor files */
  PetscBool have_ic,have_dc,have_sigma;
  PetscBool use_imex;		/* -imex: ARKIMEX with implicit local state terms (rsf_imex.c) */
  /* Task 1 outputs (all default off) */
  PetscBool cat_enable;		/* -rsf_catalog */
  PetscBool rup_enable;		/* -rsf_rupture_time */
  PetscBool budget_enable;	/* -slip_budget */
  PetscReal rupture_vth;	/* -rupture_vth [m/s]; <=0 means use vel_event */
  char geom_file[STRLEN],rsf_file[STRLEN],rsf_ic_file[STRLEN],rsf_dc_file[STRLEN],rsf_sigma_file[STRLEN];
};


PetscErrorCode rsf_get_settings(int,char **,struct interact_ctx *,struct rsf_solve_settings *);
PetscErrorCode rsf_solve_run(int,char **,struct interact_ctx *,struct rsf_solve_settings *);
int rsf_build_groups(struct med *,struct flt *,struct rsf_group_grid **);
void rsf_free_groups(struct rsf_group_grid *,int);
void rsf_write_group_geometry(const struct rsf_group_grid *,struct flt *,double,const char *);
PetscErrorCode rsf_finalize_event(struct rsf_out_ctx *,PetscReal,Vec);
PetscErrorCode rsf_init_monitor_and_event(struct rsf_out_ctx *,struct interact_ctx *,
						 PetscReal,PetscReal,PetscReal,PetscReal,
						 PetscReal,PetscReal,PetscReal,PetscBool,
						 PetscReal,Vec,PetscReal,
						 PetscBool,PetscInt,PetscReal,
						 struct rsf_group_grid *,int,double *,
						 PetscBool);
PetscErrorCode rsf_finalize_monitor_and_event(struct rsf_out_ctx *);
PetscErrorCode rsf_event_function(TS,PetscReal,Vec,PetscScalar *,void *);
PetscErrorCode rsf_post_event(TS,PetscInt,PetscInt[],PetscReal,Vec,PetscBool,void *);



/* context access for the domain check, whose callback has no user pointer */
extern struct interact_ctx *rsf_par_static;
#define __READ_RSF_HEADER__
#endif
