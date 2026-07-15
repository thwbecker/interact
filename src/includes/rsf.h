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
						 struct rsf_group_grid *,int,double *);
PetscErrorCode rsf_finalize_monitor_and_event(struct rsf_out_ctx *);
PetscErrorCode rsf_event_function(TS,PetscReal,Vec,PetscScalar *,void *);
PetscErrorCode rsf_post_event(TS,PetscInt,PetscInt[],PetscReal,Vec,PetscBool,void *);



/* context access for the domain check, whose callback has no user pointer */
extern struct interact_ctx *rsf_par_static;
