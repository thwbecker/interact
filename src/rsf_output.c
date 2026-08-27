#include "interact.h"
#include "properties.h"
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#include "rsf.h"



/*
   per fault group monitor output (-rsf_monitor_by_group).  Writes the
   same columns as RSF_MONITOR_FILE on the same cadence, but reduced over
   the patches of one group (the last column of the geometry input)
   instead of over the whole fault.  Distinct group ids are collected in
   first seen order, matching rsf_build_groups, so the numbering agrees
   with the one used by the compact field output.  The per group patch
   counts are read off the global fault array, which every rank holds, so
   the means need no extra reduction.  Returns the number of groups, or 0
   if the feature could not be set up, in which case the caller leaves it
   off.  All ranks take the same branch here, so a rank 0 file open
   failure disables only the writing, not the reductions.
*/
static int rsf_init_monitor_groups(struct rsf_out_ctx *uc)
{
  struct med *medium;
  struct flt *fault;
  int i,j,k,n,ng,id,found;
  char fname[STRLEN];

  medium = uc->par->medium;
  fault  = uc->par->fault;
  n      = medium->nrflt;
  if(n < 1)
    return 0;
  uc->mgrp_id = (int *)malloc((size_t)n*sizeof(int)); /* at most n distinct */
  if(!uc->mgrp_id)
    return 0;
  ng = 0;
  for(i=0;i < n;i++){		/* distinct group ids, first seen order */
    id = (int)fault[i].group;
    found = 0;
    for(j=0;j < ng;j++)
      if(uc->mgrp_id[j] == id){
	found = 1;
	break;
      }
    if(!found)
      uc->mgrp_id[ng++] = id;
  }
  uc->mgrp_np   = (int *)calloc((size_t)ng,sizeof(int));
  uc->mgrp_map  = (int *)malloc((size_t)medium->rn*sizeof(int));
  uc->mgrp_lsum = (PetscReal *)malloc((size_t)MGRP_NSUM*(size_t)ng*sizeof(PetscReal));
  uc->mgrp_gsum = (PetscReal *)malloc((size_t)MGRP_NSUM*(size_t)ng*sizeof(PetscReal));
  uc->mgrp_lmx  = (PetscReal *)malloc((size_t)3*(size_t)ng*sizeof(PetscReal));
  uc->mgrp_gmx  = (PetscReal *)malloc((size_t)3*(size_t)ng*sizeof(PetscReal));
  if((!uc->mgrp_np)||(!uc->mgrp_map)||(!uc->mgrp_lsum)||(!uc->mgrp_gsum)||
     (!uc->mgrp_lmx)||(!uc->mgrp_gmx)){
    fprintf(stderr,"rsf_init_monitor_groups: alloc failed for %i group(s), per group monitor off\n",ng);
    return 0;
  }
  for(i=0;i < n;i++)		/* global patch count per group */
    for(j=0;j < ng;j++)
      if(uc->mgrp_id[j] == (int)fault[i].group){
	uc->mgrp_np[j]++;
	break;
      }
  for(i=medium->rs,k=0;i < medium->re;i++,k++){ /* owned patch -> group slot */
    uc->mgrp_map[k] = 0;
    for(j=0;j < ng;j++)
      if(uc->mgrp_id[j] == (int)fault[i].group){
	uc->mgrp_map[k] = j;
	break;
      }
  }
  HEADNODE{
    uc->mgrp_fout = (FILE **)calloc((size_t)ng,sizeof(FILE *));
    if(!uc->mgrp_fout){
      fprintf(stderr,"rsf_init_monitor_groups: alloc failed for %i output file(s), no per group files written\n",ng);
    }else{
      for(j=0;j < ng;j++){
	snprintf(fname,STRLEN,RSF_MONITOR_GROUP_FORMAT,uc->mgrp_id[j]);
	uc->mgrp_fout[j] = myopen(fname,(uc->restarted)?("a"):("w"));
	if(uc->restarted)
	  fprintf(uc->mgrp_fout[j],"# restarted\n");
	fprintf(uc->mgrp_fout[j],"# fault group %i, %i of %i patches, cadence as in %s\n",
		uc->mgrp_id[j],uc->mgrp_np[j],n,RSF_MONITOR_FILE);
	fprintf(uc->mgrp_fout[j],"# step time[s] time[yr] dt[s] log10(max|v|[m/s]) mean_slip[m] mean_mu max_sigma[Pa] min_sigma[Pa]\n");
      }
      fprintf(stderr,"rsf_init_monitor_groups: per group monitor on, %i group(s), %s\n",
	      ng,RSF_MONITOR_GROUP_FORMAT);
    }
  }
  return ng;
}

/* 
   set up the monitor and event tracking environment,
   cf. init_monitor_and_event in ode_solve_test.c
*/
PetscErrorCode rsf_init_monitor_and_event(struct rsf_out_ctx *uc,struct interact_ctx *par,
					  PetscReal dt_monitor,PetscReal adx_monitor,
					  PetscReal rdx_monitor,PetscReal monitor_tmin,
					  PetscReal event_tmin,PetscReal vel_event,
					  PetscReal vel_event_hyst,PetscBool track_events,
					  PetscReal t_init,Vec X0,PetscReal vel_init,
					  PetscBool field_enable,PetscInt field_step_interval,
					  PetscReal field_tmin,struct rsf_group_grid *groups,
					  int ngroup,double *vbuf,
					  PetscBool monitor_by_group)
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
  /*
    Initialize the slipping state from the TRUE initial maximum slip rate over
    the fault, not from the scalar vel_init.  With a per-cell initial-condition
    file (e.g. the BP5 nucleation patch seeded at 3e-2 m/s) the scalar vel_init
    is only the background rate (vpl), so the old vel_init > vel_event test
    mislabeled the decaying seed as a spurious onset then arrest.  Computing the
    real field maximum here makes the seeded first event an in-progress event
    from t=0, which is what the catalog/rupture-time capture and rsf_events.dat
    both want.  This changes only the labeling of the initial seeded transient;
    the spontaneous recurrence times are unaffected.
  */
  {
    const PetscScalar *x0a;PetscInt ii,jj;PetscReal vv,d1,d2,d3,lv=0.0,gv=0.0;
    struct rsf_vars *rsfv = medium->rsf;struct flt *fltv = par->fault;
    PetscCall(VecGetArrayRead(X0,&x0a));
    for(ii=medium->rs,jj=0;ii<medium->re;ii++,jj+=rsfv->dim){
      vv = fabs(vel_from_rsf(x0a[jj+1],x0a[jj+2],x0a[jj],fltv[ii].mu_sa,rsfv->v0,&d1,&d2,&d3,medium));
      if(vv > lv)lv = vv;
    }
    PetscCall(VecRestoreArrayRead(X0,&x0a));
    PetscCallMPI(MPI_Allreduce(&lv,&gv,1,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD));
    uc->slipping = (gv > vel_event)?(PETSC_TRUE):(PETSC_FALSE);
  }
  uc->nevent = 0;
  uc->field_out = 0;
  uc->next_print_time = t_init;	/* fields at start, then every print_interval */
  /* compact per-fault slip-rate field output state */
  uc->field_enable = field_enable;
  uc->field_step_interval = field_step_interval;
  uc->field_tmin = field_tmin;
  uc->field_stress = PETSC_FALSE;	/* stress frames alongside the velocity
					   frames; same cadence, same frame ids */
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-field_stress",&uc->field_stress,NULL));
  uc->field_frame = 0;
  uc->fout_field_times = NULL;
  uc->groups = groups;
  uc->ngroup = ngroup;
  uc->vbuf = vbuf;
  uc->cat_enable = uc->rup_enable = uc->budget_enable = PETSC_FALSE;
  uc->rup_armed = uc->rup_done = PETSC_FALSE;
  uc->fout_catalog = uc->fout_budget = NULL;
  uc->snap_tau0 = uc->snap_slip0 = uc->rup_time = NULL;
  uc->cell_ruptured = NULL;
  uc->rupture_vth = vel_event;
  uc->shear_modulus = uc->vpl = uc->total_area = 0.0;
  uc->peakv_local = 0.0;
  uc->onset_time = 0.0;
  /* force the first monitor call to log */
  uc->old_time = t_init - 2.0*dt_monitor;
  uc->fout_monitor = uc->fout_event = NULL;
  /* per group monitor: all state off and unallocated unless asked for */
  uc->mgrp_n = 0;
  uc->mgrp_id = uc->mgrp_np = uc->mgrp_map = NULL;
  uc->mgrp_lsum = uc->mgrp_gsum = uc->mgrp_lmx = uc->mgrp_gmx = NULL;
  uc->mgrp_fout = NULL;
  if(monitor_by_group)
    uc->mgrp_n = rsf_init_monitor_groups(uc);
  HEADNODE{
    uc->fout_monitor = myopen(RSF_MONITOR_FILE,(uc->restarted)?("a"):("w"));
    if(uc->restarted) fprintf(uc->fout_monitor,"# restarted\n");
    fprintf(uc->fout_monitor,"# step time[s] time[yr] dt[s] log10(max|v|[m/s]) mean_slip[m] mean_mu max_sigma[Pa] min_sigma[Pa]\n");
    if(uc->track_events){
      uc->fout_event = myopen(RSF_EVENTS_FILE,(uc->restarted)?("a"):("w"));
      if(uc->restarted) fprintf(uc->fout_event,"# restarted\n");
      fprintf(uc->fout_event,"# time[s] time[yr] onset(1)/arrest(-1) log10(max|v|[m/s]) mean_slip[m] mean_mu, |v| threshold %.3e m/s\n",
	      uc->vel_event);
    }
    if(uc->field_enable){
      int ierr_dir = system("mkdir -p tmp_rsf");
      if(ierr_dir)
	fprintf(stderr,"rsf_init_monitor_and_event: WARNING: could not make tmp_rsf directory\n");
      uc->fout_field_times = myopen(RSF_VEL_TIME_FILE,(uc->restarted)?("a"):("w"));
      if(uc->fout_field_times){
	fprintf(uc->fout_field_times,"# frame step time[yr] time[s] log10(max|v|[m/s]) mean|v|[m/s] std|v|[m/s] min|v|[m/s] mean_slip[m]\n");
	fprintf(uc->fout_field_times,"# field per group in tmp_rsf/rsf_vel.gGGG.NNNNNN.bin (float32 along_strike,down_dip,log10|v| triples, xyz2grd -bi3f); geometry rsf_geom.gGGG.dat\n");
      }
    }
  }
  PetscCall(VecDuplicate(X0,&uc->Xold));
  PetscCall(VecCopy(X0,uc->Xold));
  /* gather context for the full-field output */
  PetscCall(VecScatterCreateToZero(X0,&uc->gather,&uc->gathered));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* close output files and free work space, cf. ode_solve_test.c */
PetscErrorCode rsf_finalize_monitor_and_event(struct rsf_out_ctx *uc)
{
  struct med *medium;
  int ig;
  PetscFunctionBeginUser;
  medium = uc->par->medium;
  HEADNODE{
    if(uc->fout_monitor)fclose(uc->fout_monitor);
    if(uc->fout_event){
      fclose(uc->fout_event);
      fprintf(stderr,"rsf_finalize_monitor_and_event: tracked %i events (|v| through %.3e m/s)\n",
	      uc->nevent,uc->vel_event);
    }
    if(uc->fout_field_times){
      fclose(uc->fout_field_times);
      fprintf(stderr,"rsf_finalize_monitor_and_event: wrote %i slip-rate field frame(s) across %i group(s) to tmp_rsf/\n",
	      uc->field_frame,uc->ngroup);
    }
    rsf_free_groups(uc->groups,uc->ngroup);
    uc->groups = NULL;uc->ngroup = 0;
    if(uc->vbuf){free(uc->vbuf);uc->vbuf = NULL;}
    if(uc->mgrp_fout){
      for(ig=0;ig < uc->mgrp_n;ig++)
	if(uc->mgrp_fout[ig])fclose(uc->mgrp_fout[ig]);
      free(uc->mgrp_fout);uc->mgrp_fout = NULL;
      fprintf(stderr,"rsf_finalize_monitor_and_event: wrote %i per group monitor file(s)\n",uc->mgrp_n);
    }
  }
  /* the reduction scratch lives on every rank */
  if(uc->mgrp_id){free(uc->mgrp_id);uc->mgrp_id = NULL;}
  if(uc->mgrp_np){free(uc->mgrp_np);uc->mgrp_np = NULL;}
  if(uc->mgrp_map){free(uc->mgrp_map);uc->mgrp_map = NULL;}
  if(uc->mgrp_lsum){free(uc->mgrp_lsum);uc->mgrp_lsum = NULL;}
  if(uc->mgrp_gsum){free(uc->mgrp_gsum);uc->mgrp_gsum = NULL;}
  if(uc->mgrp_lmx){free(uc->mgrp_lmx);uc->mgrp_lmx = NULL;}
  if(uc->mgrp_gmx){free(uc->mgrp_gmx);uc->mgrp_gmx = NULL;}
  uc->mgrp_n = 0;
  PetscCall(VecDestroy(&uc->Xold));
  PetscCall(VecScatterDestroy(&uc->gather));
  PetscCall(VecDestroy(&uc->gathered));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  SEAS-style event catalog, rupture-time field, and slip-budget
  diagnostic.  Called by the driver right after
  rsf_init_monitor_and_event.  Allocates the rank-local per-cell
  scratch only for the features that are on, opens the rank-0 output
  files, and stashes G, vpl, the rupture-front threshold, and the
  total patch area (needed for the moment sum and the slip budget) on
  the context.

*/

/*
   per fault group event catalogs (-rsf_catalog together with
   -rsf_monitor_by_group).  These do NOT re-detect events: they reuse the
   same global onset/arrest brackets that RSF_CATALOG_FILE is built on, and
   split each bracket by group.  That is complete rather than restrictive,
   because the bracket opens whenever the maximum |v| over the WHOLE fault
   crosses vel_event, so an episode carried by one splay alone opens a
   bracket just as a megathrust rupture does; what the split adds is which
   groups took part, how much each slipped, and when each of them went
   relative to the others.

   The group partition is the one rsf_init_monitor_groups already built, so
   this returns 0 unless that ran.  Ordering is guaranteed: rsf_solve calls
   rsf_init_monitor_and_event before rsf_init_catalog.

   Returns the number of groups, or 0 if the feature could not be set up,
   in which case the caller leaves it off.  All ranks take the same branch,
   so a rank 0 file open failure disables only the writing, not the
   reductions.
*/
static int rsf_init_catalog_groups(struct rsf_out_ctx *uc)
{
  struct med *medium;
  int ng,j;
  char fname[STRLEN];

  medium = uc->par->medium;
  ng = uc->mgrp_n;
  if(ng < 1)
    return 0;
  uc->cgrp_id   = (int *)malloc((size_t)ng*sizeof(int));
  uc->cgrp_peak = (PetscReal *)malloc((size_t)ng*sizeof(PetscReal));
  uc->cgrp_t0   = (PetscReal *)malloc((size_t)ng*sizeof(PetscReal));
  uc->cgrp_t1   = (PetscReal *)malloc((size_t)ng*sizeof(PetscReal));
  uc->cg_lsum   = (PetscReal *)malloc((size_t)7*(size_t)ng*sizeof(PetscReal));
  uc->cg_gsum   = (PetscReal *)malloc((size_t)7*(size_t)ng*sizeof(PetscReal));
  uc->cg_lmx    = (PetscReal *)malloc((size_t)4*(size_t)ng*sizeof(PetscReal));
  uc->cg_gmx    = (PetscReal *)malloc((size_t)4*(size_t)ng*sizeof(PetscReal));
  uc->cg_lmn    = (PetscReal *)malloc((size_t)ng*sizeof(PetscReal));
  uc->cg_gmn    = (PetscReal *)malloc((size_t)ng*sizeof(PetscReal));
  if((!uc->cgrp_id)||(!uc->cgrp_peak)||(!uc->cgrp_t0)||(!uc->cgrp_t1)||(!uc->cg_lsum)||
     (!uc->cg_gsum)||(!uc->cg_lmx)||(!uc->cg_gmx)||(!uc->cg_lmn)||(!uc->cg_gmn)){
    fprintf(stderr,"rsf_init_catalog_groups: alloc failed for %i group(s), per group catalog off\n",ng);
    return 0;
  }
  for(j=0;j < ng;j++){
    uc->cgrp_id[j] = uc->mgrp_id[j];
    uc->cgrp_peak[j] = 0.0;
    uc->cgrp_t0[j] = -1.0;
    uc->cgrp_t1[j] = -1.0;
  }
  HEADNODE{
    uc->cgrp_ncat = (int *)calloc((size_t)ng,sizeof(int));
    uc->cgrp_fout = (FILE **)calloc((size_t)ng,sizeof(FILE *));
    if((!uc->cgrp_ncat)||(!uc->cgrp_fout)){
      fprintf(stderr,"rsf_init_catalog_groups: alloc failed for %i output file(s), no per group files written\n",ng);
      if(uc->cgrp_ncat){free(uc->cgrp_ncat);uc->cgrp_ncat = NULL;}
      if(uc->cgrp_fout){free(uc->cgrp_fout);uc->cgrp_fout = NULL;}
    }else{
      for(j=0;j < ng;j++){
	snprintf(fname,STRLEN,RSF_CATALOG_GROUP_FORMAT,uc->cgrp_id[j]);
	uc->cgrp_fout[j] = myopen(fname,(uc->restarted)?("a"):("w"));
	if(uc->restarted)
	  fprintf(uc->cgrp_fout[j],"# restarted\n");
	fprintf(uc->cgrp_fout[j],"# fault group %i, %i of %i patches\n",
		uc->cgrp_id[j],uc->mgrp_np[j],medium->nrflt);
	fprintf(uc->cgrp_fout[j],
		"# one row per event in which this group ruptured; rows are only\n");
	fprintf(uc->cgrp_fout[j],
		"# written for groups that actually reached the rupture threshold,\n");
	fprintf(uc->cgrp_fout[j],
		"# so a quiet group simply has no row for that event.  The event\n");
	fprintf(uc->cgrp_fout[j],
		"# index and the onset/arrest times are the GLOBAL ones of %s,\n",RSF_CATALOG_FILE);
	fprintf(uc->cgrp_fout[j],
		"# so rows can be joined to it on column 1; t_first and t_last are\n");
	fprintf(uc->cgrp_fout[j],
		"# this group's own, and are resolved to the accepted time steps\n");
	fprintf(uc->cgrp_fout[j],
		"# rather than root found, unlike the global onset and arrest.\n");
	fprintf(uc->cgrp_fout[j],
		"# All quantities right of t_last are summed over THIS group only,\n");
	fprintf(uc->cgrp_fout[j],
		"# with the same definitions and units as %s.\n",RSF_CATALOG_FILE);
	fprintf(uc->cgrp_fout[j],
		"# ievent t_onset[yr] t_arrest[yr] t_first[yr] t_last[yr] ncell "
		"area[m^2] mean_slip[m] max_slip[m] mean_drop[MPa] max_drop[MPa] "
		"peak_v[m/s] M0[Nm] Mw\n");
      }
      fprintf(stderr,"rsf_init_catalog_groups: per group catalog on, %i group(s), %s\n",
	      ng,RSF_CATALOG_GROUP_FORMAT);
    }
  }
  return ng;
}

PetscErrorCode rsf_init_catalog(struct rsf_out_ctx *uc,struct interact_ctx *par,
				struct rsf_solve_settings *set,
				PetscReal shear_modulus_si,PetscReal vpl,
				Vec x0,PetscReal t_init)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  PetscInt i,rn = medium->rn;
  PetscReal larea,garea;
  PetscFunctionBeginUser;
  uc->cat_enable    = set->cat_enable;
  uc->rup_enable    = set->rup_enable;
  uc->budget_enable = set->budget_enable;
  uc->shear_modulus = shear_modulus_si;
  uc->vpl           = vpl;
  /* rupture-front threshold: default to the event threshold if not set */
  uc->rupture_vth   = (set->rupture_vth > 0.0)?(set->rupture_vth):(uc->vel_event);
  uc->peakv_local   = 0.0;
  uc->onset_time    = t_init;
  uc->ev_open       = PETSC_FALSE;
  uc->ncat          = 0;
  uc->rup_armed     = PETSC_FALSE;
  uc->rup_done      = PETSC_FALSE;
  uc->snap_tau0 = uc->snap_slip0 = uc->rup_time = NULL;
  uc->cell_ruptured = NULL;
  uc->fout_catalog = uc->fout_budget = NULL;
  uc->cgrp_n = 0;
  uc->cgrp_fout = NULL; uc->cgrp_ncat = NULL; uc->cgrp_id = NULL;
  uc->cgrp_peak = uc->cgrp_t0 = uc->cgrp_t1 = NULL;
  uc->cg_lsum = uc->cg_gsum = uc->cg_lmx = uc->cg_gmx = NULL;
  uc->cg_lmn = uc->cg_gmn = NULL;
  /*
    The catalog and the rupture-time field are built on the SAME onset/arrest
    crossings that the event tracker (rsf_events.dat) already detects, so they
    require -track_events (on by default).  They do not re-detect events; they
    only add the per-event quantities that rsf_events.dat does not carry (the
    coseismic slip and stress drop from the onset-to-arrest state difference,
    the peak slip rate over the event, the ruptured area and seismic moment,
    and the rupture-time field).  If event tracking is off, disable them.
  */
  if((uc->cat_enable || uc->rup_enable) && (!uc->track_events)){
    HEADNODE
      fprintf(stderr,"rsf_init_catalog: WARNING: -rsf_catalog/-rsf_rupture_time need -track_events; disabling them\n");
    uc->cat_enable = uc->rup_enable = PETSC_FALSE;
  }
  /* total patch area = sum over ALL patches of 4*l*w (l,w are half-lengths);
     summed locally then reduced so it is correct in parallel */
  larea = 0.0;
  for(i = medium->rs;i < medium->re;i++)
    larea += 4.0*fault[i].l*fault[i].w;
  PetscCallMPI(MPI_Allreduce(&larea,&garea,1,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD));
  uc->total_area = garea;
  /* per-cell rank-local scratch for the catalog and the rupture-time field */
  if(uc->cat_enable || uc->rup_enable){
    if(rn > 0){
      uc->cell_ruptured = (int *)      calloc((size_t)rn,sizeof(int));
      uc->rup_time      = (PetscReal *)malloc((size_t)rn*sizeof(PetscReal));
      if((!uc->cell_ruptured)||(!uc->rup_time)){
	fprintf(stderr,"rsf_init_catalog: per-cell scratch alloc failed (rn=%i)\n",(int)rn);
	exit(-1);
      }
      for(i=0;i < rn;i++)uc->rup_time[i] = -1.0; /* sentinel: never ruptured */
    }
    if(uc->cat_enable && (rn > 0)){
      uc->snap_tau0  = (PetscReal *)malloc((size_t)rn*sizeof(PetscReal));
      uc->snap_slip0 = (PetscReal *)malloc((size_t)rn*sizeof(PetscReal));
      if((!uc->snap_tau0)||(!uc->snap_slip0)){
	fprintf(stderr,"rsf_init_catalog: snapshot alloc failed (rn=%i)\n",(int)rn);
	exit(-1);
      }
    }
  }
  HEADNODE{
    if(uc->cat_enable){
      uc->fout_catalog = myopen(RSF_CATALOG_FILE,(uc->restarted)?("a"):("w"));
      if(uc->restarted) fprintf(uc->fout_catalog,"# restarted\n");
      fprintf(uc->fout_catalog,
	      "# SEAS-style event catalog from rsf_solve; complements %s\n",RSF_EVENTS_FILE);
      fprintf(uc->fout_catalog,
	      "# same onset/arrest events as %s (needs -track_events); rows are per completed event\n", RSF_EVENTS_FILE);
      fprintf(uc->fout_catalog,
	      "# onset/arrest |v| threshold = %.3e m/s ; rupture-front |v| threshold = %.3e m/s\n",
	      uc->vel_event,uc->rupture_vth);
      fprintf(uc->fout_catalog,
	      "# stress drop = |tau_onset| - |tau_arrest| (positive = drop); slip enters as |dslip|; means are over ruptured cells\n");
      fprintf(uc->fout_catalog,
	      "# for 2D geometries, M0 and Mw are equivalent values assuming along-strike length = ruptured in-plane length\n");
      fprintf(uc->fout_catalog,
	      "# ev onset[yr] arrest[yr] duration[s] n_ruptured area_ruptured[m^2] mean_slip[m] max_slip[m] mean_drop[MPa] max_drop[MPa] peak_sliprate[m/s] M0[Nm] Mw\n");
    }
    if(uc->budget_enable){
      uc->fout_budget = myopen(RSF_SLIP_BUDGET_FILE,(uc->restarted)?("a"):("w"));
      if(uc->restarted) fprintf(uc->fout_budget,"# restarted\n");
      fprintf(uc->fout_budget,
	      "# long-term slip budget: on-fault area-integrated slip vs plate-rate reference\n");
      fprintf(uc->fout_budget,
	      "# total area = %.6e m^2 ; vpl = %.3e m/s\n",uc->total_area,uc->vpl);
      fprintf(uc->fout_budget,
	      "# time[yr] slip_integral[m*m^2] vpl_t_area[m*m^2] residual[m*m^2] residual_over_ref\n");
    }
    if(uc->cat_enable)
      fprintf(stderr,"rsf_init_catalog: event catalog on (rupture-front |v| = %.3e m/s)\n",uc->rupture_vth);
    if(uc->rup_enable)
      fprintf(stderr,"rsf_init_catalog: rupture-time field on (event 1, |v| = %.3e m/s)\n",uc->rupture_vth);
    if(uc->budget_enable)
      fprintf(stderr,"rsf_init_catalog: slip-budget diagnostic on\n");
  }
  /* split the catalog by fault group when the per group monitor is on.
     Not gated on HEADNODE: the reductions below are collective */
  if(uc->cat_enable)
    uc->cgrp_n = rsf_init_catalog_groups(uc);
  /*
    Seeded-start case: if the fault is already above the event threshold at
    t_init (e.g. the BP5 nucleation patch), there is no onset crossing for the
    first event, so open it here from the initial state and arm the rupture-time
    capture, so that event 1 is the seeded rupture, matching the benchmark.
  */
  if((uc->cat_enable || uc->rup_enable) && uc->slipping){
    const PetscScalar *x0a;
    PetscInt j,k;
    PetscCall(VecGetArrayRead(x0,&x0a));
    for(i = medium->rs, j=0, k=0; i < medium->re; i++, j+=medium->rsf->dim, k++){
      if(uc->snap_tau0) uc->snap_tau0[k]  = x0a[j+1];
      if(uc->snap_slip0)uc->snap_slip0[k] = x0a[j+3];
    }
    PetscCall(VecRestoreArrayRead(x0,&x0a));
    uc->ev_open    = PETSC_TRUE;
    uc->onset_time = t_init;
    uc->peakv_local = 0.0;
    if(uc->rup_enable)uc->rup_armed = PETSC_TRUE;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  gather the event-1 rupture-time field to rank 0 in geometry order
  and write it, close the catalog/budget files, and free the per-cell
  scratch.  The gather uses a Vec with the interaction-matrix row
  layout (patch i owned by the rank whose [rs,re) contains it),
  scattered to rank 0, so the written order matches the geometry file.
*/
PetscErrorCode rsf_finalize_catalog(struct rsf_out_ctx *uc)
{
  struct med *medium = uc->par->medium;
  struct flt *fault = uc->par->fault;
  PetscFunctionBeginUser;
  if(uc->rup_enable){
    /* pack local rup_time into a row-layout Vec, gather to rank 0 */
    Vec rvec,rgath;VecScatter rsc;
    PetscInt i,k;PetscScalar *ra;
    PetscCall(VecCreate(PETSC_COMM_WORLD,&rvec));
    PetscCall(VecSetSizes(rvec,medium->rn,medium->nrflt));
    PetscCall(VecSetFromOptions(rvec));
    PetscCall(VecGetArray(rvec,&ra));
    for(k=0;k < medium->rn;k++)
      ra[k] = (uc->rup_time)?(uc->rup_time[k]):(-1.0);
    PetscCall(VecRestoreArray(rvec,&ra));
    PetscCall(VecScatterCreateToZero(rvec,&rsc,&rgath));
    PetscCall(VecScatterBegin(rsc,rvec,rgath,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(  rsc,rvec,rgath,INSERT_VALUES,SCATTER_FORWARD));
    HEADNODE{
      const PetscScalar *g;
      FILE *out = myopen(RSF_RUPTURE_TIME_FILE,"w");
      int nrup=0;
      PetscReal t0=1e300;	/* initiation = earliest crossing */
      PetscCall(VecGetArrayRead(rgath,&g));
      for(i=0;i < medium->nrflt;i++)
	if(((PetscReal)g[i] >= 0.0) && ((PetscReal)g[i] < t0))t0 = (PetscReal)g[i];
      if(t0 > 1e299)t0 = 0.0;	/* no cell ruptured */
      fprintf(out,"# rupture-time field for event 1 (first tracked event)\n");
      fprintf(out,"# rupture-front |v| threshold = %.3e m/s ; t is seconds after initiation (earliest crossing)\n",
	      uc->rupture_vth);
      fprintf(out,"# ip xc[m] yc[m] zc[m] t_rupture[s] ruptured(1/0)\n");
      for(i=0;i < medium->nrflt;i++){
	PetscReal ta = (PetscReal)g[i];
	int rup = (ta >= 0.0)?(1):(0);
	PetscReal tr = rup?(ta - t0):(-1.0);
	if(rup)nrup++;
	fprintf(out,"%6i %14.6e %14.6e %14.6e %15.7e %2i\n",
		(int)i,fault[i].x[INT_X],fault[i].x[INT_Y],fault[i].x[INT_Z],tr,rup);
      }
      PetscCall(VecRestoreArrayRead(rgath,&g));
      fclose(out);
      fprintf(stderr,"rsf_finalize_catalog: wrote rupture-time field for event 1 (%i/%i cells ruptured)\n",
	      nrup,(int)medium->nrflt);
    }
    PetscCall(VecScatterDestroy(&rsc));
    PetscCall(VecDestroy(&rgath));
    PetscCall(VecDestroy(&rvec));
  }
  HEADNODE{
    if(uc->fout_catalog){
      fclose(uc->fout_catalog);
      fprintf(stderr,"rsf_finalize_catalog: wrote event catalog (%i completed events)\n",uc->ncat);
    }
    if(uc->cgrp_fout){
      int ig;
      for(ig=0;ig < uc->cgrp_n;ig++){
	if(uc->cgrp_fout[ig]){
	  fprintf(stderr,"rsf_finalize_catalog: group %i ruptured in %i of %i events\n",
		  uc->cgrp_id[ig],(uc->cgrp_ncat)?(uc->cgrp_ncat[ig]):(0),uc->ncat);
	  fclose(uc->cgrp_fout[ig]);
	}
      }
      free(uc->cgrp_fout);uc->cgrp_fout = NULL;
    }
    if(uc->cgrp_ncat){free(uc->cgrp_ncat);uc->cgrp_ncat = NULL;}
    if(uc->fout_budget)fclose(uc->fout_budget);
  }
  if(uc->snap_tau0){free(uc->snap_tau0);uc->snap_tau0=NULL;}
  if(uc->snap_slip0){free(uc->snap_slip0);uc->snap_slip0=NULL;}
  if(uc->cell_ruptured){free(uc->cell_ruptured);uc->cell_ruptured=NULL;}
  /* the per group catalog scratch lives on every rank */
  if(uc->cgrp_id){free(uc->cgrp_id);uc->cgrp_id=NULL;}
  if(uc->cgrp_peak){free(uc->cgrp_peak);uc->cgrp_peak=NULL;}
  if(uc->cgrp_t0){free(uc->cgrp_t0);uc->cgrp_t0=NULL;}
  if(uc->cgrp_t1){free(uc->cgrp_t1);uc->cgrp_t1=NULL;}
  if(uc->cg_lsum){free(uc->cg_lsum);uc->cg_lsum=NULL;}
  if(uc->cg_gsum){free(uc->cg_gsum);uc->cg_gsum=NULL;}
  if(uc->cg_lmx){free(uc->cg_lmx);uc->cg_lmx=NULL;}
  if(uc->cg_gmx){free(uc->cg_gmx);uc->cg_gmx=NULL;}
  if(uc->cg_lmn){free(uc->cg_lmn);uc->cg_lmn=NULL;}
  if(uc->cg_gmn){free(uc->cg_gmn);uc->cg_gmn=NULL;}
  uc->cgrp_n = 0;
  if(uc->rup_time){free(uc->rup_time);uc->rup_time=NULL;}
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
#define RSF_CKPT_MAGIC 20260726.0
#define RSF_CKPT_VERSION 1.0
#define RSF_CKPT_NMETA 16

/* write a checkpoint: a 16-entry metadata vector followed by the full
   solution vector, PETSc binary, crash-safe via write-to-tmp-and-rename;
   the previous checkpoint is kept as <file>.prev */
PetscErrorCode rsf_write_checkpoint(TS ts, Vec X, struct rsf_out_ctx *uc)
{
  Vec meta;
  PetscViewer viewer;
  PetscReal t,dt;
  PetscInt step;
  PetscMPIInt rank;
  char tmpf[350],prevf[350];
  PetscFunctionBegin;
  PetscCall(TSGetTime(ts,&t));
  PetscCall(TSGetTimeStep(ts,&dt));
  PetscCall(TSGetStepNumber(ts,&step));
  PetscCall(VecCreate(PETSC_COMM_WORLD,&meta));
  PetscCall(VecSetSizes(meta,PETSC_DECIDE,RSF_CKPT_NMETA));
  PetscCall(VecSetFromOptions(meta));
  PetscCall(VecSet(meta,0.0));
  PetscCall(VecSetValue(meta,0,RSF_CKPT_MAGIC,INSERT_VALUES));
  PetscCall(VecSetValue(meta,1,RSF_CKPT_VERSION,INSERT_VALUES));
  PetscCall(VecSetValue(meta,2,(PetscScalar)uc->par->medium->nrflt,INSERT_VALUES));
  PetscCall(VecSetValue(meta,3,(PetscScalar)uc->ckpt_dim,INSERT_VALUES));
  PetscCall(VecSetValue(meta,4,(PetscScalar)uc->ckpt_slip_mode,INSERT_VALUES));
  PetscCall(VecSetValue(meta,5,(PetscScalar)uc->ckpt_law,INSERT_VALUES));
  PetscCall(VecSetValue(meta,6,(PetscScalar)step,INSERT_VALUES));
  PetscCall(VecSetValue(meta,7,(PetscScalar)t,INSERT_VALUES));
  PetscCall(VecSetValue(meta,8,(PetscScalar)dt,INSERT_VALUES));
  /* visco-elastic memory states (rsf_ve.c): slots 9/10 flag the VE
     block appended after X; 0 for elastic runs, so pre-VE
     checkpoints stay compatible */
  PetscCall(VecSetValue(meta,9,(PetscScalar)uc->par->medium->rsf->ve_np,INSERT_VALUES));
  PetscCall(VecSetValue(meta,10,(PetscScalar)uc->par->medium->rsf->ve_t0,INSERT_VALUES));
  PetscCall(VecAssemblyBegin(meta));
  PetscCall(VecAssemblyEnd(meta));
  snprintf(tmpf,sizeof(tmpf),"%s.tmp",uc->ckpt_file);
  snprintf(prevf,sizeof(prevf),"%s.prev",uc->ckpt_file);
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpf,FILE_MODE_WRITE,&viewer));
  PetscCall(VecView(meta,viewer));
  PetscCall(VecView(X,viewer));
  if(uc->par->medium->rsf->ve_np > 0){
    PetscInt p;
    for(p=0;p < uc->par->medium->rsf->ve_np;p++)
      PetscCall(VecView(uc->par->medium->rsf->ve_h[p],viewer));
    PetscCall(VecView(uc->par->medium->rsf->ve_slip_prev,viewer));
  }
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(VecDestroy(&meta));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  if(rank == 0){
    char tmpi[360],maini[360];
    rename(uc->ckpt_file,prevf);	/* ok to fail on the first write */
    snprintf(tmpi,sizeof(tmpi),"%s.info",tmpf);
    snprintf(maini,sizeof(maini),"%s.info",uc->ckpt_file);
    rename(tmpi,maini);			/* PETSc viewer companion file */
    if(rename(tmpf,uc->ckpt_file) != 0)
      fprintf(stderr,"rsf_write_checkpoint: WARNING: rename to %s failed\n",uc->ckpt_file);
    else
      fprintf(stderr,"rsf_write_checkpoint: step %012ld t %12.5e yr dt %.3e s -> %s\n",
	      (long)step,(double)t/SEC_PER_YEAR,(double)dt,uc->ckpt_file);
  }
  PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* load a checkpoint into the (already created, correctly sized)
   solution vector; validates patch count and dim, warns on law or
   slip-mode differences; returns time, dt, and step */
PetscErrorCode rsf_read_checkpoint(const char *file, Vec X, struct rsf_out_ctx *uc,
				   PetscReal *t, PetscReal *dt, PetscInt *step)
{
  Vec meta;
  PetscViewer viewer;
  const PetscScalar *mv;
  Vec mseq;
  VecScatter scat;
  PetscReal magic;
  PetscInt nrflt,dim,slip_mode,law;
  struct med *medium;
  PetscFunctionBegin;
  medium = uc->par->medium;
  PetscCall(VecCreate(PETSC_COMM_WORLD,&meta));
  PetscCall(VecSetSizes(meta,PETSC_DECIDE,RSF_CKPT_NMETA));
  PetscCall(VecSetFromOptions(meta));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&viewer));
  PetscCall(VecLoad(meta,viewer));
  /* every rank needs the metadata */
  PetscCall(VecScatterCreateToAll(meta,&scat,&mseq));
  PetscCall(VecScatterBegin(scat,meta,mseq,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecScatterEnd(scat,meta,mseq,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecGetArrayRead(mseq,&mv));
  magic = (PetscReal)mv[0];
  nrflt = (PetscInt)mv[2];
  dim = (PetscInt)mv[3];
  slip_mode = (PetscInt)mv[4];
  law = (PetscInt)mv[5];
  *step = (PetscInt)mv[6];
  *t = (PetscReal)mv[7];
  *dt = (PetscReal)mv[8];
  {				/* VE block flags, 0 in elastic checkpoints */
    PetscInt cnp = (PetscInt)mv[9];
    struct rsf_vars *rsf = medium->rsf;
    if(cnp != rsf->ve_np)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_INCOMP,
	      "rsf_read_checkpoint: checkpoint has %ld VE memory terms, run has %ld (elastic and VE runs cannot restart from each other)",
	      (long)cnp,(long)rsf->ve_np);
    if(rsf->ve_np > 0)
      rsf->ve_t0 = (PetscReal)mv[10];
  }
  PetscCall(VecRestoreArrayRead(mseq,&mv));
  PetscCall(VecScatterDestroy(&scat));
  PetscCall(VecDestroy(&mseq));
  PetscCall(VecDestroy(&meta));
  if(magic != RSF_CKPT_MAGIC)
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_UNEXPECTED,
	    "rsf_read_checkpoint: %s is not an rsf_solve checkpoint",file);
  if(nrflt != uc->par->medium->nrflt)
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_INCOMP,
	    "rsf_read_checkpoint: checkpoint has %ld patches, geometry has %i",
	    (long)nrflt,uc->par->medium->nrflt);
  if(dim != uc->ckpt_dim)
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_INCOMP,
	    "rsf_read_checkpoint: checkpoint dim %ld, run dim %ld",
	    (long)dim,(long)uc->ckpt_dim);
  HEADNODE{
    if(slip_mode != uc->ckpt_slip_mode)
      fprintf(stderr,"rsf_read_checkpoint: WARNING: checkpoint slip mode %ld, run %ld\n",
	      (long)slip_mode,(long)uc->ckpt_slip_mode);
    if(law != uc->ckpt_law)
      fprintf(stderr,"rsf_read_checkpoint: WARNING: checkpoint state law %ld, run %ld\n",
	      (long)law,(long)uc->ckpt_law);
  }
  PetscCall(VecLoad(X,viewer));
  if(medium->rsf->ve_np > 0){	/* the VE memory block after X */
    PetscInt p;
    for(p=0;p < medium->rsf->ve_np;p++)
      PetscCall(VecLoad(medium->rsf->ve_h[p],viewer));
    PetscCall(VecLoad(medium->rsf->ve_slip_prev,viewer));
  }
  PetscCall(PetscViewerDestroy(&viewer));
  HEADNODE
    fprintf(stderr,"rsf_read_checkpoint: restarting from %s: step %ld t %.8e s (%.4f yr) dt %.3e s%s\n",
	    file,(long)*step,(double)*t,(double)(*t)/SEC_PER_YEAR,(double)*dt,
	    (medium->rsf->ve_np > 0)?(" (with VE memory states)"):(""));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode rsf_TS_Monitor(TS ts,PetscInt step,PetscReal time,Vec X,void *ptr)
{
  const PetscScalar *x;
  struct med *medium;struct flt *fault;
  struct rsf_out_ctx *uc;
  Vec DX;
  PetscInt i,j;
  PetscBool bail;
  int ig;
  PetscInt kk;
  PetscReal v,lsum[MGRP_NSUM],gsum[MGRP_NSUM],lminmax[3],gminmax[3],dt,d1,d2,d3,dx_norm,x_norm;
  struct rsf_vars *rsf;
  PetscFunctionBeginUser;
  {
    struct rsf_out_ctx *uc_ck = (struct rsf_out_ctx *)ptr;
    if((uc_ck->ckpt_every > 0) && (step > uc_ck->ckpt_step0) &&
       (step % uc_ck->ckpt_every == 0))
      PetscCall(rsf_write_checkpoint(ts,X,uc_ck));
  }
  uc = (struct rsf_out_ctx *)ptr;
  medium = uc->par->medium;
  rsf = medium->rsf;
  fault = uc->par->fault;
  if(step < 0)	      /* negative indicates an interpolated solution */
    PetscFunctionReturn(PETSC_SUCCESS);
  /*
    while an event is in progress (uc->slipping), sample every
    accepted step to build the per-cell peak slip rate, the ruptured
    mask, and (for event 1) the rupture-time field.  This runs OUTSIDE
    the change-triggered monitor gate below, so it uses the solver's
    own step density, which collapses through the coseismic phase and
    therefore resolves the rupture front well.  It reads the state
    only and never touches the integrator.
  */
  if((uc->cat_enable || uc->rup_enable) && uc->slipping){
    const PetscScalar *xt;
    PetscInt k;
    PetscReal vv,d1t,d2t,d3t;
    PetscCall(VecGetArrayRead(X,&xt));
    for(i = medium->rs, j=0, k=0; i < medium->re; i++, j+=rsf->dim, k++){
      vv = fabs(vel_from_rsf(xt[j+1],xt[j+2],xt[j],fault[i].mu_sa,rsf->v0,&d1t,&d2t,&d3t,medium));
      if(vv > uc->peakv_local)uc->peakv_local = vv;
      if(uc->cgrp_n > 0){	/* same running peak, per group */
	int igc = uc->mgrp_map[k];
	if(vv > uc->cgrp_peak[igc])uc->cgrp_peak[igc] = vv;
      }
      if(vv >= uc->rupture_vth){
	if(uc->cell_ruptured)uc->cell_ruptured[k] = 1;
	if(uc->cgrp_n > 0){	/* when this group first and last ruptured */
	  int igc = uc->mgrp_map[k];
	  if(uc->cgrp_t0[igc] < 0.0)uc->cgrp_t0[igc] = time;
	  uc->cgrp_t1[igc] = time;
	}
	/* first crossing time for event 1 only (absolute; referenced to the
	   earliest crossing, i.e. the initiation time, when written) */
	if(uc->rup_enable && uc->rup_armed && (!uc->rup_done) &&
	   uc->rup_time && (uc->rup_time[k] < 0.0))
	  uc->rup_time[k] = time;
      }
    }
    PetscCall(VecRestoreArrayRead(X,&xt));
  }
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
      lsum[0]=lsum[1]=lsum[2]=0.0;
      lminmax[0] = -1e30;		/* max |v| */
      lminmax[1] = -1e30;		/* max sigma */
      lminmax[2] = -1e30;		/* -min sigma */
      for(ig=0;ig < uc->mgrp_n;ig++){ /* same, per fault group */
	uc->mgrp_lsum[MGRP_NSUM*ig]   = uc->mgrp_lsum[MGRP_NSUM*ig+1] =
	  uc->mgrp_lsum[MGRP_NSUM*ig+2] = 0.0;
	uc->mgrp_lmx[3*ig] = uc->mgrp_lmx[3*ig+1] = uc->mgrp_lmx[3*ig+2] = -1e30;
      }
      for (i = medium->rs, j=0, kk=0; i < medium->re; i++, j+=rsf->dim, kk++) {
	v = vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_sa,rsf->v0,&d1,&d2,&d3,medium);
	v = fabs(v);
	if(v > lminmax[0])lminmax[0] = v;
	if( x[j+2] > lminmax[1])lminmax[1] =  x[j+2];
	if(-x[j+2] > lminmax[2])lminmax[2] = -x[j+2];
	lsum[0] += x[j+3];		/* slip */
	lsum[1] += x[j+1]/x[j+2];	/* mu */
	lsum[2] += v;			/* |v|, the mean companion to max|v| */
	if(uc->mgrp_n > 0){		/* and into this patch's group */
	  ig = uc->mgrp_map[kk];
	  if(v > uc->mgrp_lmx[3*ig])uc->mgrp_lmx[3*ig] = v;
	  if( x[j+2] > uc->mgrp_lmx[3*ig+1])uc->mgrp_lmx[3*ig+1] =  x[j+2];
	  if(-x[j+2] > uc->mgrp_lmx[3*ig+2])uc->mgrp_lmx[3*ig+2] = -x[j+2];
	  uc->mgrp_lsum[MGRP_NSUM*ig]   += x[j+3];
	  uc->mgrp_lsum[MGRP_NSUM*ig+1] += x[j+1]/x[j+2];
	  uc->mgrp_lsum[MGRP_NSUM*ig+2] += v;
	}
      }
      PetscCall(VecRestoreArrayRead(X,&x));
      PetscCallMPI(MPI_Reduce(lsum,gsum,MGRP_NSUM,MPIU_REAL,MPI_SUM,0,PETSC_COMM_WORLD));
      PetscCallMPI(MPI_Reduce(lminmax,gminmax,3,MPIU_REAL,MPI_MAX,0,PETSC_COMM_WORLD));
      if(uc->mgrp_n > 0){
	PetscCallMPI(MPI_Reduce(uc->mgrp_lsum,uc->mgrp_gsum,MGRP_NSUM*uc->mgrp_n,MPIU_REAL,MPI_SUM,0,PETSC_COMM_WORLD));
	PetscCallMPI(MPI_Reduce(uc->mgrp_lmx, uc->mgrp_gmx, 3*uc->mgrp_n,MPIU_REAL,MPI_MAX,0,PETSC_COMM_WORLD));
      }
      if((medium->comm_rank == 0) && uc->fout_monitor){
	fprintf(uc->fout_monitor,"%9i %20.8e %17.10f %15.8e %12.7f %15.8e %10.6f %15.8e %15.8e %12.7f\n",
		(int)step,time,time/SEC_PER_YEAR,dt,
		log10(gminmax[0]),gsum[0]/(PetscReal)medium->nrflt,
		gsum[1]/(PetscReal)medium->nrflt,
		gminmax[1],-gminmax[2],
		/* column 10 is appended, so 1 to 9 keep their meaning and
		   mean|v| can be read as $NF.  the mean shear stress is not
		   written separately: it is controlled by mu, which is
		   already column 7 */
		log10((gsum[2]/(PetscReal)medium->nrflt > 1e-300)?
		      (gsum[2]/(PetscReal)medium->nrflt):(1e-300)));
	fflush(uc->fout_monitor); /* so a long run does not look dead: the file is
				     otherwise buffered and appears empty while running */
      }
      if((medium->comm_rank == 0) && (uc->mgrp_n > 0) && uc->mgrp_fout){
	for(ig=0;ig < uc->mgrp_n;ig++){
	  if(!uc->mgrp_fout[ig])
	    continue;
	  fprintf(uc->mgrp_fout[ig],"%9i %20.8e %17.10f %15.8e %12.7f %15.8e %10.6f %15.8e %15.8e %12.7f\n",
		  (int)step,time,time/SEC_PER_YEAR,dt,
		  log10((uc->mgrp_gmx[3*ig] > 1e-300)?(uc->mgrp_gmx[3*ig]):(1e-300)),
		  uc->mgrp_gsum[MGRP_NSUM*ig]/(PetscReal)uc->mgrp_np[ig],
		  uc->mgrp_gsum[MGRP_NSUM*ig+1]/(PetscReal)uc->mgrp_np[ig],
		  uc->mgrp_gmx[3*ig+1],-uc->mgrp_gmx[3*ig+2],
		  log10((uc->mgrp_gsum[MGRP_NSUM*ig+2]/(PetscReal)uc->mgrp_np[ig] > 1e-300)?
			(uc->mgrp_gsum[MGRP_NSUM*ig+2]/(PetscReal)uc->mgrp_np[ig]):(1e-300)));
	  fflush(uc->mgrp_fout[ig]);
	}
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
    HEADNODE{
      PetscScalar *values;
      PetscReal sum[3],du1,du2,du3;
      PetscReal slip_integral=0.0;	/* slip budget: sum area_i*slip_i */
      FILE *fout2;
      char vel_file[STRLEN];
      int n = medium->nrflt,ierr2;
      PetscCall(VecGetArray(uc->gathered,&values));
      for(sum[2]=0.,i=0,j=0;i < n;i++,j += rsf->dim){
	fault[i].s[STRIKE] = values[j+1]; /* shear stress */
	fault[i].s[NORMAL] = values[j+2]; /* normal stress */
	/* the rate-and-state solve carries ONE slip component (rsf->dim = 4:
	   psi,tau,sigma,slip), and rsf->slip_mode says whether it is strike or
	   dip. vel_from_rsf returns that component, so store it in the matching
	   slot rather than always in u[STRIKE], which would mislabel a dip-slip
	   run for anything downstream reading fault[].u */
	fault[i].u[rsf->slip_mode] =
	  vel_from_rsf(values[j+1],values[j+2],values[j],fault[i].mu_sa,
		       rsf->v0,&du1,&du2,&du3,medium);
	sum[2] += values[j+3];				     /* slip */
	if(uc->budget_enable)
	  slip_integral += 4.0*fault[i].l*fault[i].w*values[j+3]; /* area*slip */
      }
      if(uc->budget_enable && uc->fout_budget){
	PetscReal ref = uc->vpl*time*uc->total_area; /* cumulative slip*area at plate rate */
	PetscReal resid = slip_integral - ref;
	fprintf(uc->fout_budget,"%.8f %.8e %.8e %.8e %.6e\n",
		time/SEC_PER_YEAR,slip_integral,ref,resid,(ref!=0.0)?(resid/ref):(0.0));
      }
      if(time - medium->slip_line_time > medium->slip_line_dt){
	if(!uc->field_out){
	  ierr2 = system("mkdir -p tmp_rsf");
	  if(ierr2)
	    fprintf(stderr,"rsf_TS_Monitor: WARNING: error making tmp_rsf output directory\n");
	}
	snprintf(vel_file,STRLEN,"tmp_rsf/vel-%012.5e-%06i-gmt",time/SEC_PER_YEAR,uc->field_out);
	fout2 = myopen(vel_file,"w");
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
  /*
     compact per-fault slip-rate field frame on a fixed accepted-step
     cadence, independent of the time-based stats gather above.  The
     scatter is collective so every rank enters; only rank 0 writes.
     One binary per group holds xyz2grd-ready (along_strike, down_dip,
     log10|v|) float triples.  Reads the state only, never the
     integrator, so it cannot perturb the trajectory
  */
  if(uc->field_enable && (uc->field_step_interval > 0) &&
     ((step % uc->field_step_interval) == 0) && (time >= uc->field_tmin)){
    PetscScalar *fvals;
    PetscCall(VecScatterBegin(uc->gather,X,uc->gathered,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(  uc->gather,X,uc->gathered,INSERT_VALUES,SCATTER_FORWARD));
    HEADNODE{
      PetscReal du1,du2,du3,vmax=0.0,vmin=1e30,vsum=0.0,vsqsum=0.0,ssum=0.0,vmean,vstd;
      int ig,k,ip,nf = medium->nrflt;
      PetscCall(VecGetArray(uc->gathered,&fvals));
      for(ip=0;ip < nf;ip++){
	PetscReal vabs,vv = vel_from_rsf(fvals[ip*rsf->dim+1],fvals[ip*rsf->dim+2],
					 fvals[ip*rsf->dim+0],fault[ip].mu_sa,rsf->v0,
					 &du1,&du2,&du3,medium);
	uc->vbuf[ip] = (double)vv;
	/* statistics are on the slip SPEED |v|: vel_from_rsf is signed (the sinh
	   regularization admits v < 0), and a signed mean or max would not be a
	   meaningful summary of how fast the fault is going.  This is the single
	   solved slip component, whose direction is rsf->slip_mode; if the solve
	   is ever generalized to carry strike and dip simultaneously, the speed
	   here becomes sqrt(v_strike^2 + v_dip^2) */
	vabs = fabs(vv);
	if(vabs > vmax)vmax = vabs;
	if(vabs < vmin)vmin = vabs;
	vsum   += vabs;
	vsqsum += vabs*vabs;
	ssum   += fvals[ip*rsf->dim+3];	/* slip */
      }
      vmean = (nf > 0)?(vsum/(PetscReal)nf):(0.0);
      vstd  = (nf > 1)?(sqrt(fabs((PetscReal)nf*vsqsum - vsum*vsum)/
			     ((PetscReal)nf*(PetscReal)(nf-1)))):(0.0);
      for(ig=0;ig < uc->ngroup;ig++){
	struct rsf_group_grid *g = uc->groups+ig;
	char ffile[STRLEN];FILE *fb;
	if(g->np < 1)
	  continue;
	for(k=0;k < g->np;k++){
	  double vv = fabs(uc->vbuf[g->idx[k]]);
	  if(vv < 1e-30)vv = 1e-30;
	  g->buf[3*k+0] = (float)g->xs[k];
	  g->buf[3*k+1] = (float)g->ys[k];
	  g->buf[3*k+2] = (float)log10(vv);
	}
	snprintf(ffile,STRLEN,"tmp_rsf/rsf_vel.g%03d.%06i.bin",g->id,uc->field_frame);
	fb = myopen(ffile,"wb");
	if(fb){
	  fwrite(g->buf,sizeof(float),(size_t)3*(size_t)g->np,fb);
	  fclose(fb);
	}
	if(uc->field_stress){
	  /* shear-stress frame, (x, y, tau[MPa]) triples, same layout;
	     tau is the |tau| stress component the solve carries (state
	     slot 1), i.e. the resolved shear stress in the slip mode */
	  for(k=0;k < g->np;k++){
	    g->buf[3*k+0] = (float)g->xs[k];
	    g->buf[3*k+1] = (float)g->ys[k];
	    g->buf[3*k+2] = (float)(fvals[g->idx[k]*rsf->dim+1]/1e6);
	  }
	  snprintf(ffile,STRLEN,"tmp_rsf/rsf_tau.g%03d.%06i.bin",g->id,uc->field_frame);
	  fb = myopen(ffile,"wb");
	  if(fb){
	    fwrite(g->buf,sizeof(float),(size_t)3*(size_t)g->np,fb);
	    fclose(fb);
	  }
	}
      }
      if(uc->fout_field_times){
	fprintf(uc->fout_field_times,"%6i %9i %.8f %.8e %.6f %14.6e %14.6e %14.6e %14.6e\n",
		uc->field_frame,(int)step,time/SEC_PER_YEAR,time,
		log10((vmax>0.0)?(vmax):(1e-30)),
		vmean,vstd,(vmin<1e30)?(vmin):(0.0),
		(nf>0)?(ssum/(PetscReal)nf):(0.0));
	fflush(uc->fout_field_times);
      }
      PetscCall(VecRestoreArray(uc->gathered,&fvals));
      uc->field_frame++;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* 
   log the event, cf. myPostEventFunction in ode_solve_test.c: write
   time, onset/arrest flag, and the global state summary
*/
PetscErrorCode rsf_post_event(TS ts,PetscInt nevents,PetscInt event_list[],
			      PetscReal t,Vec X,PetscBool forwardsolve,void *ctx)
{
  const PetscScalar *x;
  struct rsf_out_ctx *uc;
  struct med *medium;struct flt *fault;
  struct rsf_vars *rsf;
  PetscInt i,j;
  PetscReal v,lsum[2],gsum[2],lvmax,gvmax,d1,d2,d3;
  PetscFunctionBeginUser;
  uc = (struct rsf_out_ctx *)ctx;
  medium = uc->par->medium;fault = uc->par->fault;
  rsf = medium->rsf;
  uc->slipping = (uc->slipping)?(PETSC_FALSE):(PETSC_TRUE); /* threshold crossed */
  if(t >= uc->event_tmin){
    lsum[0]=lsum[1]=0.0;lvmax = 0.0;
    PetscCall(VecGetArrayRead(X,&x));
    for (i = medium->rs, j=0; i < medium->re; i++, j+=rsf->dim) {
      v = fabs(vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_sa,rsf->v0,&d1,&d2,&d3,medium));
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
	      t,t/SEC_PER_YEAR,(uc->slipping)?(1):(-1),
	      log10((gvmax > 1e-300)?(gvmax):(1e-300)),
	      gsum[0]/(PetscReal)medium->nrflt,gsum[1]/(PetscReal)medium->nrflt);
      fflush(uc->fout_event);
    }
  }
  /*
    catalog capture, on the SAME crossing this event handler just
    logged to rsf_events.dat.  After the toggle above, uc->slipping ==
    TRUE means this crossing is an onset, FALSE means an arrest.
    Snapshots are captured on onset and differenced on arrest; the
    reduces are collective so all ranks enter.  This runs regardless
    of event_tmin so onset/arrest stay paired; only the written
    catalog row is gated by event_tmin.
  */
  if(uc->cat_enable || uc->rup_enable){
    const PetscScalar *xc;
    PetscInt ii,jj,kk;
    if(uc->slipping){
      /* ONSET: snapshot tau,slip; reset the per-event trackers */
      PetscCall(VecGetArrayRead(X,&xc));
      for(ii=medium->rs,jj=0,kk=0; ii<medium->re; ii++,jj+=rsf->dim,kk++){
	if(uc->snap_tau0) uc->snap_tau0[kk]  = xc[jj+1];
	if(uc->snap_slip0)uc->snap_slip0[kk] = xc[jj+3];
	if(uc->cell_ruptured)uc->cell_ruptured[kk] = 0;
      }
      PetscCall(VecRestoreArrayRead(X,&xc));
      for(kk=0;kk < (PetscInt)uc->cgrp_n;kk++){
	uc->cgrp_peak[kk] = 0.0;
	uc->cgrp_t0[kk] = uc->cgrp_t1[kk] = -1.0;
      }
      uc->ev_open     = PETSC_TRUE;
      uc->onset_time  = t;
      uc->peakv_local = 0.0;
      if(uc->rup_enable && (!uc->rup_done) && (!uc->rup_armed))
	uc->rup_armed = PETSC_TRUE; /* arm the rupture-time field for event 1 */
    }else if(uc->ev_open){
      /* ARREST: close the event out and write its catalog row */
      PetscCall(rsf_finalize_event(uc,t,X));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Close out an event and write its catalog row: coseismic slip and stress drop
  summed over the cells that reached rupture_vth between onset and arrest.  t_arr
  is the model time of the arrest that ended the event and X the solution there.
  Collective: all ranks must enter.
*/
PetscErrorCode rsf_finalize_event(struct rsf_out_ctx *uc,PetscReal t_arr,Vec X)
{
  const PetscScalar *xc;
  struct med *medium;struct flt *fault;struct rsf_vars *rsf;
  PetscInt ii,jj,kk;
  int igc,ngc;
  PetscReal lred[7],gred[7],lmx[2],gmx[2],gpeak,gcnt;
  PetscFunctionBeginUser;
  medium = uc->par->medium;fault = uc->par->fault;rsf = medium->rsf;
  if(!uc->ev_open)
    PetscFunctionReturn(PETSC_SUCCESS);
  /* lred: 0 sum |dslip|, 1 sum drop, 2 count, 3 area, 4 moment,
     5 2D per-unit-thickness moment [Nm/m], 6 2D ruptured in-plane
     length [m].  Slip enters through its absolute value and the
     stress drop through the change of |tau|, so both slip senses
     (e.g. backslip-loaded normal faulting with negative rates)
     yield positive moment and drop; without this, negative-sense
     events produced M0 < 0 and the Mw guard value of -99 */
  lred[0]=lred[1]=lred[2]=lred[3]=lred[4]=lred[5]=lred[6]=0.0;
  lmx[0]=lmx[1]=0.0;		/* max |dslip|, max drop */
  /* the same accumulators, per fault group: sums 0-6 as in lred, maxima
     0 max |dslip|, 1 max drop, 2 peak |v|, 3 last rupture time, and a
     single minimum, the first rupture time.  Groups that never reached
     rupture_vth stay at zero count and get no row */
  ngc = uc->cgrp_n;
  for(igc=0;igc < ngc;igc++){
    uc->cg_lsum[7*igc+0]=uc->cg_lsum[7*igc+1]=uc->cg_lsum[7*igc+2]=0.0;
    uc->cg_lsum[7*igc+3]=uc->cg_lsum[7*igc+4]=0.0;
    uc->cg_lsum[7*igc+5]=uc->cg_lsum[7*igc+6]=0.0;
    uc->cg_lmx[4*igc+0]=uc->cg_lmx[4*igc+1]=0.0;
    uc->cg_lmx[4*igc+2]=uc->cgrp_peak[igc];
    uc->cg_lmx[4*igc+3]=uc->cgrp_t1[igc];
    uc->cg_lmn[igc]=(uc->cgrp_t0[igc] < 0.0)?(1e300):(uc->cgrp_t0[igc]);
  }
  PetscCall(VecGetArrayRead(X,&xc));
  for(ii=medium->rs,jj=0,kk=0; ii<medium->re; ii++,jj+=rsf->dim,kk++){
    if(uc->cell_ruptured && uc->cell_ruptured[kk]){
      PetscReal dslip = fabs(xc[jj+3] - (uc->snap_slip0?uc->snap_slip0[kk]:0.0));
      PetscReal ddrop = fabs(uc->snap_tau0?uc->snap_tau0[kk]:0.0) - fabs(xc[jj+1]);
      PetscReal area  = 4.0*fault[ii].l*fault[ii].w;
      lred[0]+=dslip; lred[1]+=ddrop; lred[2]+=1.0; lred[3]+=area;
      lred[4]+=uc->shear_modulus*area*dslip;
      if((fault[ii].type == TWO_DIM_SEGMENT_PLANE_STRAIN)||
	 (fault[ii].type == TWO_DIM_SEGMENT_PLANE_STRESS)||
	 (fault[ii].type == TWO_DIM_HALFPLANE_PLANE_STRAIN)){
	lred[5]+=uc->shear_modulus*2.0*fault[ii].l*dslip;
	lred[6]+=2.0*fault[ii].l;
      }
      if(dslip > lmx[0])lmx[0]=dslip;
      if(ddrop > lmx[1])lmx[1]=ddrop;
      if(ngc > 0){
	igc = uc->mgrp_map[kk];
	uc->cg_lsum[7*igc+0]+=dslip; uc->cg_lsum[7*igc+1]+=ddrop;
	uc->cg_lsum[7*igc+2]+=1.0;   uc->cg_lsum[7*igc+3]+=area;
	uc->cg_lsum[7*igc+4]+=uc->shear_modulus*area*dslip;
	if((fault[ii].type == TWO_DIM_SEGMENT_PLANE_STRAIN)||
	   (fault[ii].type == TWO_DIM_SEGMENT_PLANE_STRESS)||
	   (fault[ii].type == TWO_DIM_HALFPLANE_PLANE_STRAIN)){
	  uc->cg_lsum[7*igc+5]+=uc->shear_modulus*2.0*fault[ii].l*dslip;
	  uc->cg_lsum[7*igc+6]+=2.0*fault[ii].l;
	}
	if(dslip > uc->cg_lmx[4*igc+0])uc->cg_lmx[4*igc+0]=dslip;
	if(ddrop > uc->cg_lmx[4*igc+1])uc->cg_lmx[4*igc+1]=ddrop;
      }
    }
  }
  PetscCall(VecRestoreArrayRead(X,&xc));
  /* global ruptured-cell count on every rank, so the real-event decision below is
     identical on all ranks (the collective reductions and the rup_done flag must
     stay consistent across ranks) */
  PetscCallMPI(MPI_Allreduce(&lred[2],&gcnt,1,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD));
  /* a bracket where no cell reached rupture_vth (gcnt == 0) is afterslip or a
     sub-threshold transient, not a rupture: skip it.  rupture_vth is the
     definition of rupture, so gcnt == 0 is not an event. */
  if(gcnt > 0.0){
    PetscCallMPI(MPI_Reduce(lred,gred,7,MPIU_REAL,MPI_SUM,0,PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Reduce(lmx,gmx,2,MPIU_REAL,MPI_MAX,0,PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Reduce(&uc->peakv_local,&gpeak,1,MPIU_REAL,MPI_MAX,0,PETSC_COMM_WORLD));
    if((medium->comm_rank == 0) && uc->fout_catalog && (t_arr >= uc->event_tmin)){
      PetscReal mean_slip = gred[0]/gcnt;
      PetscReal mean_drop = gred[1]/gcnt;
      PetscReal M0        = gred[4];
      PetscReal Mw;
      if(gred[6] > 0.0){
	/* 2D geometry: the physical quantity is the per-unit-thickness
	   moment gred[5] [Nm/m]; for an equivalent magnitude assume an
	   along-strike dimension equal to the ruptured in-plane length
	   (aspect ratio one), a common convention for 2D cycle models.
	   The M0 column then carries this equivalent moment [Nm] */
	M0 = gred[5]*gred[6];
      }
      Mw = mwfromm0(M0);
      uc->ncat++;
      fprintf(uc->fout_catalog,
	      "%5i %17.10f %17.10f %13.6e %8i %14.6e %13.6e %13.6e %12.6e %12.6e %13.6e %13.6e %8.4f\n",
	      uc->ncat,uc->onset_time/SEC_PER_YEAR,t_arr/SEC_PER_YEAR,t_arr-uc->onset_time,
	      (int)gcnt,gred[3],mean_slip,gmx[0],mean_drop/1e6,gmx[1]/1e6,gpeak,M0,Mw);
      fflush(uc->fout_catalog);
    }
    /* per group rows for the same event.  ncat was just incremented on
       rank 0 for a written row; when the row was suppressed by
       event_tmin the group rows are suppressed with it, so the two
       files stay joinable on the event index */
    if(ngc > 0){
      PetscCallMPI(MPI_Reduce(uc->cg_lsum,uc->cg_gsum,7*ngc,MPIU_REAL,MPI_SUM,0,PETSC_COMM_WORLD));
      PetscCallMPI(MPI_Reduce(uc->cg_lmx, uc->cg_gmx, 4*ngc,MPIU_REAL,MPI_MAX,0,PETSC_COMM_WORLD));
      PetscCallMPI(MPI_Reduce(uc->cg_lmn, uc->cg_gmn, ngc,  MPIU_REAL,MPI_MIN,0,PETSC_COMM_WORLD));
      if((medium->comm_rank == 0) && uc->cgrp_fout && uc->fout_catalog &&
	 (t_arr >= uc->event_tmin)){
	for(igc=0;igc < ngc;igc++){
	  PetscReal gc = uc->cg_gsum[7*igc+2];
	  PetscReal M0g,Mwg,t0g,t1g;
	  if((gc <= 0.0) || (!uc->cgrp_fout[igc]))
	    continue;		/* this group did not rupture in this event */
	  M0g = uc->cg_gsum[7*igc+4];
	  if(uc->cg_gsum[7*igc+6] > 0.0)
	    M0g = uc->cg_gsum[7*igc+5]*uc->cg_gsum[7*igc+6];
	  Mwg = mwfromm0(M0g);
	  t0g = (uc->cg_gmn[igc] > 1e299)?(uc->onset_time):(uc->cg_gmn[igc]);
	  t1g = (uc->cg_gmx[4*igc+3] < 0.0)?(t_arr):(uc->cg_gmx[4*igc+3]);
	  uc->cgrp_ncat[igc]++;
	  fprintf(uc->cgrp_fout[igc],
		  "%5i %17.10f %17.10f %17.10f %17.10f %8i %14.6e %13.6e %13.6e %12.6e %12.6e %13.6e %13.6e %8.4f\n",
		  uc->ncat,uc->onset_time/SEC_PER_YEAR,t_arr/SEC_PER_YEAR,
		  t0g/SEC_PER_YEAR,t1g/SEC_PER_YEAR,(int)gc,
		  uc->cg_gsum[7*igc+3],uc->cg_gsum[7*igc+0]/gc,uc->cg_gmx[4*igc+0],
		  uc->cg_gsum[7*igc+1]/gc/1e6,uc->cg_gmx[4*igc+1]/1e6,
		  uc->cg_gmx[4*igc+2],M0g,Mwg);
	  fflush(uc->cgrp_fout[igc]);
	}
      }
    }
    /* the rupture-time field is for the first real rupture; freeze it here */
    if(uc->rup_enable && uc->rup_armed && (!uc->rup_done))
      uc->rup_done = PETSC_TRUE;
  }
  uc->ev_open = PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*
   write one group's static geometry.  If the patches form a single
   regular grid the header carries a ready GMT -R/-I that xyz2grd can
   use directly; otherwise grid the field frames with surface or
   nearneighbor instead.  Columns are documented in the header line
*/
void rsf_write_group_geometry(const struct rsf_group_grid *g,struct flt *fault,
				     double sigma0,const char *fname)
{
  FILE *out;int k;
  if(g->np < 1)
    return;
  out=myopen(fname,"w");
  if(!out){
    fprintf(stderr,"rsf_write_group_geometry: cannot open %s for writing\n",fname);
    return;
  }
  fprintf(out,"# rsf_solve static fault geometry, group %d, %d patches, written once\n",g->id,g->np);
  if(g->regular){
    fprintf(out,"# grid regular nx %d ny %d\n",g->nx,g->ny);
    fprintf(out,"# gmt_region -R%.10g/%.10g/%.10g/%.10g\n",
	    g->xmin-0.5*g->dx,g->xmax+0.5*g->dx,g->ymin-0.5*g->dy,g->ymax+0.5*g->dy);
    fprintf(out,"# gmt_inc -I%.10g/%.10g\n",g->dx,g->dy);
  }else{
    fprintf(out,"# grid irregular (grid the field frames with surface or nearneighbor)\n");
    fprintf(out,"# gmt_region -R%.10g/%.10g/%.10g/%.10g\n",g->xmin,g->xmax,g->ymin,g->ymax);
  }
  fprintf(out,"# col: ip along_strike[m] down_dip[m] xc[m] yc[m] zc[m] strike[deg] dip[deg] half_len[m] half_wid[m] area[m^2] group a b sigma0[Pa]\n");
  for(k=0;k<g->np;k++){
    int ip=g->idx[k];
    fprintf(out,"%6d %14.6e %14.6e %14.6e %14.6e %14.6e %8.3f %8.3f %12.5e %12.5e %12.5e %5d %12.6e %12.6e %12.6e\n",
	    ip,g->xs[k],g->ys[k],
	    fault[ip].x[INT_X],fault[ip].x[INT_Y],fault[ip].x[INT_Z],
	    fault[ip].strike,fault[ip].dip,
	    fault[ip].l,fault[ip].w,fault[ip].area,
	    fault[ip].group,fault[ip].mu_sa,fault[ip].mu_db,sigma0);
  }
  fclose(out);
}
/*
   small comparison helper and a robust estimate of the smallest
   positive coordinate spacing, used to guess the GMT grid increment
*/
int rsf_dcmp(const void *a,const void *b)
{
  double d = *(const double *)a - *(const double *)b;
  return (d<0.0)?(-1):((d>0.0)?(1):(0));
}
double rsf_min_pos_spacing(const double *v,int n)
{
  double *c,m=1e30;int i;
  if(n < 2)return 0.0;
  c=(double *)malloc((size_t)n*sizeof(double));
  if(!c)return 0.0;
  for(i=0;i<n;i++)c[i]=v[i];
  qsort(c,(size_t)n,sizeof(double),rsf_dcmp);
  for(i=1;i<n;i++){
    double d=c[i]-c[i-1];
    if((d>1e-6)&&(d<m))m=d;	/* 1e-6 m tolerance ignores numerical duplicates */
  }
  free(c);
  return (m<1e30)?(m):(0.0);
}
/*
   along-strike and down-dip coordinates [m] of every patch in a group,
   projected onto the strike and dip unit vectors of the group's first
   patch and shifted so the minima are zero, plus an estimate of the
   grid increment, extent, and whether the patches form a single
   regular grid (so the GMT -R/-I in the header can be used directly)
*/
void rsf_group_coords(struct med *medium,struct flt *fault,struct rsf_group_grid *g)
{
  int k,l,ref,nlev;
  double ts[3],td[3],xmin=1e30,xmax=-1e30,arctot;
  double *pdip,*lev,*cx,*cy,*cz;
  int *lidx,*lcnt;
  (void)medium;
  if(g->np < 1)
    return;
  ref=g->idx[0];
  for(k=0;k<3;k++){
    ts[k]=(double)fault[ref].t_strike[k];
    td[k]=(double)fault[ref].t_dip[k];
  }
  pdip=(double *)malloc((size_t)g->np*sizeof(double));
  lidx=(int *)   malloc((size_t)g->np*sizeof(int));
  /* along-strike coordinate: projection onto the reference strike axis. This is
     regular as long as the strike is constant, which holds for a fault that
     bends only in dip. The dip-axis projection pdip is used only to identify and
     order the down-dip rows below, not as the coordinate itself. */
  for(k=0;k<g->np;k++){
    int ip=g->idx[k];
    g->xs[k]=fault[ip].x[INT_X]*ts[INT_X]+fault[ip].x[INT_Y]*ts[INT_Y]+fault[ip].x[INT_Z]*ts[INT_Z];
    pdip[k]=fault[ip].x[INT_X]*td[INT_X]+fault[ip].x[INT_Y]*td[INT_Y]+fault[ip].x[INT_Z]*td[INT_Z];
    if(g->xs[k]<xmin)xmin=g->xs[k];
  }
  for(k=0;k<g->np;k++){
    g->xs[k]-=xmin;
    if(g->xs[k]>xmax)xmax=g->xs[k];
  }
  g->xmin=0.0;g->xmax=xmax;
  g->dx=rsf_min_pos_spacing(g->xs,g->np);
  g->nx=(g->dx>0.0)?((int)floor(xmax/g->dx+0.5)+1):(0);
  /*
    down-dip coordinate: the old code projected every patch onto the dip axis of
    the reference patch, which is only correct for a planar (constant-dip) fault.
    On a listric fault the dip changes down dip, so that projection compresses
    the deeper rows unevenly and the grid is (wrongly) flagged irregular. Here
    the down-dip rows are identified by their dip-axis projection (constant
    within a row), ordered, and the coordinate is set to the true arc length
    along the fault (cumulative centroid-to-centroid distance), laid out on a
    uniform down-dip increment so the grid is exactly regular. For a planar fault
    this reduces to the original behavior.
  */
  {
    double tol,*sp;
    sp=(double *)malloc((size_t)g->np*sizeof(double));
    for(k=0;k<g->np;k++)sp[k]=pdip[k];
    qsort(sp,(size_t)g->np,sizeof(double),rsf_dcmp);	/* ascending */
    tol=rsf_min_pos_spacing(pdip,g->np)*0.25;		/* below one row gap */
    if(tol<=0.0)tol=1e-3;
    lev=(double *)malloc((size_t)g->np*sizeof(double));
    nlev=0;
    for(k=0;k<g->np;k++)
      if((nlev==0)||(sp[k]-lev[nlev-1]>tol))lev[nlev++]=sp[k];
    cx=(double *)calloc((size_t)nlev,sizeof(double));
    cy=(double *)calloc((size_t)nlev,sizeof(double));
    cz=(double *)calloc((size_t)nlev,sizeof(double));
    lcnt=(int *)calloc((size_t)nlev,sizeof(int));
    for(k=0;k<g->np;k++){
      int ip=g->idx[k],best=0;
      double bd=fabs(pdip[k]-lev[0]);
      for(l=1;l<nlev;l++){
	double d=fabs(pdip[k]-lev[l]);
	if(d<bd){bd=d;best=l;}
      }
      lidx[k]=best;
      cx[best]+=fault[ip].x[INT_X];cy[best]+=fault[ip].x[INT_Y];cz[best]+=fault[ip].x[INT_Z];
      lcnt[best]++;
    }
    for(l=0;l<nlev;l++)
      if(lcnt[l]>0){cx[l]/=lcnt[l];cy[l]/=lcnt[l];cz[l]/=lcnt[l];}
    arctot=0.0;
    for(l=1;l<nlev;l++){
      double ddx=cx[l]-cx[l-1],ddy=cy[l]-cy[l-1],ddz=cz[l]-cz[l-1];
      arctot+=sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
    }
    g->dy=(nlev>1)?(arctot/(double)(nlev-1)):(0.0);
    /* reference the down-dip coordinate to the shallow (surface) end: index 0
       from the sort above is the deepest row (smallest up-dip projection), so
       flip it here.  down_dip then increases with depth from 0 at the surface,
       which is what the plotting script's positive-down y axis expects, so cross
       sections come out with the surface at the top */
    for(k=0;k<g->np;k++)g->ys[k]=(double)((nlev-1)-lidx[k])*g->dy;
    g->ymin=0.0;g->ymax=(nlev>1)?(arctot):(0.0);
    g->ny=nlev;
    free(sp);free(lev);free(cx);free(cy);free(cz);free(lcnt);
  }
  g->regular=((g->nx>0)&&(g->ny>0)&&(g->nx*g->ny==g->np));
  free(pdip);free(lidx);
}
/*
  collect the distinct fault groups, gather each group's patch indices
  (geometry order), and compute its on-fault coordinates and grid.
  Returns the number of groups and sets *pg; rank 0 only
*/
int rsf_build_groups(struct med *medium,struct flt *fault,struct rsf_group_grid **pg)
{
  int i,j,n=medium->nrflt,ng=0;
  int *gid;
  struct rsf_group_grid *g;
  *pg=NULL;
  if(n < 1)
    return 0;
  gid=(int *)malloc((size_t)n*sizeof(int));
  if(!gid)
    return 0;
  for(i=0;i<n;i++){		/* distinct group ids, first-seen order */
    int id=fault[i].group,found=0;
    for(j=0;j<ng;j++)
      if(gid[j]==id){found=1;break;}
    if(!found)
      gid[ng++]=id;
  }
  g=(struct rsf_group_grid *)calloc((size_t)ng,sizeof(struct rsf_group_grid));
  if(!g){free(gid);return 0;}
  for(j=0;j<ng;j++){
    int np=0,k=0;
    g[j].id=gid[j];
    for(i=0;i<n;i++)
      if(fault[i].group==gid[j])np++;
    g[j].np=np;
    g[j].idx=(int    *)malloc((size_t)np*sizeof(int));
    g[j].xs =(double *)malloc((size_t)np*sizeof(double));
    g[j].ys =(double *)malloc((size_t)np*sizeof(double));
    g[j].buf=(float  *)malloc((size_t)3*(size_t)np*sizeof(float));
    if((!g[j].idx)||(!g[j].xs)||(!g[j].ys)||(!g[j].buf)){
      fprintf(stderr,"rsf_build_groups: alloc failed for group %d (%d patches)\n",gid[j],np);
      g[j].np=0;		/* leave empty so the writer skips it */
      continue;
    }
    for(i=0;i<n;i++)
      if(fault[i].group==gid[j])g[j].idx[k++]=i;
    rsf_group_coords(medium,fault,g+j);
  }
  free(gid);
  *pg=g;
  return ng;
}
void rsf_free_groups(struct rsf_group_grid *g,int ng)
{
  int j;
  if(!g)
    return;
  for(j=0;j<ng;j++){
    if(g[j].idx)free(g[j].idx);
    if(g[j].xs)free(g[j].xs);
    if(g[j].ys)free(g[j].ys);
    if(g[j].buf)free(g[j].buf);
  }
  free(g);
}
#endif
