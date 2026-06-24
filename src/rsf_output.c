#include "interact.h"
#include "properties.h"
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#include "rsf.h"


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
					  int ngroup,double *vbuf)
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
  /* compact per-fault slip-rate field output state */
  uc->field_enable = field_enable;
  uc->field_step_interval = field_step_interval;
  uc->field_tmin = field_tmin;
  uc->field_frame = 0;
  uc->fout_field_times = NULL;
  uc->groups = groups;
  uc->ngroup = ngroup;
  uc->vbuf = vbuf;
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
    if(uc->field_enable){
      int ierr_dir = system("mkdir -p tmp_rsf");
      if(ierr_dir)
	fprintf(stderr,"rsf_init_monitor_and_event: WARNING: could not make tmp_rsf directory\n");
      uc->fout_field_times = fopen("rsf_vel.times","w");
      if(uc->fout_field_times){
	fprintf(uc->fout_field_times,"# frame step time[yr] time[s] log10(max|v|[m/s])\n");
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
    if(uc->fout_field_times){
      fclose(uc->fout_field_times);
      fprintf(stderr,"rsf_finalize_monitor_and_event: wrote %i slip-rate field frame(s) across %i group(s) to tmp_rsf/\n",
	      uc->field_frame,uc->ngroup);
    }
    rsf_free_groups(uc->groups,uc->ngroup);
    uc->groups = NULL;uc->ngroup = 0;
    if(uc->vbuf){free(uc->vbuf);uc->vbuf = NULL;}
  }
  PetscCall(VecDestroy(&uc->Xold));
  PetscCall(VecScatterDestroy(&uc->gather));
  PetscCall(VecDestroy(&uc->gathered));
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
  struct rsf_vars *rsf;
  const PetscReal sec_per_year = 365.25*24.*60.*60.;
  PetscFunctionBeginUser;
  uc = (struct rsf_out_ctx *)ptr;
  medium = uc->par->medium;
  rsf = medium->rsf;
  fault = uc->par->fault;
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
      for (i = medium->rs, j=0; i < medium->re; i++, j+=rsf->dim) {
	v = vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,rsf->v0,&d1,&d2,&d3,medium);
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
	    i=0,j=0;i < n;i++,j += rsf->dim){
	fault[i].s[STRIKE] = values[j+1]; /* shear stress */
	fault[i].s[NORMAL] = values[j+2]; /* normal stress */
	fault[i].u[STRIKE] = vel_from_rsf(values[j+1],values[j+2],values[j],fault[i].mu_s,
					  rsf->v0,&du1,&du2,&du3,medium);
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
    if(medium->comm_rank == 0){
      PetscReal du1,du2,du3,vmax=0.0;
      int ig,k,ip,nf = medium->nrflt;
      PetscCall(VecGetArray(uc->gathered,&fvals));
      for(ip=0;ip < nf;ip++){
	PetscReal vv = vel_from_rsf(fvals[ip*rsf->dim+1],fvals[ip*rsf->dim+2],
				    fvals[ip*rsf->dim+0],fault[ip].mu_s,rsf->v0,
				    &du1,&du2,&du3,medium);
	uc->vbuf[ip] = (double)vv;
	if(fabs(vv) > vmax)vmax = fabs(vv);
      }
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
	fb = fopen(ffile,"wb");
	if(fb){
	  fwrite(g->buf,sizeof(float),(size_t)3*(size_t)g->np,fb);
	  fclose(fb);
	}
      }
      if(uc->fout_field_times)
	fprintf(uc->fout_field_times,"%6i %9i %.8f %.8e %.6f\n",
		uc->field_frame,(int)step,time/sec_per_year,time,log10((vmax>0.0)?(vmax):(1e-30)));
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
  const PetscReal sec_per_year = 365.25*24.*60.*60.;
  PetscFunctionBeginUser;
  uc = (struct rsf_out_ctx *)ctx;
  medium = uc->par->medium;fault = uc->par->fault;
  rsf = medium->rsf;
  uc->slipping = (uc->slipping)?(PETSC_FALSE):(PETSC_TRUE); /* threshold crossed */
  if(t >= uc->event_tmin){
    lsum[0]=lsum[1]=0.0;lvmax = 0.0;
    PetscCall(VecGetArrayRead(X,&x));
    for (i = medium->rs, j=0; i < medium->re; i++, j+=rsf->dim) {
      v = fabs(vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,rsf->v0,&d1,&d2,&d3,medium));
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
  out=fopen(fname,"w");
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
	    (double)fault[ip].x[INT_X],(double)fault[ip].x[INT_Y],(double)fault[ip].x[INT_Z],
	    (double)fault[ip].strike,(double)fault[ip].dip,
	    (double)fault[ip].l,(double)fault[ip].w,(double)fault[ip].area,
	    fault[ip].group,(double)fault[ip].mu_s,(double)fault[ip].mu_d,sigma0);
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
  int k,ref;
  double ts[3],td[3],xmin=1e30,ymin=1e30,xmax=-1e30,ymax=-1e30;
  (void)medium;
  if(g->np < 1)
    return;
  ref=g->idx[0];
  for(k=0;k<3;k++){
    ts[k]=(double)fault[ref].t_strike[k];
    td[k]=(double)fault[ref].t_dip[k];
  }
  for(k=0;k<g->np;k++){
    int ip=g->idx[k];
    g->xs[k]=fault[ip].x[INT_X]*ts[INT_X]+fault[ip].x[INT_Y]*ts[INT_Y]+fault[ip].x[INT_Z]*ts[INT_Z];
    g->ys[k]=fault[ip].x[INT_X]*td[INT_X]+fault[ip].x[INT_Y]*td[INT_Y]+fault[ip].x[INT_Z]*td[INT_Z];
    if(g->xs[k]<xmin)xmin=g->xs[k];
    if(g->ys[k]<ymin)ymin=g->ys[k];
  }
  for(k=0;k<g->np;k++){
    g->xs[k]-=xmin;
    g->ys[k]-=ymin;
    if(g->xs[k]>xmax)xmax=g->xs[k];
    if(g->ys[k]>ymax)ymax=g->ys[k];
  }
  g->xmin=0.0;g->xmax=xmax;g->ymin=0.0;g->ymax=ymax;
  g->dx=rsf_min_pos_spacing(g->xs,g->np);
  g->dy=rsf_min_pos_spacing(g->ys,g->np);
  g->nx=(g->dx>0.0)?((int)floor(xmax/g->dx+0.5)+1):(0);
  g->ny=(g->dy>0.0)?((int)floor(ymax/g->dy+0.5)+1):(0);
  g->regular=((g->nx>0)&&(g->ny>0)&&(g->nx*g->ny==g->np));
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
