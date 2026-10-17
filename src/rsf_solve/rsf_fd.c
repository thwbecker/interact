/*
  fully dynamic extension of rsf_solve

  M1 (2026-10-17): option parsing, rsf->fd allocation, monitor that logs
  threshold crossings, -fd_mode -1 hand-over exercise (TSRestartStep
  from inside a monitor).
  M2 (2026-10-18): segment detection on the 2-D geometry, dynamic
  correction tables per segment pair (fd_tables.c, FFT storage), gate of
  the tables' static limit against the assembled Is operator (sign and
  convention check, as rsf_ve does for external kernels). The dynamic
  marcher arrives with M3/M4; see fully_dynamic/fd_c_design.md.

  Contract: with -fd_mode absent nothing here is called from the time
  stepping (rsf_fd_setup leaves rsf->fd NULL and the monitor is not
  registered), so the quasi-dynamic code path is unchanged.

  options:
    -fd_mode <0|1|-1>   0 off (default), 1 dynamic phases, -1 hand-over test
    -fd_vsw <m/s>       switch threshold (default 1e-3)
    -fd_vback <m/s>     return threshold (default 1e-4)
    -fd_tol <f>         window tolerance (default 1e-2)
    -fd_et <f>          collocation offset (default 1.0)
    -fd_xi <f>          variable-step factor (default 0.25)
    -fd_kmax <n>        largest substep multiple (default 64)
    -fd_cfl <f>         beta dt_min / dx (default 0.5)
    -fd_lmax <n>        cap on the window in fine-grid lags (default 0: prototype bound)
    -fd_gate_tol <f>    static gate tolerance (default 0.05; abort at 5x)
    -fd_length_scale    geometry unit in metres (default 1.0)

  geometry restrictions of this version (refused with a message): 2-D
  elements (types 30 to 33) that form straight runs of equal, collinear
  elements with one common element size; full plane for in-plane, full or
  half plane for antiplane (the half-plane image is added in M4).
*/
#include "interact.h"
#include "properties.h"
#include "petsc_prototypes.h"
#include "rsf.h"
#include "fd.h"

static PetscErrorCode fd_build_tables(struct interact_ctx *par, struct fd_vars *fd);
static PetscErrorCode fd_gate_static(struct interact_ctx *par, struct fd_vars *fd, PetscReal scale);

PetscErrorCode rsf_fd_setup(struct interact_ctx *par, PetscReal shear_modulus_si, PetscReal s_wave_speed_si)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  struct rsf_vars *rsf = medium->rsf;
  struct fd_vars *fd;
  PetscInt mode = 0, i, n2d = 0, nanti = 0;
  PetscBool set = PETSC_FALSE;
  PetscMPIInt size;
  PetscFunctionBeginUser;
  rsf->fd = NULL;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_mode",&mode,&set));
  if((!set) || (mode == 0))
    PetscFunctionReturn(PETSC_SUCCESS);
  /* restrictions of the first version */
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
  PetscCheck(size == 1,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: one MPI rank only in this version");
  PetscCheck(rsf->ve_np == 0,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: not combined with -ve_mode in this version");
  for(i = 0;i < medium->nrflt;i++){
    if((fault[i].type == TWO_DIM_SEGMENT_PLANE_STRAIN) || (fault[i].type == TWO_DIM_HALFPLANE_PLANE_STRAIN))
      n2d++;
    else if(fault[i].type == TWO_DIM_ANTIPLANE)
      nanti++;
  }
  PetscCheck((n2d + nanti) == medium->nrflt,PETSC_COMM_WORLD,PETSC_ERR_SUP,
	     "-fd_mode: all elements must be 2-D segments (plane strain or antiplane) in this version");
  PetscCheck((n2d == 0) || (nanti == 0),PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: mixed in-plane and antiplane elements");
  PetscCall(PetscCalloc1(1,&fd));
  fd->mode = mode;
  fd->vsw = 1e-3; fd->vback = 1e-4; fd->tol = 1e-2; fd->et = 1.0; fd->xi = 0.25; fd->kmax = 64;
  fd->cfl = 0.5; fd->lmax = 0; fd->gate_tol = 0.05; fd->length_scale = 1.0;
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_vsw",&fd->vsw,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_vback",&fd->vback,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_tol",&fd->tol,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_et",&fd->et,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_xi",&fd->xi,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_kmax",&fd->kmax,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_cfl",&fd->cfl,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_lmax",&fd->lmax,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_gate_tol",&fd->gate_tol,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_length_scale",&fd->length_scale,NULL));
  PetscCheck(fd->vback < fd->vsw,PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"-fd_vback must be below -fd_vsw");
  PetscCheck(fd->kmax >= 1 && (fd->kmax & (fd->kmax - 1)) == 0,PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"-fd_kmax must be a power of two");
  /* medium: shear modulus and S speed from the run settings (SI), Poisson from interact */
  fd->mode2d = (nanti > 0) ? FD_MODE_ANTIPLANE : FD_MODE_INPLANE;
  fd->halfplane = PETSC_FALSE;
  for(i = 0;i < medium->nrflt;i++)
    if(fault[i].type == TWO_DIM_HALFPLANE_PLANE_STRAIN)
      fd->halfplane = PETSC_TRUE;
  if((fd->mode2d == FD_MODE_ANTIPLANE) && (!medium->full_space))
    fd->halfplane = PETSC_TRUE;
  PetscCheck(!fd->halfplane,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: half-plane geometries are M4 work; use a full-plane geometry for now");
  fd->mu = shear_modulus_si; fd->beta = s_wave_speed_si;
  fd->nu = medium->elastic.poisson;
  fd->rho = fd->mu/(fd->beta*fd->beta);
  fd->alpha = fd->beta*sqrt(2.0*(1.0 - fd->nu)/(1.0 - 2.0*fd->nu));
  fd->nphase = 0; fd->nflag = 0; fd->t_last_flag = -1.0; fd->restart_done = PETSC_FALSE;
  rsf->fd = fd;
  PetscCall(fd_build_tables(par,fd));
  PetscCall(fd_gate_static(par,fd,shear_modulus_si/medium->elastic.shear));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
		"rsf_fd_setup: fully dynamic extension, mode %d (M2: tables only, no marcher yet), %s, %ld segments, dt_min %.3e s, "
		"vsw %.1e vback %.1e m/s, tol %.1e, et %.2f, xi %.2f, kmax %d, tables %.1f MB, fft %s\n",
		(int)fd->mode,fd->mode2d == FD_MODE_ANTIPLANE ? "antiplane" : "in-plane",(long)fd->segs->nseg,(double)fd->dt_min,
		(double)fd->vsw,(double)fd->vback,(double)fd->tol,(double)fd->et,(double)fd->xi,(int)fd->kmax,
		(double)fd->nbytes/1048576.0,fd_fft_backend()));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* segments from the 2-D geometry and one pair table per (source, receiver) */
static PetscErrorCode fd_build_tables(struct interact_ctx *par, struct fd_vars *fd)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  PetscInt n = medium->nrflt, i, k, l, nm;
  double *cx,*cy,*tx,*ty,*hl,dx0 = -1.0,tn;
  char err[256];
  PetscFunctionBeginUser;
  PetscCall(PetscMalloc5(n,&cx,n,&cy,n,&tx,n,&ty,n,&hl));
  for(i = 0;i < n;i++){
    cx[i] = fault[i].x[INT_X]*fd->length_scale; cy[i] = fault[i].x[INT_Y]*fd->length_scale;
    tn = hypot(fault[i].t_strike[INT_X],fault[i].t_strike[INT_Y]);
    tx[i] = fault[i].t_strike[INT_X]/tn; ty[i] = fault[i].t_strike[INT_Y]/tn;
    hl[i] = fault[i].l*fd->length_scale;
  }
  fd->segs = fd_detect_segments((int)n,cx,cy,tx,ty,hl,1e-6,err);
  PetscCall(PetscFree5(cx,cy,tx,ty,hl));
  PetscCheck(fd->segs,PETSC_COMM_WORLD,PETSC_ERR_SUP,"rsf_fd_setup: segment detection failed: %s",err);
  for(k = 0;k < fd->segs->nseg;k++){
    if(dx0 < 0.0) dx0 = fd->segs->seg[k].dx;
    PetscCheck(fabs(fd->segs->seg[k].dx - dx0) < 1e-9*dx0,PETSC_COMM_WORLD,PETSC_ERR_SUP,
	       "rsf_fd_setup: segments with different element sizes (%g vs %g): one common size required in this version",
	       fd->segs->seg[k].dx,dx0);
    for(l = 0;l < fd->segs->nseg;l++)
      PetscCheck(fd_segments_parallel(fd->segs->seg + k,fd->segs->seg + l),PETSC_COMM_WORLD,PETSC_ERR_SUP,
		 "rsf_fd_setup: segments %ld and %ld are not parallel: non-parallel pairs are M7 work",(long)k,(long)l);
  }
  fd->dt_min = fd->cfl*dx0/fd->beta;
  /* one padded FFT length for all pairs: smallest power of two >= 2 nmax - 1 */
  nm = 0;
  for(k = 0;k < fd->segs->nseg;k++)
    if(fd->segs->seg[k].n > nm) nm = fd->segs->seg[k].n;
  fd->M = fd_fft_size_for(2*nm - 1);
  PetscCall(PetscCalloc1((size_t)fd->segs->nseg*fd->segs->nseg,&fd->pairs));
  fd->nbytes = 0; fd->L = 0;
  for(k = 0;k < fd->segs->nseg;k++)		/* source */
    for(l = 0;l < fd->segs->nseg;l++){		/* receiver */
      struct fd_pair *p = fd_pair_build(fd->segs->seg + k,fd->segs->seg + l,fd->mode2d,fd->dt_min,fd->et,
					fd->alpha,fd->beta,fd->rho,fd->tol,(int)fd->lmax,fd->lmax > 0,fd->M,err);
      PetscCheck(p,PETSC_COMM_WORLD,PETSC_ERR_SUP,"rsf_fd_setup: table for pair (%ld,%ld) failed: %s",(long)k,(long)l,err);
      fd->pairs[k*fd->segs->nseg + l] = p;
      fd->nbytes += (double)p->L*p->M*16.0*(p->has_normal ? 2 : 1);
      if(p->L > fd->L) fd->L = p->L;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: pair (%ld -> %ld): %d x %d elements, window %d lags (%.2f s)\n",
			    (long)k,(long)l,p->src->n,p->rcv->n,p->L,(double)(p->L*fd->dt_min)));
    }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   compare the static limit of the tables (stress-from-slip, SI) with the
   assembled Is (stress-drop convention, interact units times scale):
   Is*scale should equal -Ks. Near-diagonal entries of the same segment
   are skipped, as in the VE gate, because the self terms may be
   regularised differently. A pure sign error gives a deviation of 2.
*/
static PetscErrorCode fd_gate_static(struct interact_ctx *par, struct fd_vars *fd, PetscReal scale)
{
  struct med *medium = par->medium;
  PetscInt k, l, i, j, gi, gj, ncols, c;
  const PetscInt *cols;
  const PetscScalar *vals;
  double *Ks, *Kn, dev = 0.0, nrm = 0.0, dv, *row;
  int ns, nr;
  PetscFunctionBeginUser;
  PetscCall(PetscCalloc1(medium->nrflt,&row));
  for(gi = 0;gi < medium->nrflt;gi++){
    /* assembled row, scaled */
    for(c = 0;c < medium->nrflt;c++) row[c] = 0.0;
    PetscCall(MatGetRow(medium->Is,gi,&ncols,&cols,&vals));
    for(c = 0;c < ncols;c++) row[cols[c]] = PetscRealPart(vals[c])*scale;
    PetscCall(MatRestoreRow(medium->Is,gi,&ncols,&cols,&vals));
    /* our value for every source element gj: locate (segment, element) of gi and gj */
    for(k = 0;k < fd->segs->nseg;k++)
      for(i = 0;i < fd->segs->seg[k].n;i++)
	if(fd->segs->gidx[k][i] == gi)
	  for(l = 0;l < fd->segs->nseg;l++){
	    struct fd_pair *p = fd->pairs[l*fd->segs->nseg + k];	/* source l, receiver k */
	    ns = p->src->n; nr = p->rcv->n;
	    PetscCall(PetscMalloc2((size_t)ns*nr,&Ks,(size_t)ns*nr,&Kn));
	    fd_pair_static_dense(p,Ks,Kn);
	    for(j = 0;j < ns;j++){
	      gj = fd->segs->gidx[l][j];
	      if((l == k) && (labs((long)i - (long)j) <= 1))
		continue;
	      dv = row[gj] - (-Ks[i*ns + j]);
	      dev += dv*dv; nrm += row[gj]*row[gj];
	    }
	    PetscCall(PetscFree2(Ks,Kn));
	  }
  }
  PetscCall(PetscFree(row));
  dev = sqrt(dev/(nrm + 1e-300));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: static gate: rel dev of -Ks(tables) from Is*scale %.3e (tol %g)\n",dev,(double)fd->gate_tol));
  PetscCheck(dev <= 5.0*fd->gate_tol,PETSC_COMM_WORLD,PETSC_ERR_USER,
	     "rsf_fd_setup: static limit of the dynamic tables deviates %.3e from the assembled operator: convention mismatch (a sign flip gives 2)",dev);
  if(dev > fd->gate_tol)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: WARNING: static gate deviation %.3e exceeds %g\n",dev,(double)fd->gate_tol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  accepted-step monitor. rsf->vel holds the slip rates of the last RHS
  evaluation; for the switch decision that is close enough to the
  accepted state (the prototype found the threshold itself immaterial
  between 1e-3 and 1e-2 m/s)
*/
PetscErrorCode rsf_fd_monitor(TS ts, PetscInt step, PetscReal time, Vec X, void *ctx)
{
  struct interact_ctx *par = (struct interact_ctx *)ctx;
  struct rsf_vars *rsf = par->medium->rsf;
  struct fd_vars *fd = rsf->fd;
  PetscReal vmax;
  PetscInt imax;
  PetscFunctionBeginUser;
  if(!fd)
    PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(VecMax(rsf->vel,&imax,&vmax));
  if(vmax > fd->vsw){
    fd->nflag++;
    if((fd->nflag == 1) || (time - fd->t_last_flag > 86400.0))   /* log the onset, not every step */
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,
		  "rsf_fd_monitor: step %ld t %.6f yr: max slip rate %.3e m/s at patch %ld above vsw %.1e (marcher not yet implemented)\n",
		  (long)step,(double)(time/SEC_PER_YEAR),(double)vmax,(long)imax,(double)fd->vsw));
    fd->t_last_flag = time;
    if((fd->mode == -1) && (!fd->restart_done)){
      /* hand-over exercise: restart the TS from the current (unchanged)
	 state with a small step, as the marcher will do after writing
	 the post-event state into X; the run must continue unchanged
	 apart from the step sequence */
      PetscCall(TSSetTimeStep(ts,1e-3));
      PetscCall(TSRestartStep(ts));
      fd->restart_done = PETSC_TRUE;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_monitor: mode -1: TSRestartStep issued at t %.6f yr\n",(double)(time/SEC_PER_YEAR)));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode rsf_fd_finalize(struct interact_ctx *par)
{
  struct rsf_vars *rsf = par->medium->rsf;
  struct fd_vars *fd = rsf->fd;
  PetscInt k;
  PetscFunctionBeginUser;
  if(!fd)
    PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_finalize: %ld accepted steps above vsw, %ld dynamic phases\n",
		(long)fd->nflag,(long)fd->nphase));
  if(fd->pairs){
    for(k = 0;k < fd->segs->nseg*fd->segs->nseg;k++)
      fd_pair_free(fd->pairs[k]);
    PetscCall(PetscFree(fd->pairs));
  }
  fd_seglist_free(fd->segs);
  PetscCall(PetscFree(fd));
  rsf->fd = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}
