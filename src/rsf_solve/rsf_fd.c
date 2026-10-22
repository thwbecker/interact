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
    -fd_max_steps <n>   abort a dynamic phase beyond this many steps (default 2e6)
    -fd_dt_restart <s>  TS step after a phase (default 1e-2)
    -fd_conv <0|1|2>    history convolution: 0 direct sum, 1 partitioned FFT (about twice
                        the memory), 2 hybrid: partitioned for short substeps, direct for
                        long ones (about three times the memory; default)
    -fd_conv_block <n>  partition length in lags, a power of two (default 0: automatic)
    -fd_coarse <0|1>    1 (default): steps of kmax substeps start at multiples of kmax,
                        and once a pair's whole window consists of such steps its history
                        is summed over blocks of kmax with aggregated kernels (exact, kmax
                        times cheaper); 0: off
    -fd_back_eps <f>    0 (default): return to QD after the slip rate has stayed below
                        -fd_vback for one window; > 0: return once every slip rate has
                        changed by less than this fraction over one window and the
                        maximum is below -fd_back_vmax (default 0.5 vsw); see rsf_solve_fd.md
    -fd_stations i,j,.. global patch indices whose t, slip, slip rate, tau are written
                        during every phase to fd_phase_<n>_st<i>.dat; the final slip of
                        each phase goes to fd_phase_<n>_slip.dat

  geometry restrictions of this version (refused with a message): 2-D
  elements (types 30 to 33) that form straight runs of equal, collinear
  elements with one common element size; full plane for in-plane, full
  plane or half plane with vertical segments for antiplane; aging law;
  one MPI rank; no -ve_mode.
*/
#include "interact.h"
#include "properties.h"
#include "petsc_prototypes.h"
#include "rsf.h"
#include "fd.h"
#include "fd_march.h"
#include "fd3d.h"

static PetscErrorCode fd_build_tables(struct interact_ctx *par, struct fd_vars *fd);
static PetscErrorCode fd_gate_static(struct interact_ctx *par, struct fd_vars *fd);

PetscErrorCode rsf_fd_setup(struct interact_ctx *par, PetscReal shear_modulus_si, PetscReal s_wave_speed_si)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  struct rsf_vars *rsf = medium->rsf;
  struct fd_vars *fd;
  PetscInt mode = 0, i, n2d = 0, nanti = 0, n3d = 0;
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
    else if(fault[i].type == OKADA_PATCH)
      n3d++;
  }
  PetscCheck((n3d == 0) || (n3d == medium->nrflt),PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: mixed 2-D and 3-D elements");
  PetscCheck((n3d > 0) || ((n2d + nanti) == medium->nrflt),PETSC_COMM_WORLD,PETSC_ERR_SUP,
	     "-fd_mode: all elements must be 2-D segments (plane strain or antiplane) in this version");
  PetscCheck((n2d == 0) || (nanti == 0),PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: mixed in-plane and antiplane elements");
  PetscCall(PetscCalloc1(1,&fd));
  fd->mode = mode;
  fd->vsw = 1e-3; fd->vback = 1e-4; fd->tol = 1e-2; fd->et = 1.0; fd->xi = 0.25; fd->kmax = 64;
  fd->cfl = 0.5; fd->lmax = 0; fd->gate_tol = 0.05; fd->length_scale = 1.0;
  fd->max_steps = 2000000; fd->dt_restart = 1e-2;
  fd->conv = 2; fd->conv_block = 0;
  fd->back_eps = 0.0; fd->back_vmax = -1.0; fd->coarse = 1;
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
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_max_steps",&fd->max_steps,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_dt_restart",&fd->dt_restart,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_conv",&fd->conv,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_conv_block",&fd->conv_block,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_back_eps",&fd->back_eps,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_coarse",&fd->coarse,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_back_vmax",&fd->back_vmax,NULL));
  if(fd->back_vmax <= 0.0) fd->back_vmax = 0.5*fd->vsw;
  fd->nl = 0;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_nl",&fd->nl,NULL));
  fd->nstat = 16;
  PetscCall(PetscOptionsGetIntArray(NULL,NULL,"-fd_stations",fd->stat,&fd->nstat,&set));
  if(!set) fd->nstat = 0;
  PetscCheck(fd->vback < fd->vsw,PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"-fd_vback must be below -fd_vsw");
  PetscCheck(fd->kmax >= 1 && fd->kmax <= 32768 && (fd->kmax & (fd->kmax - 1)) == 0,PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"-fd_kmax must be a power of two up to 32768");
  if(n3d > 0){
    /* 3-D: square elements on one planar fault in a whole space (fd3d.c); so far
       setup and static gate only, the dynamic phases are not yet connected */
    fd->dim3 = PETSC_TRUE;
    fd->mu = shear_modulus_si; fd->beta = s_wave_speed_si;
    fd->nu = medium->elastic.poisson;
    fd->rho = fd->mu/(fd->beta*fd->beta);
    fd->alpha = fd->beta*sqrt(2.0*(1.0 - fd->nu)/(1.0 - 2.0*fd->nu));
    rsf->fd = fd;
    PetscCall(fd3d_setup_from_geometry(par,fd));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D fully dynamic extension, mode %d, %s, vsw %.1e vback %.1e m/s, fixed fine steps, %s time marching\n",
			  (int)fd->mode,fd->mirror3 ? "half space (mirror-image free surface, exact static Is)" : "whole space",(double)fd->vsw,(double)fd->vback,fd->nl ? "second-order (-fd_nl 1)" : "first-order"));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCheck(!fd->nl || !rsf->calc_sigma_dot,PETSC_COMM_WORLD,PETSC_ERR_SUP,
	     "-fd_nl 1: the second-order time marching does not include normal-stress coupling yet (drop -calc_sigma_dot or use -fd_nl 0)");
  /* medium: shear modulus and S speed from the run settings (SI), Poisson from interact */
  fd->mode2d = (nanti > 0) ? FD_MODE_ANTIPLANE : FD_MODE_INPLANE;
  fd->halfplane = PETSC_FALSE;
  for(i = 0;i < medium->nrflt;i++)
    if(fault[i].type == TWO_DIM_HALFPLANE_PLANE_STRAIN)
      fd->halfplane = PETSC_TRUE;
  if((fd->mode2d == FD_MODE_ANTIPLANE) && (!medium->full_space))
    fd->halfplane = PETSC_TRUE;
  PetscCheck(!(fd->halfplane && (fd->mode2d == FD_MODE_INPLANE)),PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: in-plane half-plane geometries are not supported (M7)");
  PetscCheck(rsf->state_law == RSF_AGING_LAW,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: aging law only (-state_law 1) in this version");
  fd->mu = shear_modulus_si; fd->beta = s_wave_speed_si;
  fd->nu = medium->elastic.poisson;
  fd->rho = fd->mu/(fd->beta*fd->beta);
  fd->alpha = fd->beta*sqrt(2.0*(1.0 - fd->nu)/(1.0 - 2.0*fd->nu));
  fd->nphase = 0; fd->nflag = 0; fd->t_last_flag = -1.0; fd->restart_done = PETSC_FALSE;
  rsf->fd = fd;
  PetscCall(fd_build_tables(par,fd));
  PetscCall(fd_gate_static(par,fd));
  /* accepted-step record for the history seeding */
  fd->nrec = 64; fd->irec = 0; fd->nrec_used = 0;
  PetscCall(PetscMalloc2(fd->nrec,&fd->rec_t,(size_t)fd->nrec*medium->nrflt,&fd->rec_v));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
		"rsf_fd_setup: fully dynamic extension, mode %d, %s%s, %ld segments, dt_min %.3e s, "
		"vsw %.1e vback %.1e m/s, tol %.1e, et %.2f, xi %.2f, kmax %d, tables %.1f MB, fft %s\n",
		(int)fd->mode,fd->mode2d == FD_MODE_ANTIPLANE ? "antiplane" : "in-plane",fd->halfplane ? " half plane (image sources)" : "",(long)fd->segs->nseg,(double)fd->dt_min,
		(double)fd->vsw,(double)fd->vback,(double)fd->tol,(double)fd->et,(double)fd->xi,(int)fd->kmax,
		(double)fd->nbytes/1048576.0,fd_fft_backend()));
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* static operator for the 3-D marcher: interact's exact Is (also in a half space),
   in the marcher's element order */
static void fd3d_static_is(void *ctx, const double *delta, double *out)
{
  struct interact_ctx *par = (struct interact_ctx *)ctx;
  struct fd_vars *fd = par->medium->rsf->fd;
  PetscScalar *a;
  const PetscScalar *b;
  int i, n = fd->f3->N;
  PetscCallAbort(PETSC_COMM_WORLD,VecGetArray(fd->v3a,&a));
  for(i = 0;i < n;i++) a[fd->g_of_s[i]] = delta[i];
  PetscCallAbort(PETSC_COMM_WORLD,VecRestoreArray(fd->v3a,&a));
  PetscCallAbort(PETSC_COMM_WORLD,MatMult(par->medium->Is,fd->v3a,fd->v3b));
  PetscCallAbort(PETSC_COMM_WORLD,VecGetArrayRead(fd->v3b,&b));
  for(i = 0;i < n;i++) out[i] = PetscRealPart(b[fd->g_of_s[i]]);
  PetscCallAbort(PETSC_COMM_WORLD,VecRestoreArrayRead(fd->v3b,&b));
}

/* 3-D dynamic phase (fd3d.c), counterpart of fd_run_phase: fixed fine steps */
static PetscErrorCode fd3d_run_phase(TS ts, struct interact_ctx *par, struct fd_vars *fd, PetscReal t, Vec X, PetscReal *t_end)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  struct rsf_vars *rsf = medium->rsf;
  struct fd3d *f = fd->f3;
  PetscScalar *x;
  int N = f->N, i, g, l, nsteps = 0, L = fd->tab3->L;
  double dt = fd->tab3->dt, eta = 0.5*fd->mu/fd->beta, vmax, vpeak = 0.0, below = 0.0, slipsum = 0.0, *vtmp, *load, t_wall0, t_wall1;
  PetscFunctionBeginUser;
  PetscCall(PetscTime(&t_wall0));
  PetscCall(PetscMalloc2(N,&vtmp,N,&load));
  fd3d_reset(f);
  PetscCall(VecGetArray(X,&x));
  vmax = 0.0;
  for(i = 0;i < N;i++){
    g = fd->g_of_s[i];
    f->sigma[i] = x[4*g+2];
    f->V[i] = 2.0*rsf->v0*exp(-x[4*g]/fault[g].mu_sa)*sinh(x[4*g+1]/(fault[g].mu_sa*x[4*g+2]));
    f->theta[i] = f->Lc[i]/f->V0*exp((x[4*g] - f->f0)/f->b[i]);       /* psi = f0 + b ln(V0 theta/L) */
    f->tau0[i] = x[4*g+1] + eta*f->V[i];                               /* elastic stress at the switch */
    if(fabs(f->V[i]) > vmax) vmax = fabs(f->V[i]);
  }
  /* history from the accepted-step record (marcher order) */
  for(l = L;l >= 1;l--){
    double tl = t - l*dt; int k, i0 = -1, i1 = -1; double w;
    for(k = 0;k < fd->nrec_used;k++){
      int idx = (fd->irec - fd->nrec_used + k + fd->nrec) % fd->nrec;
      if(fd->rec_t[idx] <= tl) i0 = idx;
      if((fd->rec_t[idx] >= tl) && (i1 < 0)) i1 = idx;
    }
    if(fd->nrec_used == 0){ for(i = 0;i < N;i++) vtmp[i] = 0.0; }
    else{
      if(i0 < 0) i0 = (fd->irec - fd->nrec_used + fd->nrec) % fd->nrec;
      if(i1 < 0) i1 = (fd->irec - 1 + fd->nrec) % fd->nrec;
      w = (fd->rec_t[i1] > fd->rec_t[i0]) ? (tl - fd->rec_t[i0])/(fd->rec_t[i1] - fd->rec_t[i0]) : 0.0;
      if(w < 0.0) w = 0.0;
      if(w > 1.0) w = 1.0;
      for(i = 0;i < N;i++) vtmp[i] = (1.0 - w)*fd->rec_v[(size_t)i0*N + i] + w*fd->rec_v[(size_t)i1*N + i];
    }
    fd3d_push(f,vtmp);
  }
  {
    int a;
    for(a = 0;a < fd->nstat;a++){
      char fn[128]; sprintf(fn,"fd_phase_%ld_st%ld.dat",(long)fd->nphase + 1,(long)fd->stat[a]);
      fd->fstat[a] = fopen(fn,"w"); fd->sstat[a] = -1;
      for(i = 0;i < N;i++) if(fd->g_of_s[i] == fd->stat[a]) fd->sstat[a] = i;
    }
  }
  for(;;){
    double trel = f->time + dt;
    int a;
    for(i = 0;i < N;i++) load[i] = fault[fd->g_of_s[i]].sinc[0]*trel;    /* backslip loading since the switch */
    fd3d_step(f,load);
    nsteps++;
    for(a = 0;a < fd->nstat;a++)
      if(fd->sstat[a] >= 0)
	fprintf(fd->fstat[a],"%.10e %.10e %.10e %.10e\n",f->time,f->delta[fd->sstat[a]],f->V[fd->sstat[a]],f->tau[fd->sstat[a]]);
    vmax = 0.0;
    for(i = 0;i < N;i++) if(fabs(f->V[i]) > vmax) vmax = fabs(f->V[i]);
    if(vmax > vpeak) vpeak = vmax;
    if(vmax < fd->vback) below += dt; else below = 0.0;
    if(below >= L*dt) break;
    PetscCheck(nsteps < fd->max_steps,PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,
	       "rsf_fd: 3-D phase exceeded -fd_max_steps = %ld steps (t - t_switch %g s, max V %g)",(long)fd->max_steps,f->time,vmax);
  }
  { int a; for(a = 0;a < fd->nstat;a++) if(fd->fstat[a]){ fclose(fd->fstat[a]); fd->fstat[a] = NULL; } }
  for(i = 0;i < N;i++){
    g = fd->g_of_s[i];
    x[4*g] = f->f0 + f->b[i]*log(f->V0*f->theta[i]/f->Lc[i]);
    x[4*g+1] = f->tau[i];
    x[4*g+3] += f->delta[i];
    slipsum += f->delta[i];
  }
  PetscCall(VecRestoreArray(X,&x));
  *t_end = t + f->time;
  PetscCall(PetscTime(&t_wall1));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd: 3-D phase %ld: t %.8f yr, %d steps over %.1f s, peak slip rate %.3f m/s, mean slip %.4f m, wall %.1f s\n",
			(long)fd->nphase + 1,(double)(t/SEC_PER_YEAR),nsteps,f->time,vpeak,slipsum/N,(double)(t_wall1 - t_wall0)));
  PetscCall(PetscFree2(vtmp,load));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  3-D: square elements of equal size forming a complete rectangular grid on
  one plane, whole space, slip along strike (-rsf_slip_mode strike) or dip.
  Builds the fd3d tables (grid axis 1 = slip direction) and checks their
  static limit against the assembled Is, as for 2-D.
*/
PetscErrorCode fd3d_setup_from_geometry(struct interact_ctx *par, struct fd_vars *fd)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  struct rsf_vars *rsf = medium->rsf;
  PetscInt n = medium->nrflt, k, c, ncols, gi;
  const PetscInt *cols;
  const PetscScalar *vals;
  double dx, e1[3], e2[3], nn[3], x0[3], lo1 = 1e300, lo2 = 1e300, hi1 = -1e300, hi2 = -1e300, *pa, *pb, dev = 0.0, nrm = 0.0, ls;
  int nx, nz, *cnt, *s_of_g, i;
  PetscFunctionBeginUser;
  PetscCheck(medium->full_space || (rsf->slip_mode == STRIKE),PETSC_COMM_WORLD,PETSC_ERR_SUP,
	     "-fd_mode: 3-D in a half space needs a vertical fault with strike slip (mirror-image free surface)");
  PetscCheck((rsf->slip_mode == STRIKE) || (rsf->slip_mode == DIP),PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: 3-D needs a fixed slip direction (-rsf_slip_mode strike or dip)");
  ls = fd->length_scale;
  dx = 2.0*fault[0].l*ls;
  for(k = 0;k < 3;k++){
    e1[k] = (rsf->slip_mode == STRIKE) ? fault[0].t_strike[k] : fault[0].t_dip[k];
    e2[k] = (rsf->slip_mode == STRIKE) ? fault[0].t_dip[k] : fault[0].t_strike[k];
    nn[k] = fault[0].normal[k]; x0[k] = fault[0].x[k]*ls;
  }
  PetscCall(PetscMalloc2(n,&pa,n,&pb));
  for(i = 0;i < n;i++){
    double d[3], dn;
    PetscCheck((fabs(2.0*fault[i].l*ls - dx) < 1e-6*dx) && (fabs(2.0*fault[i].w*ls - dx) < 1e-6*dx),PETSC_COMM_WORLD,PETSC_ERR_SUP,
	       "-fd_mode: 3-D elements must be squares of one size (patch %d: half-lengths %g, %g)",i,fault[i].l,fault[i].w);
    PetscCheck(fabs(fault[i].normal[0]*nn[0] + fault[i].normal[1]*nn[1] + fault[i].normal[2]*nn[2]) > 1.0 - 1e-9,PETSC_COMM_WORLD,PETSC_ERR_SUP,
	       "-fd_mode: 3-D elements must lie on one plane (patch %d has a different orientation)",i);
    for(k = 0;k < 3;k++) d[k] = fault[i].x[k]*ls - x0[k];
    dn = d[0]*nn[0] + d[1]*nn[1] + d[2]*nn[2];
    PetscCheck(fabs(dn) < 1e-6*dx,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: 3-D elements must lie on one plane (patch %d is %g m off)",i,dn);
    pa[i] = d[0]*e1[0] + d[1]*e1[1] + d[2]*e1[2]; pb[i] = d[0]*e2[0] + d[1]*e2[1] + d[2]*e2[2];
    if(pa[i] < lo1) lo1 = pa[i];
    if(pa[i] > hi1) hi1 = pa[i];
    if(pb[i] < lo2) lo2 = pb[i];
    if(pb[i] > hi2) hi2 = pb[i];
  }
  nx = (int)floor((hi1 - lo1)/dx + 0.5) + 1; nz = (int)floor((hi2 - lo2)/dx + 0.5) + 1;
  fd->mirror3 = PETSC_FALSE;
  if(!medium->full_space){
    /* half space: vertical fault reaching the surface z = 0, mirror image in the surface.
       Order the rows by depth (row 0 at the surface) */
    double zmax = -1e300;
    PetscCheck(fabs(nn[2]) < 1e-9,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: 3-D half space needs a vertical fault");
    for(i = 0;i < n;i++) if(fault[i].x[INT_Z]*ls > zmax) zmax = fault[i].x[INT_Z]*ls;
    PetscCheck(fabs(zmax + 0.5*dx) < 1e-6*dx,PETSC_COMM_WORLD,PETSC_ERR_SUP,
	       "-fd_mode: 3-D half space: the shallowest row must reach the surface (centre at -dx/2, have %g)",zmax);
    if(e2[2] > 0.0){		/* dip axis points up: flip so that row index grows with depth */
      for(i = 0;i < n;i++) pb[i] = -pb[i];
      { double tmp = lo2; lo2 = -hi2; hi2 = -tmp; }
      for(k = 0;k < 3;k++) e2[k] = -e2[k];
    }
    fd->mirror3 = PETSC_TRUE;
  }
  PetscCheck(nx*nz == n,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: 3-D elements must form a complete %d x %d grid (have %d)",nx,nz,(int)n);
  PetscCall(PetscCalloc1(n,&cnt)); PetscCall(PetscMalloc1(n,&fd->g_of_s)); PetscCall(PetscMalloc1(n,&s_of_g));
  for(i = 0;i < n;i++){
    double fi = (pa[i] - lo1)/dx, fj = (pb[i] - lo2)/dx;
    int ii = (int)floor(fi + 0.5), jj = (int)floor(fj + 0.5);
    PetscCheck((fabs(fi - ii) < 1e-6) && (fabs(fj - jj) < 1e-6),PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: patch %d is off the regular 3-D grid",i);
    cnt[ii*nz + jj]++;
    fd->g_of_s[ii*nz + jj] = i; s_of_g[i] = ii*nz + jj;
  }
  for(i = 0;i < n;i++) PetscCheck(cnt[i] == 1,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: duplicate or missing patch in the 3-D grid");
  PetscCall(PetscFree(cnt)); PetscCall(PetscFree2(pa,pb));
  fd->dt_min = fd->cfl*dx/fd->beta;
  fd->nx3 = nx; fd->nz3 = nz;
  fd->tab3 = fd3d_tab_create(nx,fd->mirror3 ? 2*nz : nz,dx,fd->dt_min,fd->alpha,fd->beta,fd->rho,(int)fd->lmax);
  PetscCheck(fd->tab3,PETSC_COMM_WORLD,PETSC_ERR_MEM,"rsf_fd_setup: 3-D table allocation failed");
  /* static gate: Is(receiver, source) against the offset table, skipping the
     3 x 3 self neighbourhood whose regularisation may differ */
  for(gi = 0;gi < n;gi++){
    int si = s_of_g[gi]/nz, sj = s_of_g[gi]%nz;
    PetscCall(MatGetRow(medium->Is,gi,&ncols,&cols,&vals));
    for(c = 0;c < ncols;c++){
      int ti = s_of_g[cols[c]]/nz, tj = s_of_g[cols[c]]%nz, di = si - ti, dj = sj - tj;
      double ours, dv;
      if((abs(di) <= 1) && (abs(dj) <= 1)) continue;
      if(fd->mirror3){
	int nt = 2*nz, djr = (nz + sj) - (nz + tj), dji = (nz + sj) - (nz - 1 - tj);   /* real source and its image */
	ours = fd->tab3->ks_off[(size_t)(di + nx - 1)*(2*nt - 1) + (djr + nt - 1)] + fd->tab3->ks_off[(size_t)(di + nx - 1)*(2*nt - 1) + (dji + nt - 1)];
      }else
	ours = fd->tab3->ks_off[(size_t)(di + nx - 1)*(2*nz - 1) + (dj + nz - 1)];
      dv = PetscRealPart(vals[c]) - ours;
      dev += dv*dv; nrm += PetscRealPart(vals[c])*PetscRealPart(vals[c]);
    }
    PetscCall(MatRestoreRow(medium->Is,gi,&ncols,&cols,&vals));
  }
  PetscCall(PetscFree(s_of_g));
  dev = sqrt(dev/(nrm + 1e-300));
  fd->f3 = fd3d_create2(fd->tab3,nz,fd->mirror3 ? 1 : 0);
  PetscCheck(fd->f3,PETSC_COMM_WORLD,PETSC_ERR_MEM,"rsf_fd_setup: 3-D marcher allocation failed");
  fd->f3->law = 1; fd->f3->f0 = rsf->f0; fd->f3->V0 = rsf->v0;
  fd->f3->nl = (int)fd->nl;
  PetscCall(PetscMalloc4(n,&fd->f3->a,n,&fd->f3->b,n,&fd->f3->Lc,n,&fd->f3->theta));
  for(i = 0;i < n;i++){
    int g = fd->g_of_s[i];
    fd->f3->a[i] = fault[g].mu_sa; fd->f3->b[i] = fault[g].mu_db; fd->f3->Lc[i] = rsf->dc_vec ? rsf->dc_vec[g]*1.0 : rsf->dc;
  }
  fd->f3->static_fn = fd3d_static_is; fd->f3->static_ctx = (void *)par;
  PetscCall(MatCreateVecs(medium->Is,&fd->v3a,&fd->v3b));
  fd->nrec = 64; fd->irec = 0; fd->nrec_used = 0;
  PetscCall(PetscMalloc2(fd->nrec,&fd->rec_t,(size_t)fd->nrec*n,&fd->rec_v));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D, %d x %d squares of %g m, dt_min %.3e s, window %d lags, FFT %d x %d, tables %.1f MB\n",
			nx,nz,dx,(double)fd->dt_min,fd->tab3->L,fd->tab3->Mx,fd->tab3->Mz,
			(double)fd->tab3->L*fd->tab3->Mx*fd->tab3->Mh*sizeof(fd3d_cplx)/1048576.0));
  if(fd->mirror3)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D half space by mirror image: static operator of fault plus image deviates %.3e from the exact half-space Is "
			  "(expected, a few percent); the dynamic phases use the exact Is for the static part\n",dev));
  else
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D static gate: rel dev of Ks(tables) from Is %.3e (tol %g)\n",dev,(double)fd->gate_tol));
  PetscCheck(fd->mirror3 || (dev <= 5.0*fd->gate_tol),PETSC_COMM_WORLD,PETSC_ERR_USER,
	     "rsf_fd_setup: 3-D static limit deviates %.3e from the assembled Is: convention mismatch (a sign flip gives 2)",dev);
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
  /* sources: the real segments, plus for the antiplane half plane the image of
     each (vertical) segment in the free surface y = 0, same tangent, reversed
     element order (image element i is the mirror of real element n-1-i) */
  fd->nsrc = fd->segs->nseg*(fd->halfplane ? 2 : 1);
  PetscCall(PetscMalloc3(fd->nsrc,&fd->srcseg,fd->nsrc,&fd->src_of,fd->nsrc,&fd->src_rev));
  for(k = 0;k < fd->segs->nseg;k++){
    fd->srcseg[k] = fd->segs->seg[k]; fd->src_of[k] = (int)k; fd->src_rev[k] = 0;
    if(fd->halfplane){
      struct fd_segment *im = fd->srcseg + fd->segs->nseg + k, *re = fd->segs->seg + k;
      PetscCheck(fabs(re->s[0]) < 1e-9,PETSC_COMM_WORLD,PETSC_ERR_SUP,
		 "rsf_fd_setup: antiplane half plane needs vertical segments (image of a dipping segment is not parallel, M7)");
      *im = *re;
      im->X0[0] = re->X0[0];
      im->X0[1] = -(re->X0[1] + (re->n - 1)*re->dx*re->s[1]);
      fd->src_of[fd->segs->nseg + k] = (int)k; fd->src_rev[fd->segs->nseg + k] = 1;
    }
  }
  /* one padded FFT length for all pairs: smallest power of two >= 2 nmax - 1 */
  nm = 0;
  for(k = 0;k < fd->segs->nseg;k++)
    if(fd->segs->seg[k].n > nm) nm = fd->segs->seg[k].n;
  fd->M = fd_fft_size_for(2*nm - 1);
  PetscCall(PetscCalloc1((size_t)fd->nsrc*fd->segs->nseg,&fd->pairs));
  fd->nbytes = 0; fd->L = 0;
  for(k = 0;k < fd->nsrc;k++)			/* source */
    for(l = 0;l < fd->segs->nseg;l++){		/* receiver */
      struct fd_pair *p = fd_pair_build(fd->srcseg + k,fd->segs->seg + l,fd->mode2d,fd->dt_min,fd->et,
					fd->alpha,fd->beta,fd->rho,fd->tol,(int)fd->lmax,fd->lmax > 0,fd->M,err);
      PetscCheck(p,PETSC_COMM_WORLD,PETSC_ERR_SUP,"rsf_fd_setup: table for pair (%ld,%ld) failed: %s",(long)k,(long)l,err);
      fd->pairs[k*fd->segs->nseg + l] = p;
      fd->nbytes += (double)p->L*p->M*16.0*(p->has_normal ? 2 : 1);
      if(p->L > fd->L) fd->L = p->L;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: pair (%s%ld -> %ld): %d x %d elements, window %d lags (%.2f s)\n",
			    k >= fd->segs->nseg ? "image " : "",(long)fd->src_of[k],(long)l,p->src->n,p->rcv->n,p->L,(double)(p->L*fd->dt_min)));
    }
  fd->sys = fd_sys_create(fd->segs,fd->pairs,(int)fd->nsrc,fd->src_of,fd->src_rev,fd->mode2d,fd->dt_min,fd->alpha,fd->beta,fd->rho,+1);
  PetscCheck(fd->sys,PETSC_COMM_WORLD,PETSC_ERR_MEM,"rsf_fd_setup: marcher system allocation failed");
  if(fd->coarse && (fd->conv != 1) && (fd->kmax >= 2)){
    PetscCheck(fd_sys_set_coarse(fd->sys,(int)fd->kmax) == 0,PETSC_COMM_WORLD,PETSC_ERR_LIB,"rsf_fd_setup: coarse history setup failed");
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: coarse history for aligned steps of %d substeps\n",(int)fd->kmax));
  }
  if(fd->conv){
    PetscCheck((fd->conv == 1) || (fd->conv == 2),PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"-fd_conv must be 0, 1 or 2");
    PetscCheck(fd_sys_set_conv(fd->sys,(int)fd->conv,(int)fd->conv_block,(int)fd->kmax) == 0,PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,
	       "rsf_fd_setup: -fd_conv_block must be a power of two (or 0 for automatic)");
    if(fd->sys->conv)
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: %s convolution: blocks of %d lags, lags 1..%d direct, FFT length %d\n",
			    fd->sys->conv == 2 ? "hybrid partitioned" : "partitioned",fd->sys->B,fd->sys->D,fd->sys->F));
    else
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: window (%d lags) too short for the partitioned convolution, direct sum used\n",fd->sys->L));
  }
  fd->nbytes = fd_sys_bytes(fd->sys);
  /* element map: segment order <-> global */
  PetscCall(PetscMalloc1(medium->nrflt,&fd->g_of_s));
  for(k = 0;k < fd->segs->nseg;k++)
    for(l = 0;l < fd->segs->seg[k].n;l++)
      fd->g_of_s[fd->sys->off[k] + l] = fd->segs->gidx[k][l];
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
   compare the static limit of the tables with the assembled operators.
   rsf_solve assembles Is in SI with the stress-from-slip sign (positive
   slip reduces the shear traction on the slipping patch, negative
   diagonal), so Is must equal +Ks of the tables; In is assembled
   compression positive, so In must equal -Kn (in-plane only). Near-
   diagonal entries of the same segment are skipped, as in the VE gate.
   A pure sign error gives a relative deviation of 2.
*/
static PetscErrorCode fd_gate_one(struct fd_vars *fd, Mat A, const double *K, PetscInt nrflt, double *dev2, double *nrm2)
{
  PetscInt gi, c, ncols, si, sj;
  const PetscInt *cols;
  const PetscScalar *vals;
  double *row, dev = 0.0, nrm = 0.0, dv;
  int N = fd->sys->N, *s_of_g, k, i;
  PetscFunctionBeginUser;
  PetscCall(PetscCalloc1(nrflt,&row));
  PetscCall(PetscMalloc1(nrflt,&s_of_g));
  for(k = 0;k < fd->segs->nseg;k++)
    for(i = 0;i < fd->segs->seg[k].n;i++)
      s_of_g[fd->segs->gidx[k][i]] = fd->sys->off[k] + i;
  for(gi = 0;gi < nrflt;gi++){
    for(c = 0;c < nrflt;c++) row[c] = 0.0;
    PetscCall(MatGetRow(A,gi,&ncols,&cols,&vals));
    for(c = 0;c < ncols;c++) row[cols[c]] = PetscRealPart(vals[c]);
    PetscCall(MatRestoreRow(A,gi,&ncols,&cols,&vals));
    si = s_of_g[gi];
    for(c = 0;c < nrflt;c++){
      sj = s_of_g[c];
      /* skip near-diagonal self terms within one segment */
      if(labs((long)si - (long)sj) <= 1){
	for(k = 0;k < fd->segs->nseg;k++)
	  if((si >= fd->sys->off[k]) && (si < fd->sys->off[k+1]) && (sj >= fd->sys->off[k]) && (sj < fd->sys->off[k+1]))
	    break;
	if(k < fd->segs->nseg) continue;
      }
      dv = row[c] - K[(size_t)si*N + sj];
      dev += dv*dv; nrm += row[c]*row[c];
    }
  }
  PetscCall(PetscFree(row)); PetscCall(PetscFree(s_of_g));
  *dev2 = dev; *nrm2 = nrm;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode fd_gate_static(struct interact_ctx *par, struct fd_vars *fd)
{
  struct med *medium = par->medium;
  struct rsf_vars *rsf = medium->rsf;
  double d2, n2, dev, devn = -1.0, dn2, nn2;
  PetscFunctionBeginUser;
  PetscCall(fd_gate_one(fd,medium->Is,fd->sys->Ks,medium->nrflt,&d2,&n2));
  dev = sqrt(d2/(n2 + 1e-300));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: static gate: rel dev of Ks(tables) from Is %.3e (tol %g)\n",dev,(double)fd->gate_tol));
  if((fd->mode2d == FD_MODE_INPLANE) && rsf->calc_sigma_dot){
    /* In is assembled compression positive for interact's positive (left-lateral) slip,
       which equals +Kn of the tables (tension positive for the prototype's right-lateral
       positive slip); a numerically zero In (single planar segment) is judged against
       the shear norm */
    PetscCall(fd_gate_one(fd,medium->In,fd->sys->Kn,medium->nrflt,&dn2,&nn2));
    devn = sqrt(dn2/fmax(nn2,1e-16*n2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: static gate: rel dev of Kn(tables) from In %.3e (normal norm / shear norm %.2e)\n",devn,sqrt(nn2/(n2 + 1e-300))));
  }
  PetscCheck(dev <= 5.0*fd->gate_tol,PETSC_COMM_WORLD,PETSC_ERR_USER,
	     "rsf_fd_setup: static limit of the dynamic tables deviates %.3e from the assembled Is: convention mismatch (a sign flip gives 2)",dev);
  PetscCheck(devn <= 5.0*fd->gate_tol,PETSC_COMM_WORLD,PETSC_ERR_USER,
	     "rsf_fd_setup: normal static limit deviates %.3e from the assembled In (a sign flip gives 2)",devn);
  if((dev > fd->gate_tol) || (devn > fd->gate_tol))
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: WARNING: static gate deviation exceeds %g\n",(double)fd->gate_tol));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* interpolate the accepted-step record at time t (segment order output) */
static void fd_record_interp(struct fd_vars *fd, PetscReal t, double *V)
{
  int n = fd->nrec_used, i, i0, i1, N = fd->sys->N, j;
  double w;
  PetscReal *tt = fd->rec_t;
  if(n == 0){ for(j = 0;j < N;j++) V[j] = 0.0; return; }
  /* records are stored cyclically; find the two bracketing entries (times increase) */
  i0 = -1; i1 = -1;
  for(i = 0;i < n;i++){
    int idx = (fd->irec - n + i + fd->nrec) % fd->nrec;
    if(tt[idx] <= t) i0 = idx;
    if((tt[idx] >= t) && (i1 < 0)) i1 = idx;
  }
  if(i0 < 0) i0 = (fd->irec - n + fd->nrec) % fd->nrec;
  if(i1 < 0) i1 = (fd->irec - 1 + fd->nrec) % fd->nrec;
  w = (tt[i1] > tt[i0]) ? (t - tt[i0])/(tt[i1] - tt[i0]) : 0.0;
  if(w < 0.0) w = 0.0;
  if(w > 1.0) w = 1.0;
  for(j = 0;j < N;j++)
    V[j] = (1.0 - w)*fd->rec_v[(size_t)i0*N + j] + w*fd->rec_v[(size_t)i1*N + j];
}

/* one dynamic phase: from X at time t to X at the time the slip rate has
   stayed below vback for one window; returns the new time */
static PetscErrorCode fd_run_phase(TS ts, struct interact_ctx *par, struct fd_vars *fd, PetscReal t, Vec X, PetscReal *t_end)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  struct rsf_vars *rsf = medium->rsf;
  struct fd_sys *sy = fd->sys;
  struct fd_state *st = fd_state_create(sy->N);
  struct fd_rsf fr;
  PetscScalar *x;
  int N = sy->N, i, g, l, k, nsteps = 0, kk;
  double vmax = 0.0, below = 0.0, dcmin = 1e30, *vtmp, *vref, vpeak = 0.0, slipsum = 0.0, t_wall0, t_wall1;
  PetscFunctionBeginUser;
  PetscCall(PetscTime(&t_wall0));
  fr.N = N; fr.v0 = rsf->v0; fr.f0 = rsf->f0; fr.tol = 1e-10; fr.max_iter = 50;
  PetscCall(PetscMalloc3(N,&fr.a,N,&fr.b,N,&fr.dc));
  PetscCall(PetscMalloc1(N,&vtmp));
  PetscCall(PetscMalloc1(N,&vref));
  st->normal_coupling = (fd->mode2d == FD_MODE_INPLANE) && rsf->calc_sigma_dot;
  st->nl = (int)fd->nl;
  PetscCall(VecGetArray(X,&x));
  for(i = 0;i < N;i++){
    g = fd->g_of_s[i];
    fr.a[i] = fault[g].mu_sa; fr.b[i] = fault[g].mu_db; fr.dc[i] = rsf->dc_vec ? rsf->dc_vec[g] : rsf->dc;
    if(fr.dc[i] < dcmin) dcmin = fr.dc[i];
    st->psi[i] = x[4*g]; st->sigma0[i] = x[4*g+2];
    st->V[i] = 2.0*rsf->v0*exp(-x[4*g]/fault[g].mu_sa)*sinh(x[4*g+1]/(fault[g].mu_sa*x[4*g+2]));
    st->tau0[i] = x[4*g+1] + sy->eta*st->V[i];	/* elastic stress at the switch */
    st->taudot[i] = fault[g].sinc[0];
    st->Hn_dot[i] = rsf->calc_sigma_dot ? fault[g].sinc[1] : 0.0;
    st->sigma[i] = st->sigma0[i]; st->tau[i] = x[4*g+1];
    if(st->V[i] > vmax) vmax = st->V[i];
  }
  /* seed the fine-grid history from the accepted-step record */
  fd_sys_reset_history(sy);
  for(l = sy->L;l >= 1;l--){
    fd_record_interp(fd,t - l*sy->dt,vtmp);
    fd_sys_push(sy,vtmp);
  }
  st->t = 0.0;			/* marcher time is relative to the switch (loading = taudot * (t - t_switch)) */
  if(getenv("FD_DEBUG")){
    /* consistency at the switch: elastic stress tau0 against strength plus damping */
    double rmax = 0.0, vm = 0.0; int im = 0;
    for(i = 0;i < N;i++){
      double r = st->tau0[i] - sy->eta*st->V[i] - fd_rsf_tau(&fr,i,st->V[i],st->psi[i],st->sigma0[i]);
      if(fabs(r) > rmax) rmax = fabs(r);
      if(st->V[i] > vm){ vm = st->V[i]; im = i; }
    }
    fprintf(stderr,"FD_DEBUG: switch at t %.6f yr: vmax %.3e at s-index %d (z %.0f m), max |tau0 - eta V - tau_f| %.3e Pa, eta %.3e, L %d, nsrc %d, N %d\n",
	    t/SEC_PER_YEAR,vm,im,-(double)(fd->segs->seg[0].X0[1] + im*fd->segs->seg[0].dx*fd->segs->seg[0].s[1]),rmax,sy->eta,sy->L,sy->nsrc,N);
  }
  {
    int a;
    for(a = 0;a < fd->nstat;a++){
      char fn[128]; sprintf(fn,"fd_phase_%ld_st%ld.dat",(long)fd->nphase + 1,(long)fd->stat[a]);
      fd->fstat[a] = fopen(fn,"w");
      fd->sstat[a] = -1;
      for(i = 0;i < N;i++) if(fd->g_of_s[i] == fd->stat[a]) fd->sstat[a] = i;
    }
  }
  for(;;){
    int a;
    double kd = fd->xi*dcmin/(vmax*sy->dt);
    k = 1;
    /* powers of two up to kmax and xi dc / (V dt); with the coarse history the
       step must also start at a multiple of its length (aligned blocks) */
    while((2*k <= fd->kmax) && (2*k <= kd) && (!sy->Kc || (sy->hk[0] % (2*k) == 0))) k *= 2;
    if(fd_step_rsf(sy,st,&fr,k) != 0)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_LIB,"rsf_fd: marcher step failed");
    nsteps++;
    vmax = 0.0;
    for(i = 0;i < N;i++) if(st->V[i] > vmax) vmax = st->V[i];
    if(getenv("FD_DEBUG") && (nsteps <= 3)){
      int im = 0; double hs = 0.0;
      for(i = 0;i < N;i++){ if(st->V[i] == vmax) im = i; if(fabs(st->hs[i]) > hs) hs = fabs(st->hs[i]); }
      fprintf(stderr,"FD_DEBUG: step %d k %d: vmax %.3e at s-index %d, V[0] %.2e V[N/2] %.2e V[N-1] %.2e, max|hist| %.3e, tau at im %.4e, H at im %.4e\n",
	      nsteps,k,vmax,im,st->V[0],st->V[N/2],st->V[N-1],hs,st->tau[im],st->H[im]);
    }
    for(a = 0;a < fd->nstat;a++)
      if(fd->sstat[a] >= 0)
	fprintf(fd->fstat[a],"%.10e %.10e %.10e %.10e %.4e %d\n",st->t,st->delta[fd->sstat[a]],st->V[fd->sstat[a]],st->tau[fd->sstat[a]],vmax,k);
    if(vmax > vpeak) vpeak = vmax;
    if(fd->back_eps > 0.0){
      /* quasi-steady hand-back: the hereditary term of a slip rate that is
	 constant over one window equals the QD radiation damping, so the
	 phase can end once every element's slip rate has changed by less
	 than back_eps (relative, floored at vback) over a full window and the
	 fault-wide maximum is below back_vmax */
      int quiet = (vmax < fd->back_vmax);
      if(quiet && (below > 0.0))
	for(i = 0;i < N;i++)
	  if(fabs(st->V[i] - vref[i]) > fd->back_eps*(fabs(vref[i]) + fd->vback)){ quiet = 0; break; }
      if(!quiet || (below == 0.0)){
	memcpy(vref,st->V,sizeof(double)*N);
	below = quiet ? 1e-30 : 0.0;
      }else
	below += k*sy->dt;
    }else{
      if(vmax < fd->vback) below += k*sy->dt; else below = 0.0;
    }
    if(below >= sy->L*sy->dt) break;
    if(nsteps >= fd->max_steps)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,"rsf_fd: dynamic phase exceeded -fd_max_steps = %ld steps (t - t_switch = %g s, max V %g)",(long)fd->max_steps,(double)st->t,vmax);
  }
  {
    int a; char fn[128]; FILE *fp;
    for(a = 0;a < fd->nstat;a++) if(fd->fstat[a]) fclose(fd->fstat[a]);
    sprintf(fn,"fd_phase_%ld_slip.dat",(long)fd->nphase + 1); fp = fopen(fn,"w");
    for(i = 0;i < N;i++) fprintf(fp,"%d %ld %.10e %.6e %.6e\n",i,(long)fd->g_of_s[i],st->delta[i],st->V[i],st->psi[i]);
    fclose(fp);
  }
  /* hand back: psi, tau (= strength at the final slip rate = elastic stress at the end of the last step), sigma, slip */
  for(i = 0;i < N;i++){
    g = fd->g_of_s[i];
    x[4*g] = st->psi[i]; x[4*g+1] = st->tau[i]; x[4*g+2] = st->sigma[i]; x[4*g+3] += st->delta[i];
    slipsum += st->delta[i];
  }
  PetscCall(VecRestoreArray(X,&x));
  *t_end = t + st->t;
  PetscCall(PetscTime(&t_wall1));
  if(getenv("FD_PROFILE")){
    extern double fd_prof[6];
    fprintf(stderr,"FD_PROFILE phase: history %.1f s, static %.1f s, solve %.1f s, push %.1f s, other %.1f s\n",fd_prof[0],fd_prof[1],fd_prof[2],fd_prof[3],fd_prof[5]);
    for(i = 0;i < 6;i++) fd_prof[i] = 0.0;
  }
  for(kk = 0,i = 0;i < sy->nops;i++) if(sy->ops[i].bw > kk) kk = sy->ops[i].bw;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
		"rsf_fd: phase %ld: t %.8f yr, %d steps over %.1f s, peak slip rate %.3f m/s, mean slip %.4f m, max band %d, wall %.1f s\n",
		(long)fd->nphase + 1,(double)(t/SEC_PER_YEAR),nsteps,(double)st->t,vpeak,slipsum/N,kk,(double)(t_wall1 - t_wall0)));
  PetscCall(PetscFree3(fr.a,fr.b,fr.dc)); PetscCall(PetscFree(vtmp)); PetscCall(PetscFree(vref));
  fd_state_free(st);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  accepted-step monitor: records (t, slip rate) for the history seeding.
  rsf->vel holds the slip rates of the last RHS evaluation, close enough
  to the accepted state for this purpose. The switch itself is a TS event
  (slot uc->fd_event_index in rsf_event_function) whose post-event
  callback rsf_fd_postevent runs the phase: PETSc locks the solution during
  monitors and its event machinery must be reinitialised after a
  discontinuous change, both of which the event route handles.
*/
PetscErrorCode rsf_fd_monitor(TS ts, PetscInt step, PetscReal time, Vec X, void *ctx)
{
  struct interact_ctx *par = (struct interact_ctx *)ctx;
  struct rsf_vars *rsf = par->medium->rsf;
  struct fd_vars *fd = rsf->fd;
  const PetscScalar *v;
  PetscInt i;
  PetscFunctionBeginUser;
  if(!fd)
    PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(VecGetArrayRead(rsf->vel,&v));
  fd->rec_t[fd->irec] = time;
  {
    int N = fd->dim3 ? fd->f3->N : fd->sys->N;
    for(i = 0;i < N;i++) fd->rec_v[(size_t)fd->irec*N + i] = PetscRealPart(v[fd->g_of_s[i]]);
  }
  PetscCall(VecRestoreArrayRead(rsf->vel,&v));
  fd->irec = (fd->irec + 1) % fd->nrec; if(fd->nrec_used < fd->nrec) fd->nrec_used++;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* post-event callback for the switch event: dynamic phase from (t, X),
   solution replaced in place, TS clock advanced to the end of the phase */
PetscErrorCode rsf_fd_postevent(TS ts, struct interact_ctx *par, PetscReal t, Vec X)
{
  struct rsf_vars *rsf = par->medium->rsf;
  struct fd_vars *fd = rsf->fd;
  PetscReal t_end = t;
  PetscFunctionBeginUser;
  if(!fd)
    PetscFunctionReturn(PETSC_SUCCESS);
  fd->nflag++;
  if(fd->mode == 1){
    if(fd->dim3) PetscCall(fd3d_run_phase(ts,par,fd,t,X,&t_end));
    else PetscCall(fd_run_phase(ts,par,fd,t,X,&t_end));
    fd->nphase++;
    fd->nrec_used = 0; fd->irec = 0;	/* the record is stale after the phase */
    PetscCall(TSSetTime(ts,t_end));
    PetscCall(TSSetTimeStep(ts,fd->dt_restart));
  }else
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_postevent: t %.6f yr: switch threshold crossed (mode %d: no phase)\n",
			  (double)(t/SEC_PER_YEAR),(int)fd->mode));
  PetscFunctionReturn(PETSC_SUCCESS);
}


/* before the integration starts: if the initial state already slips faster than
   the switch threshold (e.g. BP5's seeded nucleation patch), run a dynamic
   phase right away, since no upward crossing will be detected */
PetscErrorCode rsf_fd_initial(TS ts, struct interact_ctx *par, Vec X)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  struct rsf_vars *rsf = medium->rsf;
  struct fd_vars *fd = rsf->fd;
  const PetscScalar *x;
  PetscReal t, t_end, vmax = 0.0, v;
  PetscInt i;
  PetscFunctionBeginUser;
  if(!fd || (fd->mode != 1))
    PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(VecGetArrayRead(X,&x));
  for(i = 0;i < medium->nrflt;i++){
    v = 2.0*rsf->v0*exp(-x[4*i]/fault[i].mu_sa)*sinh(x[4*i+1]/(fault[i].mu_sa*x[4*i+2]));
    if(fabs(v) > vmax) vmax = fabs(v);
  }
  PetscCall(VecRestoreArrayRead(X,&x));
  if(vmax <= fd->vsw)
    PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(TSGetTime(ts,&t));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_initial: initial max slip rate %.3e m/s above vsw: dynamic phase from t = %g s\n",(double)vmax,(double)t));
  if(fd->dim3) PetscCall(fd3d_run_phase(ts,par,fd,t,X,&t_end));
  else PetscCall(fd_run_phase(ts,par,fd,t,X,&t_end));
  fd->nphase++; fd->nflag++;
  fd->nrec_used = 0; fd->irec = 0;
  PetscCall(TSSetTime(ts,t_end));
  PetscCall(TSSetTimeStep(ts,fd->dt_restart));
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
  if(fd->f3){
    PetscCall(PetscFree4(fd->f3->a,fd->f3->b,fd->f3->Lc,fd->f3->theta));
    fd3d_free(fd->f3);
    PetscCall(VecDestroy(&fd->v3a)); PetscCall(VecDestroy(&fd->v3b));
  }
  if(fd->tab3) fd3d_tab_free(fd->tab3);
  if(fd->sys) fd_sys_free(fd->sys);
  if(fd->pairs){
    for(k = 0;k < fd->nsrc*fd->segs->nseg;k++)
      fd_pair_free(fd->pairs[k]);
    PetscCall(PetscFree(fd->pairs));
  }
  if(fd->srcseg) PetscCall(PetscFree3(fd->srcseg,fd->src_of,fd->src_rev));
  if(fd->g_of_s) PetscCall(PetscFree(fd->g_of_s));
  if(fd->rec_t) PetscCall(PetscFree2(fd->rec_t,fd->rec_v));
  if(fd->segs) fd_seglist_free(fd->segs);
  PetscCall(PetscFree(fd));
  rsf->fd = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}
