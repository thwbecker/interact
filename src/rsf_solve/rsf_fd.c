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
extern double fd3d_layer_prof[4];

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
  /* MPI: the 3-D marcher splits its tables, histories and local solves over the ranks
     (fd_plan.md S2); the 2-D code runs on one rank */
  PetscCheck((size == 1) || (n3d > 0),PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: 2-D faults on one MPI rank only in this version");
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
  fd->back_dyn = 0.0; fd->back_nwin = 10; fd->reentry = 1.0;
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_back_dyn",&fd->back_dyn,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_back_nwin",&fd->back_nwin,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_reentry",&fd->reentry,NULL));
  if(fd->back_vmax <= 0.0) fd->back_vmax = 0.5*fd->vsw;
  fd->vsw_eff = fd->vsw;
  fd->nl = 1;			/* second-order time marching (Noda et al. 2020) by default */
  fd->layer_w = 0; fd->layer_m = 2; fd->layer_margin = -1; fd->surface_opt = 0; fd->surf_ext = -1;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_surface",&fd->surface_opt,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_surface_ext",&fd->surf_ext,NULL));
  fd->surf_m = 4; fd->surf_grade = 1; fd->surf_extend = -1;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_surface_extend",&fd->surf_extend,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_surface_grade",&fd->surf_grade,NULL));
  {
    PetscBool set = PETSC_FALSE;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_surface_m",&fd->surf_m,&set));
    fd->surf_m_default = set ? PETSC_FALSE : PETSC_TRUE;
  }
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_layer",&fd->layer_w,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_layer_m",&fd->layer_m,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_layer_margin",&fd->layer_margin,NULL));
  fd->kmax3 = 0;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_kmax3",&fd->kmax3,NULL));
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
    {
      /* 3-D defaults for speed (fd_plan.md L2, L3), unless given: quasi-steady hand-back,
	 hand-back by traction difference (10 kPa), re-entry only above 10 vsw until the
	 peak rate has dropped below the hand-back limit */
      PetscBool set;
      PetscCall(PetscOptionsHasName(NULL,NULL,"-fd_back_eps",&set));
      if(!set) fd->back_eps = 0.05;
      PetscCall(PetscOptionsHasName(NULL,NULL,"-fd_back_dyn",&set));
      if(!set) fd->back_dyn = 1e4;
      PetscCall(PetscOptionsHasName(NULL,NULL,"-fd_reentry",&set));
      if(!set) fd->reentry = 10.0;
    }
    fd->mu = shear_modulus_si; fd->beta = s_wave_speed_si;
    fd->nu = medium->elastic.poisson;
    fd->rho = fd->mu/(fd->beta*fd->beta);
    fd->alpha = fd->beta*sqrt(2.0*(1.0 - fd->nu)/(1.0 - 2.0*fd->nu));
    rsf->fd = fd;
    PetscCall(fd3d_setup_from_geometry(par,fd));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D fully dynamic extension, mode %d, %s, vsw %.1e vback %.1e m/s, fixed fine steps, %s time marching\n",
			  (int)fd->mode,fd->mirror3 ? "half space (mirror-image free surface, exact static Is)" : (fd->surf3 ? "half space (virtual free surface, exact static Is)" : "whole space"),(double)fd->vsw,(double)fd->vback,fd->nl ? "second-order (-fd_nl 1)" : "first-order"));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
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
static PetscErrorCode fd_full_scatter(Vec X, VecScatter *sc, Vec *Xall, PetscScalar **x);
/* full copy of a distributed vector on every rank (MPI) or the array itself (one rank) */
static PetscErrorCode fd_full_get(Vec X, VecScatter *sc, Vec *Xall, PetscScalar **x)
{
  PetscMPIInt size;
  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)X),&size));
  if(size == 1){ PetscCall(VecGetArray(X,x)); PetscFunctionReturn(PETSC_SUCCESS); }
  return fd_full_scatter(X,sc,Xall,x);
}
/* read-only variant: on one rank the vector's state is not touched (VecGetArrayRead), so the
   QD integrator sees exactly what it saw before */
static PetscErrorCode fd_full_get_read(Vec X, VecScatter *sc, Vec *Xall, PetscScalar **x)
{
  PetscMPIInt size;
  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)X),&size));
  if(size == 1){ PetscCall(VecGetArrayRead(X,(const PetscScalar **)x)); PetscFunctionReturn(PETSC_SUCCESS); }
  return fd_full_scatter(X,sc,Xall,x);
}
static PetscErrorCode fd_full_put_read(Vec X, VecScatter *sc, Vec *Xall, PetscScalar **x)
{
  PetscMPIInt size;
  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)X),&size));
  if(size == 1){ PetscCall(VecRestoreArrayRead(X,(const PetscScalar **)x)); PetscFunctionReturn(PETSC_SUCCESS); }
  (void)sc;
  PetscCall(VecRestoreArray(*Xall,x));
  PetscFunctionReturn(PETSC_SUCCESS);
}
static PetscErrorCode fd_full_scatter(Vec X, VecScatter *sc, Vec *Xall, PetscScalar **x)
{
  PetscFunctionBeginUser;
  if(!*sc) PetscCall(VecScatterCreateToAll(X,sc,Xall));
  PetscCall(VecScatterBegin(*sc,X,*Xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecScatterEnd(*sc,X,*Xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecGetArray(*Xall,x));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* release; write = 1 copies the (replicated) full array back into the owned parts of X */
static PetscErrorCode fd_full_put(Vec X, VecScatter *sc, Vec *Xall, PetscScalar **x, int write)
{
  PetscMPIInt size;
  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)X),&size));
  if(size == 1){ PetscCall(VecRestoreArray(X,x)); PetscFunctionReturn(PETSC_SUCCESS); }
  PetscCall(VecRestoreArray(*Xall,x));
  if(write){
    PetscCall(VecScatterBegin(*sc,*Xall,X,INSERT_VALUES,SCATTER_REVERSE));
    PetscCall(VecScatterEnd(*sc,*Xall,X,INSERT_VALUES,SCATTER_REVERSE));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static void fd3d_static_is(void *ctx, const double *delta, double *out)
{
  struct interact_ctx *par = (struct interact_ctx *)ctx;
  struct fd_vars *fd = par->medium->rsf->fd;
  PetscScalar *a, *b;
  PetscInt lo, hi;
  int i, n = fd->f3->N;
  PetscCallAbort(PETSC_COMM_WORLD,VecGetOwnershipRange(fd->v3a,&lo,&hi));
  PetscCallAbort(PETSC_COMM_WORLD,VecGetArray(fd->v3a,&a));
  for(i = 0;i < n;i++) if((fd->g_of_s[i] >= lo) && (fd->g_of_s[i] < hi)) a[fd->g_of_s[i] - lo] = delta[i];
  PetscCallAbort(PETSC_COMM_WORLD,VecRestoreArray(fd->v3a,&a));
  PetscCallAbort(PETSC_COMM_WORLD,MatMult(par->medium->Is,fd->v3a,fd->v3b));
  PetscCallAbort(PETSC_COMM_WORLD,fd_full_get_read(fd->v3b,&fd->sc3,&fd->v3all,&b));
  for(i = 0;i < n;i++) out[i] = PetscRealPart(b[fd->g_of_s[i]]);
  PetscCallAbort(PETSC_COMM_WORLD,fd_full_put_read(fd->v3b,&fd->sc3,&fd->v3all,&b));
}

/* exact static normal traction (tension positive) of the phase slip: -In delta (In: compression positive) */
static void fd3d_static_in(void *ctx, const double *delta, double *out)
{
  struct interact_ctx *par = (struct interact_ctx *)ctx;
  struct fd_vars *fd = par->medium->rsf->fd;
  PetscScalar *a, *b;
  PetscInt lo, hi;
  int i, n = fd->f3->N;
  PetscCallAbort(PETSC_COMM_WORLD,VecGetOwnershipRange(fd->v3a,&lo,&hi));
  PetscCallAbort(PETSC_COMM_WORLD,VecGetArray(fd->v3a,&a));
  for(i = 0;i < n;i++) if((fd->g_of_s[i] >= lo) && (fd->g_of_s[i] < hi)) a[fd->g_of_s[i] - lo] = delta[i];
  PetscCallAbort(PETSC_COMM_WORLD,VecRestoreArray(fd->v3a,&a));
  PetscCallAbort(PETSC_COMM_WORLD,MatMult(par->medium->In,fd->v3a,fd->v3b));
  PetscCallAbort(PETSC_COMM_WORLD,fd_full_get_read(fd->v3b,&fd->sc3,&fd->v3all,&b));
  for(i = 0;i < n;i++) out[i] = -PetscRealPart(b[fd->g_of_s[i]]);
  PetscCallAbort(PETSC_COMM_WORLD,fd_full_put_read(fd->v3b,&fd->sc3,&fd->v3all,&b));
}

/* 3-D dynamic phase (fd3d.c), counterpart of fd_run_phase: fixed fine steps */

/* 3-D phase for two parallel planes: as fd3d_run_phase, both marchers stepped together by
   fd3d_pair_step (fixed fine steps), combined index p*N + i, locked cells skipped */
static PetscErrorCode fd3d_run_phase_pair(TS ts, struct interact_ctx *par, struct fd_vars *fd, PetscReal t, Vec X, PetscReal *t_end)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  struct rsf_vars *rsf = medium->rsf;
  struct fd3d *fm[2] = {fd->f3, fd->f3b};
  struct fd3d_pair *pp = fd->pair3;
  PetscScalar *x;
  int N = fd->f3->N, NT = 2*N, i, p, g, l, a, nsteps = 0, L = fd->tab3->L, ndyn_quiet = 0;
  double dt = fd->tab3->dt, eta = 0.5*fd->mu/fd->beta, vmax, vpeak = 0.0, below = 0.0, slipsum = 0.0, *vtmp, *load[2], *vref, t_wall0, t_wall1;
  PetscFunctionBeginUser;
  (void)ts;
  PetscCall(PetscTime(&t_wall0));
  PetscCall(PetscMalloc4(NT,&vtmp,N,&load[0],N,&load[1],NT,&vref));
  for(p = 0;p < 2;p++) fd3d_reset(fm[p]);
  PetscCall(fd_full_get(X,&fd->scX,&fd->Xall,&x));
  vmax = 0.0;
  for(p = 0;p < 2;p++) for(i = 0;i < N;i++){
    struct fd3d *f = fm[p];
    g = fd->g_of_s[p*N + i];
    if(g < 0){ f->sigma[i] = 0.0; f->V[i] = 0.0; f->tau0[i] = 0.0; pp->sig0[p][i] = 0.0; continue; }
    f->sigma[i] = x[4*g+2]; pp->sig0[p][i] = x[4*g+2];
    f->V[i] = 2.0*rsf->v0*exp(-x[4*g]/fault[g].mu_sa)*sinh(x[4*g+1]/(fault[g].mu_sa*x[4*g+2]));
    f->theta[i] = f->Lc[i]/f->V0*exp((x[4*g] - f->f0)/f->b[i]);
    f->tau0[i] = x[4*g+1] + eta*f->V[i];
    if(fabs(f->V[i]) > vmax) vmax = fabs(f->V[i]);
  }
  /* history from the accepted-step record (combined order) */
  for(l = L;l >= 1;l--){
    double tl = t - l*dt, w; int k, i0 = -1, i1 = -1;
    for(k = 0;k < fd->nrec_used;k++){
      int idx = (fd->irec - fd->nrec_used + k + fd->nrec) % fd->nrec;
      if(fd->rec_t[idx] <= tl) i0 = idx;
      if((fd->rec_t[idx] >= tl) && (i1 < 0)) i1 = idx;
    }
    if(fd->nrec_used == 0){ for(i = 0;i < NT;i++) vtmp[i] = 0.0; }
    else{
      if(i0 < 0) i0 = (fd->irec - fd->nrec_used + fd->nrec) % fd->nrec;
      if(i1 < 0) i1 = (fd->irec - 1 + fd->nrec) % fd->nrec;
      w = (fd->rec_t[i1] > fd->rec_t[i0]) ? (tl - fd->rec_t[i0])/(fd->rec_t[i1] - fd->rec_t[i0]) : 0.0;
      if(w < 0.0) w = 0.0;
      if(w > 1.0) w = 1.0;
      for(i = 0;i < NT;i++) vtmp[i] = (1.0 - w)*fd->rec_v[(size_t)i0*NT + i] + w*fd->rec_v[(size_t)i1*NT + i];
    }
    for(p = 0;p < 2;p++) fd3d_push(fm[p],vtmp + p*N);
  }
  for(a = 0;a < fd->nstat;a++){
    char fn[128]; PetscMPIInt rk;
    sprintf(fn,"fd_phase_%ld_st%ld.dat",(long)fd->nphase + 1,(long)fd->stat[a]);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rk); fd->fstat[a] = (rk == 0) ? fopen(fn,"w") : NULL;
    fd->sstat[a] = -1;
    for(i = 0;i < NT;i++) if(fd->g_of_s[i] == fd->stat[a]) fd->sstat[a] = i;
  }
  for(;;){
    double trel = fm[0]->time + dt;
    for(p = 0;p < 2;p++) for(i = 0;i < N;i++){ g = fd->g_of_s[p*N + i]; load[p][i] = (g >= 0) ? fault[g].sinc[0]*trel : 0.0; }
    fd3d_pair_step(pp,load[0],load[1]);
    nsteps++;
    for(a = 0;a < fd->nstat;a++)
      if(fd->sstat[a] >= 0 && fd->fstat[a]){
        struct fd3d *f = fm[fd->sstat[a]/N]; int s = fd->sstat[a]%N;
        fprintf(fd->fstat[a],"%.10e %.10e %.10e %.10e %.6e\n",f->time,f->delta[s],f->V[s],f->tau[s],f->sigma[s]);
      }
    vmax = 0.0;
    for(p = 0;p < 2;p++) for(i = 0;i < N;i++) if(fabs(fm[p]->V[i]) > vmax) vmax = fabs(fm[p]->V[i]);
    if(vmax > vpeak) vpeak = vmax;
    if(fd->back_dyn > 0.0 && vmax < fd->vsw*(fd->reentry > 1.0 ? fd->reentry : 1.0) && pp->static_fn){
      /* as in fd3d_run_phase: FD traction against the QD traction of both planes */
      double dmax = 0.0, *s0 = vtmp, *s1 = vtmp + N;
      pp->static_fn(pp->static_ctx, fm[0]->delta, fm[1]->delta, s0, s1);
      for(p = 0;p < 2;p++) for(i = 0;i < N;i++){
	double qd, dd;
	if(fm[p]->locked[i]) continue;
	qd = fm[p]->tau0[i] + (p ? s1[i] : s0[i]) + load[p][i] - eta*fm[p]->V[i]; dd = fabs(fm[p]->tau[i] - qd);
	if(dd > dmax) dmax = dd;
      }
      if(dmax < fd->back_dyn) ndyn_quiet++; else ndyn_quiet = 0;
      if(ndyn_quiet >= fd->back_nwin) break;
    }
    if(fd->back_eps > 0.0){
      int quiet = (vmax < fd->back_vmax);
      if(quiet && (below > 0.0))
	for(p = 0;p < 2 && quiet;p++) for(i = 0;i < N;i++)
	  if(fabs(fm[p]->V[i] - vref[p*N + i]) > fd->back_eps*(fabs(vref[p*N + i]) + fd->vback)){ quiet = 0; break; }
      if(!quiet || (below == 0.0)){
	for(p = 0;p < 2;p++) memcpy(vref + p*N,fm[p]->V,sizeof(double)*N);
	below = quiet ? 1e-30 : 0.0;
      }else
	below += dt;
    }else{
      if(vmax < fd->vback) below += dt; else below = 0.0;
    }
    if(below >= L*dt) break;
    PetscCheck(nsteps < fd->max_steps,PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,
	       "rsf_fd: 3-D phase exceeded -fd_max_steps = %ld steps (t - t_switch %g s, max V %g)",(long)fd->max_steps,fm[0]->time,vmax);
  }
  for(a = 0;a < fd->nstat;a++) if(fd->fstat[a]){ fclose(fd->fstat[a]); fd->fstat[a] = NULL; }
  {
    double slip[2] = {0, 0}; int cnt[2] = {0, 0};
    for(p = 0;p < 2;p++) for(i = 0;i < N;i++){
      struct fd3d *f = fm[p];
      g = fd->g_of_s[p*N + i];
      if(g < 0) continue;
      x[4*g] = f->f0 + f->b[i]*log(f->V0*f->theta[i]/f->Lc[i]);
      x[4*g+1] = f->tau[i];
      x[4*g+3] += f->delta[i];
      if(pp->normal_coupling) x[4*g+2] = f->sigma[i];
      if(getenv("FD_DEBUG_REENTRY") && (!isfinite(x[4*g]) || !isfinite(x[4*g+1]) || !(x[4*g+2] > 0.0)))
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"hand-back: patch %d plane %d cell %d: psi %g tau %g sigma %g theta %g V %g\n",g,p,i,(double)x[4*g],(double)x[4*g+1],(double)x[4*g+2],f->theta[i],f->V[i]));
      slip[p] += f->delta[i]; cnt[p]++; slipsum += f->delta[i];
    }
    PetscCall(fd_full_put(X,&fd->scX,&fd->Xall,&x,1));
    *t_end = t + fm[0]->time;
    PetscCall(PetscTime(&t_wall1));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd: 3-D phase %ld (two planes): t %.8f yr, %d steps over %.1f s, peak slip rate %.3f m/s, mean slip %.4f m (plane 1) and %.4f m (plane 2), wall %.1f s\n",
			  (long)fd->nphase + 1,(double)(t/SEC_PER_YEAR),nsteps,fm[0]->time,vpeak,slip[0]/(cnt[0] > 0 ? cnt[0] : 1),slip[1]/(cnt[1] > 0 ? cnt[1] : 1),t_wall1 - t_wall0));
  }
  PetscCall(PetscFree4(vtmp,load[0],load[1],vref));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode fd3d_run_phase(TS ts, struct interact_ctx *par, struct fd_vars *fd, PetscReal t, Vec X, PetscReal *t_end)
{
  if(fd->pair3) return fd3d_run_phase_pair(ts,par,fd,t,X,t_end);
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
  PetscCall(fd_full_get(X,&fd->scX,&fd->Xall,&x));
  vmax = 0.0;
  for(i = 0;i < N;i++){
    g = fd->g_of_s[i];
    f->sigma[i] = x[4*g+2];
    f->V[i] = 2.0*rsf->v0*exp(-x[4*g]/fault[g].mu_sa)*sinh(x[4*g+1]/(fault[g].mu_sa*x[4*g+2]));
    f->theta[i] = f->Lc[i]/f->V0*exp((x[4*g] - f->f0)/f->b[i]);       /* psi = f0 + b ln(V0 theta/L) */
    f->tau0[i] = x[4*g+1] + eta*f->V[i];                               /* elastic stress at the switch */
    if(f->layer && f->layer->surf) f->layer->sig0[i] = x[4*g+2];
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
      { PetscMPIInt rk; MPI_Comm_rank(PETSC_COMM_WORLD,&rk); fd->fstat[a] = (rk == 0) ? fopen(fn,"w") : NULL; }   /* rank 0 writes */
      fd->sstat[a] = -1;
      for(i = 0;i < N;i++) if(fd->g_of_s[i] == fd->stat[a]) fd->sstat[a] = i;
    }
  }
  double vprev = 0.0, *vref; int kprev = 0, nretry = 0, khist[7] = {0, 0, 0, 0, 0, 0, 0}, ndyn_quiet = 0;
  PetscCall(PetscMalloc1(N,&vref));
  for(;;){
    double trel, kd, dcmin = 1e300;
    int a, k = 1;
    /* variable substeps as in 2-D: powers of two up to kmax and xi dc/(V dt) */
    vmax = 0.0;
    for(i = 0;i < N;i++){ if(fabs(f->V[i]) > vmax) vmax = fabs(f->V[i]); if(f->Lc[i] < dcmin) dcmin = f->Lc[i]; }
    kd = fd->xi*dcmin/((vmax > 1e-30 ? vmax : 1e-30)*dt);
    /* also limit the step by how fast the peak slip rate changes (relative change of
       at most 10 percent per step, from the previous step), and grow k by at most a
       factor of two per step: in 3-D the fine step can be long in absolute time
       (dt = 0.29 s at 2 km), and the slip-distance criterion alone lets a slowly
       starting rupture be stepped over */
    if(kprev > 0 && vprev > 0.0){
      double rate = fabs(vmax - vprev)/(kprev*dt), ka = (rate > 0.0) ? 0.1*vmax/(rate*dt) : 1e300;
      if(ka < kd) kd = ka;
      if(kd > 2*kprev) kd = 2*kprev;
    }
    if(!fd->kmax3) kd = 1.0;              /* long steps in 3-D only on request (-fd_kmax3 1) */
    while((2*k <= fd->kmax) && (2*k <= kd) && (2*k <= 64)) k *= 2;
    vprev = vmax; kprev = k;
    trel = f->time + k*dt;
    for(i = 0;i < N;i++) load[i] = fault[fd->g_of_s[i]].sinc[0]*trel;    /* backslip loading since the switch */
    for(;;){
      int e = fd3d_step_k(f,load,k);
      if(e >= 0) break;
      PetscCheck((e == -2) && (k > 1),PETSC_COMM_WORLD,PETSC_ERR_LIB,"rsf_fd: 3-D step of %d substeps failed (%d)",k,e);
      k /= 2; kprev = k; nretry++;       /* not converged: halve the step and retry */
      trel = f->time + k*dt;
      for(i = 0;i < N;i++) load[i] = fault[fd->g_of_s[i]].sinc[0]*trel;
    }
    nsteps++;
    { int q = 0; while((1 << q) < k) q++; if(q < 7) khist[q]++; }
    for(a = 0;a < fd->nstat;a++)
      if(fd->sstat[a] >= 0 && fd->fstat[a])
	fprintf(fd->fstat[a],"%.10e %.10e %.10e %.10e %d\n",f->time,f->delta[fd->sstat[a]],f->V[fd->sstat[a]],f->tau[fd->sstat[a]],k);
    vmax = 0.0;
    for(i = 0;i < N;i++) if(fabs(f->V[i]) > vmax) vmax = fabs(f->V[i]);
    if(vmax > vpeak) vpeak = vmax;
    if(fd->back_dyn > 0.0 && vmax < fd->vsw*(fd->reentry > 1.0 ? fd->reentry : 1.0) && f->static_fn){
      /* hand back once the fully dynamic traction equals the quasi-dynamic one (elastostatic
	 interaction plus radiation damping) to within back_dyn everywhere, for back_nwin steps
	 in a row (fd_plan.md L2): from then on the QD integrator reproduces the phase */
      double dmax = 0.0;
      f->static_fn(f->static_ctx, f->delta, vtmp);
      for(i = 0;i < N;i++){
	double qd = f->tau0[i] + vtmp[i] + load[i] - eta*f->V[i], dd = fabs(f->tau[i] - qd);
	if(dd > dmax) dmax = dd;
      }
      if(dmax < fd->back_dyn) ndyn_quiet++; else ndyn_quiet = 0;
      if(ndyn_quiet >= fd->back_nwin) break;
    }
    if(fd->back_eps > 0.0){
      /* quasi-steady hand-back as in 2-D: every slip rate changed by less than back_eps
	 (relative, floored at vback) over a full window, peak below back_vmax */
      int quiet = (vmax < fd->back_vmax);
      if(quiet && (below > 0.0))
	for(i = 0;i < N;i++)
	  if(fabs(f->V[i] - vref[i]) > fd->back_eps*(fabs(vref[i]) + fd->vback)){ quiet = 0; break; }
      if(!quiet || (below == 0.0)){
	memcpy(vref,f->V,sizeof(double)*N);
	below = quiet ? 1e-30 : 0.0;
      }else
	below += k*dt;
    }else{
      if(vmax < fd->vback) below += k*dt; else below = 0.0;
    }
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
    if(f->layer && f->layer->surf && rsf->calc_sigma_dot) x[4*g+2] = f->sigma[i];   /* normal stress from the phase */
    slipsum += f->delta[i];
  }
  PetscCall(fd_full_put(X,&fd->scX,&fd->Xall,&x,1));
  *t_end = t + f->time;
  PetscCall(PetscTime(&t_wall1));
  if(f->layer && f->layer->n_steps > 0)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd: free-surface layer: %.1f outer and %.1f inner iterations per step; time in known part %.1f s, fault solves %.1f s, layer solves %.1f s, finish %.1f s\n",
			  (double)f->layer->n_outer/f->layer->n_steps,(double)f->layer->n_inner/f->layer->n_steps,
			  fd3d_layer_prof[0],fd3d_layer_prof[1],fd3d_layer_prof[2],fd3d_layer_prof[3]));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd: 3-D phase steps by substeps k = 1..64: %d %d %d %d %d %d %d, retries %d\n",
			khist[0],khist[1],khist[2],khist[3],khist[4],khist[5],khist[6],nretry));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd: 3-D phase %ld: t %.8f yr, %d steps over %.1f s, peak slip rate %.3f m/s, mean slip %.4f m, wall %.1f s\n",
			(long)fd->nphase + 1,(double)(t/SEC_PER_YEAR),nsteps,f->time,vpeak,slipsum/N,(double)(t_wall1 - t_wall0)));
  PetscCall(PetscFree2(vtmp,load));
  PetscCall(PetscFree(vref));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  3-D: square elements of equal size forming a complete rectangular grid on
  one plane, whole space, slip along strike (-rsf_slip_mode strike) or dip.
  Builds the fd3d tables (grid axis 1 = slip direction) and checks their
  static limit against the assembled Is, as for 2-D.
*/

/* exact statics for two planes: shear traction (Is) and normal traction (-In, tension positive)
   on both planes from both planes' slips; combined marcher index p*N + i, locked cells -1 */
static void fd3d_pair_apply(struct interact_ctx *par, Mat A, double sgn, const double *d0, const double *d1, double *o0, double *o1)
{
  struct fd_vars *fd = par->medium->rsf->fd;
  PetscScalar *a, *b;
  PetscInt lo, hi;
  int i, p, N = fd->f3->N;
  const double *din[2] = {d0, d1};
  double *out[2] = {o0, o1};
  PetscCallAbort(PETSC_COMM_WORLD,VecGetOwnershipRange(fd->v3a,&lo,&hi));
  PetscCallAbort(PETSC_COMM_WORLD,VecGetArray(fd->v3a,&a));
  for(p = 0;p < 2;p++) for(i = 0;i < N;i++){
    int g = fd->g_of_s[p*N + i];
    if((g >= lo) && (g < hi)) a[g - lo] = din[p][i];
  }
  PetscCallAbort(PETSC_COMM_WORLD,VecRestoreArray(fd->v3a,&a));
  PetscCallAbort(PETSC_COMM_WORLD,MatMult(A,fd->v3a,fd->v3b));
  PetscCallAbort(PETSC_COMM_WORLD,fd_full_get_read(fd->v3b,&fd->sc3,&fd->v3all,&b));
  for(p = 0;p < 2;p++) for(i = 0;i < N;i++){
    int g = fd->g_of_s[p*N + i];
    out[p][i] = (g >= 0) ? sgn*PetscRealPart(b[g]) : 0.0;
  }
  PetscCallAbort(PETSC_COMM_WORLD,fd_full_put_read(fd->v3b,&fd->sc3,&fd->v3all,&b));
}
static void fd3d_pair_static_is(void *ctx, const double *d0, const double *d1, double *o0, double *o1)
{
  struct interact_ctx *par = (struct interact_ctx *)ctx;
  fd3d_pair_apply(par,par->medium->Is,1.0,d0,d1,o0,o1);
}
static void fd3d_pair_static_in(void *ctx, int p, const double *d_other, const double *d_own, double *out)
{
  struct interact_ctx *par = (struct interact_ctx *)ctx;
  struct fd_vars *fd = par->medium->rsf->fd;
  int N = fd->f3->N;
  double *o0 = malloc(sizeof(double)*N), *o1 = malloc(sizeof(double)*N);
  if(p == 0) fd3d_pair_apply(par,par->medium->In,-1.0,d_own,d_other,o0,o1);
  else fd3d_pair_apply(par,par->medium->In,-1.0,d_other,d_own,o0,o1);
  memcpy(out,p == 0 ? o0 : o1,sizeof(double)*N);
  free(o0); free(o1);
}

/* setup for two parallel planes of equal squares (stepovers): strike slip, whole space or
   (vertical planes) half space by the mirror image; one lattice for both (along strike and
   by depth), cells of a plane without a patch are locked */
static PetscErrorCode fd3d_setup_pair(struct interact_ctx *par, struct fd_vars *fd)
{
  struct med *medium = par->medium;
  struct flt *fault = par->fault;
  struct rsf_vars *rsf = medium->rsf;
  PetscInt n = medium->nrflt, gi, c, ncols;
  const PetscInt *cols;
  const PetscScalar *vals;
  double dx, ls = fd->length_scale, e1[3], e2[3], nn[3], x0[3], e3t[3], dpl[2] = {0, 0}, lo1 = 1e300, lo2 = 1e300, hi1 = -1e300, hi2 = -1e300, *pa, *pb, sgn, Wt;
  double dev[2] = {0, 0}, nrm[2] = {0, 0}, devn = 0.0, nrmn = 0.0;
  int *plane, nx, nz, N, i, k, p;
  PetscFunctionBeginUser;
  PetscCheck(rsf->slip_mode == STRIKE,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: two 3-D planes need strike slip (-rsf_slip_mode 0)");
  dx = 2.0*fault[0].l*ls;
  for(k = 0;k < 3;k++){ e1[k] = fault[0].t_strike[k]; e2[k] = fault[0].t_dip[k]; nn[k] = fault[0].normal[k]; x0[k] = fault[0].x[k]*ls; }
  if(e2[2] > 0.0) for(k = 0;k < 3;k++) e2[k] = -e2[k];          /* rows by depth */
  PetscCall(PetscMalloc3(n,&pa,n,&pb,n,&plane));
  dpl[1] = 1e300;
  for(i = 0;i < n;i++){
    double d[3], dn;
    PetscCheck((fabs(2.0*fault[i].l*ls - dx) < 1e-6*dx) && (fabs(2.0*fault[i].w*ls - dx) < 1e-6*dx),PETSC_COMM_WORLD,PETSC_ERR_SUP,
	       "-fd_mode: 3-D elements must be squares of one size (patch %d)",i);
    PetscCheck(fabs(fault[i].normal[0]*nn[0] + fault[i].normal[1]*nn[1] + fault[i].normal[2]*nn[2]) > 1.0 - 1e-9,PETSC_COMM_WORLD,PETSC_ERR_SUP,
	       "-fd_mode: the 3-D planes must be parallel (patch %d)",i);
    for(k = 0;k < 3;k++) d[k] = fault[i].x[k]*ls - x0[k];
    dn = d[0]*nn[0] + d[1]*nn[1] + d[2]*nn[2];
    if(fabs(dn) < 1e-6*dx) plane[i] = 0;
    else{
      if(dpl[1] > 1e299) dpl[1] = dn;
      PetscCheck(fabs(dn - dpl[1]) < 1e-6*dx,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: at most two parallel 3-D planes (patch %d at %g m)",i,dn);
      plane[i] = 1;
    }
    pa[i] = d[0]*e1[0] + d[1]*e1[1] + d[2]*e1[2]; pb[i] = d[0]*e2[0] + d[1]*e2[1] + d[2]*e2[2];
    if(pa[i] < lo1) lo1 = pa[i];
    if(pa[i] > hi1) hi1 = pa[i];
    if(pb[i] < lo2) lo2 = pb[i];
    if(pb[i] > hi2) hi2 = pb[i];
  }
  nx = (int)floor((hi1 - lo1)/dx + 0.5) + 1; nz = (int)floor((hi2 - lo2)/dx + 0.5) + 1; N = nx*nz;
  fd->mirror3 = PETSC_FALSE; fd->surf3 = PETSC_FALSE;
  if(!medium->full_space){
    double zmax = -1e300;
    PetscCheck(fabs(nn[2]) < 1e-9,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: two 3-D planes in a half space must be vertical (mirror image)");
    for(i = 0;i < n;i++) if(fault[i].x[INT_Z]*ls > zmax) zmax = fault[i].x[INT_Z]*ls;
    PetscCheck(fabs(zmax + 0.5*dx) < 1e-6*dx,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: 3-D half space: the shallowest row must reach the surface");
    fd->mirror3 = PETSC_TRUE;
  }
  /* the receiver plane's offset along the kernel frame's x3: the kernels take the slip as the
     motion of the +x3 side relative to the -x3 side, interact as that of the side its normal
     points to, so x3 is measured along interact's normal (the in-plane axes then follow; the
     strike-slip kernels are even in x2). Checked by the normal-traction gate against In. */
  (void)e3t; sgn = 1.0;
  Wt = sgn*dpl[1];
  PetscCall(PetscMalloc1(2*N,&fd->g_of_s));
  for(i = 0;i < 2*N;i++) fd->g_of_s[i] = -1;
  for(i = 0;i < n;i++){
    double fi = (pa[i] - lo1)/dx, fj = (pb[i] - lo2)/dx;
    int ii = (int)floor(fi + 0.5), jj = (int)floor(fj + 0.5), s = plane[i]*N + ii*nz + jj;
    PetscCheck((fabs(fi - ii) < 1e-6) && (fabs(fj - jj) < 1e-6),PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: patch %d is off the common 3-D lattice",i);
    PetscCheck(fd->g_of_s[s] < 0,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: duplicate patch in the 3-D lattice (patch %d)",i);
    fd->g_of_s[s] = i;
  }
  fd->dt_min = fd->cfl*dx/fd->beta;
  fd->nx3 = nx; fd->nz3 = nz; fd->n3tot = 2*N;
#ifdef FD3D_MPI
  { PetscMPIInt sz; PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&sz)); if(sz > 1) fd3d_set_comm(PETSC_COMM_WORLD); }
#endif
  fd->tab3 = fd3d_tab_create(nx,fd->mirror3 ? 2*nz : nz,dx,fd->dt_min,fd->alpha,fd->beta,fd->rho,(int)fd->lmax);
  PetscCheck(fd->tab3,PETSC_COMM_WORLD,PETSC_ERR_MEM,"rsf_fd_setup: 3-D table allocation failed");
  fd->f3 = fd3d_create2(fd->tab3,nz,fd->mirror3 ? 1 : 0); fd->f3b = fd3d_create2(fd->tab3,nz,fd->mirror3 ? 1 : 0);
  PetscCheck(fd->f3 && fd->f3b,PETSC_COMM_WORLD,PETSC_ERR_MEM,"rsf_fd_setup: 3-D marcher allocation failed");
  for(p = 0;p < 2;p++){
    struct fd3d *f = p ? fd->f3b : fd->f3;
    f->law = 1; f->f0 = rsf->f0; f->V0 = rsf->v0; f->nl = (int)fd->nl;
    PetscCall(PetscMalloc4(N,&f->a,N,&f->b,N,&f->Lc,N,&f->theta));
    f->locked = (int *)calloc(N,sizeof(int));
    for(i = 0;i < N;i++){
      int g = fd->g_of_s[p*N + i];
      if(g >= 0){ f->a[i] = fault[g].mu_sa; f->b[i] = fault[g].mu_db; f->Lc[i] = rsf->dc_vec ? rsf->dc_vec[g]*1.0 : rsf->dc; }
      else{ f->a[i] = 0.01; f->b[i] = 0.01; f->Lc[i] = 1.0; f->locked[i] = 1; f->theta[i] = 1.0; }
    }
  }
  fd->pair3 = fd3d_pair_create(fd->f3,fd->f3b,Wt,1);
  PetscCheck(fd->pair3,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: the two 3-D planes (%g m apart) interact within one step; use smaller elements",fabs(dpl[1]));
  fd->pair3->static_fn = fd3d_pair_static_is; fd->pair3->static_ctx = (void *)par;
  fd->pair3->normal_coupling = rsf->calc_sigma_dot ? 1 : 0;
  if(rsf->calc_sigma_dot){ fd->pair3->static_n_fn = fd3d_pair_static_in; fd->pair3->static_n_ctx = (void *)par; }
  /* static gate: Is against the self (p = q) and cross (p != q) tables, off the 3 x 3 self
     neighbourhood; the cross tables' static slopes are in pair3->ks_off */
  {
    PetscInt rlo, rhi;
    int *pl_of = malloc(sizeof(int)*n), *s_of = malloc(sizeof(int)*n);
    for(i = 0;i < 2*N;i++) if(fd->g_of_s[i] >= 0){ pl_of[fd->g_of_s[i]] = i/N; s_of[fd->g_of_s[i]] = i%N; }
    PetscCall(MatGetOwnershipRange(medium->Is,&rlo,&rhi));
    for(gi = rlo;gi < rhi && gi < n;gi++){
      int si = s_of[gi]/nz, sj = s_of[gi]%nz, pr = pl_of[gi];
      PetscCall(MatGetRow(medium->Is,gi,&ncols,&cols,&vals));
      for(c = 0;c < ncols;c++){
        int ti = s_of[cols[c]]/nz, tj = s_of[cols[c]]%nz, di = si - ti, dj = sj - tj, ps = pl_of[cols[c]], self = (ps == pr), nt = fd->mirror3 ? 2*nz : nz;
        double ours = 0.0, dv;
        const double *ko = self ? fd->tab3->ks_off : fd->pair3->ks_off[ps][0];
        if(self && (abs(di) <= 1) && (abs(dj) <= 1)) continue;
        if(fd->mirror3){
          int djr = (nz + sj) - (nz + tj), dji = (nz + sj) - (nz - 1 - tj);
          ours = ko[(size_t)(di + nx - 1)*(2*nt - 1) + (djr + nt - 1)] + ko[(size_t)(di + nx - 1)*(2*nt - 1) + (dji + nt - 1)];
        }else
          ours = ko[(size_t)(di + nx - 1)*(2*nt - 1) + (dj + nt - 1)];
        dv = PetscRealPart(vals[c]) - ours;
        dev[!self] += dv*dv; nrm[!self] += PetscRealPart(vals[c])*PetscRealPart(vals[c]);
      }
      PetscCall(MatRestoreRow(medium->Is,gi,&ncols,&cols,&vals));
    }
    if(medium->In){                                 /* normal traction across planes: tables (tension +) against -In */
      for(gi = rlo;gi < rhi && gi < n;gi++){
        int si = s_of[gi]/nz, sj = s_of[gi]%nz, pr = pl_of[gi];
        PetscCall(MatGetRow(medium->In,gi,&ncols,&cols,&vals));
        for(c = 0;c < ncols;c++){
          int ti = s_of[cols[c]]/nz, tj = s_of[cols[c]]%nz, di = si - ti, dj = sj - tj, ps = pl_of[cols[c]], nt = fd->mirror3 ? 2*nz : nz;
          double ours, dv;
          const double *ko = fd->pair3->ks_off[ps][1];
          if(ps == pr) continue;
          if(fd->mirror3){
            int djr = (nz + sj) - (nz + tj), dji = (nz + sj) - (nz - 1 - tj);
            ours = ko[(size_t)(di + nx - 1)*(2*nt - 1) + (djr + nt - 1)] + ko[(size_t)(di + nx - 1)*(2*nt - 1) + (dji + nt - 1)];
          }else
            ours = ko[(size_t)(di + nx - 1)*(2*nt - 1) + (dj + nt - 1)];
          dv = -PetscRealPart(vals[c]) - ours;
          devn += dv*dv; nrmn += PetscRealPart(vals[c])*PetscRealPart(vals[c]);
        }
        PetscCall(MatRestoreRow(medium->In,gi,&ncols,&cols,&vals));
      }
    }
    { double sums[6] = {dev[0], dev[1], nrm[0], nrm[1], devn, nrmn}; PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE,sums,6,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD));
      dev[0] = sqrt(sums[0]/(sums[2] + 1e-300)); dev[1] = sqrt(sums[1]/(sums[3] + 1e-300)); devn = (sums[5] > 0.0) ? sqrt(sums[4]/sums[5]) : -1.0; }
    free(pl_of); free(s_of);
  }
  PetscCall(PetscFree3(pa,pb,plane));
  PetscCall(MatCreateVecs(medium->Is,&fd->v3a,&fd->v3b));
  fd->nrec = 64; fd->irec = 0; fd->nrec_used = 0;
  PetscCall(PetscMalloc2(fd->nrec,&fd->rec_t,(size_t)fd->nrec*2*N,&fd->rec_v));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D, two parallel planes %g m apart, one lattice of %d x %d squares of %g m, dt_min %.3e s, window %d lags, FFT %d x %d, %s\n",
			fabs(dpl[1]),nx,nz,dx,(double)fd->dt_min,fd->tab3->L,fd->tab3->Mx,fd->tab3->Mz,fd->mirror3 ? "half space (mirror image)" : "whole space"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D static gate, two planes: rel dev of the tables from Is %.3e within planes, %.3e across, normal traction across from -In %.3e%s; normal-stress coupling %s\n",
			dev[0],dev[1],devn,fd->mirror3 ? " (mirror: a few percent expected)" : "",fd->pair3->normal_coupling ? "on (exact In)" : "off"));
  PetscCheck((devn < 0.0) || fd->mirror3 || (devn <= 5.0*fd->gate_tol),PETSC_COMM_WORLD,PETSC_ERR_USER,"rsf_fd_setup: 3-D normal-traction gate failed for the two planes (%g)",devn);
  PetscCheck(fd->mirror3 || ((dev[0] <= 5.0*fd->gate_tol) && (dev[1] <= 5.0*fd->gate_tol)),PETSC_COMM_WORLD,PETSC_ERR_USER,
	     "rsf_fd_setup: 3-D static gate failed for the two planes (%g, %g)",dev[0],dev[1]);
  PetscFunctionReturn(PETSC_SUCCESS);
}

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
  PetscCheck((rsf->slip_mode == STRIKE) || (rsf->slip_mode == DIP),PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: 3-D needs a fixed slip direction (-rsf_slip_mode strike or dip)");
  ls = fd->length_scale;
  dx = 2.0*fault[0].l*ls;
  {
    /* two parallel planes (stepovers, fd_plan.md P2): elements off the plane of patch 0 */
    int off = 0;
    for(i = 0;i < n;i++){
      double dn = 0.0;
      for(k = 0;k < 3;k++) dn += (fault[i].x[k] - fault[0].x[k])*ls*fault[0].normal[k];
      if(fabs(dn) > 1e-6*dx) off = 1;
    }
    if(off) return fd3d_setup_pair(par,fd);
  }
  /* half space: a vertical strike-slip fault reaching the surface uses the mirror image (plus
     -fd_layer); any other planar fault (or -fd_surface 1) uses the full virtual surface, with
     grid axis 1 along strike and slip along axis 1 (strike) or axis 2 (dip) */
  fd->surf3 = PETSC_FALSE;
  if(!medium->full_space){
    int vertical = fabs(fault[0].normal[2]) < 1e-9;
    if(fd->surface_opt || !vertical || (rsf->slip_mode == DIP)) fd->surf3 = PETSC_TRUE;
  }
  for(k = 0;k < 3;k++){
    if(fd->surf3){ e1[k] = fault[0].t_strike[k]; e2[k] = fault[0].t_dip[k]; }
    else{
      e1[k] = (rsf->slip_mode == STRIKE) ? fault[0].t_strike[k] : fault[0].t_dip[k];
      e2[k] = (rsf->slip_mode == STRIKE) ? fault[0].t_dip[k] : fault[0].t_strike[k];
    }
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
  if(fd->surf3){
    /* rows ordered by depth (row 0 shallowest): flip the dip axis if it points up */
    if(e2[2] > 0.0){
      for(i = 0;i < n;i++) pb[i] = -pb[i];
      { double tmp = lo2; lo2 = -hi2; hi2 = -tmp; }
      for(k = 0;k < 3;k++) e2[k] = -e2[k];
    }
  }else if(!medium->full_space){
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
#ifdef FD3D_MPI
  { PetscMPIInt sz; PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&sz)); if(sz > 1) fd3d_set_comm(PETSC_COMM_WORLD); }
#endif
  if(fd->surf3){
    /* a fault whose top edge lies in the surface continues above it by next rows tied to the
       top row (fd_plan.md B7): the table has nz + next rows */
    double sind0 = fabs(e2[2]), zc0 = fault[fd->g_of_s[0]].x[INT_Z]*ls, ztop0 = zc0 + 0.5*dx*sind0;
    fd->next3 = (fd->surf_extend >= 0) ? (int)fd->surf_extend : ((fabs(ztop0) < 1e-6*dx) ? 1 : 0);
    fd->tab3 = fd3d_tab_create_c(nx,nz + fd->next3,dx,fd->dt_min,fd->alpha,fd->beta,fd->rho,(int)fd->lmax,(rsf->slip_mode == STRIKE) ? 1 : 2);
  }else
    fd->tab3 = fd3d_tab_create(nx,fd->mirror3 ? 2*nz : nz,dx,fd->dt_min,fd->alpha,fd->beta,fd->rho,(int)fd->lmax);
  PetscCheck(fd->tab3,PETSC_COMM_WORLD,PETSC_ERR_MEM,"rsf_fd_setup: 3-D table allocation failed");
  /* static gate: Is(receiver, source) against the offset table, skipping the
     3 x 3 self neighbourhood whose regularisation may differ */
  {
    PetscInt rlo, rhi;                             /* MPI: own rows of Is, sums reduced below */
    PetscCall(MatGetOwnershipRange(medium->Is,&rlo,&rhi));
  for(gi = rlo;gi < rhi && gi < n;gi++){
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
    { double sums[2] = {dev, nrm}; PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE,sums,2,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD)); dev = sums[0]; nrm = sums[1]; }
  }
  PetscCall(PetscFree(s_of_g));
  dev = sqrt(dev/(nrm + 1e-300));
  fd->f3 = fd->surf3 ? fd3d_create3(fd->tab3,nz,0,fd->next3) : fd3d_create2(fd->tab3,nz,fd->mirror3 ? 1 : 0);
  PetscCheck(fd->f3,PETSC_COMM_WORLD,PETSC_ERR_MEM,"rsf_fd_setup: 3-D marcher allocation failed");
  fd->f3->law = 1; fd->f3->f0 = rsf->f0; fd->f3->V0 = rsf->v0;
  fd->f3->nl = (int)fd->nl;
  PetscCall(PetscMalloc4(n,&fd->f3->a,n,&fd->f3->b,n,&fd->f3->Lc,n,&fd->f3->theta));
  for(i = 0;i < n;i++){
    int g = fd->g_of_s[i];
    fd->f3->a[i] = fault[g].mu_sa; fd->f3->b[i] = fault[g].mu_db; fd->f3->Lc[i] = rsf->dc_vec ? rsf->dc_vec[g]*1.0 : rsf->dc;
  }
  fd->f3->static_fn = fd3d_static_is; fd->f3->static_ctx = (void *)par;
  {
    PetscInt nk = 0;                              /* -fd_nk 1: Newton-Krylov for long steps (fd_plan.md L1) */
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_nk",&nk,NULL));
    if(nk) fd->f3->nk = (int)nk;
  }
  {
    PetscInt ts = 0;                              /* -fd_table_statics 1: the tables' own statics instead of Is (experiments) */
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_table_statics",&ts,NULL));
    if(ts){ fd->f3->static_fn = NULL; fd->f3->static_ctx = NULL; PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D static part from the tables (-fd_table_statics 1), not from Is\n")); }
  }
  if(fd->surf3){
    /* full virtual surface: the marcher's frame has x' along strike (column 0 at x' = 0),
       y' horizontal towards the up-dip side, z up; element (0, 0) centred at y' = 0 */
    double sind = -e2[2], cosd = sqrt(1.0 - sind*sind), dipdeg = asin(sind)*180.0/M_PI, zc0, ytop, ztop, ys0, Wd, fsign;
    /* surface elements of dx/m: m = 4 by default. With the NL scheme the surface's shear
       components grew in a random-perturbation test at cS dt/hs = 1 (m = 2) and decayed at 0.5
       (m = 1) and 2 (m = 4); m = 4 also meets the static accuracy rule (a quarter of the fault
       element near the trace) */
    int m = (int)fd->surf_m, mg = (int)(fd->surf_ext > 0 ? fd->surf_ext : (nz > 2 ? nz : 2)), nys, nxs, Lsf;
    /* the unstable band (cS dt/hs about 0.85 to 1.05, fd_plan.md B3): with the default, keep
       cS dt/hs = m cfl at 2 or more (m = 4 at cfl 0.5, more for smaller steps) */
    if(fd->surf_m_default) while(m*fd->cfl < 1.999) m++;
    if(m*fd->cfl > 0.8 && m*fd->cfl < 1.2)
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: WARNING: surface elements of dx/%d give cS dt/hs = %.2f, inside the band where the surface's shear modes grow (0.85 to 1.05)\n",m,(double)(m*fd->cfl)));
    zc0 = fault[fd->g_of_s[0]].x[INT_Z]*ls;
    ytop = 0.5*dx*cosd; ztop = zc0 + 0.5*dx*sind;
    PetscCheck(ztop <= 1e-6*dx,PETSC_COMM_WORLD,PETSC_ERR_SUP,"-fd_mode: the fault must lie below the surface (top edge at z = %g m)",ztop);
    Wd = nz*dx*cosd;
    /* rows: the fault's surface projection plus mg fault elements either side; the top edge line
       (y' = ytop) on a row boundary when the fault reaches the surface */
    ys0 = ytop - Wd - mg*dx;
    ys0 = ytop - ceil((ytop - ys0)/(dx/m))*(dx/m);
    nys = (int)ceil((ytop + mg*dx - ys0)/(dx/m));
    nxs = (nx + 2*mg)*m;
    Lsf = (int)ceil(sqrt(pow((nx + 2*mg)*dx, 2) + pow(Wd + 2*mg*dx, 2) + pow(nz*dx*sind - ztop, 2))/(fd->beta*fd->dt_min)) + 2;
    {
      /* sign of the QD slip direction in the surface code's frame (x' = strike, e2 = up dip,
         y' = horizontal part of up dip): dip slip along interact's t_dip, strike slip along
         t_strike in a frame that may be left-handed */
      double up[3] = {-e2[0], -e2[1], -e2[2]}, yh[3] = {up[0], up[1], 0.0}, yn = sqrt(up[0]*up[0] + up[1]*up[1]), cr;
      if(yn > 1e-12){ yh[0] /= yn; yh[1] /= yn; } else { yh[0] = -e1[1]; yh[1] = e1[0]; }
      cr = e1[0]*yh[1] - e1[1]*yh[0];                                  /* (x' x y') . z */
      fsign = (rsf->slip_mode == STRIKE) ? (cr > 0 ? 1.0 : -1.0) : (fault[0].t_dip[2] > 0.0 ? 1.0 : -1.0);
      {
        /* side: the kernels move the side their normal e3 = e1 x e2 (frame of the surface code,
           in global coordinates (x' along e1, y' horizontal up dip, z)) points to, interact the side
           its normal points to; a mismatch reverses the slip variable as well (fd_plan.md P2) */
        /* in the surface code's own frame (x' = e1, y' = horizontal up dip, z up), where the
           kernels live: its e3 = (0, -sin(dip), cos(dip)); interact's normal there */
        double nloc[3] = {nn[0]*e1[0] + nn[1]*e1[1] + nn[2]*e1[2], nn[0]*yh[0] + nn[1]*yh[1], nn[2]};
        double ss = (-sind*nloc[1] + cosd*nloc[2]) > 0.0 ? 1.0 : -1.0;
        if(ss < 0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: virtual surface: interact's fault normal points to the kernel frame's -e3 side\n"));
        fsign *= ss;
      }
    }
    if(fd->surf_grade){
      /* graded rows: elements of dx/m within one fault element of the trace (y' = ytop), then
         rows doubling in width (up to dx) out to the margins; all boundaries on the dx/m lattice
         through the trace */
      /* row width at y': no more than half the depth of the fault below that point (within the
         fault's surface projection d = top depth + (ytop - y') tan(dip); beyond the bottom edge
         the bottom depth; above the trace the distance from it), a power of two times dx/m
         between dx/m and dx, growing at most by two from row to row */
      double hx = dx/m, *yc, *wy, y, w, ylo = ytop - Wd - mg*dx, yhi = ytop + mg*dx;
      int nr = 0, cap = 4096;
      yc = malloc(sizeof(double)*cap); wy = malloc(sizeof(double)*cap);
#define ROW_W(dep, wprev) do{ double wa = 0.5*(dep); if(wa > dx) wa = dx; if(wa > 2*(wprev)) wa = 2*(wprev); w = hx; while(2*w <= wa + 1e-9) w *= 2.0; }while(0)
      /* dep: distance from the surface point to the nearest part of the fault (vertical depth over
         the projection, distance to the top or bottom edge line outside it) */
      /* dep: distance from the surface point (y, 0) to the fault's cross-section, the segment from
         the top edge (ytop, ztop) to the bottom edge (ytop - Wd, ztop - nz dx sin(dip)) */
#define SEG_DIST(yy) ({ double ay = ytop, az = ztop, by = ytop - Wd, bz = ztop - nz*dx*sind, vy = by - ay, vz = bz - az, \
        tt = ((yy - ay)*vy + (0.0 - az)*vz)/(vy*vy + vz*vz + 1e-30); if(tt < 0) tt = 0; if(tt > 1) tt = 1; \
        sqrt(pow(yy - (ay + tt*vy), 2) + pow(0.0 - (az + tt*vz), 2)); })
      for(y = ytop, w = hx; y < yhi - 1e-9; ){
        double dep = SEG_DIST(y);
        ROW_W(dep, w); if(y + w > yhi + 1e-9) break; yc[nr] = y + 0.5*w; wy[nr++] = w; y += w;
      }
      for(y = ytop, w = hx; y > ylo + 1e-9; ){
        double dep = SEG_DIST(y);
        ROW_W(dep, w); if(y - w < ylo - 1e-9) break; yc[nr] = y - 0.5*w; wy[nr++] = w; y -= w;
      }
#undef SEG_DIST
#undef ROW_W
      { int a,b; for(a = 0; a < nr; a++) for(b = a + 1; b < nr; b++) if(yc[b] < yc[a]){ double t1 = yc[a], t2 = wy[a]; yc[a] = yc[b]; wy[a] = wy[b]; yc[b] = t1; wy[b] = t2; } }
      nys = nr;
      fd->f3->layer = fd3d_surf_create_rows(fd->f3,-0.5*dx,ytop,ztop,dipdeg,(rsf->slip_mode == STRIKE) ? 1 : 2,m,nxs,nys,-0.5*dx - mg*dx,yc,wy,Lsf,fsign);
      free(yc); free(wy);
    }else
      fd->f3->layer = fd3d_surf_create(fd->f3,-0.5*dx,ytop,ztop,dipdeg,(rsf->slip_mode == STRIKE) ? 1 : 2,m,nxs,nys,-0.5*dx - mg*dx,ys0,Lsf,fsign);
    PetscCheck(fd->f3->layer,PETSC_COMM_WORLD,PETSC_ERR_LIB,"rsf_fd_setup: virtual surface setup failed");
    fd->f3->layer->normal_coupling = rsf->calc_sigma_dot ? 1 : 0;
    if(rsf->calc_sigma_dot){ fd->f3->layer->static_n_fn = fd3d_static_in; fd->f3->layer->static_n_ctx = (void *)par; }
    {
      struct fd3d_layer *ly = fd->f3->layer;
      double mb = 16.0*((double)5*(ly->L + 2)*ly->Mh*ly->nys*ly->nys + 9.0*(ly->L + 2)*ly->Mh*ly->nys*ly->nz + (double)ly->L*ly->Mh*(3*ly->nys + ly->nz))/1048576.0;
      double wmin = 1e300, wmax = 0.0; int q; for(q = 0; q < ly->nys; q++){ if(ly->wy[q] < wmin) wmin = ly->wy[q]; if(ly->wy[q] > wmax) wmax = ly->wy[q]; }
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: virtual free surface: dip %.1f deg, %s slip, %d columns of %g m x %d rows of %g to %g m (3 unknowns each), %d lags, about %.0f MB; normal-stress changes %s\n",
			    dipdeg,(rsf->slip_mode == STRIKE) ? "strike" : "dip",ly->nxs,ly->hx,ly->nys,wmin,wmax,ly->L,mb,rsf->calc_sigma_dot ? "on (exact static In)" : "off (no -calc_sigma_dot)"));
      if(fd->next3) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: the fault reaches the surface: its plane continues %d row(s) above it, tied to the top row\n",(int)fd->next3));
    }
  }
  if(fd->mirror3 && (fd->layer_w > 0)){
    /* free-surface layer completing the mirror (fd_plan.md, stage A): a strip of opening
       elements of dx/m along the trace, +-layer_w fault elements wide, layer_margin fault
       elements beyond the fault ends */
    int m = (int)fd->layer_m, mg = (int)fd->layer_margin, Lly;
    double x0f = -0.5*dx, rmax;                     /* the marcher's along-strike origin: element 0 centred at 0 */
    if(mg < 0) mg = (int)fd->layer_w;
    rmax = sqrt(pow((nx + 2*mg)*dx, 2) + pow((nz + fd->layer_w)*dx, 2));
    Lly = (int)ceil(rmax/(fd->beta*fd->dt_min)) + 2;
    fd->f3->layer = fd3d_layer_create_fft(fd->f3,x0f,m,(nx + 2*mg)*m,(int)(2*fd->layer_w*m),x0f - mg*dx,Lly);
    PetscCheck(fd->f3->layer,PETSC_COMM_WORLD,PETSC_ERR_LIB,"rsf_fd_setup: free-surface layer setup failed");
    {
      struct fd3d_layer *ly = fd->f3->layer;
      double mb = 16.0*((double)ly->L*ly->Mx2*ly->Mh2 + 2.0*ly->L*ly->Mh*ly->nys*ly->nz + (double)ly->L*ly->Mh*(ly->nys + ly->nz) + (double)ly->L*ly->Mx2*ly->Mh2)/1048576.0;
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: free-surface layer: %d x %d opening elements of %g m (+-%d fault elements across the trace, %d beyond the ends), %d lags, about %.0f MB\n",
			    ly->nxs,ly->nys,ly->hx,(int)fd->layer_w,mg,ly->L,mb));
    }
  }
  PetscCall(MatCreateVecs(medium->Is,&fd->v3a,&fd->v3b));
  fd->nrec = 64; fd->irec = 0; fd->nrec_used = 0;
  PetscCall(PetscMalloc2(fd->nrec,&fd->rec_t,(size_t)fd->nrec*n,&fd->rec_v));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D, %d x %d squares of %g m, dt_min %.3e s, window %d lags, FFT %d x %d, tables %.1f MB\n",
			nx,nz,dx,(double)fd->dt_min,fd->tab3->L,fd->tab3->Mx,fd->tab3->Mz,
			(double)fd->tab3->L*fd->tab3->Mx*fd->tab3->Mh*sizeof(fd3d_cplx)/1048576.0));
  if(fd->mirror3)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D half space by mirror image: static operator of fault plus image deviates %.3e from the exact half-space Is "
			  "(expected, a few percent); the dynamic phases use the exact Is for the static part\n",dev));
  else if(!fd->surf3)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D static gate: rel dev of Ks(tables) from Is %.3e (tol %g)\n",dev,(double)fd->gate_tol));
  if(fd->surf3)
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_setup: 3-D half space by a virtual surface: the full-space tables deviate %.3e from the exact half-space Is "
			  "(the free-surface part); the dynamic phases use the exact Is (and In) for the static part\n",dev));
  PetscCheck(fd->mirror3 || fd->surf3 || (dev <= 5.0*fd->gate_tol),PETSC_COMM_WORLD,PETSC_ERR_USER,
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
  PetscScalar *v;
  PetscInt i;
  PetscFunctionBeginUser;
  if(!fd)
    PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(fd_full_get_read(rsf->vel,&fd->scV,&fd->Vall,&v));
  fd->rec_t[fd->irec] = time;
  {
    int N = fd->dim3 ? (fd->pair3 ? fd->n3tot : fd->f3->N) : fd->sys->N;
    double vm = 0.0;
    for(i = 0;i < N;i++){
      fd->rec_v[(size_t)fd->irec*N + i] = (fd->g_of_s[i] >= 0) ? PetscRealPart(v[fd->g_of_s[i]]) : 0.0;
      if(fabs(fd->rec_v[(size_t)fd->irec*N + i]) > vm) vm = fabs(fd->rec_v[(size_t)fd->irec*N + i]);
    }
    if((fd->vsw_eff > fd->vsw) && (vm < fd->vback)) fd->vsw_eff = fd->vsw;   /* end of the re-entry hysteresis: afterslip below vback */
  }
  PetscCall(fd_full_put_read(rsf->vel,&fd->scV,&fd->Vall,&v));
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
  if(getenv("FD_DEBUG_REENTRY")) PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_postevent: t %.9f yr, switching threshold %g\n",(double)(t/SEC_PER_YEAR),(double)fd->vsw_eff));
  if(fd->mode == 1){
    if(fd->dim3) PetscCall(fd3d_run_phase(ts,par,fd,t,X,&t_end));
    else PetscCall(fd_run_phase(ts,par,fd,t,X,&t_end));
    fd->nphase++;
    fd->nrec_used = 0; fd->irec = 0;	/* the record is stale after the phase */
    /* re-entry hysteresis (fd_plan.md L3): right after a phase, the afterslip may exceed vsw
       again; QD reproduces it (the hand-back tested that), so a new phase opens only above
       reentry vsw until the peak rate has fallen below vback (reset in rsf_fd_monitor) */
    if(fd->reentry > 1.0) fd->vsw_eff = fd->reentry*fd->vsw;
    PetscCall(TSSetTime(ts,t_end));
    PetscCall(TSSetTimeStep(ts,fd->dt_restart));
    /* the event handler otherwise resumes with its own step (the pre-event one), which after an
       early hand-back (-fd_back_dyn) can be far too long for the still-fast afterslip (the
       error estimate turned NaN and the integration stalled, fd_plan.md L3) */
    if(fd->dim3) PetscCall(TSSetPostEventStep(ts,fd->dt_restart));   /* 3-D only: 2-D runs unchanged */
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
  PetscScalar *x;
  PetscReal t, t_end, vmax = 0.0, v;
  PetscInt i;
  PetscFunctionBeginUser;
  if(!fd || (fd->mode != 1))
    PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(fd_full_get_read(X,&fd->scX,&fd->Xall,&x));
  for(i = 0;i < medium->nrflt;i++){
    v = 2.0*rsf->v0*exp(-x[4*i]/fault[i].mu_sa)*sinh(x[4*i+1]/(fault[i].mu_sa*x[4*i+2]));
    if(fabs(v) > vmax) vmax = fabs(v);
  }
  PetscCall(fd_full_put_read(X,&fd->scX,&fd->Xall,&x));
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
    if(fd->sc3){ PetscCall(VecScatterDestroy(&fd->sc3)); PetscCall(VecDestroy(&fd->v3all)); }
    if(fd->scX){ PetscCall(VecScatterDestroy(&fd->scX)); PetscCall(VecDestroy(&fd->Xall)); }
    if(fd->scV){ PetscCall(VecScatterDestroy(&fd->scV)); PetscCall(VecDestroy(&fd->Vall)); }
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
