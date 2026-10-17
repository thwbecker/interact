/*
  fully dynamic extension of rsf_solve

  Milestone M1 skeleton (2026-10-17): option parsing, allocation of
  rsf->fd, a monitor that logs when the maximum slip rate exceeds the
  switch threshold, and (mode -1) an exercise of the hand-over
  mechanism, TSRestartStep from inside a monitor without changing the
  state. The dynamic marcher (tables, history, banded friction solve)
  arrives with M2 to M4; see fully_dynamic/fd_c_design.md.

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
*/
#include "interact.h"
#include "properties.h"
#include "petsc_prototypes.h"
#include "rsf.h"
#include "fd.h"

PetscErrorCode rsf_fd_setup(struct interact_ctx *par)
{
  struct med *medium = par->medium;
  struct rsf_vars *rsf = medium->rsf;
  struct fd_vars *fd;
  PetscInt mode = 0;
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
  PetscCall(PetscCalloc1(1,&fd));
  fd->mode = mode;
  fd->vsw = 1e-3; fd->vback = 1e-4; fd->tol = 1e-2; fd->et = 1.0; fd->xi = 0.25; fd->kmax = 64;
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_vsw",&fd->vsw,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_vback",&fd->vback,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_tol",&fd->tol,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_et",&fd->et,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-fd_xi",&fd->xi,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-fd_kmax",&fd->kmax,NULL));
  PetscCheck(fd->vback < fd->vsw,PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"-fd_vback must be below -fd_vsw");
  PetscCheck(fd->kmax >= 1 && (fd->kmax & (fd->kmax - 1)) == 0,PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"-fd_kmax must be a power of two");
  fd->nphase = 0; fd->nflag = 0; fd->t_last_flag = -1.0; fd->restart_done = PETSC_FALSE;
  rsf->fd = fd;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
		"rsf_fd_setup: fully dynamic extension, mode %d (M1 skeleton: logging only), vsw %.1e vback %.1e m/s, tol %.1e, et %.2f, xi %.2f, kmax %d, fft %s\n",
		(int)fd->mode,(double)fd->vsw,(double)fd->vback,(double)fd->tol,(double)fd->et,(double)fd->xi,(int)fd->kmax,fd_fft_backend()));
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
  PetscFunctionBeginUser;
  if(!rsf->fd)
    PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"rsf_fd_finalize: %ld accepted steps above vsw, %ld dynamic phases\n",
		(long)rsf->fd->nflag,(long)rsf->fd->nphase));
  PetscCall(PetscFree(rsf->fd));
  rsf->fd = NULL;
  PetscFunctionReturn(PETSC_SUCCESS);
}
