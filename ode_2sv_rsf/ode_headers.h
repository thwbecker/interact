/*

  solving ordinary differential equations originally based on
  https://www.cs.usask.ca/~spiteri/CMPT851/notes/PETSc.pdf

  this examples solves a 3-D rate state friction example, see Becker
  (2000)

  this holds RHS versions

*/
#include "interact.h"

#define LHEADNODE if(par->medium->comm_rank==0)
#ifdef USE_PETSC
#include "petscts.h"



#define BASIN_NSTATE 3
#define BASIN_NFULL 12		/* state + 3x3 tangent */
#define BASIN_MAX_CLUSTER 4096

/*
  parameters needed by the derivative function, and the monitor/event
  functions
*/
struct AppCtx{
  int n,imode;
  PetscReal b1,b2,r,k,knd;
  struct med medium[1];
  PetscReal old_time,event_tmin,t_final,monitor_tmin,dt_monitor,adx_monitor,rdx_monitor;
  PetscBool track_events,log_state;
  FILE *fout_monitor,*fout_event;
  char fname_monitor[PETSC_MAX_PATH_LEN],fname_event[PETSC_MAX_PATH_LEN];
  Vec Xold;
  /* Lyapunov bookkeeping */
  PetscBool do_lyap;
  PetscReal logsum[3];
  PetscReal last_renorm_t,dt_renorm,renorm_max;
  /* events (Poincare section x=0), detected in PostStep between
     accepted steps to avoid TSEvent (unstable for this stiff system
     in PETSc <= 3.20) */
  PetscInt nevent,mevent;
  PetscReal *ey,*ez;
  PetscInt *edir;		/* +1 upward, -1 downward crossing */
  PetscBool have_prev;
  PetscReal prev_t,prev_x[3],prev_f[3];
  PetscBool overflowed;		/* domain errors encountered */
};

/*
User-defined routines
*/
PetscErrorCode RHSFunction3D(TS,PetscReal,Vec,Vec,void*);
PetscErrorCode RHSFunction4D(TS,PetscReal,Vec,Vec,void*);
PetscErrorCode RHSFunction3DTangent(TS ,PetscReal ,Vec ,Vec ,void *);
void eval_rhs(const struct AppCtx *,const PetscReal [3],PetscReal [3]);
PetscErrorCode  init_stiff_par(struct AppCtx *, PetscReal );
PetscErrorCode unstable_plane(const struct AppCtx *,PetscReal [3],PetscReal [3],
			      PetscReal *,PetscReal *);

#endif
