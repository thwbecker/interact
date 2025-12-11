#include "../interact.h"
#define LHEADNODE if(par->medium->comm_rank==0)
/*

  solving ordinary differential equations originally based on
  https://www.cs.usask.ca/~spiteri/CMPT851/notes/PETSc.pdf

  this examples solves a 3-D rate state friction example
  
*/
#ifdef USE_PETSC

#include "petscts.h"
/*
  parameters needed by the derivative function, and the monitor
  function (could be separate, but why not (?))
*/
struct AppCtx{
  int n;
  PetscReal b1,b2,r,k,knd;
  struct med medium[1];
  /*  */
  PetscReal old_time;
  PetscReal dt_monitor,adx_monitor,rdx_monitor;
  FILE *fout_monitor;
  Vec Xold;
};
/*
User-defined routines
*/
static PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);
static PetscErrorCode myMonitor(TS , PetscInt , PetscReal , Vec , void *);
static PetscErrorCode init_monitor(void *, PetscReal , PetscReal ,  PetscReal , PetscReal , Vec );
static PetscErrorCode finalize_monitor(void *);
#endif

int main(int argc,char **argv)
{
#ifdef USE_PETSC
  TS ts; /* timestepping context */
  Vec x; /* solution, residual vectors */
  PetscReal t_final,time,t_init,atol,rtol;
  PetscReal kcr1,kcr2,adx_monitor,rdx_monitor,dt_monitor;
  PetscBool flag_set;
  TSAdapt        adapt;
  struct AppCtx par[1]; /* user-defined work context */

  PetscFunctionBegin;
  PetscInitialize(&argc,&argv,NULL,NULL);
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &par->medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &par->medium->comm_rank));
  /*  
      control parmater: non dimensional stiffness
  */
  PetscOptionsGetReal(NULL, NULL, "-knd", &par->knd, &flag_set);
  if(!flag_set){
    //par->knd = 0.88; % 2 orbit
    par->knd = 0.86;// % 4 
    //par->knd = 0.856; //% 8 
    //par->knd = 0.8552; //% 16
    //par->knd = 0.8525 ; // attractor
    LHEADNODE
      fprintf(stderr,"%s: using default nd stiffness of %g\n",argv[0],par->knd);
  }else{
    LHEADNODE
      fprintf(stderr,"%s: command line override nd stiffness of %g\n",argv[0],par->knd);
  }
  /*  */
  PetscOptionsGetReal(NULL, NULL, "-t_finel", &t_final, &flag_set);
  if(!flag_set)
    t_final=1e5;			/* stop time */
  /* 
     solver control 
  */
  atol = 1e-4;
  rtol = 1e-12;
  /* 
     output control 
  */
  dt_monitor = 5;	
  adx_monitor = 1e-2;		/* absolute */
  rdx_monitor = 1e-4;		/* relative */

  /* 
     parameters
  */
  /* dimensions */
  par->n = 3;
  /*  */
  par->b1=1.0;
  par->b2=0.84;
  par->r=0.048;
  /* non dim stiffnesses */
  kcr1 = par->b1 - 1.0;
  kcr2=(kcr1+par->r*(2.*par->b1+(par->b2-1.)*(2.+par->r))+
	sqrt(4.*par->r*par->r*(kcr1+par->b2)+pow(kcr1+par->r*par->r*(par->b2-1.),2)))/(2.+2.*par->r);
  par->k = par->knd * kcr2;	/* set actual stiffness */

  /*  */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
  PetscCall(VecSetSizes(x,PETSC_DECIDE,par->n));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecAssemblyBegin(x));
  PetscCall(VecAssemblyEnd(x));
  /* initial condition */
  PetscCall(VecSet(x, 0.001));
  fprintf(stderr,"%s: initializing vector with\n",argv[0]);
  PetscCall(VecView(x, PETSC_VIEWER_STDERR_SELF));
  
  /*
    Create timestepper context
  */
  PetscCall(TSCreate(PETSC_COMM_WORLD,&ts));
  PetscCall(TSSetProblemType(ts,TS_NONLINEAR)); /* can be linear? */
  PetscCall(TSSetSolution(ts,x));
  PetscCall(TSSetRHSFunction(ts, NULL, RHSFunction,&par));
  /*  */
  t_init = 0;
  
  /* init my control environment for output */
  init_monitor(par,dt_monitor,adx_monitor,rdx_monitor,t_init,x);
  /*
    use runge kutta
  */
  //PetscCall(TSSetType(ts,TSRK));
  PetscCall(TSSetType(ts,TSARKIMEX));
  /*  */
  PetscCall(TSSetTolerances(ts, atol, NULL, rtol, NULL));
  PetscCall(TSGetAdapt(ts,&adapt));
  //PetscCall(TSAdaptSetStepLimits(adapt,1e-15, dt_monitor));
  /*
    Set the initial time and the initial timestep given above.
  */
  PetscCall(TSSetTime(ts,t_init));	/* initial time */
  PetscCall(TSSetTimeStep(ts,1e-2)); /* initial timestep */
  PetscCall(TSSetFromOptions(ts));
  
  /* allow for rough final time (else might have hard time finding
     exaxt state */
  PetscCall(TSSetExactFinalTime(ts,  TS_EXACTFINALTIME_STEPOVER ));

  /* set the monitor and output function */
  PetscCall(TSMonitorSet(ts, myMonitor, (void *)par, NULL));
  
  PetscCall(TSSetMaxTime(ts,t_final));
  PetscCall(TSSolve(ts, x));
  PetscCall(TSGetTime(ts,&time)); /* get the current time */
  LHEADNODE
    fprintf(stderr,"done at t %g x ",time);
  finalize_monitor(par);
  /*View information about the time-stepping method and the solutionat the end time.*/
  PetscCall(TSView(ts, PETSC_VIEWER_STDOUT_SELF));

  /*Free the data structures constructed above*/
  PetscCall(VecDestroy(&x));
  PetscCall(TSDestroy(&ts));
  PetscCall(PetscFinalize());
  exit(0); 
#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
  exit(-1); 
#endif

}

#ifdef USE_PETSC

/* set output checks */
static PetscErrorCode init_monitor(void *ctx, PetscReal dt_monitor, PetscReal adx_monitor,
				   PetscReal rdx_monitor,PetscReal t_init, Vec X0)
{
  struct AppCtx *par;
  char fname[PETSC_MAX_PATH_LEN];
  PetscFunctionBeginUser;
  par = (struct AppCtx *)ctx;
  par->dt_monitor = dt_monitor;			
  par->adx_monitor = adx_monitor;
  par->rdx_monitor = rdx_monitor;
  par->old_time = t_init;
  
  LHEADNODE{
    snprintf(fname,PETSC_MAX_PATH_LEN,"state.%020.16f.dat",par->knd);
    par->fout_monitor = fopen(fname,"w");
    fprintf(stderr,"writing state to %s\n",fname);
  }
  PetscCall(VecDuplicate(X0, &par->Xold));
  PetscCall(VecCopy(X0,par->Xold));
  
  PetscFunctionReturn(PETSC_SUCCESS);
}
static PetscErrorCode finalize_monitor(void *ctx)
{
  struct AppCtx *par;
  PetscFunctionBeginUser;
  par = (struct AppCtx *)ctx;
  LHEADNODE{
    fclose(par->fout_monitor);	/* what else? */
  }
  PetscCall(VecDestroy(&par->Xold));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode myMonitor(TS ts, PetscInt step, PetscReal t, Vec X, void *ctx)
{
  Vec DX;
  PetscReal dx_norm,x_norm;
  int i;
  PetscBool rdx_bail = PETSC_FALSE;
  struct AppCtx *par;
  const PetscScalar *x;
  PetscFunctionBeginUser;
  par = (struct AppCtx *)ctx;
  if (step < 0)
    PetscFunctionReturn(PETSC_SUCCESS); /* negative one is used to indicate an interpolated solution */
  PetscCall(VecDuplicate(X, &DX));	/* make room (do it every time?) */

  PetscCall(VecWAXPY(DX, -1.0, X, par->Xold));
  PetscCall(VecNorm(X, NORM_2, &x_norm));
  PetscCall(VecNorm(DX, NORM_2, &dx_norm));
  if(x_norm > 1e-15){
    if(dx_norm/x_norm > par->rdx_monitor)
      rdx_bail = PETSC_TRUE;
  }
  if(rdx_bail || (dx_norm > par->adx_monitor) || (fabs(t-par->old_time) > par->dt_monitor)){
    PetscCall(VecGetArrayRead(X,&x));
    LHEADNODE{
      fprintf(par->fout_monitor,"%20.15e\t ",t);
      for(i=0;i<par->n;i++)
	fprintf(par->fout_monitor,"%20.10e ",x[i]);
      fprintf(par->fout_monitor,"\n");
    }
    PetscCall(VecRestoreArrayRead(X,&x));
  }
  /* store old state */
  par->old_time = t;
  PetscCall(VecCopy(X,par->Xold));
  /* free space */
  PetscCall(VecDestroy(&DX));
  
  
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RHSFunction(TS ts,PetscReal time,Vec X,
				  Vec F,void *ptr)
{
  PetscScalar *f,expx;
  const PetscScalar *x;
  struct AppCtx *par;
  PetscFunctionBeginUser;
  par = (struct AppCtx *)ptr;

  PetscCall(VecGetArrayRead(X,&x));
  PetscCall(VecGetArray(F,&f));
  
  expx = PetscExpReal(x[0]);
  /* 
     rate state friction with two state variables from Becker (2000) 
  */

  f[1] =  (1.0 - expx) * par->k;
  f[2] = -expx * par->r * (par->b2 * x[0] + x[2]);
  f[0] = expx * ((par->b1 - 1.0) * x[0] + x[1] - x[2]) + f[1] -  f[2];

  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(F,&f));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
