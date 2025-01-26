#include "../interact.h"

/*

  solving ordinary differential equations originally based on
  https://www.cs.usask.ca/~spiteri/CMPT851/notes/PETSc.pdf

  this examples solves a 3-D rate state friction example
  
*/
#ifdef USE_PETSC
#include "petscts.h"
/*
  parameters needed by the derivative function
*/
struct AppCtx{
  PetscReal b1,b2,r,k;
  struct med medium[1];
} ;
/*
User-defined routines
*/
static PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);

#endif

int main(int argc,char **argv)
{
#ifdef USE_PETSC
  TS ts; /* timestepping context */
  Vec x; /* solution, residual vectors */
  PetscReal end_time,dt_out,time;
  PetscReal kcr1,kcr2;
  const PetscScalar *xprint;
  FILE *out;
  PetscInt m,i;
  struct AppCtx par[1]; /* user-defined work context */

  PetscFunctionBegin;
  PetscInitialize(&argc,&argv,NULL,NULL);
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &par->medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &par->medium->comm_rank));

  end_time=1e5;			/* stop time */
  dt_out = 1;			/* for output, not actual timestep*/

  /* 
     parameters
  */
  par->b1=1.0;
  par->b2=0.84;
  par->r=0.048;
  /* non dim stiffnesses */
  kcr1 = par->b1 - 1.0;
  kcr2=(kcr1+par->r*(2.*par->b1+(par->b2-1.)*(2.+par->r))+
	sqrt(4.*par->r*par->r*(kcr1+par->b2)+pow(kcr1+par->r*par->r*(par->b2-1.),2)))/(2.+2.*par->r);
  //par->k = 0.9 * kcr2; % 2 orbit
  //par->k = 0.86 * kcr2; % 4 
  //par->k = 0.856 * kcr2; % 8 
  //par->k = 0.8552 * kcr2; % 16
  par->k = 0.8525 * kcr2; // attractor

  //par->k = 0.86 * kcr2; // 4 orbit
  
  
  m=3;				/* 3-D vector */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
  PetscCall(VecSetSizes(x,PETSC_DECIDE,m));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecAssemblyBegin(x));
  PetscCall(VecAssemblyEnd(x));
  /* initial condition */
  PetscCall(VecSet(x, 0.001));
  PetscCall(VecView(x, PETSC_VIEWER_STDOUT_SELF));
  
  /*
    Create timestepper context
  */
  PetscCall(TSCreate(PETSC_COMM_WORLD,&ts));
  PetscCall(TSSetProblemType(ts,TS_NONLINEAR)); /* can be linear? */
  PetscCall(TSSetSolution(ts,x));
  PetscCall(TSSetRHSFunction(ts, NULL, RHSFunction,&par));
  /*
    use runge kutta
  */
  PetscCall(TSSetType(ts,TSRK));
  /*
    Set the initial time and the initial timestep given above.
  */
  PetscCall(TSSetTime(ts,0.0));	/* initial time */
  PetscCall(TSSetTimeStep(ts,0.01)); /* initial timestep */
  
  PetscCall(TSSetFromOptions(ts));
  /* allow for rough final time (else might have hard time finding
     exaxt state */
  PetscCall(TSSetExactFinalTime(ts,  TS_EXACTFINALTIME_STEPOVER ));
  out = fopen("tmp.dat","w");
  while(time < end_time){	/* advance until next output */
    PetscCall(TSSetMaxTime(ts,time+dt_out));
    PetscCall(TSSolve(ts, x));
    PetscCall(TSGetTime(ts,&time)); /* get the current time */
    /* prep for output */
    PetscCall(VecGetArrayRead(x,&xprint));
    fprintf(out,"%g\t",time);
    for(i=0;i<m;i++)
      fprintf(out,"%20.10e ",xprint[i]);
    fprintf(out,"\n");
    PetscCall(VecRestoreArrayRead(x,&xprint));
  }
  fclose(out);
  
  /*View information about the time-stepping method and the solutionat the end time.*/
  PetscCall(TSView(ts, PETSC_VIEWER_STDOUT_SELF));

  /*Free the data structures constructed above*/
  PetscCall(VecDestroy(&x));
  PetscCall(TSDestroy(&ts));
  PetscCall(PetscFinalize());
#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);
}

#ifdef USE_PETSC
static PetscErrorCode RHSFunction(TS ts,PetscReal time,Vec X,Vec F,void *ptr)
{
  PetscScalar *f,expx;
  const PetscScalar *x;
  struct AppCtx *par;
  PetscFunctionBeginUser;
  par = (struct AppCtx *)ptr;

  PetscCall(VecGetArrayRead(X,&x));
  PetscCall(VecGetArray(F,&f));
  
  expx = PetscExpReal(x[0]);
  /* rate state friction with two state variables from Becker (2000) */

  f[1] =  (1.0 - expx) * par->k;
  f[2] = -expx * par->r * (par->b2 * x[0] + x[2]);
  f[0] = expx * ((par->b1 - 1.0) * x[0] + x[1] - x[2]) + f[1] -  f[2];

  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(F,&f));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
