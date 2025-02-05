#include "interact.h"

/*
  
*/
#ifdef USE_PETSC

#include "petsc_prototypes.h"
/*
User-defined routines
*/


#endif

int main(int argc,char **argv)
{
#ifdef USE_PETSC
  TS ts; /* timestepping context */
  const PetscScalar *xprint,*values=NULL;
  const PetscReal shear_modulus_si = 32.04e9,s_wave_speed_si = 3.464e3;
  VecScatter ctx,xout;
  Vec slip_rate,stress_rate;
  FILE *out;
  struct med *medium;struct flt *fault;
  PetscInt m,i,n,lm,ln,j;
  struct interact_ctx par[1]; /* user-defined work context */
  char geom_file[STRLEN]="geom.in",rsf_file[STRLEN]="rsf.dat";
  PetscBool read_value,use_h = PETSC_TRUE;
  /*  */
  par->medium=(struct med *)calloc(1,sizeof(struct med));
  fault = par->fault;
  medium = par->medium;
  /*  */
  PetscFunctionBegin;
  PetscInitialize(&argc,&argv,NULL,NULL);
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  if(medium->comm_size>1)
    if(medium->comm_size%3!=0){
      fprintf(stderr,"%s: number of cores %i not divisible by 3\n",argv[0],medium->comm_size);
      PetscCall(PetscFinalize());
 
    }
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-rsf_file", rsf_file, STRLEN,&read_value));
  HEADNODE{
    if(read_value)
      fprintf(stderr,"%s: reading geometry from %s as set by -geom_file\n",argv[0],geom_file);
    else
      fprintf(stderr,"%s: reading geometry from default, %s\n",argv[0],geom_file);
  }
  /* get the geometry */
  read_geometry(geom_file,&medium,&fault,TRUE,FALSE,FALSE,FALSE);
  n = medium->nrflt;
  /* 
     define general parmaters 
  */
  medium->f0  = 0.6; /* f_0 reference friction */
  medium->shear_mod_over_2cs_si = shear_modulus_si /(2.0*s_wave_speed_si);
  medium->dc = 0.02;
  medium->b = 0.02;
  medium->vpl = 1e-9;
  medium->v0 = medium->vpl;
  /* now, read in variations or assign based on geometry */
  read_rsf(rsf_file,&medium,&par->fault,TRUE,FALSE,FALSE,FALSE);
  /* 
     create and calculate interaction matrices 
  */
  calc_petsc_Isn_matrices(medium,par->fault,use_h,shear_modulus_si/SHEAR_MODULUS);
  /* make sure we have the ownership rante */
  PetscCall(MatGetOwnershipRange(medium->Is, &medium->rs, &medium->re));
  medium->rn = medium->re  - medium->rs; /* number of local elements */
  fprintf(stderr,"%s: core %03i computing for patches %06i to %06i\n",argv[0],medium->comm_rank,medium->rs, medium->re);
  /* make slip and shear and normal stressing vectors */
  PetscCall(MatCreateVecs(medium->Is, &stress_rate, &slip_rate)); /* For A x = b: x -> left, b -> right */
  /* 
     compute backslip stressing 
  */
  PetscCall(VecSet(slip_rate,-medium->vpl));
  PetscCall(VecAssemblyBegin(slip_rate));PetscCall(VecAssemblyEnd(slip_rate));
  /* assign to stored vec */
  for(i=0;i<2;i++){
    if(i==0)
      PetscCall(MatMult(medium->Is,slip_rate,stress_rate)); /* this is shear */
    else
      PetscCall(MatMult(medium->In,slip_rate,stress_rate)); /* this is normal */
    /* get to the faults, in case */
    PetscCall(VecScatterCreateToZero(stress_rate,&ctx,&xout));
    PetscCall(VecScatterBegin(ctx,stress_rate,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx,stress_rate,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGetArray(xout,&values));
    for(j=0;j < medium->nrflt;j++) /* store the backslip loading rates */
      fault[j].sinc[i] = values[j];
    PetscCall(VecRestoreArray(xout,&values));
  }
  /* timestepping */
  medium->stop_time=1e5;			/* stop time */
  medium->slip_line_dt = 1;			/* for output, not actual timestep*/

  /* solution vector */
  m=3*n;			
  PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
  PetscCall(VecSetSizes(x,PETSC_DECIDE,m));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecAssemblyBegin(x));
  PetscCall(VecAssemblyEnd(x));
  /* initial condition */
  PetscCall(VecSet(x, 0.0001));
  
  /*
    Create timestepper context
  */
  PetscCall(TSCreate(PETSC_COMM_WORLD,&ts));
  PetscCall(TSSetProblemType(ts,TS_NONLINEAR)); /* can be linear? */
  PetscCall(TSSetSolution(ts,x));
  PetscCall(TSSetRHSFunction(ts, NULL, rsf_ODE_RHSFunction,&par));
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
  while(medium->time < medium->stop_time){	/* advance until next output */
    PetscCall(TSSetMaxTime(ts,(medium->time + medium->slip_line_dt)));
    PetscCall(TSSolve(ts, x));
    PetscCall(TSGetTime(ts,&(medium->time))); /* get the current time */
    /* prep for output */
    PetscCall(VecGetArrayRead(x,&xprint));
    fprintf(out,"%g\t",medium->time);
    fprintf(out,"\n");
    PetscCall(VecRestoreArrayRead(x,&xprint));
  }
  fclose(out);
  
  /*View information about the time-stepping method and the solutionat the end time.*/
  PetscCall(TSView(ts, PETSC_VIEWER_STDOUT_SELF));

  PetscCall(MatDestroy(&medium->Is));
  PetscCall(MatDestroy(&medium->In));
  PetscCall(VecDestroy(&x));
  /*  */
  PetscCall(VecDestroy(&slip_rate));
  PetscCall(VecDestroy(&stress_rate));
  /*  */
  PetscCall(VecDestroy(&xout));
  PetscCall(VecScatterDestroy(&ctx));
  /*  */
  PetscCall(TSDestroy(&ts)); 
  PetscCall(PetscFinalize());

#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);
}

#ifdef USE_PETSC
/* compute time-derivative for state vector, state, shear, normal stress for each n patch
   x = (ψ1,τ1,σ1,ψ2,τ2,σ2,ψ3,τ3,σ3,....) 3*n
*/
static PetscErrorCode rsf_ODE_RHSFunction(TS ts,PetscReal time,Vec X,Vec F,void *ptr)
{
  PetscScalar *f,cosh_fac,*two_v0,*sinh_fac,*mu,*scaled_tau,b_over_dc,two_v0;
  PetscScalar dvdtau,dvdsigma,dvdstate,loc_slip_rate;
  const PetscScalar *x;struct med *medium;struct flt *fault;
  PetscInt n,i,ln,j,k;
  Vec tau_dot,sigma_dot;
  struct interact_ctx *par;
  PetscFunctionBeginUser;
  par = (struct interact_ctx *)ptr;medium = par->medium;fault = par->fault;
  n = medium->nrflt;
  /* unpack */
  PetscCall(VecGetLocalSize(X, &ln));
  PetscCall(VecGetArrayRead(X,&x));
  PetscCall(VecGetArray(F,&f));
  /* make room for factors */
  mu =            (PetscScalar *)malloc(sizeof(PetscScalar)*medium->rn);
  scaled_tau =    (PetscScalar *)malloc(sizeof(PetscScalar)*medium->rn);
  loc_slip_rate = (PetscScalar *)malloc(sizeof(PetscScalar)*medium->rn);
  exp_fac =       (PetscScalar *)malloc(sizeof(PetscScalar)*medium->rn);
  sinh_fac =      (PetscScalar *)malloc(sizeof(PetscScalar)*medium->rn);

  PetscCall(MatCreateVecs(medium->Is, &tau_dot, &par->slip_rate)); 
  /* i global, j global*3, k local */
  for (i = medium->rs, j=i*3, k=0; i < medium->re; i++,j+=3,k++) { /* compute
								      prefactors
								      and
								      slip
								      rate
								      vector */
    //state = x[j];
    //tau =   x[j+1];
    //sigma = x[j+2];
    /* v = 2 v_0 exp(-psi/a) sinh(tau/sigma/a) */
    loc_slip_rate[k] = vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].a,medium->v0,
				    (mu+k), (scaled_tau+k), (exp_fac+k), (sinh_fac+k));    
    PetscCall(VecSetValue(par->slip_rate, i, loc_slip_rate[k], INSERT_VALUES));
  }
  PetscCall(VecAssemblyBegin(par->slip_rate));PetscCall(VecAssemblyEnd(par->slip_rate));
  
  /* make room for loading rate */
  PetscCall(VecDuplicate(tau_dot, &sigma_dot));
  /* forward solve effect of all slipping faults */
  PetscCall(MatMult(medium->Is, par->slip_rate, tau_dot));
  PetscCall(MatMult(medium->In, par->slip_rate, sigma_dot));

  /* compute derivatives */
  b_over_dc = medium->b/medium->dc;
  two_v0 = 2.0*medium->v0;
  /* i global, j global*3, k local */
  for (i = medium->rs,j=i*3,k=0; i < medium->re; i++,j+=3,k++) {
    two_v0_div_a = two_v0/fault[i].a;
    /* d state/ dt 
     */
    f[j] = b_over_dc *(medium->v0 * PetscExp((medium->f0 - x[j])/medium->b) - loc_slip_rate[k]);
    /* d sigma/dt */
    PetscCall(VecGetValues(sigma_dot,1,&i,(f+j+2)));//f[j+2] = sigma_dot[i];
    f[j+2] +=  fault[i].sinc[1];		    /* backgroudn */
    /* d tau/dt */
    /* dv/dstate */
    dvdstate = -two_v0/fault[i].a * exp_fac[k]*sinh_fac[k];
    /* dv/dtau */
    pre_fac =  two_v0_div_a/x[j+2];
    cosh_fac  = PetscCosh(scaled_tau[k]); /*  */
    cosh_fac *= exp_fac[k];
    dvdtau   =   pre_fac       * cosh_fac;
    /* dv/dsigma */
    dvdsigma =  -pre_fac*mu[k] * cosh_fac;
    /* d tau/dt */
    f[j+1]  = 1./(1.0+medium->shear_mod_over_2cs * dvdtau);
    PetscCall(VecGetValues(tau_dot,1,&i,&loc_tau_dot));
    f[j+1] *= (loc_tau_dot  + fault[i].sinc[0] - medium->shear_mod_over_2cs * (dvdsigma * f[j+2] + dvdstate*f[j]));
  }
  /*  */
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(F,&f));
  /*  */
  PetscCall(VecDestroy(sigma_dot));
  PetscCall(VecDestroy(tau_dot));
  /*  */
  free(mu);free(scaled_tau);free(loc_slip_rate);free(exp_fac);free(sinh_fac);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* 
   compute the sliding velocity and a bunch of helper factor to be
   reused 
   exp_fac = exp(-state/a) 
   scaled_tau = tau/mu/a
 */
PetscReal vel_from_rsf(PetscReal tau, PetscReal sigma, PetscReal state, PetscReal a,PetscReal v0,
		       PetscReal *mu, PetscReal *scaled_tau, PetscReal *exp_fac, PetscReal *sinh_fac)
{
  *mu = tau/sigma;
  *scaled_tau = (*mu)/a;
  *exp_fac = PetscExp(-state/a);
  *sinh_fac = PetscSinh(*scaled_tau);
  /* v = 2 v_0 exp(-psi/a) sinh(tau/sigma/a) */
  return 2.*v0 * (*sinh_fac) * (*exp_fac);
}

#endif
