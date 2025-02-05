#include "interact.h"
#include "properties.h"
/*
  
*/
#ifdef USE_PETSC
#include "petsc_prototypes.h"


int main(int argc,char **argv)
{

  TS ts; /* timestepping context */
  const PetscScalar *xprint;
  PetscReal *values = NULL;
  /* material parameters */
  const PetscReal shear_modulus_si = 32.04e9,s_wave_speed_si = 3.464e3;
  /*  */
  VecScatter ctx;
  Vec stress_rate,xout,x,F;
  FILE *out;
  struct med *medium;struct flt *fault;
  PetscInt m,i,n,j;
  PetscReal state,tau,sigma,tmp1,tmp2,tmp3,tmp4;
  struct interact_ctx par[1]; /* user-defined work context */
  char geom_file[STRLEN]="geom.in",rsf_file[STRLEN]="rsf.dat";
  PetscBool read_value,use_h = PETSC_TRUE;
  char *home_dir = getenv("HOME");char par_file[STRLEN];
  sprintf(par_file,"%s/progs/src/interact/petsc_settings.yaml",home_dir);
  /* set up structure */
  par->medium=(struct med *)calloc(1,sizeof(struct med));
  medium = par->medium;
  /*  */
  PetscFunctionBegin;
  PetscInitialize(&argc,&argv,par_file,NULL); /* read defaults */
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  if(medium->comm_size>1)
    if(medium->comm_size%3!=0){
      fprintf(stderr,"%s: number of cores %i not divisible by 3\n",argv[0],medium->comm_size);
      PetscCall(PetscFinalize());
 
    }
  /* options for this code */
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-use_h", &use_h,&read_value)); /* H matrix or dense */
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-rsf_file", rsf_file, STRLEN,&read_value));
  /* input parameter */
  /* get the geometry */
  read_geometry(geom_file,&medium,&par->fault,TRUE,FALSE,FALSE,FALSE);
  fault = par->fault;
  n = medium->nrflt;
  HEADNODE{
    if(read_value)
      fprintf(stderr,"%s: reading geometry from %s as set by -geom_file, %i patches\n",argv[0],geom_file,n);
    else
      fprintf(stderr,"%s: reading geometry from default, %s, %i patches\n",argv[0],geom_file,n);
  }
  /* 
     define general frictional paramaters
  */
  medium->f0  = 0.6; /* f_0 reference friction */
  medium->shear_mod_over_2cs_si = shear_modulus_si /(2.0*s_wave_speed_si);
  medium->dc = 0.02;
  medium->b = 0.02;
  medium->vpl = 1e-9;
  medium->v0 = medium->vpl;
  medium->vmin = 1e-20;
  /* 
     now, read in variations or assign based on geometry 
  */
  read_rsf(rsf_file,medium,fault);
  /* 
     create and calculate interaction matrices 
  */
  calc_petsc_Isn_matrices(medium,fault,use_h,shear_modulus_si/SHEAR_MODULUS,0,&medium->Is); /* shear stress */
  calc_petsc_Isn_matrices(medium,fault,use_h,shear_modulus_si/SHEAR_MODULUS,1,&medium->In); /* normal stress */

  
  /* make slip and shear and normal stressing vectors */
  PetscCall(MatCreateVecs(medium->Is,&par->slip_rate,&stress_rate)); /* For A x = b: x -> left, b -> right */
  /* 
     compute backslip stressing 
  */
  PetscCall(VecSet(par->slip_rate,-medium->vpl));
  PetscCall(VecAssemblyBegin(par->slip_rate));
  PetscCall(VecAssemblyEnd(par->slip_rate));
  //VecView(par->slip_rate,PETSC_VIEWER_STDOUT_WORLD);

  /* assign to stored vec */
  for(i=0;i < 2;i++){
    if(i==0)
      PetscCall(MatMult(medium->Is,par->slip_rate,stress_rate)); /* this is shear */
    else
      PetscCall(MatMult(medium->In,par->slip_rate,stress_rate)); /* this is normal */
    //VecView(stress_rate,PETSC_VIEWER_STDOUT_WORLD);
    /* assign to the faults */
    PetscCall(VecScatterCreateToAll(stress_rate,&ctx,&xout));
    PetscCall(VecScatterBegin(ctx,stress_rate,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(  ctx,stress_rate,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGetArray(xout,&values));
    for(j=0;j < medium->nrflt;j++) /* store the backslip loading rates */
      fault[j].sinc[i] = values[j];
    PetscCall(VecRestoreArray(xout,&values));

  }
  PetscCall(VecDestroy(&stress_rate));
  /* 
     timestepping 
  */
  medium->stop_time=1e5;			/* stop time */
  medium->slip_line_dt = 1;			/* for output, not actual timestep*/
  /* 
     solution vector 
  */
  m=3*n;			
  PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
  PetscCall(VecSetSizes(x,PETSC_DECIDE,m));
  PetscCall(VecSetFromOptions(x));
  /* initialize */
  for (i = medium->rs; i < medium->re; i++){
    sigma = 0.001 * shear_modulus_si;tau = medium->f0 * sigma;state = medium->dc/medium->vpl;
    PetscCall(VecSetValue(x, (PetscInt)(i*3+0),state, INSERT_VALUES)); /* steady state state = D_c/v */
    PetscCall(VecSetValue(x, (PetscInt)(i*3+1),tau, INSERT_VALUES)); /* shear stress */
    PetscCall(VecSetValue(x, (PetscInt)(i*3+2),sigma, INSERT_VALUES)); /* normal stress */

    fprintf(stderr,"patch %i a %g b %g dc %e v0 %e tau %e sigma %e state %e - vel %e\n",
	    i,fault[i].mu_s,fault[i].mu_d,medium->dc,medium->v0,tau,sigma,state,
	    vel_from_rsf(tau, sigma, state, fault[i].mu_s,medium->v0,&tmp1,&tmp2,&tmp3,&tmp4,medium));
  }
  PetscCall(VecAssemblyBegin(x));PetscCall(VecAssemblyEnd(x));
  //PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(VecDuplicate(x, &F));

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

  //PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));
  rsf_ODE_RHSFunction(ts,medium->time,x,F,par);
  //PetscCall(VecView(F,PETSC_VIEWER_STDOUT_WORLD));
  exit(-1);
  
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
  PetscCall(VecDestroy(&par->slip_rate));

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
/* 
   compute time-derivative for state vector, state, shear, normal stress for each n patch
   x = (ψ1,τ1,σ1,ψ2,τ2,σ2,ψ3,τ3,σ3,....) 3*n

   (see hmatrix_test/ode_solve_test.c for a working example)

   a = fault[].mu_s
   b = fault[].mu_d 
*/
PetscErrorCode rsf_ODE_RHSFunction(TS ts,PetscReal time,Vec X,Vec F,void *ptr)
{
  PetscScalar *f,cosh_fac,*sinh_fac,*mu,*scaled_tau,b_over_dc,two_v0,two_v0_div_a,pre_fac,loc_tau_dot ;
  const PetscScalar *x;
  PetscScalar dvdtau,dvdsigma,dvdstate,*loc_slip_rate,*exp_fac;
  struct med *medium;struct flt *fault;
  PetscInt i,j,k,ln;
  Vec tau_dot,sigma_dot;
  struct interact_ctx *par;
  PetscFunctionBeginUser;
  par = (struct interact_ctx *)ptr;medium = par->medium;fault = par->fault;
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
  /* parallel assembly of slip_rate vector */
  /* i global, j global*3, k local */
  for (i = medium->rs, j=i*3, k=0; i < medium->re; i++,j+=3,k++) { /* compute
								      prefactors
								      and
								      slip
								      rate
								      vector */
    //state = x[i*3  ];
    //tau =   x[i*3+1];
    //sigma = x[i*3+2];
    /* v = 2 v_0 exp(-psi/a) sinh(tau/sigma/a) */
    loc_slip_rate[k] = vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,medium->v0,
				    (mu+k), (scaled_tau+k), (exp_fac+k), (sinh_fac+k),medium);    
    PetscCall(VecSetValue(par->slip_rate, i, loc_slip_rate[k], INSERT_VALUES));
  }
  PetscCall(VecAssemblyBegin(par->slip_rate));
  PetscCall(VecAssemblyEnd(par->slip_rate));
  //VecView(par->slip_rate,PETSC_VIEWER_STDOUT_WORLD);
  
  /* make room for loading rate */
  PetscCall(VecDuplicate(tau_dot, &sigma_dot));
  /* forward solve effect of all slipping faults */
  PetscCall(MatMult(medium->Is, par->slip_rate, tau_dot));
  PetscCall(MatMult(medium->In, par->slip_rate, sigma_dot));
  /* 
     compute derivatives 
  */
  b_over_dc = medium->b/medium->dc; /* for now, leave b constant */
  two_v0 = 2.0*medium->v0;
  /* i global, j global*3, k local */
  for (i = medium->rs,j=i*3,k=0; i < medium->re; i++,j+=3,k++) {
    /* 2v0/a */
    two_v0_div_a = two_v0/fault[i].mu_s;
    /* 
       d state/ dt = b/dc (v0 exp((f0-state)/b)-v) 
    */
    f[j] = b_over_dc *(medium->v0 * PetscExpReal((medium->f0 - x[j])/medium->b) - loc_slip_rate[k]);
    /* 
       d sigma/dt  = 
    */
    PetscCall(VecGetValues(sigma_dot,1,&i,(f+j+2)));//f[j+2] = sigma_dot[i];
    f[j+2] +=  fault[i].sinc[1];		    /* background */
    /* 
       d tau/dt 
    */
    /* dv/dstate */
    dvdstate = -two_v0_div_a * exp_fac[k] * sinh_fac[k];
    /* dv/dtau */
    pre_fac =  two_v0_div_a/x[j+2];
    cosh_fac  = PetscCoshReal(scaled_tau[k]); /*  */
    cosh_fac *= exp_fac[k];
    dvdtau   =   pre_fac       * cosh_fac;
    /* dv/dsigma */
    dvdsigma =  -pre_fac*mu[k] * cosh_fac;
    /* d tau/dt */
    f[j+1]  = 1./(1.0+medium->shear_mod_over_2cs_si * dvdtau);
    PetscCall(VecGetValues(tau_dot,1,&i,&loc_tau_dot));
    f[j+1] *= (loc_tau_dot  + fault[i].sinc[0] - medium->shear_mod_over_2cs_si * (dvdsigma * f[j+2] + dvdstate * f[j]));
  }
  /*  */
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(F,&f));
  /*  */
  PetscCall(VecDestroy(&sigma_dot));
  PetscCall(VecDestroy(&tau_dot));
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
		       PetscReal *mu, PetscReal *scaled_tau, PetscReal *exp_fac, PetscReal *sinh_fac,
		       struct med *medium)
{
  PetscReal vel;
  *mu = tau/sigma;
  *scaled_tau = (*mu)/a;
  *sinh_fac = PetscSinhReal(*scaled_tau);
  *exp_fac = PetscExpReal(-state/a);
  //fprintf(stderr,"%g %g %g %g %g\n",state,*mu,*scaled_tau,*sinh_fac,*exp_fac);
  /* v = 2 v_0 exp(-psi/a) sinh(tau/sigma/a) */
  vel =  2.*v0 * (*sinh_fac) * (*exp_fac);
  if(vel < medium->vmin)
    vel = medium->vmin;
  return vel;
}

#endif
