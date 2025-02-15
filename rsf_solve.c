#include "interact.h"
#include "properties.h"
/*
  
*/
#ifdef USE_PETSC
#include "petsc_prototypes.h"

#define INT_RSF_DIM 4		/* ψ1,τ1,σ1,u(slip) */

int main(int argc,char **argv)
{

  TS ts; /* timestepping context */
  PetscReal *values = NULL;
  /* material parameters */
  const PetscReal shear_modulus_si = 32.04e9,s_wave_speed_si = 3.464e3, sec_per_year = 365.25*24.*60.*60.;
  /*  */
  VecScatter ctx;
  Vec stress_rate,xout,x,F,islip_rate_vec;
  FILE *fout1,*fout2;
  struct med *medium;struct flt *fault;
  PetscInt m,i,n,j,field_out,ierr;
  PetscRandom prand;
  PetscReal state,slip,sum[3],patch_l,stable_l,rand_fac,dummy[3],vmean,vstd,vmin,vmax,smean;
  struct interact_ctx par[1]; /* user-defined work context */
  char geom_file[STRLEN]="geom.in",rsf_file[STRLEN]="rsf.dat";
  PetscBool read_value,use_h = PETSC_TRUE,warned = PETSC_FALSE;
  char *home_dir = getenv("HOME");char par_file[STRLEN],vel_file[STRLEN];
  snprintf(par_file,STRLEN,"%s/progs/src/interact/petsc_settings.yaml",home_dir);
  /* set up structure */
  par->medium=(struct med *)calloc(1,sizeof(struct med));
  medium = par->medium;
  /* 
     define general frictional paramaters
  */
  medium->f0  = 0.6; /* f_0 reference friction */
  medium->shear_mod_over_2cs_si = shear_modulus_si /(2.0*s_wave_speed_si); /* G/(2v_s) */
  medium->dc = 0.02;		/* D_c */
  medium->vpl = 1.5e-9;		/* plate motion */
  medium->v0 = medium->vpl;	/* reference speed */
  /* 
     timestepping 
  */
  medium->time = medium->slip_line_time = 0.;
  medium->stop_time = 5000*sec_per_year;			/* stop time */
  /* output */
  medium->print_interval = 0.01*sec_per_year;			/* for average property output */
  medium->slip_line_dt = sec_per_year;			/* for output of velocity grid */
  
  /* start up MPI */
  PetscFunctionBegin;
  PetscInitialize(&argc,&argv,par_file,NULL); /* read defaults */
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  if(medium->comm_size > 1){
    if((medium->comm_size%INT_RSF_DIM) != 0){
      fprintf(stderr,"%s: number of cores %i not divisible by %i\n",argv[0],medium->comm_size,INT_RSF_DIM);
      PetscCall(PetscFinalize());
 
    }
  }
  //fprintf(stderr,"%s: running on %i cores, I am core %i\n",argv[0],medium->comm_size,medium->comm_rank);

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
      fprintf(stderr,"%s: reading geometry from %s as set by -geom_file, %i patches\n",
	      argv[0],geom_file,n);
    else
      fprintf(stderr,"%s: reading geometry from default, %s, %i patches\n",
	      argv[0],geom_file,n);
  }
 
  /* 
     now, read in variations or assign based on geometry 
  */
  read_rsf(rsf_file,medium,fault);
  /* 
     create and calculate interaction matrices 
  */
  calc_petsc_Isn_matrices(medium,fault,use_h,shear_modulus_si/SHEAR_MODULUS,0,&medium->Is); /* shear stress */
  calc_petsc_Isn_matrices(medium,fault,use_h,shear_modulus_si/SHEAR_MODULUS,1,&medium->In); /* normal stress */
  if(medium->use_h)
    PetscCall(MatView(medium->Is,PETSC_VIEWER_STDOUT_WORLD));
  
  /* make slip and shear and normal stressing vectors */
  PetscCall(MatCreateVecs(medium->Is,&islip_rate_vec,&stress_rate)); /* For A x = b: x -> left, b -> right */
  /* 
     compute backslip stressing 
  */
  PetscCall(VecSet(islip_rate_vec,-medium->vpl));
  PetscCall(VecAssemblyBegin(islip_rate_vec));
  PetscCall(VecAssemblyEnd(islip_rate_vec));
  //VecView(islip_rate,PETSC_VIEWER_STDOUT_WORLD);

  /* assign to stored vec */
  for(i=0;i < 2;i++){
    if(i==0)
      PetscCall(MatMult(medium->Is,islip_rate_vec,stress_rate)); /* this is shear */
    else
      PetscCall(MatMult(medium->In,islip_rate_vec,stress_rate)); /* this is normal */
    //VecView(stress_rate,PETSC_VIEWER_STDOUT_WORLD);
    /* assign to the faults for every node */
    PetscCall(VecScatterCreateToAll(stress_rate,&ctx,&xout));
    PetscCall(VecScatterBegin(ctx,stress_rate,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(  ctx,stress_rate,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecGetArray(xout,&values));
    for(j=0;j < medium->nrflt;j++) /* store the backslip loading rates */
      fault[j].sinc[i] = values[j];
    PetscCall(VecRestoreArray(xout,&values));
    PetscCall(VecScatterDestroy(&ctx));
  }
  PetscCall(VecDestroy(&stress_rate));
  PetscCall(VecDestroy(&islip_rate_vec));

  /* 
     solution vector 
  */
  m=INT_RSF_DIM*n;
  /*  */
  PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
  PetscCall(VecSetSizes(x,PETSC_DECIDE,m));
  PetscCall(VecSetFromOptions(x));
  /* 
     initialize x vector
  */
  /* random numbers */
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&prand));
  PetscCall(PetscRandomSetType(prand,PETSCRAND48));
  PetscCall(PetscRandomSetInterval(prand,0.95,1.05)); 
  for (i = medium->rs,j=i*INT_RSF_DIM; i < medium->re; i++,j+=INT_RSF_DIM){
    PetscCall(PetscRandomGetValue(prand,&rand_fac));
    fault[i].s[NORMAL] = 0.02 * shear_modulus_si; /* pick normal stress */
    fault[i].s[STRIKE] = medium->f0 * fault[i].s[NORMAL];
    /* define state such that fault starts at v_p */
    fault[i].u[0] = medium->vpl;
    /* have fault move at plate rate v = vp, compute state from that */
    state = -fault[i].mu_s*log(fault[i].u[0]/(2.0*medium->v0)/sinh(fault[i].s[STRIKE]/fault[i].s[NORMAL]/fault[i].mu_s))* rand_fac ;
    PetscCall(VecSetValue(x, (PetscInt)(j+0),state, INSERT_VALUES)); /* steady state state = D_c/v */
    PetscCall(VecSetValue(x, (PetscInt)(j+1),fault[i].s[STRIKE], INSERT_VALUES)); /* shear stress */
    PetscCall(VecSetValue(x, (PetscInt)(j+2),fault[i].s[NORMAL], INSERT_VALUES)); /* normal stress */
    dummy[0] = 0.;
    PetscCall(VecSetValue(x, (PetscInt)(j+3),dummy[0], INSERT_VALUES)); /* total slip */
    /* check dimensions */
    patch_l = sqrt(fault[i].l*fault[i].w*4.);
    /* L crierion */
    stable_l = shear_modulus_si*medium->dc/fault[i].mu_d/fault[i].s[NORMAL];
    if(patch_l>stable_l){			/* check */
      if(!warned){
	fprintf(stderr,"patch %05i a %7.3f b %7.3f a-b %9.5f dc %8.2e v0 %8.2e tau %8.3e sigma %8.3e state %8.3e - vel %8.3e - Lp %.3e Lps %.3e - %s\n",
		i,fault[i].mu_s,fault[i].mu_d,fault[i].mu_s-fault[i].mu_d,medium->dc,medium->v0,
		fault[i].s[STRIKE],fault[i].s[NORMAL],state,
		vel_from_rsf(fault[i].s[STRIKE], fault[i].s[NORMAL], state,
			     fault[i].mu_s,medium->v0,dummy,(dummy+1),(dummy+2),medium),
		patch_l,stable_l,(patch_l>stable_l)?("exceeds"):("OK"));
	fprintf(stderr,"suppressing further warnings\n");
      }
      warned = PETSC_TRUE;

    }
  }
  PetscCall(VecAssemblyBegin(x));PetscCall(VecAssemblyEnd(x));
  //PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));
  //  exit(-1);
  PetscCall(VecDuplicate(x, &F)); /* make room for F */
  PetscCall(PetscRandomDestroy(&prand));
  /*
    Create timestepper context
  */
  PetscCall(TSCreate(PETSC_COMM_WORLD,&ts));
  if(0){
    /* do some simple testing  */
    medium->dt = 0.001*sec_per_year;
    /* test Euler steps */
    while(medium->time < medium->stop_time){	
      rsf_ODE_RHSFunction(ts,medium->time,x,F,par);
      VecView(x,PETSC_VIEWER_STDOUT_WORLD);
      //VecView(F,PETSC_VIEWER_STDOUT_WORLD);
      /* update */
      PetscCall(VecAXPY(x, medium->dt, F));
      medium->time += medium->dt;
    }
    exit(-1);
  }


  
  PetscCall(TSSetProblemType(ts,TS_NONLINEAR)); /* can be linear? */
  PetscCall(TSSetSolution(ts,x));
  PetscCall(TSSetRHSFunction(ts, NULL, rsf_ODE_RHSFunction,&par));
  /*
    use runge kutta (for now)
  */
  PetscCall(TSSetType(ts,TSRK));
  /*
    Set the initial time and the initial timestep given above.
  */
  PetscCall(TSSetTime(ts,medium->time));	/* initial time, set above */
  PetscCall(TSSetTimeStep(ts,0.001*sec_per_year)); /* initial timestep */
  PetscCall(TSSetTolerances(ts,1e-3,NULL,1e-7,NULL));
  PetscCall(TSSetFromOptions(ts));

  field_out = 0;
  HEADNODE{
    fout1 = fopen("rsf_stats.dat","w");
  }
  /* for sending to head node */
  PetscCall(VecScatterCreateToZero(x,&ctx,&xout));		/* solution vector */
  while(medium->time < medium->stop_time){	/* advance until next output */
    PetscCall(TSSetMaxTime(ts,(medium->time + medium->print_interval))); /* advance solution to next output time */
    PetscCall(TSSolve(ts, x));
    PetscCall(TSGetTime(ts,&(medium->time))); /* get the current time */
    HEADNODE{
      /* 
	 distribute to zero node
      */
      
      // scatter as many times as you need
      PetscCall(VecScatterBegin(ctx,x,xout,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecScatterEnd(ctx,x,xout,INSERT_VALUES,SCATTER_FORWARD));
      //PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));
      /* assign to x solution vector */
      HEADNODE{
	PetscCall(VecGetArray(xout,&values));
	/* all patches loop */
	for(sum[0]=sum[1]=sum[2]=0.,vmin=1e20,vmax=1e-20,
	      i=0,j=i*INT_RSF_DIM;i < n;i++,j += INT_RSF_DIM){
	  state = values[j];
	  fault[i].s[STRIKE] = values[j+1]; /* shear stress */
	  fault[i].s[NORMAL] = values[j+2]; /* normal stress */
	  slip = values[j+3];
	  /* velocity */
	  fault[i].u[STRIKE] = vel_from_rsf(fault[i].s[STRIKE],fault[i].s[NORMAL],
					    state,fault[i].mu_s,medium->v0,dummy,(dummy+1),(dummy+2),medium);
	  /* compute average slip date */
	  if(fault[i].u[STRIKE] < vmin)
	    vmin = fault[i].u[STRIKE];
	  if(fault[i].u[STRIKE] > vmax)
	    vmax = fault[i].u[STRIKE];
	  sum[0] += fault[i].u[STRIKE];
	  sum[1] += fault[i].u[STRIKE] * fault[i].u[STRIKE]; /* mean and std of log10(vel) */
	  sum[2] += slip;				     /* mean slip */
	}
	/* output */
	vstd = sqrt(((COMP_PRECISION)n * sum[1] - sum[0] * sum[0]) / ((COMP_PRECISION)(n*(n-1))));
	vmean = sum[0]/(COMP_PRECISION)n;
	smean = sum[2]/(COMP_PRECISION)n;
	/* print some stats

	   time[yr] mean_vel std_vel min_vel max_vel mean_slip 
	*/
	fprintf(fout1,"%.5f %e %e %e %e %e\n",medium->time/sec_per_year,vmean,vstd,vmin,vmax,smean);
	/* some more output */
	if(medium->time - medium->slip_line_time > medium->slip_line_dt){
	  if(!field_out){
	    ierr = system("mkdir -p tmp_rsf");
	    if(ierr){
	      fprintf(stderr,"%s: error making output directory\n",argv[0]);
	      exit(-1);
	    }
	  }
	  /* 
	     print log10(slip rate) to file 
	  */
	  snprintf(vel_file,STRLEN,"tmp_rsf/vel-%012.5e-%06i-gmt",medium->time/sec_per_year,field_out);
	  fout2 = fopen(vel_file,"w");
	  for(i=0;i < n;i++)
	    print_patch_geometry_and_bc(i,fault,PSXYZ_STRIKE_DISP_OUT_MODE,medium->time,FALSE,fout2,FALSE,dummy);
	  fclose(fout2);
	  field_out++;
	  medium->slip_line_time = medium->time;
	}
	PetscCall(VecRestoreArray(xout,&values));
      }
    }
  }
  HEADNODE{
    fclose(fout1);
  }
  PetscCall(VecScatterDestroy(&ctx));
  
  /*View information about the time-stepping method and the solutionat the end time.*/
  PetscCall(TSView(ts, PETSC_VIEWER_STDOUT_SELF));

  PetscCall(MatDestroy(&medium->Is));
  PetscCall(MatDestroy(&medium->In));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&xout));

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
   x = (ψ1,τ1,σ1,s1,ψ2,τ2,σ2,s2,ψ3,τ3,σ3,s3,....) INT_RSF_DIM*n

   (see hmatrix_test/ode_solve_test.c for a working example)

   this will serve to update x and F (but we cannot use the slip_rate
   vector globally because it might be off state)

*/
PetscErrorCode rsf_ODE_RHSFunction(TS ts,PetscReal time,Vec X,Vec F,void *ptr)
{
  PetscScalar *f,cosh_fac,*mu,*scaled_tau,two_v0,two_v0_div_a,pre_fac,loc_tau_dot ;
  PetscScalar dvdtau,dvdsigma,dvdstate,loc_slip_rate,*exp_fac;
  const PetscScalar *x;
  struct med *medium;struct flt *fault;
  PetscInt i,j,k,ln;
  Vec tau_dot,sigma_dot,slip_rate;
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
  exp_fac =       (PetscScalar *)malloc(sizeof(PetscScalar)*medium->rn);
  /* make room to compute shear and normal stressing rate */
  PetscCall(MatCreateVecs(medium->Is, &tau_dot, &slip_rate));
  PetscCall(VecDuplicate(tau_dot, &sigma_dot));

  /* parallel assembly of slip_rate vector */
  /* i global, j global*INT_RSF_DIM, k local */
  for (i = medium->rs, j=i*INT_RSF_DIM, k=0; i < medium->re; i++,j+=INT_RSF_DIM,k++) { /* compute
								      prefactors
								      and
								      slip
								      rate
								      vector */
    //state = x[i*INT_RSF_DIM  ];
    //tau =   x[i*INT_RSF_DIM+1];
    //sigma = x[i*INT_RSF_DIM+2];
    //slip =  x[i*INT_RSF_DIM+3];
    
    /* compute velocity */
    /* v = 2 v_0 exp(-psi/a) sinh(tau/sigma/a) */
    /* and store the mu=tau/sigma, mu/a, exp(-psi/a) factors */
    PetscCall(VecSetValue(slip_rate, i, vel_from_rsf(x[j+1],x[j+2],x[j],fault[i].mu_s,medium->v0,
						     (mu+k), (scaled_tau+k), (exp_fac+k), medium),INSERT_VALUES));
  }
  /* assemble the slip rate vector */
  PetscCall(VecAssemblyBegin(slip_rate));
  PetscCall(VecAssemblyEnd(slip_rate));
  //VecView(slip_rate,PETSC_VIEWER_STDOUT_WORLD);

  /* forward solve effect of all slipping faults, this works for dense
     and H matrix */
  PetscCall(MatMult(medium->Is, slip_rate, tau_dot));	/* background loading rate will be added below */
  PetscCall(MatMult(medium->In, slip_rate, sigma_dot)); /*  */
  //PetscCall(VecView(tau_dot,PETSC_VIEWER_STDOUT_WORLD));
  //PetscCall(VecView(sigma_dot,PETSC_VIEWER_STDOUT_WORLD));
  /* 
     compute derivatives 
  */
  two_v0 = 2.0*medium->v0;
  /* loop through the node local indices */
  /* i global, j global*INT_RSF_DIM, k local, as above for the assembly */
  for (i = medium->rs,j=i*INT_RSF_DIM,k=0; i < medium->re; i++,j+=INT_RSF_DIM,k++) {

    /* state = x[j], tau = x[j+1], sigma = x[j+2], slip = x[j+3] */
    /* mu_s = a, mu_d = b */
    
    /* get path local shear stress rate and slip rate */
    PetscCall(VecGetValues(tau_dot,1,&i,&loc_tau_dot)); /* shear stressing rate on this patch from multiplication outcome */
    PetscCall(VecGetValues(slip_rate,1,&i,&loc_slip_rate)); /* local slip rate */
    /* 2v0/a */
    two_v0_div_a = two_v0/fault[i].mu_s;
    /* 
       d state/ dt = b/dc (v0 exp((f0-state)/b) - v) 
       
       this is the ageing law
    */
    f[j] = (fault[i].mu_d/medium->dc) * (medium->v0 * PetscExpReal((medium->f0 - x[j])/fault[i].mu_d) - loc_slip_rate);
    /* this is how it's implemented in HBI, using |v| instead of v */
    //f[j] = (fault[i].mu_d/medium->dc) * (medium->v0 * PetscExpReal((medium->f0 - x[j])/fault[i].mu_d) - fabs(loc_slip_rate));
    /* 
       d sigma/dt  = 
    */
    /* from interactions, based on multiplication outcome */
    PetscCall(VecGetValues(sigma_dot,1,&i,(f+j+2)));//f[j+2] = sigma_dot[i];
    //fprintf(stderr,"%g %g %g \n",loc_slip_rate,loc_tau_dot,f[j+2]);
    f[j+2] +=  fault[i].sinc[1];		    /* background normal stressing rate*/
    /* 
       d tau/dt 
    */
    /* dv/dstate */
    //dvdstate = -two_v0_div_a * exp(-psi/a) * sinh(tau/sigma/a) = -v/a
    dvdstate = -loc_slip_rate/fault[i].mu_s;
    /* dv/dtau */
    pre_fac =  two_v0_div_a/x[j+2];	      /* 2v0/(a*sigma) */
    cosh_fac  = PetscCoshReal(scaled_tau[k]); /* cosh(tau/sigma/a) */
    cosh_fac *= exp_fac[k];		      /* exp(-state/a) */
    dvdtau   =   pre_fac * cosh_fac;
    /* dv/dsigma */
    dvdsigma =  -pre_fac * cosh_fac * mu[k]; /* mu = tau/sigma */
    /* d tau/dt */
    f[j+1]  = 1./(1.0+medium->shear_mod_over_2cs_si * dvdtau);
    f[j+1] *= (loc_tau_dot  + fault[i].sinc[0] - medium->shear_mod_over_2cs_si * (dvdsigma * f[j+2] + dvdstate * f[j]));
    /* d slip/dt */
    f[j+3] = loc_slip_rate;
  }
  /*  */
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(F,&f));
  /*  */
  PetscCall(VecDestroy(&sigma_dot));
  PetscCall(VecDestroy(&tau_dot));
  PetscCall(VecDestroy(&slip_rate));
  /*  */
  free(mu);free(scaled_tau);free(exp_fac);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* 
   compute the sliding velocity from tau, sigma, and state

   and a bunch of helper factors to be reused
   mu = tau/sigma
   scaled_faut = mu/a
   exp_fac = exp(-state/a) 
 */
PetscReal vel_from_rsf(PetscReal tau, PetscReal sigma, PetscReal state, PetscReal a,PetscReal v0,
		       PetscReal *mu, PetscReal *scaled_tau, PetscReal *exp_fac, 
		       struct med *medium)
{
  PetscReal vel;
  *mu = tau/sigma;
  *scaled_tau = (*mu)/a;
  *exp_fac = PetscExpReal(-state/a);
  vel =  2.*v0 * (*exp_fac) * PetscSinhReal(*scaled_tau);
  return vel;
}

#endif
