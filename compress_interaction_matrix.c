#include "interact.h"
#ifdef USE_PETSC
#include <petscksp.h>
/* this function is in interact.c */
PetscErrorCode GenEntries(PetscInt , PetscInt , PetscInt ,
				 const PetscInt *, const PetscInt *, PetscScalar *, void *);
#endif

/*
  reads in geometry file and calculates the interaction matrix, and
  then compresses it, testing forward and inverse computations

  can use HTOOLS Petsc options 

  also

  -geom_file geom.in - for the fault geometry
  -nrandom number (default: 0) - to run a number of solves with random vectors for timing)
  -test_forward true/false (default: true) - run matrix multiplication or inversion
  
  
*/

int main(int argc, char **argv)
{
#ifdef USE_PETSC
  struct med *medium;
  struct flt *fault;
  struct interact_ctx ictx[1];
  clock_t start_time,stop_time;
  double *bglobal,cpu_time_used;
  KSP               ksp,ksph;
  PC                pc,pch;
  Vec         x, xh, b, bh, bout,d;
  Mat         Ah,Ah_dense;
  PetscReal   *coords,*avalues=NULL,*bvalues=NULL,norm[3];
  PetscInt    ndim, n, m, lm,ln,i,j,k,dn,on, *col_idx=NULL,rs,re;
  PetscInt nrandom = 0;	/* for timing tests */
  VecScatter ctx;
  PetscRandom rand_str;
  PetscBool read_value,flg,test_forward=PETSC_TRUE;
  /* defined in interact.c */
  MatHtoolKernelFn *kernel = GenEntries;
  char geom_file[STRLEN]="geom.in";
  /* generate frameworks */
  medium=(struct med *)calloc(1,sizeof(struct med)); /* make one zero medium structure */
  ictx->medium = medium;
  ictx->src_slip_mode = STRIKE;
  ictx->rec_stress_mode = STRIKE;
  ndim = 3;
  /* 
     start up petsc 
  */
  PetscFunctionBegin;
  /* set defaults, can always override */
  PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_eta","100")); /* not sure  */
  PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_epsilon","1e-3"));
  PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_compressor","SVD"));
  /* start up Petsc proper */
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rand_str));
  PetscCall(PetscRandomSetFromOptions(rand_str));
  /* 
     
     read in geometry

  */
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nrandom", &nrandom,&read_value));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_forward", &test_forward,&read_value));
  HEADNODE{
    if(read_value)
      fprintf(stderr,"%s: reading geometry from %s as set by -geom_file\n",argv[0],geom_file);
    else
      fprintf(stderr,"%s: reading geometry from default, %s\n",argv[0],geom_file);
  }
  read_geometry(geom_file,&medium,&fault,TRUE,FALSE,FALSE,FALSE);
  
  ictx->fault = fault;
  m = n = medium->nrflt;
  coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*n);
  bglobal = (double *)malloc(sizeof(double)*m);
  HEADNODE{
    fprintf(stderr,"%s: computing %i by %i matrix\n",argv[0], m,n);
    
  }
  for(i=0;i<m;i++)		/* all sources or receiveer coordinates  */
    for(k=0;k<3;k++)
      coords[i*ndim+k] = fault[i].x[k];
 
  /* dense matrix setup */
  PetscCall(MatCreate(PETSC_COMM_WORLD, &medium->pA));
  PetscCall(MatSetSizes(medium->pA, PETSC_DECIDE, PETSC_DECIDE, m, n));
  PetscCall(MatSetType(medium->pA, MATDENSE));
  PetscCall(MatSetUp(medium->pA));
  PetscCall(MatGetLocalSize(medium->pA, &lm, &ln));
  dn = ln;on = n - ln;
  PetscCall(MatSeqAIJSetPreallocation(medium->pA, n, NULL));
  PetscCall(MatMPIAIJSetPreallocation(medium->pA, dn, NULL, on, NULL));
  PetscCall(MatGetOwnershipRange(medium->pA, &medium->rs, &medium->re));

  /*  */
  medium->rn = medium->re  - medium->rs; /* number of local elements */
  //fprintf(stderr,"%s: core %i: dn %i on %i n %i rs %i re %i \n",argv[0],medium->comm_rank,dn,on,n,medium->rs,medium->re);
  /*  */
  PetscCall(PetscCalloc(m*sizeof(PetscScalar), &avalues));
  PetscCall(PetscCalloc(n*sizeof(PetscInt), &col_idx));
  for (i=0; i < n; i++) 
    col_idx[i] = i;

  /* 
     assemble dense matrix 
  */
  fprintf(stderr,"%s: core %03i/%03i: assigning dense row %5i to %5i\n",
	  argv[0],medium->comm_rank,medium->comm_size,medium->rs,medium->re);
  for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
    GenEntries(ndim,1,n,&j, col_idx, avalues,ictx);
    PetscCall(MatSetValues(medium->pA, 1, &j, n, col_idx,avalues, INSERT_VALUES));
  }
  PetscCall(PetscFree(avalues));
  PetscCall(PetscFree(col_idx));
  PetscCall(MatAssemblyBegin(medium->pA, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(medium->pA, MAT_FINAL_ASSEMBLY));
  if(n<20){
    HEADNODE
      fprintf(stderr,"%s: dense matrix:\n",argv[0]);
    
    PetscCall(MatView(medium->pA,PETSC_VIEWER_STDOUT_WORLD));
  }
 
  /* 
     hirarchical version 
  */
  PetscCall(MatCreate(PETSC_COMM_WORLD, &Ah));
  PetscCall(MatSetSizes(Ah, PETSC_DECIDE, PETSC_DECIDE, m, n));  
  PetscCall(MatSetType(Ah,MATHTOOL));
  PetscCall(MatSetUp(Ah));
  PetscCall(MatGetLocalSize(Ah, &lm, &ln));
  dn = ln;on = n - ln;
  PetscCall(MatSeqAIJSetPreallocation(Ah, n, NULL));
  PetscCall(MatMPIAIJSetPreallocation(Ah, dn, NULL, on, NULL));
  PetscCall(MatGetOwnershipRange(Ah, &rs, &re));
  
  fprintf(stderr,"%s: core %03i/%03i: assigning htool row %5i to %5i\n",argv[0],medium->comm_rank,medium->comm_size,rs,re);
    
  //PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE, m, n, ndim,(coords+rs), coords, kernel,ictx, &Ah));
  PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,lm,ln, m, n, ndim,(coords+rs), (coords+re), kernel,ictx, &Ah));
    
  PetscCall(MatSetOption(Ah, MAT_SYMMETRIC, PETSC_FALSE));
  PetscCall(MatSetFromOptions(Ah));
    
  PetscCall(MatAssemblyBegin(Ah, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Ah, MAT_FINAL_ASSEMBLY));
  /* get info on H matrix */
  MatView(Ah,PETSC_VIEWER_STDOUT_WORLD);
  if(n<20){
    HEADNODE
      fprintf(stderr,"%s: H matrix converted back to dense:\n",argv[0]);
    
    PetscCall(MatConvert(Ah,MATDENSE,MAT_INITIAL_MATRIX,&Ah_dense));
    MatView(Ah_dense,PETSC_VIEWER_STDOUT_WORLD);
  }
  PetscCall(MatCreateVecs(medium->pA, &x, &b));/* For A x = b: x -> left, b -> right */
  PetscCall(MatCreateVecs(Ah, &xh, &bh));/* For A x = b: x -> left, b -> right */
  if(test_forward){
    
    /* 
       test matrix multiplication 
       b = A x
    */
    PetscCall(VecSet(x, 1.0));
    PetscCall(VecSet(xh,1.0));
    /* do we need those? */
    PetscCall(VecAssemblyBegin(x));PetscCall(VecAssemblyEnd(x));
    PetscCall(VecAssemblyBegin(xh));PetscCall(VecAssemblyEnd(xh));
    
    
    start_time = clock();
    /* dense solver */
    PetscCall(MatMult(medium->pA, x, b));
    for(i=0;i<nrandom;i++){
      PetscCall(VecSetRandom(x,rand_str));
      PetscCall(MatMult(medium->pA, x, b));
    }
    stop_time = clock();  
    cpu_time_used = ((double)stop_time-start_time)/CLOCKS_PER_SEC;
    HEADNODE
      fprintf(stderr,"%s: it took %20.3fs for %05i dense solves\n",argv[0],cpu_time_used,nrandom+1);
    if((m<20)&&(nrandom==0))
      VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    start_time = clock();
    /* H matrix solve */
    PetscCall(MatMult(Ah, xh, bh));
    for(i=0;i<nrandom;i++){
      PetscCall(VecSetRandom(xh,rand_str));
      PetscCall(MatMult(Ah, xh, bh));
    }
    stop_time = clock();
    cpu_time_used = ((double)stop_time-start_time)/CLOCKS_PER_SEC;
    HEADNODE
      fprintf(stderr,"%s: it took %20.3fs for %05i htool solves\n",argv[0],cpu_time_used,nrandom+1);
    
    if((m<20)&&(nrandom==0))
      VecView(bh,PETSC_VIEWER_STDOUT_WORLD);

    if(nrandom==0){
      /* compute difference */
      PetscCall(VecDuplicate(b, &d));PetscCall(VecCopy(b, d));
  
      PetscCall(VecAXPY(d,-1.0,bh));
      PetscCall(VecNorm(b,NORM_2,norm));
      PetscCall(VecNorm(bh,NORM_2,(norm+1)));
      PetscCall(VecNorm(d,NORM_2,(norm+2)));
      HEADNODE
	fprintf(stdout,"%s: |b| = %20.10e |b_h| = %20.10e |b-b_h|/|b| = %20.10e\n",argv[0],norm[0],norm[1],norm[2]/norm[0]);
      
      /* get b values */
      PetscCall(VecScatterCreateToZero(b,&ctx,&bout));
      PetscCall(VecScatterBegin(ctx,b,bout,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecScatterEnd(ctx,b,bout,INSERT_VALUES,SCATTER_FORWARD));
      HEADNODE{
	PetscCall(VecGetArray(bout,&bvalues));
	for(i=0;i<m;i++)
	  bglobal[i] = bvalues[i];
	PetscCall(VecRestoreArray(bout,&bvalues));
      }
      MPI_Bcast(bglobal,m,MPI_DOUBLE,0, MPI_COMM_WORLD);
      PetscCall(VecScatterDestroy(&ctx));
      PetscCall(VecDestroy(&bout));
      HEADNODE{
	//  for(i=0;i<m;i++)
	//fprintf(stdout,"%20.10e\n",bglobal[i]);
      }
    }
  }else{
    /* 
       test inverse  x = A^-1 b
    */
    /* make constext for solver */
    /* dense */
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, medium->pA, medium->pA));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCLU));
    /* H */
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksph));
    PetscCall(KSPSetOperators(ksph, Ah, Ah));
    PetscCall(KSPSetFromOptions(ksph));
    PetscCall(KSPGetPC(ksph, &pch));
    PetscCall(PetscObjectTypeCompare((PetscObject)pch, PCHPDDM, &flg));
    /*  */
    
    PetscCall(VecSet(b, 1.0));
    PetscCall(VecSet(bh,1.0));
    /* do we need those? */
    PetscCall(VecAssemblyBegin(b));PetscCall(VecAssemblyEnd(b));
    PetscCall(VecAssemblyBegin(bh));PetscCall(VecAssemblyEnd(bh));
    /* 

       dense solver 

    */
    start_time = clock();
    PetscCall(KSPSolve(ksp, b, x));
    for(i=0;i<nrandom;i++){	/* random trials */
      PetscCall(VecSetRandom(b,rand_str));
      PetscCall(KSPSolve(ksp, b, x));
    }
    stop_time = clock();  
    cpu_time_used = ((double)stop_time-start_time)/CLOCKS_PER_SEC;
    HEADNODE
      fprintf(stderr,"%s: it took %20.3fs for %05i dense inverse solves\n",argv[0],cpu_time_used,nrandom+1);
    if((m<20)&&(nrandom==0))
      VecView(x,PETSC_VIEWER_STDOUT_WORLD);
    /* 

       H matrix solve 
    */
    start_time = clock();
    PetscCall(KSPSolve(ksph, bh, xh));
    for(i=0;i<nrandom;i++){
      PetscCall(VecSetRandom(bh,rand_str));
      PetscCall(KSPSolve(ksph, bh, xh));
    }
    stop_time = clock();
    cpu_time_used = ((double)stop_time-start_time)/CLOCKS_PER_SEC;
    HEADNODE
      fprintf(stderr,"%s: it took %20.3fs for %05i htool inverse solves\n",argv[0],cpu_time_used,nrandom+1);
    
    if((m<20)&&(nrandom==0))
      VecView(xh,PETSC_VIEWER_STDOUT_WORLD);

    if(nrandom==0){
      /* compute difference */
      PetscCall(VecDuplicate(x, &d));PetscCall(VecCopy(x, d));
  
      PetscCall(VecAXPY(d,-1.0,xh));
      PetscCall(VecNorm(x,NORM_2,norm));
      PetscCall(VecNorm(xh,NORM_2,(norm+1)));
      PetscCall(VecNorm(d,NORM_2,(norm+2)));
      HEADNODE
	fprintf(stdout,"%s: |x| = %20.10e |x_h| = %20.10e |x-x_h|/|x| = %20.10e\n",argv[0],norm[0],norm[1],norm[2]/norm[0]);
    }

  }
  
 
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(VecDestroy(&xh));
  PetscCall(VecDestroy(&bh));
  if(nrandom==0){
    PetscCall(VecDestroy(&d));
  }
  if(n<20)
    PetscCall(MatDestroy(&Ah_dense));
  free(coords);free(bglobal);
  PetscCall(MatDestroy(&medium->pA));
  PetscCall(MatDestroy(&Ah));
  PetscCall(PetscFinalize());
#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);

}

