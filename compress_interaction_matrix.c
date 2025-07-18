#include "interact.h"
#ifdef USE_PETSC_HMAT

/* this function is in interact.c */
#include "petsc_prototypes.h"


#include "kdtree.h"		/* for H2OPUS */

#endif

/*
  
  this function serves to test dense and H matrix forward and inverse operations
  
  reads in geometry file and calculates the interaction matrix, and
  then compresses it, testing forward and inverse computations

  can use HTOOLS and H2OPUS Petsc options 

  also

  -geom_file geom.in - for the fault geometry
  -nrandom number (default: 0) - to run a number of solves with random vectors for timing)
  -test_forward true/false (default: true) - run matrix multiplication or inversion
  -use_h2opus if set, use that package, else, use htools
  
  
*/

int main(int argc, char **argv)
{
#ifdef USE_PETSC_HMAT
  struct med *medium;
  struct flt *fault;
  struct interact_ctx ictx[1];
  /* timing */
  clock_t start_time,stop_time;
  PetscLogDouble t0,t1;

  double *bglobal,cpu_time_used;
  MatHtoolKernelFn *htools_kernel = GenKEntries_htools;
  MatH2OpusKernelFn *h2opus_kernel = GenKEntries_h2opus;
  KSP               ksp,ksph;
  PC                pc,pch;
  Vec         x, xh, b, bh, bout,d;
  Mat         Ah_dense;
  PetscReal   *coords,*avalues=NULL,*bvalues=NULL,norm[3];
  PetscInt    ndim, n, m, lm,ln,i,j,k,dn,on, *col_idx=NULL,rs,re;
  PetscInt nrandom = 0;	/* for timing tests */
  VecScatter ctx;
  PetscRandom rand_str;
  PetscBool read_value,flg,test_forward=PETSC_TRUE,use_h2opus=PETSC_FALSE;
  PetscBool make_matrix_externally=PETSC_FALSE; /* make matrices here on in external routine (for testing) */
  
  char geom_file[STRLEN]="geom.in";
  /* IMPORTANT */
  PetscFunctionBeginUser;

  /* start up Petsc proper */
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));

  
  /* generate frameworks */
  medium=(struct med *)calloc(1,sizeof(struct med)); /* make one zero medium structure */
  ictx->medium = medium;
  ictx->src_slip_mode = STRIKE;	/* slip mode */
  ictx->rec_stress_mode = STRIKE; /* recording stress mode */
  
  ndim = 3;
  /* 
     start up petsc 
  */

  
  /* set defaults, can always override */
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-use_h2opus", &use_h2opus,&read_value));
  if(use_h2opus){
    /*  */
    fprintf(stderr,"%s: setting up for H2OPUS\n",argv[0]);
    medium->use_hmatrix = 2;
    
    medium->h2opus_eta = 0.6;	/*  */
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-eta", &medium->h2opus_eta, NULL));
    medium->h2opus_leafsize = 32;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-leafsize", &medium->h2opus_leafsize, NULL));
    medium->h2opus_basisord = 8;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-basisord", &medium->h2opus_basisord, NULL));

  }else{

    fprintf(stderr,"%s: setting up for HTOOLS\n",argv[0]);
	
    medium->use_hmatrix = 1;
    /* how do I get those assigned internally? */
    /* HTOOLS */
    PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_eta","100")); /* not sure  */
    PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_epsilon","1e-3"));
    PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_compressor","SVD"));
    PetscCall(PetscOptionsSetValue(NULL,"-pc_type","none"));
  }

  

  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rand_str));
  PetscCall(PetscRandomSetFromOptions(rand_str));
  /* 
     
     read in geometry

  */
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  /* options */
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
  /*  */
  m = n = medium->nrflt;
  /*  */
  bglobal = (double *)malloc(sizeof(double)*m);

  
  if(!make_matrix_externally){

    HEADNODE{
      fprintf(stderr,"%s: computing %i by %i matrix LOCALLY\n",argv[0], m,n);
      
    }
    /* 
       dense matrix setup, using medium->Is
    */
    PetscCall(MatCreate(PETSC_COMM_WORLD, &medium->Is));
    PetscCall(MatSetSizes(medium->Is, PETSC_DECIDE, PETSC_DECIDE, m, n));
    
    PetscCall(MatSetType(medium->Is, MATDENSE));
    PetscCall(MatSetFromOptions(medium->Is));
    
    PetscCall(MatSetUp(medium->Is));
    PetscCall(MatGetLocalSize(medium->Is, &lm, &ln));
    dn = ln;on = n - ln;
    PetscCall(MatSeqAIJSetPreallocation(medium->Is, n, NULL));
    PetscCall(MatMPIAIJSetPreallocation(medium->Is, dn, NULL, on, NULL));
    PetscCall(MatGetOwnershipRange(medium->Is, &medium->rs, &medium->re));
    
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
      GenKEntries_htools(ndim,1,n,&j, col_idx, avalues,ictx);
      PetscCall(MatSetValues(medium->Is, 1, &j, n, col_idx,avalues, INSERT_VALUES));
    }
    PetscCall(PetscFree(avalues));
    PetscCall(PetscFree(col_idx));
    
    PetscCall(MatAssemblyBegin(medium->Is, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(medium->Is, MAT_FINAL_ASSEMBLY));

    if(medium->use_hmatrix){	/* for all, get coordinates */
      /* 
	 hirarchical version, using medium->In
      */
      coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*medium->nrflt);
      for(i=0;i < medium->nrflt;i++)		/* all sources or receiveer coordinates  */
	for(k=0;k < ndim;k++)
	  coords[i*ndim+k] = fault[i].x[k];
    }
    if(medium->use_hmatrix == 2){
      /* H2OPUS setup */
      PetscPrintf(PETSC_COMM_WORLD,"%s: kdtree setup for H2OPUS %i nodes\n",argv[0],m);
      PetscTime(&t0);
      /*  */
      KDTreeCreate(ndim,&medium->kdtree);
      KDTreeSetPoints(medium->kdtree,medium->nrflt);
    
      KDTreeGetPoints(medium->kdtree,&m,&medium->kd_nodes);
      for (i=0; i < medium->nrflt; i++) {
	medium->kd_nodes[i].index = i;
	for (k=0; k < ndim; k++) { /* speed up later, could use fault
				      structure itself for sort */
	  medium->kd_nodes[i].x[k] = (double)fault[i].x[k];
	}
      }
      KDTreeSetup(medium->kdtree);
      PetscTime(&t1);
      PetscPrintf(PETSC_COMM_WORLD,"%s: kdtree setup for H2OPUS done: %1.4e sec (%i)\n",argv[0],t1-t0,medium->nrflt);
    }
    /* here, we're using the normal stress interaction matrix for the
       H version of Is */
    
    PetscCall(MatCreate(PETSC_COMM_WORLD, &medium->In));
    PetscCall(MatSetSizes(medium->In, PETSC_DECIDE, PETSC_DECIDE, m, n));  
    
    if(medium->use_hmatrix==1)
      PetscCall(MatSetType(medium->In,MATHTOOL));
    else
      PetscCall(MatSetType(medium->In,MATH2OPUS));
    /*  */
    PetscCall(MatSetUp(medium->In));
    PetscCall(MatGetLocalSize(medium->In, &lm, &ln));
    dn = ln;on = n - ln;
    PetscCall(MatSeqAIJSetPreallocation(medium->In, n, NULL));
    PetscCall(MatMPIAIJSetPreallocation(medium->In, dn, NULL, on, NULL));
    PetscCall(MatGetOwnershipRange(medium->In, &rs, &re));
    
    fprintf(stderr,"%s: core %03i/%03i: assigning %s row %5i to %5i, lm %i ln %i m %i n %i\n",
	    argv[0],medium->comm_rank,medium->comm_size,(medium->use_hmatrix==1)?"HTOOLS":"H2OPUS",
	    rs,re,lm,ln,m,n);
    if(medium->use_hmatrix==1){
      PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,lm,ln, m, n,
					 ndim,(coords+rs), (coords+re), htools_kernel, ictx, &medium->In));
    }else{
      PetscCall(MatCreateH2OpusFromKernel(PETSC_COMM_WORLD, lm, ln, m, n,
					  ndim, (coords+rs), PETSC_FALSE, h2opus_kernel,ictx,
					  medium->h2opus_eta, medium->h2opus_leafsize,
					  medium->h2opus_basisord, &medium->In));
    }
    PetscCall(MatSetOption(medium->In, MAT_SYMMETRIC, PETSC_FALSE));      
    

    PetscCall(MatSetFromOptions(medium->In));
    
    PetscCall(MatAssemblyBegin(medium->In, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(medium->In, MAT_FINAL_ASSEMBLY));
    free(coords);
  }else{
    /* use external routines */
    calc_petsc_Isn_matrices(medium, fault,0,1.0,0,&medium->Is); /* dense */
    calc_petsc_Isn_matrices(medium, fault,1,1.0,0,&medium->In); /* Htools */
  }
  /* dense */
  if(n<20){
    HEADNODE
      fprintf(stderr,"%s: dense matrix:\n",argv[0]);
    PetscCall(MatView(medium->Is,PETSC_VIEWER_STDOUT_WORLD));
  }
  /* get info on H matrix */
  MatView(medium->In,PETSC_VIEWER_STDOUT_WORLD);
  if(n<20){
    HEADNODE
      fprintf(stderr,"%s: H matrix converted back to dense:\n",argv[0]);
    PetscCall(MatConvert(medium->In,MATDENSE,MAT_INITIAL_MATRIX,&Ah_dense));
    PetscCall(MatView(Ah_dense,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(MatDestroy(&Ah_dense));
  }
  /* 
   */
  PetscCall(MatCreateVecs(medium->Is, &x, &b));/* For A x = b: x -> left, b -> right */
  PetscCall(MatCreateVecs(medium->In, &xh, &bh));/* For A x = b: x -> left, b -> right */
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
    PetscCall(MatMult(medium->Is, x, b));
    for(i=0;i<nrandom;i++){
      PetscCall(VecSetRandom(x,rand_str));
      PetscCall(MatMult(medium->Is, x, b));
    }
    stop_time = clock();  
    cpu_time_used = ((double)stop_time-start_time)/CLOCKS_PER_SEC;
    HEADNODE
      fprintf(stderr,"%s: it took %20.3fs for %05i dense solves\n",argv[0],cpu_time_used,nrandom+1);
    if((m<20)&&(nrandom==0))
      VecView(b,PETSC_VIEWER_STDOUT_WORLD);

    start_time = clock();
    /* 
       H matrix solve 
    */
    PetscCall(MatMult(medium->In, xh, bh));
    for(i=0;i<nrandom;i++){
      PetscCall(VecSetRandom(xh,rand_str));
      PetscCall(MatMult(medium->In, xh, bh));
    }
    stop_time = clock();
    cpu_time_used = ((double)stop_time-start_time)/CLOCKS_PER_SEC;
    HEADNODE
      fprintf(stderr,"%s: it took %20.3fs for %05i H-matrix(%i) solves\n",argv[0],cpu_time_used,nrandom+1,medium->use_hmatrix);
    
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
      PetscCall(VecDestroy(&xh));
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
    /* make context for solver */
    /* dense */
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, medium->Is, medium->Is));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCLU));
    PetscCall(KSPSetFromOptions(ksp)); 
    PetscCall(VecSet(b, 1.0));
    /* do we need those? */
    PetscCall(VecAssemblyBegin(b));PetscCall(VecAssemblyEnd(b));
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
    
    if(1){
      /* 
	 H 
      */
      PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksph));
      PetscCall(KSPSetOptionsPrefix(ksph,"htool_"));
      PetscCall(KSPSetOperators(ksph, medium->In, medium->In));
      PetscCall(KSPSetFromOptions(ksph));
      PetscCall(KSPGetPC(ksph, &pch));
      PetscCall(PetscObjectTypeCompare((PetscObject)pch, PCHPDDM, &flg));
      PetscCall(VecSet(bh,1.0));
      PetscCall(VecAssemblyBegin(bh));PetscCall(VecAssemblyEnd(bh));
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
	fprintf(stderr,"%s: it took %20.3fs for %05i Hmatrix(%i) inverse solves\n",
		argv[0],cpu_time_used,nrandom+1,medium->use_hmatrix);
      
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
	PetscCall(VecDestroy(&d));
      }
      PetscCall(VecDestroy(&xh));
    }
  }
  
 
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
 
  PetscCall(VecDestroy(&bh));
  free(bglobal);
  if(medium->use_hmatrix==2)	/* free the KD tree */
    KDTreeDestroy(&medium->kdtree);

  PetscCall(MatDestroy(&medium->Is));
  PetscCall(MatDestroy(&medium->In));
  PetscCall(PetscFinalize());
#else
  fprintf(stderr,"%s only petsc H matrix version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);

}

