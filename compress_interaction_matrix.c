#include "interact.h"
#ifdef USE_PETSC
#include "petsc_prototypes.h"

#ifdef USE_PETSC_HMAT
#include "kdtree.h"		/* for H2OPUS */
#endif

/*
  
  this function serves to test dense and H matrix forward and inverse
  operations
  
  reads in geometry file and calculates the interaction matrix, and
  then compresses it, testing forward and inverse computations
  
  see compress_interaction_matrix.md

*/



int main(int argc, char **argv)
{
  struct med *medium;
  struct flt *fault;
  struct interact_ctx ictx[1];
  /* timing */
  clock_t start_time,stop_time;
  PetscLogDouble t0,t1;

  double *bglobal,cpu_time_used;
  MatHtoolKernelFn *htools_kernel = GenKEntries_htools;
  KSP               ksp,ksph;
  PC                pc,pch;
  Vec         x, xh, b, bh, bout,d;
  Mat         Ah_dense;
  PetscReal   *coords=NULL,*avalues=NULL,*bvalues=NULL,norm[3];
  PetscInt    ndim, n, m, lm,ln,i,j,k,dn,on, *col_idx=NULL,rs,re;
  PetscInt nrandom = 0;	/* for timing tests */
  VecScatter ctx;
  PetscRandom rand_str;
  PetscBool read_value,flg,test_forward=PETSC_TRUE;
  PetscBool make_matrix_externally=PETSC_FALSE; /* make matrices here on in external routine (for testing) */
  double     target_x[3], sep_x;
  kd_node    nearest_x;
  PetscMPIInt h2opus_csize,h2opus_crank;
  
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
  medium->use_hmatrix = 1;	/* default is HTOOLS: 0,1,2,3,4 are options */
  
  /* set defaults, can always override */
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-use_hmatrix", &medium->use_hmatrix,&read_value));
  if(medium->use_hmatrix == 2){
#ifdef USE_PETSC_HMAT    
    /*  */
    fprintf(stderr,"%s: setting up for H2OPUS\n",argv[0]);
    fprintf(stderr,"%s: WARNING: H2OPUS construction ASSUMES A SYMMETRIC OPERATOR (not mutable:\n",
	    argv[0]);
    fprintf(stderr,"%s: WARNING: this h2opus only implements sampling-based construction for symmetric\n",
	    argv[0]);
    fprintf(stderr,"%s: WARNING: matrices); the H matrix approximates (K+K^T)/2, and the operator\n",argv[0]);
    fprintf(stderr,"%s: WARNING: asymmetry printed below sets an irreducible error floor\n",argv[0]);
    /* this runs before medium->comm_size is set */
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&h2opus_csize));
    if(h2opus_csize > 1){
      PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&h2opus_crank));
      if(h2opus_crank == 0)
	fprintf(stderr,"%s: H2OPUS sampling-based construction (MatCreateH2OpusFromMat) is not\n%s: supported in parallel in this PETSc/h2opus version - run serially\n",
		argv[0],argv[0]);
      exit(-1);
    }
    /* defaults for H2OPUS */
    medium->h2opus_eta = 0.6;	/*  */
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-eta", &medium->h2opus_eta, NULL));
    medium->h2opus_leafsize = 32;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-leafsize", &medium->h2opus_leafsize, NULL));
    medium->h2opus_basisord = 8;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-basisord", &medium->h2opus_basisord, NULL));
#else
    fprintf(stderr,"%s: H2OPUS requested but not compiled in (compile and add USE_PETSC_HMAT in makefile.petsc)\n",argv[0]);
    exit(-1);
#endif
  }else if(medium->use_hmatrix == 1){
#ifdef USE_PETSC_HMAT    
    fprintf(stderr,"%s: setting up for HTOOLS\n",argv[0]);
    /* HTOOLS */
    /* defaults, only applied if not given on the command line */
    PetscCall(PetscOptionsHasName(NULL,NULL,"-mat_htool_eta",&flg));
    if(!flg)PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_eta","100")); /* not sure  */
    PetscCall(PetscOptionsHasName(NULL,NULL,"-mat_htool_epsilon",&flg));
    if(!flg)PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_epsilon","1e-3"));
    PetscCall(PetscOptionsHasName(NULL,NULL,"-mat_htool_compressor",&flg));
    if(!flg)PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_compressor","SVD"));
    PetscCall(PetscOptionsHasName(NULL,NULL,"-pc_type",&flg));
    if(!flg)PetscCall(PetscOptionsSetValue(NULL,"-pc_type","none"));
#else
    fprintf(stderr,"%s: HTOOLS requested but not compiled in (compile and add USE_PETSC_HMAT in makefile.petsc)\n",argv[0]);
    exit(-1);
#endif
  }else{
    fprintf(stderr,"%s: use_hmatrix %i\n",argv[0],medium->use_hmatrix);
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
    fprintf(stderr,"%s: core %03i/%03i: assigning dense  row %5i to %5i\n",
	    argv[0],medium->comm_rank,medium->comm_size,medium->rs,medium->re);
    PetscTime(&t0);
    for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
      GenKEntries_htools(ndim,1,n,&j, col_idx, avalues,ictx);
      PetscCall(MatSetValues(medium->Is, 1, &j, n, col_idx,avalues, INSERT_VALUES));
    }
    PetscCall(PetscFree(avalues));
    PetscCall(PetscFree(col_idx));
    
    PetscCall(MatAssemblyBegin(medium->Is, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(medium->Is, MAT_FINAL_ASSEMBLY));
    PetscTime(&t1);
    HEADNODE
      fprintf(stderr,"%s: dense assembly took %12.4f s\n",argv[0],t1-t0);

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
#ifdef USE_PETSC_HMAT      
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
	  medium->kd_nodes[i].x[k] = coords[i*ndim+k];
	}
      }
      KDTreeSetup(medium->kdtree);
      PetscTime(&t1);
      PetscPrintf(PETSC_COMM_WORLD,"%s: kdtree setup for H2OPUS done: %1.4e sec (%i)\n",
		  argv[0],t1-t0,medium->nrflt);

      /* check if nodal locations are found */
      for (i=0; i < medium->nrflt; i++){
	for(j=0;j<ndim;j++)
	  target_x[j] = fault[i].x[j];
	KDTreeFindNearest(medium->kdtree,target_x,&nearest_x,&sep_x);
	if((nearest_x->index != i)||(sep_x > EPS_COMP_PREC)){
	  fprintf(stderr,"KD error: orig index %i - found index %i - sep: %g\n",i,nearest_x->index,sep_x);
	  exit(-1);
	}
      }
#else
      fprintf(stderr,"%s: H2OPUS requested but not compiled in (compile and add USE_PETSC_HMAT in makefile.petsc)\n",argv[0]);
      exit(-1);
#endif
    }
    /* here, we're using the normal stress interaction matrix for the
       H version of Is */
    PetscTime(&t0);
    if(medium->use_hmatrix == 4){	/* hmmvp via MATSHELL */
#ifdef USE_HMMVP
      double *xc,*yc,*zc;
      void *hmmvp_handle;
      PetscReal hmmvp_tol = 1.0e-5, hmmvp_eta = 3.0;
      int ic,hmmvp_nthreads = 1;
      long hmmvp_nnz;
      int hmm,hmn;
      PetscCall(PetscOptionsGetReal(NULL,NULL,"-hmmvp_tol",&hmmvp_tol,NULL));
      PetscCall(PetscOptionsGetReal(NULL,NULL,"-hmmvp_eta",&hmmvp_eta,NULL));
      PetscCall(PetscOptionsGetInt(NULL,NULL,"-hmmvp_nthreads",&hmmvp_nthreads,NULL));
      xc = (double *)malloc(sizeof(double)*m);
      yc = (double *)malloc(sizeof(double)*m);
      zc = (double *)malloc(sizeof(double)*m);
      for(ic=0;ic < m;ic++){
	xc[ic] = (double)ictx->fault[ic].x[INT_X];
	yc[ic] = (double)ictx->fault[ic].x[INT_Y];
	zc[ic] = (double)ictx->fault[ic].x[INT_Z];
      }
      fprintf(stderr,"%s: core %03i/%03i: assigning hmmvp m %i n %i tol %g (whole-matrix rel Frobenius) eta %g nthreads %i\n",
	      argv[0],medium->comm_rank,medium->comm_size,m,n,(double)hmmvp_tol,
	      (double)hmmvp_eta,hmmvp_nthreads);
      if((hmmvp_nthreads > 1) && (dc3dts() == 0)){
	fprintf(stderr,"%s: -hmmvp_nthreads %i requested but dc3d.F was compiled WITHOUT\n%s: -fopenmp: the THREADPRIVATE directives are inactive and threaded kernel\n%s: calls would corrupt the matrix - rebuild with -fopenmp in FFLAGS/LDFLAGS\n",
		argv[0],hmmvp_nthreads,argv[0],argv[0]);
	exit(-1);
      }
      hmmvp_handle = chmmvp_compress_in_memory((int)m,xc,yc,zc,(double)hmmvp_tol,
					       (double)hmmvp_eta,hmmvp_nthreads,
					       (void *)ictx);
      if(!hmmvp_handle){
	fprintf(stderr,"%s: hmmvp compression failed\n",argv[0]);
	exit(-1);
      }
      chmmvp_get_info(hmmvp_handle,&hmm,&hmn,&hmmvp_nnz);
      HEADNODE
	fprintf(stderr,"%s: hmmvp %i by %i, %ld stored scalars, compression ratio %.5g\n",
		argv[0],hmm,hmn,hmmvp_nnz,
		(double)((double)m*(double)n/(double)hmmvp_nnz));
      {
	hacapk_shell_ctx *hctx;
	Vec xd;
	hctx = (hacapk_shell_ctx *)malloc(sizeof(hacapk_shell_ctx));
	hctx->handle = hmmvp_handle;
	hctx->ball = (double *)malloc(sizeof(double)*m);
	PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,
				 (void *)hctx,&medium->In));
	PetscCall(MatShellSetOperation(medium->In,MATOP_MULT,(void (*)(void))MatMult_hmmvp));
	PetscCall(MatCreateVecs(medium->In,&xd,NULL));
	PetscCall(VecScatterCreateToAll(xd,&hctx->scat,&hctx->xall));
	PetscCall(VecGetOwnershipRange(xd,&hctx->rs,&hctx->re));
	PetscCall(VecDestroy(&xd));
      }
      free(xc);free(yc);free(zc);
#else
      fprintf(stderr,"%s: hmmvp requested but not compiled in (set HMMVP_DEFINES/HMMVP_LIBS/HMMVP_DIR in makefile.petsc)\n",argv[0]);
      exit(-1);
#endif
    }else if(medium->use_hmatrix == 3){	/* HACApK via MATSHELL */
#ifdef USE_HACAPK
      double *xc,*yc,*zc;
      void *hacapk_handle;
      PetscReal hacapk_ztol = 1.0e-4;
      int ic;
      PetscCall(PetscOptionsGetReal(NULL,NULL,"-hacapk_ztol",&hacapk_ztol,NULL));
      xc = (double *)malloc(sizeof(double)*m);
      yc = (double *)malloc(sizeof(double)*m);
      zc = (double *)malloc(sizeof(double)*m);
      for(ic=0;ic < m;ic++){
	xc[ic] = (double)ictx->fault[ic].x[INT_X];
	yc[ic] = (double)ictx->fault[ic].x[INT_Y];
	zc[ic] = (double)ictx->fault[ic].x[INT_Z];
      }
      hacapk_handle = cinit_hacapk_struct((int)m,(void *)ictx);
      cset_hacapk_struct_coord(hacapk_handle,xc,yc,zc);
      fprintf(stderr,"%s: core %03i/%03i: assigning HACApK m %i n %i ztol %g\n",
	      argv[0],medium->comm_rank,medium->comm_size,m,n,(double)hacapk_ztol);
      cmake_hacapk_struct_hmat(hacapk_handle,(double)hacapk_ztol);
      {
	hacapk_shell_ctx *hctx;
	Vec xd;
	hctx = (hacapk_shell_ctx *)malloc(sizeof(hacapk_shell_ctx));
	hctx->handle = hacapk_handle;
	hctx->ball = (double *)malloc(sizeof(double)*m);
	PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,
				 (void *)hctx,&medium->In));
	PetscCall(MatShellSetOperation(medium->In,MATOP_MULT,(void (*)(void))MatMult_HACApK));
	PetscCall(MatCreateVecs(medium->In,&xd,NULL));
	PetscCall(VecScatterCreateToAll(xd,&hctx->scat,&hctx->xall));
	PetscCall(VecGetOwnershipRange(xd,&hctx->rs,&hctx->re));
	PetscCall(VecDestroy(&xd));
      }
      free(xc);free(yc);free(zc);
#else
      fprintf(stderr,"%s: HACApK requested but not compiled in (set HACAPK_DEFINES/HACAPK_LIBS in makefile.petsc)\n",argv[0]);
      exit(-1);
#endif
    }else{
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
      if(medium->use_hmatrix == 1){
	PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,lm,ln, m, n,
					   ndim,(coords+rs*ndim), (coords+rs*ndim), htools_kernel, ictx, &medium->In));
      }else{
#if defined(PETSC_HAVE_H2OPUS)
      /* 
	 construct the H2 matrix by hierarchical randomized sampling of
	 the assembled dense operator (HARA), rather than from the
	 kernel callback: the FromKernel interface evaluates the kernel
	 by Chebyshev interpolation at arbitrary points within cluster
	 bounding boxes, which is incompatible with patch-pair Green's
	 functions that are only defined at element centers (the nearest
	 neighbor mapping yields a piecewise constant surrogate whose
	 polynomial interpolation has O(1) errors in all admissible
	 blocks); sampling-based construction only requires matrix
	 vector products and reproduces the true BEM operator to
	 -mat_h2opus_rtol (default 1e-4)
      */
	PetscCall(MatCreateH2OpusFromMat(medium->Is, ndim, coords, PETSC_FALSE,
					 medium->h2opus_eta, medium->h2opus_leafsize,
					 PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
					 &medium->In));
#else
	fprintf(stderr,"%s: H2OPUS requested but PETSc was built without h2opus\n",argv[0]);
	exit(-1);
#endif
      }
    }				/* end htool/h2opus (non-HACApK) construction */
    if(medium->use_hmatrix == 2){
      /* 
	 this version of h2opus only implements sampling-based
	 construction for symmetric matrices (hlru_sym assertion); the
	 result therefore approximates the symmetrized operator
	 (K+K^T)/2 - check the printed operator asymmetry to judge the
	 error this introduces
      */
      Mat KT;
      PetscReal nrmK,nrmD;
      PetscCall(MatTranspose(medium->Is,MAT_INITIAL_MATRIX,&KT));
      PetscCall(MatNorm(medium->Is,NORM_FROBENIUS,&nrmK));
      PetscCall(MatAXPY(KT,-1.0,medium->Is,SAME_NONZERO_PATTERN));
      PetscCall(MatNorm(KT,NORM_FROBENIUS,&nrmD));
      HEADNODE
	fprintf(stderr,"%s: operator asymmetry |K-K^T|_F/|K|_F = %.6e (h2opus approximates the symmetrized operator)\n",
		argv[0],(double)(nrmD/nrmK));
      PetscCall(MatDestroy(&KT));
      PetscCall(MatSetOption(medium->In, MAT_SYMMETRIC, PETSC_TRUE));
    }else{
      PetscCall(MatSetOption(medium->In, MAT_SYMMETRIC, PETSC_FALSE));
    }

    PetscCall(MatSetFromOptions(medium->In));
    
    PetscCall(MatAssemblyBegin(medium->In, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(medium->In, MAT_FINAL_ASSEMBLY));
    PetscTime(&t1);
    HEADNODE
      fprintf(stderr,"%s: H matrix assembly took %12.4f s\n",argv[0],t1-t0);
    free(coords);
  }else{
    /* use external routines */
    calc_petsc_Isn_matrices(medium, fault,0,                  1.0,0,&medium->Is); /* dense */
    calc_petsc_Isn_matrices(medium, fault,medium->use_hmatrix,1.0,0,&medium->In); /* Htools or HACApK */
  }
  /* dense */
  if(n < 20){
    HEADNODE
      fprintf(stderr,"%s: dense matrix:\n",argv[0]);
    PetscCall(MatView(medium->Is,PETSC_VIEWER_STDOUT_WORLD));
  }
  /* get info on H matrix */
  MatView(medium->In,PETSC_VIEWER_STDOUT_WORLD);
  if(n < 20){
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
      fprintf(stderr,"%s: it took %20.3fs for %05i H-matrix(%i) solves\n",
	      argv[0],cpu_time_used,nrandom+1,medium->use_hmatrix);
    
    if((m<20)&&(nrandom==0))
      VecView(bh,PETSC_VIEWER_STDOUT_WORLD);

    if(1){			/* accuracy of the x=1 product (random
				   loops above are pure timing and leave
				   b, bh from their last random x, so
				   recompute with x=1 here) */
      PetscCall(VecSet(x, 1.0));
      PetscCall(VecSet(xh,1.0));
      PetscCall(MatMult(medium->Is, x, b));
      PetscCall(MatMult(medium->In, xh, bh));
      /* compute difference */
      PetscCall(VecDuplicate(b, &d));PetscCall(VecCopy(b, d));
  
      PetscCall(VecAXPY(d,-1.0,bh));
      PetscCall(VecNorm(b,NORM_2,norm));
      PetscCall(VecNorm(bh,NORM_2,(norm+1)));
      PetscCall(VecNorm(d,NORM_2,(norm+2)));
      HEADNODE
	fprintf(stdout,"%s: |b| = %20.10e |b_h| = %20.10e |b-b_h|/|b| = %20.10e\n",
		argv[0],norm[0],norm[1],norm[2]/norm[0]);
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
      PetscCallMPI(MPI_Bcast(bglobal,m,MPI_DOUBLE,0, MPI_COMM_WORLD));
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
      fprintf(stderr,"%s: it took %20.3fs for %05i dense inverse solves\n",
	      argv[0],cpu_time_used,nrandom+1);
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
	  fprintf(stdout,"%s: |x| = %20.10e |x_h| = %20.10e |x-x_h|/|x| = %20.10e\n",
		  argv[0],norm[0],norm[1],norm[2]/norm[0]);
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
  exit(0);

}

#endif
