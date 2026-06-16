#include "interact.h"
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#endif

/*
  
  this function serves to test dense and H matrix forward and inverse
  operations
  
  reads in geometry file and calculates the interaction matrix, and
  then compresses it, testing forward and inverse computations
  
  see compress_interaction_matrix.md for full discussion

  --------------------------------------------------------------------
  BACKEND SCALING SUMMARY (N=14400, single 48-core node, 100 matvecs,
  matched ~1e-6 error band; see hmat_scaling_test.sh and the .md):

    settings: HTOOL  -mat_htool_epsilon 3e-5 -mat_htool_eta 10  (6.6e-7)
              HACApK -hacapk_ztol 1e-1                          (2.2e-7)
              hmmvp  -hmmvp_tol 1e-7                            (1.6e-6)
    (HACApK ztol is very conservative for this smooth Okada kernel, so
     ztol 1e-1 - much looser than nominal - is what makes it comparable
     rather than near-dense; do not run it at ztol 1e-4 for speed tests.)

    assembly [s], 1->48 cores: HACApK 25.4->1.1 (fastest, ~23x)
                               hmmvp  27.3->1.4 (np1~np2: master/worker)
                               HTOOL  567.6->6.7 (most costly, ~85x super-linear)
    matvec   [s], 1->48 cores: hmmvp  2.651->0.061 (fastest from np>=4, eff ~0.91)
                               HTOOL  3.849->0.078 (eff ~1.0, 2nd at high core)
                               HACApK 2.402->0.117 (fastest at np<=2, saturates, eff ~0.43)

  RECOMMENDATION (this geometry/band/machine; may shift with N, tol,
  library version, hardware - rerun hmat_scaling_test.sh):
    - matvec-bound / iterative-solver use (e.g. rsf_solve), np>=4:
      prefer hmmvp (-use_hmatrix 4, -hmmvp_tol 1e-7) - fastest matvec,
      best scaling, near-cheapest assembly. hmmvp's error is floored
      ~1e-6 at this N (Frobenius-estimate limited); use HACApK/HTOOL if
      a tighter operator is needed.
    - assembly-heavy / low core count / deterministic operator:
      HACApK (-use_hmatrix 3, -hacapk_ztol 1e-1).
    - build-once amortized over very many matvecs: HTOOL (-use_hmatrix 1).
  --------------------------------------------------------------------

*/

int main(int argc, char **argv)
{
#ifdef USE_PETSC
  struct med *medium;
  struct flt *fault;
  struct interact_ctx ictx[1];
  /* timing */
  clock_t start_time,stop_time;
  PetscLogDouble t0,t1;

  double *bglobal,cpu_time_used;
#ifdef USE_PETSC_HMAT
  MatHtoolKernelFn *htools_kernel = GenKEntries_htools;
#endif
  KSP               ksp,ksph;
  PC                pc,pch;
  Vec         x, xh, b, bh, bout,d;
  Mat         Adense,AH,AH_dense;
  PetscReal   *coords=NULL,*avalues=NULL,*bvalues=NULL,norm[3];
  PetscInt    ndim, n, m, lm,ln,i,j,k,dn,on, *col_idx=NULL,rs,re;
  PetscInt nrandom = 0;	/* for timing tests */
  VecScatter ctx;
  PetscRandom rand_str;
  PetscBool read_value,flg,test_forward=PETSC_TRUE;
  PetscBool make_matrix_externally=PETSC_FALSE; /* make matrices here on in external routine (for testing) */
  char geom_file[STRLEN]="geom.in";
#ifdef USE_PETSC_HMAT
  Mat KT;
  PetscReal nrmK,nrmD;
#endif
  hacapk_shell_ctx *hsc_dense,*hsc_h;
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
  double *xc,*yc,*zc;
  Vec xd;
#endif
#ifdef USE_HMMVP
#ifdef USE_HMMVP_MPI
  char hmmvp_fn[STRLEN], hmmvp_tmp[STRLEN];
  int cret;
#endif
  void *hmmvp_handle;
  long hmmvp_nnz;
  int hmm,hmn;
#endif 
#ifdef USE_HACAPK
  void *hacapk_handle;
#endif
  /* IMPORTANT */
  PetscFunctionBeginUser;

  /* start up Petsc proper */
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
  
  /* generate frameworks */
  medium=(struct med *)calloc(1,sizeof(struct med)); /* make one zero medium structure */
  ictx->medium = medium;
  /*  */
  ictx->src_slip_mode = STRIKE;	/* slip mode */
  ictx->rec_stress_mode = STRIKE; /* recording stress mode */
  
  ndim = 3;
  /* 
     start up petsc 
  */
  medium->use_hmatrix = 1;	/* default is HTOOLS: 0,1,2,3,4 are options */
  
  /* set defaults, can always override */
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-use_hmatrix", &medium->use_hmatrix,&read_value));
  /* options */
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nrandom", &nrandom,&read_value));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_forward", &test_forward,&read_value));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-make_matrix_externally", &make_matrix_externally,NULL));
  
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  
  /* set up defaults */
  switch(medium->use_hmatrix){
  case 1:
#ifdef USE_PETSC_HMAT
    HEADNODE
      fprintf(stderr,"%s: setting up for HTOOLS\n",argv[0]);
    /* HTOOLS */
    set_htools_defaults_and_options(medium);
#else
    fprintf(stderr,"%s: HTOOLS requested but not compiled in (compile and add USE_PETSC_HMAT in makefile.petsc)\n",argv[0]);
    exit(-1);
#endif
    break;
  case 2:
#ifdef USE_PETSC_HMAT
    /*  */
    fprintf(stderr,"%s: setting up for H2OPUS\n",argv[0]);
    set_h2opus_defaults_and_options(medium);
    fprintf(stderr,"%s: WARNING: H2OPUS construction ASSUMES A SYMMETRIC OPERATOR (not mutable:\n",
	    argv[0]);
    fprintf(stderr,"%s: WARNING: this h2opus only implements sampling-based construction for symmetric\n",
	    argv[0]);
    fprintf(stderr,"%s: WARNING: matrices); the H matrix approximates (K+K^T)/2, and the operator\n",argv[0]);
    fprintf(stderr,"%s: WARNING: asymmetry printed below sets an irreducible error floor\n",argv[0]);
    if(medium->comm_size > 1){
      if(medium->comm_rank == 0)
	fprintf(stderr,"%s: H2OPUS sampling-based construction (MatCreateH2OpusFromMat) is not\n%s: supported in parallel in this PETSc/h2opus version - run serially\n",
		argv[0],argv[0]);
      exit(-1);
    }
#else
    HEADNODE
      fprintf(stderr,"%s: H2OPUS requested but not compiled in (compile and add USE_PETSC_HMAT in makefile.petsc)\n",argv[0]);
    exit(-1);
#endif
    break;
  case 3:
#ifdef USE_HACAPK
    HEADNODE
      fprintf(stderr,"%s: setting up for HACApK\n",argv[0]);
    set_hacapk_defaults_and_options(medium);
#else
    fprintf(stderr,"%s: HACAPK requested but not compiled in (see USE_HACAPK and makefile.petc)\n",argv[0]);
    exit(-1);
#endif
    break;
  case 4:
#ifdef USE_HMMVP
    HEADNODE
      fprintf(stderr,"%s: setting up for HMMVP\n",argv[0]);
    set_hmmvp_defaults_and_options(medium);
#else
    fprintf(stderr,"%s: HMMVP requested but not compiled in (see USE_HMMVP, hmmvp subdirectory, and makefile.petc)\n",argv[0]);
    exit(-1);
#endif
    break;
  case 0:			/* dense */
    break;
  }

  
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rand_str));
  PetscCall(PetscRandomSetFromOptions(rand_str));
  /* 
     
     read in geometry

  */
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));

 

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
       dense matrix setup, using Adense
    */
    PetscCall(MatCreate(PETSC_COMM_WORLD, &Adense));
    PetscCall(MatSetSizes(Adense, PETSC_DECIDE, PETSC_DECIDE, m, n));
    PetscCall(MatSetType(Adense, MATDENSE));
    PetscCall(MatSetFromOptions(Adense));
    
    PetscCall(MatSetUp(Adense));
    PetscCall(MatGetLocalSize(Adense, &lm, &ln));
    dn = ln;on = n - ln;
    PetscCall(MatSeqAIJSetPreallocation(Adense, n, NULL));
    PetscCall(MatMPIAIJSetPreallocation(Adense, dn, NULL, on, NULL));
    PetscCall(MatGetOwnershipRange(Adense, &medium->rs, &medium->re));
    
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
      PetscCall(MatSetValues(Adense, 1, &j, n, col_idx,avalues, INSERT_VALUES));
    }
    PetscCall(PetscFree(avalues));
    PetscCall(PetscFree(col_idx));
    PetscCall(MatAssemblyBegin(Adense, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Adense, MAT_FINAL_ASSEMBLY));
    PetscTime(&t1);
    HEADNODE
      fprintf(stderr,"%s: dense assembly took %12.4f s\n",argv[0],t1-t0);
    /* dense done */

    /* 
       
       ASSEMBLE DIFFERENT KINDS OF H MATRICES
       
    */
    /* 
       prepare with coordinates and such  
    */
    if((medium->use_hmatrix==1)||(medium->use_hmatrix==2)){	
      coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*medium->nrflt);
      for(i=0;i < medium->nrflt;i++)		/* all sources or receiveer coordinates  */
	for(k=0;k < ndim;k++)
	  coords[i*ndim+k] = fault[i].x[k];
    }
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
    if((medium->use_hmatrix==3)||(medium->use_hmatrix==4)){	
      xc = (double *)malloc(sizeof(double)*m);
      yc = (double *)malloc(sizeof(double)*m);
      zc = (double *)malloc(sizeof(double)*m);
      for(i=0;i < m;i++){
	xc[i] = (double)ictx->fault[i].x[INT_X];
	yc[i] = (double)ictx->fault[i].x[INT_Y];
	zc[i] = (double)ictx->fault[i].x[INT_Z];
      }
    }
#endif
    PetscTime(&t0);
    switch(medium->use_hmatrix){
    case 0:
      /*
	DENSE REFERENCE (use_hmatrix=0): build In as a dense copy of Is
	so that the forward (b = A x) and inverse (x = A\b) tests below
	compare the dense operator against itself. This is a
	self-consistency check; |b-b_h|/|b| should be at the level of
	machine precision and |x-x_h|/|x| at the level of the KSP solver
	tolerance. Note that the dense baseline used for the comparison
	is always Adense - use_hmatrix=0 simply makes the second
	("H") operator dense as well. 
	
	Previously this value fell through to the MATH2OPUS branch
	without the (use_hmatrix==2 gated) coordinate/kdtree setup and
	aborted in MatAssemblyEnd.
      */
      HEADNODE
	fprintf(stderr,"%s: core %03i/%03i: assigning dense (use_hmatrix=0 self-consistency reference) m %i n %i\n",
		argv[0],medium->comm_rank,medium->comm_size,m,n);
      PetscCall(MatDuplicate(Adense, MAT_COPY_VALUES, &AH));
      break;
    case 1:
    case 2:
#ifdef USE_PETSC_HMAT
      /* HTOOLS or H2OPUS */
      PetscCall(MatCreate(PETSC_COMM_WORLD, &AH));
      PetscCall(MatSetSizes(AH, PETSC_DECIDE, PETSC_DECIDE, m, n));  
      if(medium->use_hmatrix==1)
	PetscCall(MatSetType(AH,MATHTOOL));
      else
	PetscCall(MatSetType(AH,MATH2OPUS));
      /*  */
      PetscCall(MatSetUp(AH));
      PetscCall(MatGetLocalSize(AH, &lm, &ln));
      dn = ln;on = n - ln;
      PetscCall(MatSeqAIJSetPreallocation(AH, n, NULL));
      PetscCall(MatMPIAIJSetPreallocation(AH, dn, NULL, on, NULL));
      PetscCall(MatGetOwnershipRange(AH, &rs, &re));
      
      fprintf(stderr,"%s: core %03i/%03i: assigning %s row %5i to %5i, lm %i ln %i m %i n %i\n",
	      argv[0],medium->comm_rank,medium->comm_size,(medium->use_hmatrix==1)?"HTOOLS":"H2OPUS",
	      rs,re,lm,ln,m,n);
      if(medium->use_hmatrix == 1){
	/* 
	   HTOOLS 
	*/
	PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,lm,ln, m, n,
					   ndim,(coords+rs*ndim),
					   (coords+rs*ndim), htools_kernel, ictx, &AH));
	PetscCall(MatSetOption(AH, MAT_SYMMETRIC, PETSC_FALSE));
      }else{
	/* H2opUS */
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
	PetscCall(MatCreateH2OpusFromMat(Adense, ndim, coords, PETSC_FALSE,
					 medium->h2opus_eta, medium->h2opus_leafsize,
					 PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
					 &AH));
      
	/* 
	   this version of h2opus only implements sampling-based
	   construction for symmetric matrices (hlru_sym assertion); the
	   result therefore approximates the symmetrized operator
	   (K+K^T)/2 - check the printed operator asymmetry to judge the
	   error this introduces
	*/
	PetscCall(MatTranspose(Adense,MAT_INITIAL_MATRIX,&KT));
	PetscCall(MatNorm(Adense,NORM_FROBENIUS,&nrmK));
	PetscCall(MatAXPY(KT,-1.0,Adense,SAME_NONZERO_PATTERN));
	PetscCall(MatNorm(KT,NORM_FROBENIUS,&nrmD));
	HEADNODE
	  fprintf(stderr,"%s: operator asymmetry |K-K^T|_F/|K|_F = %.6e (h2opus approximates the symmetrized operator)\n",
		  argv[0],(double)(nrmD/nrmK));
	PetscCall(MatDestroy(&KT));
	PetscCall(MatSetOption(AH, MAT_SYMMETRIC, PETSC_TRUE));
#else
	fprintf(stderr,"%s: H2OPUS requested but PETSc was built without h2opus\n",argv[0]);
	exit(-1);
#endif
#else
	fprintf(stderr,"%s: H2OPUS or HTOOLS requested by not compiled in\n",argv[0]);
	exit(-1);
      }
      break;
    case 3:
      /* 
	 HACApK via MATSHELL 
      */
#ifdef USE_HACAPK
      hacapk_handle = cinit_hacapk_struct((int)m,(void *)ictx);
      cset_hacapk_struct_coord(hacapk_handle,xc,yc,zc);
      fprintf(stderr,"%s: core %03i/%03i: assigning HACApK m %i n %i ztol %g\n",
	      argv[0],medium->comm_rank,medium->comm_size,m,n,(double)medium->hacapk_ztol);
      cmake_hacapk_struct_hmat(hacapk_handle,(double)medium->hacapk_ztol);
      hsc_h = (hacapk_shell_ctx *)malloc(sizeof(hacapk_shell_ctx));
      hsc_h->handle = hacapk_handle;
      hsc_h->ball = (double *)malloc(sizeof(double)*m);
      PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,
			       (void *)hsc_h,&AH));
      PetscCall(MatShellSetOperation(AH,MATOP_MULT,(void (*)(void))MatMult_HACApK));
      PetscCall(MatCreateVecs(AH,&xd,NULL));
      PetscCall(VecScatterCreateToAll(xd,&hsc_h->scat,&hsc_h->xall));
      PetscCall(VecGetOwnershipRange(xd,&hsc_h->rs,&hsc_h->re));
      PetscCall(VecDestroy(&xd));
      PetscCall(MatSetOption(AH, MAT_SYMMETRIC, PETSC_FALSE));
#endif
      break;
    case 4:
      /* 
	 hmmvp via MATSHELL 
      */
#ifdef USE_HMMVP
      
#ifdef USE_HMMVP_MPI
      /* MPI: compress with distributed assembly to a temporary file
	 (collective), load as a distributed MpiHmat; MatMult_hmmvp
	 gathers x to the root and scatters y back */
      HEADNODE  		/* make /tmp filename with PID */
	snprintf(hmmvp_fn,STRLEN,"/tmp/hmmvp_interact_%d.hm",(int)getpid());
      PetscCallMPI(MPI_Bcast(hmmvp_fn,STRLEN,MPI_CHAR,0,PETSC_COMM_WORLD));
      HEADNODE
	fprintf(stderr,"%s: hmmvp MPI compress m %i n %i tol %g eta %g -> %s\n",
		argv[0],m,n,(double)medium->hmmvp_tol,(double)medium->hmmvp_eta,hmmvp_fn);
      cret = chmmvp_compress_to_file((int)m,xc,yc,zc,(double)medium->hmmvp_tol,
				     (double)medium->hmmvp_eta,(void *)ictx,hmmvp_fn);
      if(cret != 0){
	HEADNODE
	  fprintf(stderr,"%s: hmmvp MPI compression failed\n",argv[0]);
	exit(-1);
      }
      hmmvp_handle = chmmvp_mpi_load(hmmvp_fn,medium->hmmvp_nthreads);
      if(!hmmvp_handle){
	HEADNODE
	  fprintf(stderr,"%s: hmmvp MPI load failed\n",argv[0]);
	exit(-1);
      }
      chmmvp_mpi_get_info(hmmvp_handle,&hmm,&hmn,&hmmvp_nnz);
      /* each rank removes its own hmmvp scratch file "<hmmvp_fn>_<rank>"
	 left by the parallel compressor (root concatenates these over MPI but
	 never unlinks them); /tmp is typically node-local, so root cannot do it */
      if(medium->comm_rank > 0){
	snprintf(hmmvp_tmp,STRLEN,"%s_%d",hmmvp_fn,medium->comm_rank);
	remove(hmmvp_tmp);
      }
      HEADNODE{			/* clean up */
	remove(hmmvp_fn);
	fprintf(stderr,"%s: hmmvp(MPI) %i by %i, %ld stored scalars, compression ratio %.5g\n",
		argv[0],hmm,hmn,hmmvp_nnz,(double)((double)m*(double)n/(double)hmmvp_nnz));
      }
      hsc_h = (hacapk_shell_ctx *)malloc(sizeof(hacapk_shell_ctx));
      hsc_h->handle = hmmvp_handle;
      hsc_h->ball = (double *)malloc(sizeof(double)*m); /* full y on every rank */
      PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,
			       (void *)hsc_h,&AH));
      PetscCall(MatShellSetOperation(AH,MATOP_MULT,(void (*)(void))MatMult_hmmvp));
      PetscCall(MatCreateVecs(AH,&xd,NULL));
      PetscCall(VecScatterCreateToZero(xd,&hsc_h->scat,&hsc_h->xall)); /* gather-to-root */
      PetscCall(VecGetOwnershipRange(xd,&hsc_h->rs,&hsc_h->re));
      PetscCall(VecDestroy(&xd));
      PetscCall(MatSetOption(AH, MAT_SYMMETRIC, PETSC_FALSE));
#else
      if((medium->hmmvp_nthreads > 1) && (dc3dts() == 0)){
	fprintf(stderr,"%s: -hmmvp_nthreads %i requested but dc3d.F was compiled WITHOUT\n%s: -fopenmp: the THREADPRIVATE directives are inactive and threaded kernel\n%s: calls would corrupt the matrix - rebuild with -fopenmp in FFLAGS/LDFLAGS\n",
		argv[0],medium->hmmvp_nthreads,argv[0],argv[0]);
	exit(-1);
      }
      /* in-memory OpenMP/serial path */
      fprintf(stderr,"%s: core %03i/%03i: assigning hmmvp m %i n %i tol %g (whole-matrix rel Frobenius) eta %g nthreads %i\n",
	      argv[0],medium->comm_rank,medium->comm_size,m,n,(double)medium->hmmvp_tol,
	      (double)medium->hmmvp_eta,medium->hmmvp_nthreads);
      hmmvp_handle = chmmvp_compress_in_memory((int)m,xc,yc,zc,(double)medium->hmmvp_tol,
					       (double)medium->hmmvp_eta,medium->hmmvp_nthreads,
					       (void *)ictx);
      if(!hmmvp_handle){
	HEADNODE
	  fprintf(stderr,"%s: hmmvp compression failed\n",argv[0]);
	exit(-1);
      }
      chmmvp_get_info(hmmvp_handle,&hmm,&hmn,&hmmvp_nnz);
      HEADNODE
	fprintf(stderr,"%s: hmmvp %i by %i, %ld stored scalars, compression ratio %.5g\n",
		argv[0],hmm,hmn,hmmvp_nnz,
		(double)((double)m*(double)n/(double)hmmvp_nnz));
      hsc_h = (hacapk_shell_ctx *)malloc(sizeof(hacapk_shell_ctx));
      hsc_h->handle = hmmvp_handle;
      hsc_h->ball = (double *)malloc(sizeof(double)*m);
      PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,
			       (void *)hsc_h,&AH));
      PetscCall(MatShellSetOperation(AH,MATOP_MULT,(void (*)(void))MatMult_hmmvp));
      PetscCall(MatCreateVecs(AH,&xd,NULL));
      PetscCall(VecScatterCreateToAll(xd,&hsc_h->scat,&hsc_h->xall));
      PetscCall(VecGetOwnershipRange(xd,&hsc_h->rs,&hsc_h->re));
      PetscCall(VecDestroy(&xd));
      PetscCall(MatSetOption(AH, MAT_SYMMETRIC, PETSC_FALSE));
#endif /* USE_HMMVP_MPI */
#endif /* USE_HMMVP */
      break;
    }
    
    if((medium->use_hmatrix==1)||(medium->use_hmatrix==2))
      free(coords);
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
    if((medium->use_hmatrix==3)||(medium->use_hmatrix==4)){
      free(xc);free(yc);free(zc);
    }
#endif    
    PetscCall(MatSetFromOptions(AH));
    
    PetscCall(MatAssemblyBegin(AH, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(AH, MAT_FINAL_ASSEMBLY));
    PetscTime(&t1);
    HEADNODE
      fprintf(stderr,"%s: H matrix assembly took %12.4f s\n",argv[0],t1-t0);

  }else{
    /* 
       use external routines 
    */
    calc_petsc_Isn_matrices(medium, fault,0,                  1.0,0,&Adense,hsc_dense); /* dense */
    calc_petsc_Isn_matrices(medium, fault,medium->use_hmatrix,1.0,0,&AH,hsc_h); /* Htools, H2opus, HACApK, or HMVVP */
  }

  /* 


     matrix assembly done, from now on, do not have to specify what
     type of matrix we are dealing with
     
  */

  
  /* dense */
  if(n < 20){
    HEADNODE
      fprintf(stderr,"%s: dense matrix:\n",argv[0]);
    PetscCall(MatView(Adense,PETSC_VIEWER_STDOUT_WORLD));
  }
  /* get info on H matrix (skip for the use_hmatrix=0 dense reference,
     where In is a full dense copy of Is and MatView would print every
     entry rather than a compact compression summary) */
  if(medium->use_hmatrix)
    MatView(AH,PETSC_VIEWER_STDOUT_WORLD);
  if(n < 20){
    HEADNODE
      fprintf(stderr,"%s: H matrix converted back to dense:\n",argv[0]);
    PetscCall(MatConvert(AH,MATDENSE,MAT_INITIAL_MATRIX,&AH_dense));
    PetscCall(MatView(AH_dense,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(MatDestroy(&AH_dense));
  }
  /* 
   */
  PetscCall(MatCreateVecs(Adense, &x, &b));/* For A x = b: x -> left, b -> right */
  PetscCall(MatCreateVecs(AH, &xh, &bh));/* For A x = b: x -> left, b -> right */
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
    PetscCall(MatMult(Adense, x, b));
    for(i=0;i<nrandom;i++){
      PetscCall(VecSetRandom(x,rand_str));
      PetscCall(MatMult(Adense, x, b));
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
    PetscCall(MatMult(AH, xh, bh));
    for(i=0;i<nrandom;i++){
      PetscCall(VecSetRandom(xh,rand_str));
      PetscCall(MatMult(AH, xh, bh));
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
      PetscCall(MatMult(Adense, x, b));
      PetscCall(MatMult(AH, xh, bh));
      /* compute difference */
      PetscCall(VecDuplicate(b, &d));
      PetscCall(VecCopy(b, d));
  
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
    PetscCall(KSPSetOperators(ksp, Adense, Adense));
    PetscCall(KSPGetPC(ksp, &pc));
    if(medium->comm_size > 1){
      /*
	 the dense reference matrix is mpidense under MPI, and sequential
	 PCLU cannot factor it; unless PETSc was built with a parallel
	 direct solver (MUMPS or SuperLU_DIST) MatGetFactor() fails with
	 "Could not locate a solver type for factorization type LU and
	 matrix type mpidense". fall back to an iterative reference solve
	 (unpreconditioned GMRES to a tight tolerance) so the inverse
	 accuracy comparison still runs in parallel. configure PETSc with
	 --download-mumps (then this branch can use PCLU/MATSOLVERMUMPS)
	 for a parallel DIRECT reference instead.
      */
      PetscCall(PCSetType(pc, PCNONE));
      PetscCall(KSPSetType(ksp, KSPGMRES));
      PetscCall(KSPSetTolerances(ksp, 1.0e-10, PETSC_DEFAULT, PETSC_DEFAULT, 100000));
    }else{
      PetscCall(PCSetType(pc, PCLU)); /* exact direct reference solve in serial */
    }
    PetscCall(KSPSetFromOptions(ksp)); /* command line still overrides the above */
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
      PetscCall(KSPSetOperators(ksph, AH, AH));
      PetscCall(KSPSetFromOptions(ksph));
      PetscCall(KSPGetPC(ksph, &pch));
      PetscCall(PetscObjectTypeCompare((PetscObject)pch, PCHPDDM, &flg));
      PetscCall(VecSet(bh,1.0));
      PetscCall(VecAssemblyBegin(bh));
      PetscCall(VecAssemblyEnd(bh));
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
	PetscCall(VecDuplicate(x, &d));
	PetscCall(VecCopy(x, d));
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

  PetscCall(MatDestroy(&Adense));
  PetscCall(MatDestroy(&AH));
  PetscCall(PetscFinalize());
#endif
  exit(0);

}
