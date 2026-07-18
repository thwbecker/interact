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
  MatHtoolKernelFn *htools_kernel = GenKEntries_petsc;
  Mat KT;
  PetscReal nrmK,nrmD;
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
  PetscBool read_value,flg,test_forward=PETSC_TRUE,use_full_space=PETSC_FALSE;
  PetscBool make_matrix_externally=PETSC_FALSE; /* make matrices here
						   on in external
						   routine (for
						   testing) */
  /* -skip_dense: build only the H-matrix, time it, and exit before the
     dense reference and the error/solve check. The dense reference is
     m*n entries, which exceeds 32-bit PetscInt and a single rank's memory
     at large N, so this is what lets the assembly be timed on one rank. */
  PetscBool skip_dense=PETSC_FALSE;
  /* -dense_reference_only: build ONLY the dense operator (external path),
     time nrandom matvecs on it, print one machine-readable line, and
     exit. For establishing the dense baseline once at large N (needs
     enough ranks that local_rows*n stays below PetscInt); the H-matrix
     sweep then reads the result instead of rebuilding dense each time. */
  PetscBool dense_reference_only=PETSC_FALSE;
  PetscReal dt0,dt1;
  /* for the unified H-operator matvec timing (runs also with -skip_dense) */
  Vec mvx,mvy;
  PetscRandom mrnd;
  PetscReal mt0,mt1;
  PetscInt ir;
  char geom_file[STRLEN]="geom.in";
  /* optional external dumps (see -dump_matrix / -dump_coords below) */
  char dump_matrix_file[STRLEN]="",dump_coords_file[STRLEN]="";
  PetscBool do_dump_matrix=PETSC_FALSE,do_dump_coords=PETSC_FALSE;
  hmat_helper_shell_ctx *hsc_dense,*hsc_h;
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
  PetscCall(interact_petsc_initialize(&argc, &argv));
  
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
  medium->use_hmatrix =  IHMAT_TYPE_HTOOLS;	/* default is HTOOLS: 0,1,2,3,4,5 are possibly options, 
						   depending on compile, see petsc_prototypes.h */
  
  /* set defaults, can always override */
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-use_hmatrix", &medium->use_hmatrix,&read_value));
  /* forward or inverse test? */
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_forward", &test_forward,&read_value));
  /* default: standard Okada half-space */
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-full_space", &use_full_space,&read_value)); /* use read in or default */
  medium->full_space = (my_boolean)use_full_space;

  /* options */
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nrandom", &nrandom,&read_value));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_forward", &test_forward,&read_value));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-make_matrix_externally", &make_matrix_externally,NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-skip_dense", &skip_dense,NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-dense_reference_only", &dense_reference_only,NULL));
  if(dense_reference_only && skip_dense){
    HEADNODE
      fprintf(stderr,"%s: -dense_reference_only and -skip_dense are mutually exclusive\n",argv[0]);
    exit(-1);
  }
  
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));

  /* propagate the half-space / full-space choice into the Okada kernel
     (dc3d / dc3d0) once, before any assembly. with full_space != 0 the
     Okada-based dense and H-matrix paths return the infinite-medium
     (full-space, real-source) term only. this is read-only during the
     subsequent (possibly threaded) assembly, so the shared flag is safe.
     a native full-space backend (e.g. BigWham, use_hmatrix=5) can then be
     compared against this on the same operator. */
  HEADNODE
    fprintf(stderr,"%s: Okada kernel mode: %s\n",argv[0],
	    medium->full_space ?
	    "FULL-SPACE (infinite medium, free surface suppressed)" :
	    "half-space (standard Okada)");
  
  /* set up defaults */
  if(medium->use_hmatrix)
    set_hmat_defaults_and_options(medium,medium->use_hmatrix);

  
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rand_str));
  PetscCall(PetscRandomSetFromOptions(rand_str));
  /* 
     
     read in geometry

  */
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  /* optional: dump the dense interaction matrix and/or the patch source
     geometry to file for external analysis (e.g. testing alternative
     cluster trees for H-matrix compression). both are off unless a file
     name is given. */
  PetscCall(PetscOptionsGetString(NULL, NULL, "-dump_matrix", dump_matrix_file, STRLEN,&do_dump_matrix));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-dump_coords", dump_coords_file, STRLEN,&do_dump_coords));

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
    /* 
       make all the matrices here in this program, only for testing purposes

     */
    HEADNODE{
      fprintf(stderr,"%s: computing %i by %i matrix LOCALLY\n",argv[0], m,n);
      
    }
    /* 
       dense matrix setup, using Adense
    */
    /* The dense reference is a MATDENSE of m by n. With 32-bit PetscInt a
       local block of local_rows*n entries overflows the index once it
       exceeds PETSC_MAX_INT, which segfaults rather than erroring cleanly.
       Detect that up front from the row split and fall back to the
       H-matrix-only path (as if -skip_dense) with a warning, so the run
       still produces the assembly instead of crashing. Adding MPI ranks
       lowers local_rows. A 64-bit PetscInt build removes the limit but is
       not an option when use_hmatrix selects HTOOL, which does not support
       64-bit indices, so adding ranks (or -skip_dense) is the way out. */
    if(!skip_dense){
      PetscInt local_rows=PETSC_DECIDE,glob=m;
      PetscCall(PetscSplitOwnership(PETSC_COMM_WORLD,&local_rows,&glob));
      if((PetscInt64)local_rows*(PetscInt64)n > (PetscInt64)PETSC_MAX_INT){
	HEADNODE
	  fprintf(stderr,"%s: WARNING: dense reference local block %ld x %i = %lld entries exceeds the PetscInt limit (%lld); skipping the dense reference and the error check. Add MPI ranks so each rank's block stays under the limit, or pass -skip_dense. (A 64-bit PetscInt build removes the limit but is incompatible with HTOOL.)\n",
		  argv[0],(long)local_rows,n,
		  (long long)((PetscInt64)local_rows*(PetscInt64)n),(long long)PETSC_MAX_INT);
	skip_dense = PETSC_TRUE;
      }
    }

    if(!skip_dense){
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
      GenKEntries_petsc(ndim,1,n,&j, col_idx, avalues,ictx);
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
    }else{
      /* -skip_dense: do not build the dense reference. Create a tiny empty
         AIJ matrix with the same global size only to reproduce the row
         distribution (lm, ln, rs, re, rn) that the H-matrix build below
         relies on. The error/solve check is skipped (we exit right after
         the H-matrix assembly is timed). */
      PetscCall(MatCreate(PETSC_COMM_WORLD, &Adense));
      PetscCall(MatSetSizes(Adense, PETSC_DECIDE, PETSC_DECIDE, m, n));
      PetscCall(MatSetType(Adense, MATAIJ));
      PetscCall(MatSeqAIJSetPreallocation(Adense, 0, NULL));
      PetscCall(MatMPIAIJSetPreallocation(Adense, 0, NULL, 0, NULL));
      PetscCall(MatSetUp(Adense));
      PetscCall(MatGetLocalSize(Adense, &lm, &ln));
      PetscCall(MatGetOwnershipRange(Adense, &medium->rs, &medium->re));
      medium->rn = medium->re - medium->rs;
      HEADNODE
	fprintf(stderr,"%s: -skip_dense: no dense reference (H-matrix assembly timing only)\n",argv[0]);
    }

    /*
       optional external dumps of the operator and its point cloud, so
       that alternative cluster trees / admissibility choices can be
       explored outside interact (e.g. fault-split vs joint geometric
       clustering for H-matrix compression).

       -dump_matrix <file>: the dense interaction matrix as raw
          row-major float64, A[i*n+j] = stress at receiver i from unit
          slip at source j (the same operator the H-matrix backends
          approximate). a companion <file>.info records "m n" and the
          layout. needs a single MPI rank so the full matrix is local.

       -dump_coords <file>: one ASCII row per patch with centroid,
          orientation, half-sizes, area, unit normal and group id, i.e.
          everything a clustering routine needs, including the normal so
          that orientation (not just centroid position) can be used.
    */
    if(do_dump_matrix){
      if(medium->comm_size != 1){
	HEADNODE
	  fprintf(stderr,"%s: -dump_matrix needs a single MPI rank (rerun with -np 1); skipping matrix dump\n",argv[0]);
      }else{
	print_petsc_matrix(Adense,n,m,dump_matrix_file);
      }
    }
    if(do_dump_coords)
      HEADNODE
	print_fault_geometry_and_normals(fault,medium->nrflt,dump_coords_file);
    
    /* 
       
       ASSEMBLE DIFFERENT KINDS OF H MATRICES
       
    */
    /* 
       prepare with coordinates and such  
    */
#ifdef USE_PETSC_HMAT
    if((medium->use_hmatrix==IHMAT_TYPE_HTOOLS)||(medium->use_hmatrix==IHMAT_TYPE_H2OPUS)){	
      coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*medium->nrflt);
      for(i=0;i < medium->nrflt;i++)		/* all sources or receiveer coordinates  */
	for(k=0;k < ndim;k++)
	  coords[i*ndim+k] = fault[i].x[k];
    }
#endif
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
    if((medium->use_hmatrix==IHMAT_TYPE_HACAPK)||(medium->use_hmatrix==IHMAT_TYPE_HMMVP)){	
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
    case IHMAT_TYPE_DENSE:
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
    case IHMAT_TYPE_HTOOLS:
    case IHMAT_TYPE_H2OPUS:
#ifdef USE_PETSC_HMAT
      /* HTOOLS or H2OPUS */
      PetscCall(MatCreate(PETSC_COMM_WORLD, &AH));
      PetscCall(MatSetSizes(AH, PETSC_DECIDE, PETSC_DECIDE, m, n));  
      if(medium->use_hmatrix==IHMAT_TYPE_HTOOLS)
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
	      argv[0],medium->comm_rank,medium->comm_size,(medium->use_hmatrix==IHMAT_TYPE_HTOOLS)?"HTOOLS":"H2OPUS",
	      rs,re,lm,ln,m,n);
      if(medium->use_hmatrix == IHMAT_TYPE_HTOOLS){
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
      }
#endif
      break;
    case IHMAT_TYPE_HACAPK:
      /* 
	 HACApK via MATSHELL 
      */
#ifdef USE_HACAPK
      hacapk_handle = cinit_hacapk_struct((int)m,(void *)ictx);
      cset_hacapk_struct_coord(hacapk_handle,xc,yc,zc);
      cset_hacapk_eta(hacapk_handle,(double)medium->hacapk_eta); /* override param(51) before the build, eta */
      cset_hacapk_inorm(hacapk_handle,medium->hacapk_inorm);     /* error norm mode */
      
      fprintf(stderr,"%s: core %03i/%03i: assigning HACApK m %i n %i ztol %g eta %g inorm: %i\n",
	      argv[0],medium->comm_rank,medium->comm_size,m,n,(double)medium->hacapk_ztol,
	      (double)medium->hacapk_eta,medium->hacapk_inorm);

      cmake_hacapk_struct_hmat(hacapk_handle,(double)medium->hacapk_ztol);

      hsc_h = (hmat_helper_shell_ctx *)malloc(sizeof(hmat_helper_shell_ctx));
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
    case IHMAT_TYPE_HMMVP:
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
				     (double)medium->hmmvp_eta,medium->hmmvp_inorm,(void *)ictx,hmmvp_fn);
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
      hsc_h = (hmat_helper_shell_ctx *)malloc(sizeof(hmat_helper_shell_ctx));
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
					       (double)medium->hmmvp_eta,medium->hmmvp_inorm,medium->hmmvp_nthreads,
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
      hsc_h = (hmat_helper_shell_ctx *)malloc(sizeof(hmat_helper_shell_ctx));
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
    case IHMAT_TYPE_BIGWHAM:
#ifdef USE_BIGWHAM
      /* BigWham full-space H matrix as a MATSHELL; builds its own mesh from
	 the patch geometry and applies the strike-slip -> strike-shear
	 sub-block of the 3N x 3N operator (see setup_bigwham_matshell). */
      PetscCall(setup_bigwham_matshell(medium,fault,1.0,0,&AH,&hsc_h));
      PetscCall(MatSetOption(AH, MAT_SYMMETRIC, PETSC_FALSE));
#else
      fprintf(stderr,"%s: BigWham requested but not compiled in\n",argv[0]);
      exit(-1);
#endif
      break;
    }
    
    if((medium->use_hmatrix==IHMAT_TYPE_HTOOLS)||(medium->use_hmatrix==IHMAT_TYPE_H2OPUS))
      free(coords);
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
    if((medium->use_hmatrix==IHMAT_TYPE_HACAPK )||(medium->use_hmatrix==IHMAT_TYPE_HMMVP)){
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
       use external routines. honor -skip_dense here as well, and apply
       the same 32-bit dense-block overflow guard as the local branch:
       the dense reference preallocates local_rows*n AIJ entries inside
       calc_petsc_Isn_matrices, which overflows PetscInt at large N
       (e.g. 265k patches on 24 ranks) and previously aborted the run
       before the H matrix was ever built.
    */
    if(!skip_dense){
      PetscInt local_rows=PETSC_DECIDE,glob=m;
      PetscCall(PetscSplitOwnership(PETSC_COMM_WORLD,&local_rows,&glob));
      if((PetscInt64)local_rows*(PetscInt64)n > (PetscInt64)PETSC_MAX_INT){
	if(dense_reference_only){
	  /* the user explicitly asked for dense; do not silently skip */
	  HEADNODE
	    fprintf(stderr,"%s: ERROR: -dense_reference_only, but the dense local block %ld x %i = %lld entries exceeds the PetscInt limit (%lld). Increase the MPI rank count so each rank's block stays under the limit.\n",
		    argv[0],(long)local_rows,n,
		    (long long)((PetscInt64)local_rows*(PetscInt64)n),(long long)PETSC_MAX_INT);
	  exit(-1);
	}
	HEADNODE
	  fprintf(stderr,"%s: WARNING: dense reference local block %ld x %i = %lld entries exceeds the PetscInt limit (%lld); proceeding as if -skip_dense was set. Add MPI ranks so each rank's block stays under the limit if the dense reference and error check are needed.\n",
		  argv[0],(long)local_rows,n,
		  (long long)((PetscInt64)local_rows*(PetscInt64)n),(long long)PETSC_MAX_INT);
	skip_dense = PETSC_TRUE;
      }
    }
    if(!skip_dense){
      PetscCall(PetscTime(&dt0));
      calc_petsc_Isn_matrices(medium, fault,0,                  1.0,0,STRIKE,&Adense,hsc_dense); /* dense */
      PetscCall(PetscTime(&dt1));
      if(dense_reference_only){
	/* time matvecs on the dense operator and leave; one line holds
	   everything the sweep script needs */
	if(nrandom > 0){
	  PetscCall(MatCreateVecs(Adense,&mvx,&mvy));
	  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&mrnd));
	  PetscCall(VecSetRandom(mvx,mrnd));
	  PetscCall(MatMult(Adense,mvx,mvy)); /* warm up */
	  PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
	  PetscCall(PetscTime(&mt0));
	  for(ir=0;ir < nrandom;ir++)
	    PetscCall(MatMult(Adense,mvx,mvy));
	  PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
	  PetscCall(PetscTime(&mt1));
	  HEADNODE
	    fprintf(stderr,"%s: dense_reference m %i n %i assembly_s %.3f applies %i per_matvec_ms %.4f\n",
		    argv[0],m,n,(double)(dt1-dt0),(int)nrandom,
		    1000.0*(double)(mt1-mt0)/(double)nrandom);
	  PetscCall(PetscRandomDestroy(&mrnd));
	  PetscCall(VecDestroy(&mvx));
	  PetscCall(VecDestroy(&mvy));
	}else{
	  HEADNODE
	    fprintf(stderr,"%s: dense_reference m %i n %i assembly_s %.3f applies 0 per_matvec_ms 0\n",
		    argv[0],m,n,(double)(dt1-dt0));
	}
	HEADNODE
	  fprintf(stderr,"%s: -dense_reference_only: exiting after dense assembly and matvec timing\n",argv[0]);
	PetscCall(MatDestroy(&Adense));
	PetscCall(PetscFinalize());
	return 0;
      }
    }else{
      /* row-distribution placeholder instead of the dense reference,
         as in the local -skip_dense path */
      PetscCall(MatCreate(PETSC_COMM_WORLD, &Adense));
      PetscCall(MatSetSizes(Adense, PETSC_DECIDE, PETSC_DECIDE, m, n));
      PetscCall(MatSetType(Adense, MATAIJ));
      PetscCall(MatSeqAIJSetPreallocation(Adense, 0, NULL));
      PetscCall(MatMPIAIJSetPreallocation(Adense, 0, NULL, 0, NULL));
      PetscCall(MatSetUp(Adense));
      PetscCall(MatGetLocalSize(Adense, &lm, &ln));
      PetscCall(MatGetOwnershipRange(Adense, &medium->rs, &medium->re));
      medium->rn = medium->re - medium->rs;
      HEADNODE
	fprintf(stderr,"%s: -skip_dense: no dense reference (H-matrix assembly and matvec timing only)\n",argv[0]);
    }
    calc_petsc_Isn_matrices(medium, fault,medium->use_hmatrix,1.0,0,STRIKE,&AH,hsc_h); /* Htools, H2opus, HACApK, or HMVVP */
  }

  /*
     unified H-operator report: for HTOOL print the operator info (which
     includes the compression ratio; the other backends print their
     storage on their own assembly lines), then time nrandom matvecs.
     This block runs also under -skip_dense, which exits right after,
     since everything below needs the dense reference.
  */
  if(medium->use_hmatrix==IHMAT_TYPE_HTOOLS)
    PetscCall(MatView(AH,PETSC_VIEWER_STDOUT_WORLD));
  {
    if(nrandom > 0){
      PetscCall(MatCreateVecs(AH,&mvx,&mvy));
      PetscCall(PetscRandomCreate(PETSC_COMM_WORLD,&mrnd));
      PetscCall(VecSetRandom(mvx,mrnd));
      PetscCall(MatMult(AH,mvx,mvy)); /* warm up */
      PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
      PetscCall(PetscTime(&mt0));
      for(ir=0;ir < nrandom;ir++)
	PetscCall(MatMult(AH,mvx,mvy));
      PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
      PetscCall(PetscTime(&mt1));
      HEADNODE
	fprintf(stderr,"%s: hmat_matvec backend %i m %i n %i applies %i total %.4f s per_matvec %.4f ms\n",
		argv[0],(int)medium->use_hmatrix,m,n,(int)nrandom,
		(double)(mt1-mt0),1000.0*(double)(mt1-mt0)/(double)nrandom);
      PetscCall(PetscRandomDestroy(&mrnd));
      PetscCall(VecDestroy(&mvx));
      PetscCall(VecDestroy(&mvy));
    }
  }
  if(skip_dense){
    HEADNODE
      fprintf(stderr,"%s: -skip_dense: exiting after H-matrix assembly and matvec timing (no error check)\n",argv[0]);
    PetscCall(MatDestroy(&AH));
    PetscCall(MatDestroy(&Adense));
    PetscCall(PetscFinalize());
    return 0;
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
      /* generic accuracy: relative matvec error on a random x, an estimate
         of the operator error ||A_h-A||/||A|| in a typical direction, the
         fundamental b = A x quality decoupled from the x=1 coherent/loading
         direction reported just below. one shared random x feeds both
         operators (x to dense, xh to H, mirroring the x=1 path), so the
         difference is meaningful; this assumes the dense and H column layouts
         match, which they do here (default PetscLayout over n on one comm). */
      {
        Vec        dr;
        PetscReal  nb, nd, re, sum_re = 0.0, max_re = 0.0;
        PetscInt   k, nacc;
        nacc = (nrandom < 1) ? 1 : ((nrandom > 20) ? 20 : nrandom);
        PetscCall(VecDuplicate(b, &dr));
        for(k=0;k<nacc;k++){
          PetscCall(VecSetRandom(x, rand_str));
          PetscCall(VecCopy(x, xh));
          PetscCall(MatMult(Adense, x,  b));
          PetscCall(MatMult(AH,     xh, bh));
          PetscCall(VecCopy(b, dr));
          PetscCall(VecAXPY(dr, -1.0, bh));
          PetscCall(VecNorm(b,  NORM_2, &nb));
          PetscCall(VecNorm(dr, NORM_2, &nd));
          re = (nb > 0.0) ? nd/nb : 0.0;
          sum_re += re; if(re > max_re) max_re = re;
        }
        PetscCall(VecDestroy(&dr));
        HEADNODE
          fprintf(stdout,"%s: random-x rel err: mean = %20.10e max = %20.10e over %d vectors\n",
                  argv[0], (double)(sum_re/nacc), (double)max_re, (int)nacc);
      }
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

