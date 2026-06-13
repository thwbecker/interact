/*
  interact: model fault interactions using dislocations in a 
            halfspace

	    (C) Thorsten Becker, thwbecker@post.harvard.edu

  Okada and other dislocation routines based elastic half-space
  fault-patch interaction program

  calculates the stresses on fault patches due to slip on other
  faults. inversion of the interaction matrix is done using various
  matrix solves, see 'interact -h' for description




*/
#include "interact.h"
#ifdef USE_PETSC
#include "petsc_prototypes.h"


#ifdef USE_HMMVP
/* 
   hmmvp (A.M. Bradley) support via hmmvp_c_shim.cpp: index-based
   block kernel (same ckernel_func as HACApK), in-memory compression
   with a WHOLE-MATRIX relative Frobenius tolerance -hmmvp_tol
   (||B-A||_F <= tol ||B||_F, i.e. tol bounds exactly the error this
   tool measures), OpenMP-threaded construction and matvec. wrapped
   as a MATSHELL like HACApK; the Mvp needs global vectors, so the
   same scatter-to-all context is used (at np>1 every rank holds the
   full H matrix and computes the full product - correct but
   redundant; intended for serial/threaded use)
*/
PetscErrorCode MatMult_hmmvp(Mat A, Vec x, Vec y)
{
  hacapk_shell_ctx *hctx;	/* same context layout as HACApK */
  const PetscScalar *xa;
  PetscScalar *ya;
  PetscInt i;
  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&hctx));
  PetscCall(VecScatterBegin(hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecScatterEnd(hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecGetArrayRead(hctx->xall,&xa));
  chmmvp_mvp(hctx->handle,(double *)xa,hctx->ball);
  PetscCall(VecRestoreArrayRead(hctx->xall,&xa));
  PetscCall(VecGetArray(y,&ya));
  for(i=hctx->rs;i < hctx->re;i++)
    ya[i-hctx->rs] = (PetscScalar)hctx->ball[i];
  PetscCall(VecRestoreArray(y,&ya));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif	/* end HMM */
#ifdef USE_HACAPK
/* MATSHELL multiply: y = A x through the HACApK H matrix */
/* 
   HACApK's adot expects the GLOBAL x vector on each rank and leaves
   the GLOBAL result on each rank (ring exchange + permutation inside
   HACApK_adot_pmt_lfmtx_p), while PETSc vectors are distributed:
   scatter x to a full-length local copy first, then copy back the
   locally owned slice of the result
*/
PetscErrorCode MatMult_HACApK(Mat A, Vec x, Vec y)
{
  hacapk_shell_ctx *hctx;
  const PetscScalar *xa;
  PetscScalar *ya;
  PetscInt i;
  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&hctx));
  PetscCall(VecScatterBegin(hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecScatterEnd(hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecGetArrayRead(hctx->xall,&xa));
  chacapk_mult_Ax_H(hctx->handle,(double *)xa,hctx->ball);
  PetscCall(VecRestoreArrayRead(hctx->xall,&xa));
  PetscCall(VecGetArray(y,&ya));
  for(i=hctx->rs;i < hctx->re;i++)
    ya[i-hctx->rs] = (PetscScalar)hctx->ball[i];
  PetscCall(VecRestoreArray(y,&ya));
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif	/* end Haca */

/* 
   
   compute fixed slip interaction matrices - this is an isolated
   assembly, compare compress_interaction_matrix for the debugging type version

   mode 0: shear stress
   mode 1: normal stress

   scale: scale factor 

   hmatrix can be 0 (dense), 1 HTOOLS or 2 HOPUS, 3 or 4 
*/

PetscErrorCode calc_petsc_Isn_matrices(struct med *medium, struct flt *fault,
				       PetscInt use_hmatrix,PetscReal scale, int mode,
				       Mat *this_mat, hacapk_shell_ctx *hctx)
{
  /* context */
  struct interact_ctx ictx[1];
  PetscReal   *avalues=NULL;
  PetscInt    n, m, lm,ln,i,j,dn,on, *col_idx=NULL;
  /* kernel function */
#ifdef USE_PETSC_HMAT		/* htools and H2opus  */
  PetscReal   *coords=NULL,*av=NULL;
  PetscInt k,jj,*ci=NULL;
  MatHtoolKernelFn *kernel = GenKEntries_htools;
  Mat dtmp;
#endif
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
  double *xc,*yc,*zc;
  Vec xd;
#endif
#ifdef USE_HMMVP
  void *hmmvp_handle;
  long hmmvp_nnz;
  int hmm,hmn;
#endif 
#ifdef USE_HACAPK
  void *hacapk_handle;
#endif
  const PetscInt ndim = 3;
  ictx->medium = medium;
  ictx->fault = fault;
  /* defines how to slip */
  ictx->src_slip_mode = 0;
  /*  */
  m = n = medium->nrflt;

  if(mode==0){
    ictx->rec_stress_mode = STRIKE;
  }else{
    ictx->rec_stress_mode = NORMAL;
  }
  switch(use_hmatrix){
  case 0:
    /* dense */
    PetscCall(PetscCalloc(m*sizeof(PetscScalar), &avalues));
    PetscCall(PetscCalloc(n*sizeof(PetscInt), &col_idx));
    for (i=0; i < n; i++) 
      col_idx[i] = i;
    break;
  case 1:	/* HTOOLS */
#ifdef USE_PETSC_HMAT
    set_htools_defaults_and_options(medium);
#else
    fprintf(stderr,"calc_petsc_Isn_matrices: HTOOLS not compiled in - check makefile.petsc\n");exit(-1);
#endif
    break;
  case 2:	/* H2OPUS */
#ifdef USE_PETSC_HMAT
    set_h2opus_defaults_and_options(medium);
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: WARNING: construction assumes a SYMMETRIC operator, see compress_interaction_matrix)\n");
#else
    fprintf(stderr,"calc_petsc_Isn_matrices: H2OPUS not compiled in - check makefile.petsc\n");exit(-1);
#endif
    break;
  case 3:
#ifdef USE_HACAPK
    set_hacapk_defaults_and_options(medium);
#else
    fprintf(stderr,"HACAPK requested but not compiled in (see USE_HACAPK and makefile.petc)\n");
    exit(-1);
#endif
    break;
  case 4:
#ifdef USE_HMMVP
    set_hmmvp_defaults_and_options(medium);
#else
    fprintf(stderr,"HMMVP requested but not compiled in (see USE_HMMVP and makefile.petc)\n");
    exit(-1);
#endif
    break;
  default:
    fprintf(stderr,"calc_petsc_Isn_matrices: hmat mode %i is undefined\n",use_hmatrix);
    exit(-1);
    break;
  }
#ifdef USE_PETSC_HMAT
  if((use_hmatrix==1)||(use_hmatrix==2)){	
    coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*medium->nrflt);
    for(i=0;i < medium->nrflt;i++)		/* all sources or receiveer coordinates  */
      for(k=0;k < ndim;k++)
	coords[i*ndim+k] = fault[i].x[k];
  }
#endif
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
  
  if((use_hmatrix==3)||(use_hmatrix==4)){	
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
 

  PetscCall(MatCreate(PETSC_COMM_WORLD, this_mat));
  PetscCall(MatSetSizes(*this_mat, PETSC_DECIDE, PETSC_DECIDE, m, n));
 

  PetscCall(MatSetUp(*this_mat));
  PetscCall(MatGetLocalSize(*this_mat, &lm, &ln));
  dn = ln;on = n - ln;
  
  PetscCall(MatSeqAIJSetPreallocation(*this_mat, n, NULL));
  PetscCall(MatMPIAIJSetPreallocation(*this_mat, dn, NULL, on, NULL));
  PetscCall(MatGetOwnershipRange(*this_mat, &medium->rs, &medium->re));
  medium->rn = medium->re  - medium->rs; /* number of local elements */

  switch(use_hmatrix){
  case 0:
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: using dense matrix for stress mode %i stress type %i\n",
	      mode,ictx->rec_stress_mode);
    PetscCall(MatSetType(*this_mat, MATDENSE));
    /* 
       assemble dense matrix 
    */
    fprintf(stderr,"calc_petsc_Isn_matrices: core %03i/%03i: assigning dense row %5i to %5i\n",
	    (int)medium->comm_rank,(int)medium->comm_size,(int)medium->rs,(int)medium->re);
    for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
      GenKEntries_htools(ndim,1,n,&j, col_idx, avalues,ictx);
      PetscCall(MatSetValues(*this_mat, 1, &j, n, col_idx,avalues, INSERT_VALUES));
    }
    break;
#ifdef USE_PETSC_HMAT
  case 1:
    /* make the matrix */
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: using HTOOL matrix for stress mode %i stress type %i\n",
	      mode,ictx->rec_stress_mode);
    PetscCall(MatSetType(*this_mat, MATHTOOL));
    /* 
       
       H matrix setup HTOOLS
       
    */
    fprintf(stderr,"calc_petsc_Isn_matrices: core %03i/%03i: assigning HTOOL  row %5i to %5i, lm %i ln %i m %i n %i\n",
	    medium->comm_rank,medium->comm_size,medium->rs,medium->re,lm,ln,m,n);
    /* target (row) and source (column) coordinates are the LOCAL
       portions, each ndim entries per point (for our square matrices,
       row and column layouts coincide) */
    PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,lm,ln, m, n,ndim,
				       (coords+medium->rs*ndim), (coords+medium->rs*ndim),
				       kernel,(void *)ictx, this_mat));
    break;
  case 2:
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: using H2OPUS matrix for stress mode %i stress type %i\n",
	      mode,ictx->rec_stress_mode);
    PetscCall(MatSetType(*this_mat, MATH2OPUS));
    /* 
       
       H matrix setup H2OPUS
       
    */
    fprintf(stderr,"calc_petsc_Isn_matrices: core %03i/%03i: assigning H2OPUS row %5i to %5i, lm %i ln %i m %i n %i\n",
	    medium->comm_rank,medium->comm_size,medium->rs,medium->re,lm,ln,m,n);

#if defined(PETSC_HAVE_H2OPUS)
    /*
       construct the H2 matrix by hierarchical randomized sampling of an
       assembled dense operator (MatCreateH2OpusFromMat / HARA), NOT from
       the kernel callback. MatCreateH2OpusFromKernel interpolates the
       kernel by Chebyshev polynomials at arbitrary points inside cluster
       bounding boxes, which is incompatible with patch-pair Green's
       functions defined only at element centers: the adaptive
       construction never reaches tolerance and hangs. This mirrors the
       inline path in compress_interaction_matrix.c - a temporary dense
       operator (same local layout) is built, sampled, and discarded.
    */
      
    PetscCall(MatDestroy(this_mat)); /* drop the placeholder created above */

    PetscCall(PetscCalloc(n*sizeof(PetscScalar), &av));
    PetscCall(PetscCalloc(n*sizeof(PetscInt), &ci));
    for(jj=0;jj < n;jj++) ci[jj] = jj;
    /* make a dens matrix */
    PetscCall(MatCreate(PETSC_COMM_WORLD,&dtmp));
    PetscCall(MatSetSizes(dtmp,lm,ln,m,n));
    PetscCall(MatSetType(dtmp,MATDENSE));
    PetscCall(MatSetUp(dtmp));
    for(jj=medium->rs;jj < medium->re;jj++){
      GenKEntries_htools(ndim,1,n,&jj,ci,av,ictx);
      PetscCall(MatSetValues(dtmp,1,&jj,n,ci,av,INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(dtmp,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(dtmp,MAT_FINAL_ASSEMBLY));
    /* dense finished, now convert */
    PetscCall(MatCreateH2OpusFromMat(dtmp, ndim, coords, PETSC_FALSE,
				     medium->h2opus_eta, medium->h2opus_leafsize,
				     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, this_mat));
    PetscCall(MatDestroy(&dtmp));
    PetscCall(PetscFree(av));
    PetscCall(PetscFree(ci));

#else
    fprintf(stderr,"calc_petsc_Isn_matrices: H2OPUS requested but PETSc was built without h2opus\n");
    exit(-1);
#endif
    break;
#endif
  case 3:
#ifdef USE_HACAPK
    hacapk_handle = cinit_hacapk_struct((int)m,(void *)ictx);
    cset_hacapk_struct_coord(hacapk_handle,xc,yc,zc);
    fprintf(stderr,"core %03i/%03i: assigning HACApK m %i n %i ztol %g\n",
	    medium->comm_rank,medium->comm_size,m,n,(double)medium->hacapk_ztol);
    cmake_hacapk_struct_hmat(hacapk_handle,(double)medium->hacapk_ztol);
    hctx = (hacapk_shell_ctx *)malloc(sizeof(hacapk_shell_ctx));
    hctx->handle = hacapk_handle;
    hctx->ball = (double *)malloc(sizeof(double)*m);
    PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,(void *)hctx,this_mat));
    PetscCall(MatShellSetOperation(*this_mat,MATOP_MULT,(void (*)(void))MatMult_HACApK));
    PetscCall(MatCreateVecs(*this_mat,&xd,NULL));
    PetscCall(VecScatterCreateToAll(xd,&hctx->scat,&hctx->xall));
    PetscCall(VecGetOwnershipRange(xd,&hctx->rs,&hctx->re));
    PetscCall(VecDestroy(&xd));
#endif
    break;
  case 4:
#ifdef USE_HMMVP
    hmmvp_handle = chmmvp_compress_in_memory((int)m,xc,yc,zc,(double)medium->hmmvp_tol,
					     (double)medium->hmmvp_eta,medium->hmmvp_nthreads,
					       (void *)ictx);
    if(!hmmvp_handle){
      fprintf(stderr,"hmmvp compression failed\n");
      exit(-1);
    }
    chmmvp_get_info(hmmvp_handle,&hmm,&hmn,&hmmvp_nnz);
    HEADNODE
      fprintf(stderr,"hmmvp %i by %i, %ld stored scalars, compression ratio %.5g\n",
	      hmm,hmn,hmmvp_nnz,
	      (double)((double)m*(double)n/(double)hmmvp_nnz));
    
    hctx = (hacapk_shell_ctx *)malloc(sizeof(hacapk_shell_ctx));
    hctx->handle = hmmvp_handle;
    hctx->ball = (double *)malloc(sizeof(double)*m);
    PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,
			     (void *)hctx,this_mat));
    PetscCall(MatShellSetOperation(*this_mat,MATOP_MULT,(void (*)(void))MatMult_hmmvp));
    PetscCall(MatCreateVecs(*this_mat,&xd,NULL));
    PetscCall(VecScatterCreateToAll(xd,&hctx->scat,&hctx->xall));
    PetscCall(VecGetOwnershipRange(xd,&hctx->rs,&hctx->re));
    PetscCall(VecDestroy(&xd));
#endif
    break;
  }
  if(use_hmatrix == 2){
    /* 
       this version of h2opus only implements sampling-based
       construction for symmetric matrices (hlru_sym assertion); the
       result therefore approximates the symmetrized operator
    */
    HEADNODE
      fprintf(stderr,"WARNING: H2OPUS only uses symmetric matrices approximation\n");
    PetscCall(MatSetOption(*this_mat, MAT_SYMMETRIC, PETSC_TRUE));
  }else{
    PetscCall(MatSetOption(*this_mat, MAT_SYMMETRIC, PETSC_FALSE));
  }

  
  PetscCall(MatSetFromOptions(*this_mat)); /* needed? */
  PetscCall(MatAssemblyBegin(*this_mat, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(*this_mat, MAT_FINAL_ASSEMBLY));
  /* 
     done, now scale (CHECK IF IMPLEMENTED)
  */
  PetscCall(MatScale(*this_mat,scale));

  /* free things */
  if(use_hmatrix==0){
    PetscCall(PetscFree(avalues));
    PetscCall(PetscFree(col_idx));
  }
  if((use_hmatrix==1) || (use_hmatrix==2))
    free(coords);
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
  if((use_hmatrix==3)||(use_hmatrix==4)){
    free(xc);free(yc);free(zc);
  }
  if(use_hmatrix==3)
    PetscCall(VecDestroy(&xd));
#endif
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}

#ifdef USE_PETSC_HMAT		/* htools and H2opus  */

/* 
   
   the approach with the KDtree never worked for H2OPUS, leave here for now

 */
void setup_kdtree(int ndim,int m,PetscReal *coords,struct flt *fault,struct med *medium)
{
  PetscLogDouble t0,t1;
  int i,j,k;
  kd_node    nearest_x;
  double     target_x[3], sep_x;
  /* H2OPUS setup */
  PetscPrintf(PETSC_COMM_WORLD,"kdtree setup for H2OPUS %i nodes\n",m);
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
  PetscPrintf(PETSC_COMM_WORLD,"kdtree setup for H2OPUS done: %1.4e sec (%i)\n",t1-t0,medium->nrflt);

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

}

PetscErrorCode set_h2opus_defaults_and_options(struct med *medium)
{
  static my_boolean init=FALSE;
  if(!init){
    /* defaults for H2OPUS */
    medium->h2opus_eta = 0.6;	/*  */
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-eta", &medium->h2opus_eta, NULL));
    medium->h2opus_leafsize = 32;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-leafsize", &medium->h2opus_leafsize, NULL));
    medium->h2opus_basisord = 8;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-basisord", &medium->h2opus_basisord, NULL));
  }
  init = TRUE;
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}

PetscErrorCode set_htools_defaults_and_options(struct med *medium)
{
  static my_boolean init = FALSE;
  PetscBool flg;
  if(!init){
    /* defaults, only applied if not given on the command line */
    PetscCall(PetscOptionsHasName(NULL,NULL,"-mat_htool_eta",&flg));
    if(!flg)
      PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_eta","100")); 
    PetscCall(PetscOptionsHasName(NULL,NULL,"-mat_htool_epsilon",&flg));
    if(!flg)
      PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_epsilon","1e-6"));
    PetscCall(PetscOptionsHasName(NULL,NULL,"-mat_htool_compressor",&flg));
    if(!flg)
      PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_compressor","SVD"));
    PetscCall(PetscOptionsHasName(NULL,NULL,"-pc_type",&flg));
    if(!flg)
      PetscCall(PetscOptionsSetValue(NULL,"-pc_type","none"));
  }
  init = TRUE;
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}
#endif
#ifdef USE_HMMVP
PetscErrorCode set_hmmvp_defaults_and_options(struct med *medium)
{
  static my_boolean init=FALSE;
  if(!init){
    medium->hmmvp_tol = 1.0e-5;
    medium->hmmvp_eta = 3.0;
    medium->hmmvp_nthreads = 1;
    
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-hmmvp_tol",&medium->hmmvp_tol,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-hmmvp_eta",&medium->hmmvp_eta,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-hmmvp_nthreads",&medium->hmmvp_nthreads,NULL));
  }
  init = TRUE;
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}
#endif
#ifdef USE_HACAPK
PetscErrorCode set_hacapk_defaults_and_options(struct med *medium)
{
  static my_boolean init=FALSE;
  if(!init){
    medium->hacapk_ztol = 1.0e-4;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-hacapk_ztol",&medium->hacapk_ztol,NULL));
  }
  init=TRUE;
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}
#endif



#endif	/* end use Petsc */
