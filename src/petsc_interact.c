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
#ifdef USE_BIGWHAM
#include "properties.h"		/* YOUNG_MODULUS, POISSON_NU for the BigWham elastic constants */
#endif

#ifdef USE_PETSC
#include "petsc_prototypes.h"

/*
   interact_petsc_initialize: single, shared PETSc startup for the interact
   tools (interact, compress_interaction_matrix, rsf_solve).

   Locates the optional PETSc options file
   $HOME/progs/src/interact/petsc_settings.yaml and, if it exists, passes it to
   PetscInitialize so its entries act as defaults (command-line options still
   override). On rank 0 it notes which file was used. This centralizes what
   rsf_solve previously did inline so every PETSc-using interact program starts
   up identically. medium is intentionally not touched here: callers still set
   comm_size / comm_rank themselves after this returns, exactly as before.
*/
PetscErrorCode interact_petsc_initialize(int *argc, char ***argv)
{
  char par_file[STRLEN];
  char *home_dir = getenv("HOME");
  FILE *tst;
  PetscBool have_file = PETSC_FALSE;
  /* only read the YAML defaults file if it exists */
  snprintf(par_file,STRLEN,"%s/progs/src/interact/petsc_settings.yaml",
	   (home_dir)?(home_dir):("."));
  tst = fopen(par_file,"r");
  if(tst){
    have_file = PETSC_TRUE;
    fclose(tst);
  }
  PetscCall(PetscInitialize(argc,argv,(have_file)?(par_file):(NULL),NULL));
  if(have_file)
    PetscCall(PetscFPrintf(PETSC_COMM_WORLD,stderr,
			   "%s: found and using Petsc options in %s\n",
			   (argv && (*argv) && (*argv)[0])?((*argv)[0]):("interact"),
			   par_file));
  return PETSC_SUCCESS;
}


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
#ifdef USE_HMMVP_MPI
  PetscMPIInt rank;
  PetscInt N;
  /*
    MPI path: hmmvp::MpiHmat::Mvp needs the full x on the root and
    leaves the full y valid on the root, with the compute distributed
    across ranks. hctx->scat is a VecScatterToZero (xall is the
    full-length sequential vector on rank 0, empty elsewhere); ball is
    the full-length y work buffer that hmmvp requires on every rank.
  */
  
  PetscCallMPI(MPI_Comm_rank(PetscObjectComm((PetscObject)A),&rank));
  /* gather distributed x -> full x on the root */
  PetscCall(VecScatterBegin(hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecScatterEnd(hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  /* collective matvec: root supplies x, every rank computes its
     blocks, y (=ball) valid on the root afterwards */
  PetscCall(VecGetArrayRead(hctx->xall,&xa));
  chmmvp_mpi_mvp(hctx->handle,(rank==0)?(double *)xa:NULL,hctx->ball);
  PetscCall(VecRestoreArrayRead(hctx->xall,&xa));
  /* copy full y back into the root sequential vector, scatter out */
  if(rank == 0){
    PetscCall(VecGetLocalSize(hctx->xall,&N));
    PetscCall(VecGetArray(hctx->xall,&ya));
    for(i=0;i < N;i++)
      ya[i] = (PetscScalar)hctx->ball[i];
    PetscCall(VecRestoreArray(hctx->xall,&ya));
  }
  PetscCall(VecScatterBegin(hctx->scat,hctx->xall,y,INSERT_VALUES,SCATTER_REVERSE));
  PetscCall(VecScatterEnd(hctx->scat,hctx->xall,y,INSERT_VALUES,SCATTER_REVERSE));
#else
  /* in-memory (serial/OpenMP) path: every rank holds the full H matrix
     and computes the full product on the gathered x (xall is a
     scatter-to-all copy) */
  PetscCall(VecScatterBegin(hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecScatterEnd(hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecGetArrayRead(hctx->xall,&xa));
  chmmvp_mvp(hctx->handle,(double *)xa,hctx->ball);
  PetscCall(VecRestoreArrayRead(hctx->xall,&xa));
  PetscCall(VecGetArray(y,&ya));
  for(i=hctx->rs;i < hctx->re;i++)
    ya[i-hctx->rs] = (PetscScalar)hctx->ball[i];
  PetscCall(VecRestoreArray(y,&ya));
#endif
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

#ifdef USE_BIGWHAM
/*
   BigWham (full-space) driver. BigWham owns its kernel and operates on a
   3N x 3N traction-from-slip operator in each element's local frame; interact's
   compress / rsf operator is the strike-slip -> strike-shear N x N sub-block.
   We assemble BigWham from interact's patch geometry, then in the matvec expand
   the N strike-slip vector to 3N, apply, and extract the strike traction.

   Sign / scale: BigWham returns traction with its own displacement-discontinuity
   sign convention; interact's resolved shear stress has the opposite sign for the
   self term (slip relieves shear). BIGWHAM_STRESS_SIGN folds that in; it is
   verified against the native full-space Okada dense operator (compress with
   -full_space 1) and may need revisiting if BigWham's convention changes.
*/
#ifndef BIGWHAM_STRESS_SIGN
#define BIGWHAM_STRESS_SIGN (-1.0)
#endif


/* y = A x through BigWham, strike-slip -> strike-shear projection */
PetscErrorCode MatMult_bigwham(Mat A, Vec x, Vec y)
{
  hacapk_shell_ctx *hctx;
  const PetscScalar *xa;
  PetscScalar *ya;
  PetscInt i,n,lo,hi;
  PetscFunctionBeginUser;
  PetscCall(MatShellGetContext(A,&hctx));
  /* gather the full x onto every rank (BigWham is OpenMP, global vectors) */
  PetscCall(VecScatterBegin(hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecScatterEnd  (hctx->scat,x,hctx->xall,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecGetArrayRead(hctx->xall,&xa));
  n = hctx->nelt;
  for(i=0;i < n;i++){		/* expand N strike-slip -> 3N (e1=strike) */
    hctx->x3[3*i+0] = (double)xa[i];
    hctx->x3[3*i+1] = 0.0;
    hctx->x3[3*i+2] = 0.0;
  }
  PetscCall(VecRestoreArrayRead(hctx->xall,&xa));
  cbigwham_mvp(hctx->handle,hctx->x3,hctx->y3);	      /* 3N -> 3N */
  for(i=0;i < n;i++)		/* extract strike traction (e1) -> N */
    hctx->ball[i] = hctx->bscale * hctx->y3[3*i+0];
  PetscCall(VecGetOwnershipRange(y,&lo,&hi));
  PetscCall(VecGetArray(y,&ya));
  for(i=lo;i < hi;i++)
    ya[i-lo] = (PetscScalar)hctx->ball[i];
  PetscCall(VecRestoreArray(y,&ya));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*
  Build (coor,conn) for BigWham 3DR0 rectangles from interact patch geometry,
  assemble the hierarchical matrix, and wrap it in a MATSHELL.

  Corner ordering matches BigWham's local-frame convention (polygon.h:
  tangent1 = v1-v0, tangent2 = v3-v0, normal = t1 x t2), so that with
    v0 = x - l*ts - w*td,  v1 = x + l*ts - w*td,
    v2 = x + l*ts + w*td,  v3 = x - l*ts + w*td,
  BigWham's local e1 = t_strike, e2 = t_dip (l,w are interact HALF lengths).
  Elastic constants come from the same macros the Okada kernel uses, so the
  comparison is consistent.
*/
PetscErrorCode setup_bigwham_matshell(struct med *medium, struct flt *fault,
				      PetscReal scale, int mode,
				      Mat *this_mat, hacapk_shell_ctx **hctx_out)
{
  hacapk_shell_ctx *hctx;
  Vec xd;
  PetscInt i,k,m,n;
  double *coor; int *conn;
  double E,nu,l,w,ts,td,c;
  void *handle;
  int bm,bn; double comp;
  PetscFunctionBeginUser;
  m = n = medium->nrflt;
  coor = (double *)malloc(sizeof(double)*3*4*(size_t)n); /* 4 corners/patch, not shared */
  conn = (int    *)malloc(sizeof(int)   *4*(size_t)n);
  for(i=0;i < n;i++){
    for(k=0;k < 3;k++){
      ts = (double)fault[i].t_strike[k];
      td = (double)fault[i].t_dip[k];
      c  = (double)fault[i].x[k];
      l  = (double)fault[i].l;
      w  = (double)fault[i].w;
      coor[3*(4*i+0)+k] = c - l*ts - w*td;
      coor[3*(4*i+1)+k] = c + l*ts - w*td;
      coor[3*(4*i+2)+k] = c + l*ts + w*td;
      coor[3*(4*i+3)+k] = c - l*ts + w*td;
    }
    conn[4*i+0] = 4*i+0; conn[4*i+1] = 4*i+1;
    conn[4*i+2] = 4*i+2; conn[4*i+3] = 4*i+3;
  }
  E  = (double)YOUNG_MODULUS;	/* = 2 G (1+nu), same constants as Okada */
  nu = (double)POISSON_NU;
  if(!medium->full_space)	/* BigWham is an infinite-medium kernel */
    HEADNODE
      fprintf(stderr,"setup_bigwham_matshell: WARNING: BigWham is full-space only; without -full_space 1 the dense reference is the Okada half-space operator and the comparison is not apples-to-apples\n");
  handle = cbigwham_create(coor,3*4*(int)n,conn,4*(int)n,"3DR0-H",E,nu,
			   (int)medium->bigwham_nthreads);
  cbigwham_build(handle,(int)medium->bigwham_max_leaf,
		 (double)medium->bigwham_eta,(double)medium->bigwham_eps_aca);
  free(coor); free(conn);
  cbigwham_get_info(handle,&bm,&bn,&comp);
  HEADNODE
    fprintf(stderr,"bigwham %i by %i (3 x %i patches), compression ratio %.5g\n",
	    bm,bn,(int)n,comp);

  hctx = (hacapk_shell_ctx *)malloc(sizeof(hacapk_shell_ctx));
  hctx->handle = handle;
  hctx->nelt   = n;
  hctx->x3     = (double *)malloc(sizeof(double)*3*(size_t)n);
  hctx->y3     = (double *)malloc(sizeof(double)*3*(size_t)n);
  hctx->ball   = (double *)malloc(sizeof(double)*m);
  hctx->bscale = (double)scale * BIGWHAM_STRESS_SIGN;
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,(void *)hctx,this_mat));
  PetscCall(MatShellSetOperation(*this_mat,MATOP_MULT,(void (*)(void))MatMult_bigwham));
  PetscCall(MatCreateVecs(*this_mat,&xd,NULL));
  PetscCall(VecScatterCreateToAll(xd,&hctx->scat,&hctx->xall));
  PetscCall(VecGetOwnershipRange(xd,&hctx->rs,&hctx->re));
  PetscCall(VecDestroy(&xd));
  *hctx_out = hctx;
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif /* USE_BIGWHAM */

/* 
   unified, backend-independent report of the assembled interaction-matrix
   storage, so every H-matrix backend emits one identically formatted line
   that can be grepped from a run log. "stored" is the GLOBAL number of
   stored scalars (-1 when the backend does not expose a count); dense_ratio
   is m*n/stored and mbytes assumes 8-byte reals. The per-backend source of
   the count differs (dense: m*n exactly; hacapk: summed over leaf blocks;
   hmmvp: its own info call; htool/h2opus: MatGetInfo when populated), so the
   numbers are only as accurate as each library's own accounting and are
   specific to the tested build. 
*/
const char *hmat_backend_name(int t)
{
  switch(t){
  case IHMAT_TYPE_DENSE:   return "dense";
  case IHMAT_TYPE_HTOOLS:  return "HTOOL";
  case IHMAT_TYPE_H2OPUS:  return "H2OPUS";
  case IHMAT_TYPE_HACAPK:  return "HACApk";
  case IHMAT_TYPE_HMMVP:   return "HMMVP";
  case IHMAT_TYPE_BIGWHAM: return "BigWham";
  default:                 return "unknown";
  }
}
void report_hmat_storage(struct med *medium, const char *backend,
			 PetscInt m, PetscInt n, long stored)
{
  HEADNODE{
    double dense = (double)m * (double)n;
    if(stored > 0)
      fprintf(stderr,"calc_petsc_Isn_matrices: hmat_storage backend %s m %ld n %ld stored_scalars %ld dense_ratio %.6g mbytes %.6g\n",
	      backend,(long)m,(long)n,stored,dense/(double)stored,
	      (double)stored*(double)sizeof(double)/1048576.0);
    else
      fprintf(stderr,"calc_petsc_Isn_matrices: hmat_storage backend %s m %ld n %ld stored_scalars NA dense_ratio NA mbytes NA\n",
	      backend,(long)m,(long)n);
  }
}

PetscErrorCode calc_petsc_Isn_matrices(struct med *medium, struct flt *fault,
				       PetscInt use_hmatrix,PetscReal scale, int mode,
				       Mat *this_mat, hacapk_shell_ctx *hctx)
{
  /* context */
  struct interact_ctx ictx[1];
  PetscReal   *avalues=NULL;
  PetscInt    n, m, lm,ln,i,j,dn,on, *col_idx=NULL;
  long hmat_stored_global = -1;
  /* kernel function */
#ifdef USE_PETSC_HMAT		/* htools and H2opus  */
  PetscReal   *coords=NULL,*av=NULL;
  PetscInt k,jj,*ci=NULL;
  MatHtoolKernelFn *kernel = GenKEntries_petsc;
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
#ifdef USE_HMMVP_MPI
  char hmmvp_tmp[STRLEN],hmmvp_fn[STRLEN];
  PetscMPIInt mrank;
  int cret;
#endif
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
  if(!use_hmatrix){
    /* dense */
    PetscCall(PetscCalloc(m*sizeof(PetscScalar), &avalues));
    PetscCall(PetscCalloc(n*sizeof(PetscInt), &col_idx));
    for (i=0; i < n; i++) 
      col_idx[i] = i;
  }else{
    set_hmat_defaults_and_options(medium,use_hmatrix);
  }

#ifdef USE_PETSC_HMAT
  if((use_hmatrix==IHMAT_TYPE_HTOOLS)||(use_hmatrix==IHMAT_TYPE_H2OPUS)){	
    coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*medium->nrflt);
    for(i=0;i < medium->nrflt;i++)		/* all sources or receiveer coordinates  */
      for(k=0;k < ndim;k++)
	coords[i*ndim+k] = fault[i].x[k];
  }
#endif
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
  if((use_hmatrix==IHMAT_TYPE_HACAPK)||(use_hmatrix==IHMAT_TYPE_HMMVP)){	
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
  
  if(use_hmatrix == IHMAT_TYPE_DENSE){
    PetscCall(MatSeqAIJSetPreallocation(*this_mat, n, NULL));
    PetscCall(MatMPIAIJSetPreallocation(*this_mat, dn, NULL, on, NULL));
  }else{
    /* the H-matrix backends below replace *this_mat with their own type
       (MATHTOOL / MATSHELL), so a full dense-row AIJ preallocation here
       (lm*n entries per rank) is allocated and then discarded. At large N
       that throwaway allocation is O(N^2/np) in work and memory and can
       dominate the build (and exhaust memory); preallocate nothing, since
       only the row distribution (lm/ln/rs/re) computed above is needed. */
    PetscCall(MatSeqAIJSetPreallocation(*this_mat, 0, NULL));
    PetscCall(MatMPIAIJSetPreallocation(*this_mat, 0, NULL, 0, NULL));
  }
  PetscCall(MatGetOwnershipRange(*this_mat, &medium->rs, &medium->re));
  medium->rn = medium->re  - medium->rs; /* number of local elements */

  switch(use_hmatrix){
  case IHMAT_TYPE_DENSE :
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: creating dense matrix for stress mode %i stress type %i\n",
	      mode,ictx->rec_stress_mode);
    PetscCall(MatSetType(*this_mat, MATDENSE));
    /* 
       assemble dense matrix 
    */
    fprintf(stderr,"calc_petsc_Isn_matrices: core %03i/%03i: assigning dense row %5i to %5i\n",
	    (int)medium->comm_rank,(int)medium->comm_size,(int)medium->rs,(int)medium->re);
    for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
      GenKEntries_petsc(ndim,1,n,&j, col_idx, avalues,ictx);
      PetscCall(MatSetValues(*this_mat, 1, &j, n, col_idx,avalues, INSERT_VALUES));
    }
    break;
#ifdef USE_PETSC_HMAT
  case IHMAT_TYPE_HTOOLS:
    /* make the matrix */
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: creating HTOOL matrix for stress mode %i stress type %i\n",
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
  case IHMAT_TYPE_H2OPUS:
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: creating H2OPUS matrix for stress mode %i stress type %i\n",
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
      GenKEntries_petsc(ndim,1,n,&jj,ci,av,ictx);
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
  case IHMAT_TYPE_HACAPK:
#ifdef USE_HACAPK
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: creating HACAPK matrix for stress mode %i stress type %i\n",
	      mode,ictx->rec_stress_mode);
    hacapk_handle = cinit_hacapk_struct((int)m,(void *)ictx);
    cset_hacapk_struct_coord(hacapk_handle,xc,yc,zc);
    cset_hacapk_eta(hacapk_handle,(double)medium->hacapk_eta); /* override param(51) before the build, eta */
    cset_hacapk_inorm(hacapk_handle,medium->hacapk_inorm);     /* error norm mode */
    
    fprintf(stderr,"core %03i/%03i: assigning HACApK m %i n %i ztol %g eta %g inorm: %i\n",
	    medium->comm_rank,medium->comm_size,m,n,(double)medium->hacapk_ztol,
	    (double)medium->hacapk_eta,medium->hacapk_inorm);
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
  case IHMAT_TYPE_HMMVP:
    
#ifdef USE_HMMVP
#ifdef USE_HMMVP_MPI
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: creating HMMVP matrix for stress mode %i stress type %i MPI mode\n",
	      mode,ictx->rec_stress_mode);
    /*
      MPI path: compress with distributed assembly to a temporary file
      (collective over all ranks), then load it as a distributed
      MpiHmat. The matvec (MatMult_hmmvp) gathers x to the root and
      scatters y back. The file name must be identical on all ranks; we
      derive it from the root pid and broadcast it.
    */
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&mrank));
    if(mrank == 0)
      snprintf(hmmvp_fn,STRLEN,"/tmp/hmmvp_interact_%d.hm",(int)getpid());
    PetscCallMPI(MPI_Bcast(hmmvp_fn,STRLEN,MPI_CHAR,0,PETSC_COMM_WORLD));
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: hmmvp MPI compress to %s (tol %g eta %g)\n",
	      hmmvp_fn,(double)medium->hmmvp_tol,(double)medium->hmmvp_eta);
    cret = chmmvp_compress_to_file((int)m,xc,yc,zc,(double)medium->hmmvp_tol,
				   (double)medium->hmmvp_eta,medium->hmmvp_inorm,(void *)ictx,hmmvp_fn);
    if(cret != 0){
      fprintf(stderr,"hmmvp MPI compression failed\n");
      exit(-1);
    }
    hmmvp_handle = chmmvp_mpi_load(hmmvp_fn,medium->hmmvp_nthreads);
    if(!hmmvp_handle){
      fprintf(stderr,"hmmvp MPI load failed\n");
      exit(-1);
    }
    chmmvp_mpi_get_info(hmmvp_handle,&hmm,&hmn,&hmmvp_nnz);
    /* each rank removes its own scratch file "<hmmvp_fn>_<rank>" left by the
       parallel compressor (root concatenates these over MPI but never
       unlinks them); /tmp is typically node-local, so root cannot do it */
    if(mrank > 0){
      snprintf(hmmvp_tmp,STRLEN,"%s_%d",hmmvp_fn,mrank);
      remove(hmmvp_tmp);
    }
    /* the file is fully read into the distributed MpiHmat at load, so
       the root can remove it now */
    if(mrank == 0)
      remove(hmmvp_fn);
  
    HEADNODE
      fprintf(stderr,"hmmvp(MPI) %i by %i, %ld stored scalars, compression ratio %.5g\n",
	      hmm,hmn,hmmvp_nnz,(double)((double)m*(double)n/(double)hmmvp_nnz));
    hctx = (hacapk_shell_ctx *)malloc(sizeof(hacapk_shell_ctx));
    hctx->handle = hmmvp_handle;
    hctx->ball = (double *)malloc(sizeof(double)*m); /* full y on every rank */
    PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,
			     (void *)hctx,this_mat));
    PetscCall(MatShellSetOperation(*this_mat,MATOP_MULT,(void (*)(void))MatMult_hmmvp));
    PetscCall(MatCreateVecs(*this_mat,&xd,NULL));
    /* gather-to-root scatter: xall is full length on rank 0, empty elsewhere */
    PetscCall(VecScatterCreateToZero(xd,&hctx->scat,&hctx->xall));
    PetscCall(VecGetOwnershipRange(xd,&hctx->rs,&hctx->re));
    PetscCall(VecDestroy(&xd));
#else  /* in-memory OpenMP/serial path */
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: creating HMMVP matrix for stress mode %i stress type %i OpenMP\n",
	      mode,ictx->rec_stress_mode);
    hmmvp_handle = chmmvp_compress_in_memory((int)m,xc,yc,zc,(double)medium->hmmvp_tol,
					     (double)medium->hmmvp_eta,medium->hmmvp_inorm,medium->hmmvp_nthreads,
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
#endif /* USE_HMMVP_MPI */
#endif /* USE_HMMVP */
    break;
  case IHMAT_TYPE_BIGWHAM:
#ifdef USE_BIGWHAM
    /* BigWham builds its own mesh from the patch geometry, so it needs none of
       the htool coords or the hacapk/hmmvp xc/yc/zc arrays. It overwrites the
       AIJ *this_mat created above with its MATSHELL, like cases 3 and 4. */
    PetscCall(setup_bigwham_matshell(medium,fault,scale,mode,this_mat,&hctx));
#endif
    break;
  }
  if(use_hmatrix == IHMAT_TYPE_H2OPUS){
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

  /* 
     unified H-matrix storage report: collect a GLOBAL stored-scalar count
     from whichever source the active backend exposes, then emit one common
     line (see report_hmat_storage above). This is post-assembly so the
     htool/h2opus structure is finalized. Collective calls (Allreduce,
     MatGetInfo, MatView) run on all ranks; the print itself is head-node only.
  */
  {

    if(use_hmatrix == IHMAT_TYPE_DENSE)
      hmat_stored_global = (long)m * (long)n;
#ifdef USE_HACAPK
    if(use_hmatrix == IHMAT_TYPE_HACAPK){
      /* cget_hacapk_nnz returns this rank's leaf-block storage; sum for global */
      long hloc = cget_hacapk_nnz(hacapk_handle), hglob = 0;
      PetscCallMPI(MPI_Allreduce(&hloc,&hglob,1,MPI_LONG,MPI_SUM,PETSC_COMM_WORLD));
      hmat_stored_global = hglob;
    }
#endif
#ifdef USE_HMMVP
    if(use_hmatrix == IHMAT_TYPE_HMMVP)
      hmat_stored_global = hmmvp_nnz; /* already the global count on the head node */
#endif
#ifdef USE_PETSC_HMAT
    if((use_hmatrix==IHMAT_TYPE_HTOOLS)||(use_hmatrix==IHMAT_TYPE_H2OPUS)){
      /* htool/h2opus implement neither MatGetInfo (calling it errors here) nor
	 a C accessor for their stored size, and their detailed info does not
	 survive capture to a string viewer, so the unified line below reports
	 stored_scalars NA for these backends. Emit PETSc's own
	 hierarchical-matrix info view to stderr instead (compression ratio,
	 space saving, epsilon, eta, compressor, ...), which otherwise appears
	 only with -mat_view on the command line; read the ratio from that
	 adjacent line. The exact fields depend on the PETSc build and version. */
      PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDERR_WORLD, PETSC_VIEWER_ASCII_INFO));
      PetscCall(MatView(*this_mat, PETSC_VIEWER_STDERR_WORLD));
      PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDERR_WORLD));
    }
#endif
    report_hmat_storage(medium, hmat_backend_name(use_hmatrix), m, n, hmat_stored_global);
  }

  /* free things */
  if(use_hmatrix==IHMAT_TYPE_DENSE ){
    PetscCall(PetscFree(avalues));
    PetscCall(PetscFree(col_idx));
  }
#ifdef USE_PETSC_HMAT
  if((use_hmatrix==IHMAT_TYPE_HTOOLS) || (use_hmatrix==IHMAT_TYPE_H2OPUS))
    free(coords);
#endif
#if ( defined(USE_HMMVP) || defined(USE_HACAPK) )
  if((use_hmatrix==IHMAT_TYPE_HACAPK)||(use_hmatrix==IHMAT_TYPE_HMMVP)){
    free(xc);free(yc);free(zc);
  }
  if(use_hmatrix==IHMAT_TYPE_HACAPK)
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

#endif


PetscErrorCode set_hmat_defaults_and_options(struct med *medium, int hmat) /*  the hmat might be
									       different from medium->hmat */
{
  PetscBool flg;
  HEADNODE
    fprintf(stderr,"set_hmat_defaults_and_options: initializing for %s\n",hmat_backend_name(hmat));
  switch(hmat){
  case IHMAT_TYPE_HTOOLS:
#ifdef USE_PETSC_HMAT
    /* htools */
    PetscCall(PetscOptionsHasName(NULL,NULL,"-mat_htool_eta",&flg));
    if(!flg)
      PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_eta","100")); 
    /* epsilon */
    PetscCall(PetscOptionsHasName(NULL,NULL,"-mat_htool_epsilon",&flg));
    if(!flg)
      PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_epsilon","1e-5"));
    /* compressor */
    PetscCall(PetscOptionsHasName(NULL,NULL,"-mat_htool_compressor",&flg));
    if(!flg)			/* this is a symmetric compressor, a
				   mismatch, but fast? */
      PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_compressor","sympartialACA"));
    /*  */
    PetscCall(PetscOptionsHasName(NULL,NULL,"-pc_type",&flg));
    if(!flg)
      PetscCall(PetscOptionsSetValue(NULL,"-pc_type","none"));
#else
    HEADNODE
      fprintf(stderr,"set_hmat_defaults_and_options: HTOOLS not compiled in - see USE_PETSC_HMAT check makefile.petsc\n");
    exit(-1);
#endif     
    break;
  case IHMAT_TYPE_H2OPUS:
#ifdef USE_PETSC_HMAT
    /* defaults for H2OPUS */
    medium->h2opus_eta = 0.6;	/*  */
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-eta", &medium->h2opus_eta, NULL));
    medium->h2opus_leafsize = 32;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-leafsize", &medium->h2opus_leafsize, NULL));
    medium->h2opus_basisord = 8;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-basisord", &medium->h2opus_basisord, NULL));
#else
    HEADNODE
      fprintf(stderr,"set_hmat_defaults_and_options: H2OPUS not compiled in - see USE_PETSC_HMAT check makefile.petsc\n");
    exit(-1);
    
#endif     
    break;
  case IHMAT_TYPE_HACAPK:
#ifdef USE_HACAPK
    /* hacapl */
    medium->hacapk_ztol = 1.0e-4; /* seems like a solid choice */
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-hacapk_ztol",&medium->hacapk_ztol,NULL));
    medium->hacapk_eta = 2.0; /* admissibility distance param(51); larger admits more far-field, HACApk default 2 */
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-hacapk_eta",&medium->hacapk_eta,NULL));
    /* error norm, 1 (absolute) or 3 (relative) (HBI default) */
    medium->hacapk_inorm = 3;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-hacapk_inorm",&medium->hacapk_inorm,NULL));
#else
    HEADNODE
      fprintf(stderr,"set_hmat_defaults_and_options: HACApk requested but not compiled in (see USE_HACAPK and makefile.petc)\n");
    exit(-1);
#endif
    break;
  case IHMAT_TYPE_HMMVP:
#ifdef USE_HMMVP
    /* hmmvp */
    medium->hmmvp_tol = 1.0e-6;	/* 1e-5 OK for low res, high res might need 1-6*/
    medium->hmmvp_eta = 3.0;
    medium->hmmvp_nthreads = 1;
    medium->hmmvp_inorm = 3;	/* tolerance norm mode, kept comparable with -hacapk_inorm:
				   1 = block-local (tm_brem_fro), else matrix-global (tm_mrem_fro,
				   the hmmvp default and the prior hardcoded behavior) */
    
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-hmmvp_tol",&medium->hmmvp_tol,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-hmmvp_eta",&medium->hmmvp_eta,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-hmmvp_nthreads",&medium->hmmvp_nthreads,NULL));
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-hmmvp_inorm",&medium->hmmvp_inorm,NULL));
#else
    HEADNODE
      fprintf(stderr,"set_hmat_defaults_and_options: HMMVP requested but not compiled in (see USE_HMMVP and makefile.petc)\n");
    exit(-1);
#endif
    break;
  case IHMAT_TYPE_BIGWHAM:
#ifdef USE_BIGWHAM
    medium->bigwham_eta       = 3.0;
    medium->bigwham_eps_aca   = 1.0e-4;
    medium->bigwham_max_leaf  = 32;
    medium->bigwham_nthreads  = 1;
    PetscOptionsGetReal(NULL,NULL,"-bigwham_eta",     &medium->bigwham_eta,      &flg);
    PetscOptionsGetReal(NULL,NULL,"-bigwham_eps_aca", &medium->bigwham_eps_aca,  &flg);
    PetscOptionsGetInt (NULL,NULL,"-bigwham_leaf",    &medium->bigwham_max_leaf, &flg);
    PetscOptionsGetInt (NULL,NULL,"-bigwham_nthreads",&medium->bigwham_nthreads, &flg);
#else
    HEADNODE
      fprintf(stderr,"set_hmat_defaults_and_options: BigWham requested but not compiled in (see USE_BIGWHAM and makefile.petsc)\n");
    exit(-1);
#endif
    break;
  default:
    HEADNODE
      fprintf(stderr,"set_hmat_defaults_and_options: not set up for %s\n",hmat_backend_name(hmat));
    exit(-1);
    break;
  }
  
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}




#endif	/* end use Petsc */
