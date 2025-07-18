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

/* 
   
   compute fixed slip interaction matrices

   mode 0: shear stress
   mode 1: normal stress

   scale: scale factor 

   hmatrix can be 0 (dense), 1 HTOOLS or 2 HOPUS
*/

PetscErrorCode calc_petsc_Isn_matrices(struct med *medium, struct flt *fault,
				       PetscInt use_hmatrix,PetscReal scale, int mode, Mat *this_mat)
{
  /* context */
  struct interact_ctx ictx[1];
  PetscReal   *avalues=NULL;
  PetscInt    n, m, lm,ln,i,j,dn,on, *col_idx=NULL;
  /* kernel function */
#ifdef USE_PETSC_HMAT
  PetscReal   *coords=NULL;
  PetscInt k;
  MatHtoolKernelFn *kernel = GenKEntries_htools;
  MatH2OpusKernelFn *h2opus_kernel = GenKEntries_h2opus;
#endif
  const PetscInt ndim = 3;
  ictx->medium = medium;
  ictx->fault = fault;
  /* defines how to slip */
  ictx->src_slip_mode = 0;
  /*  */
  medium->use_hmatrix = use_hmatrix;
  /*  */
  m = n = medium->nrflt;
#ifdef USE_PETSC_HMAT
  if(medium->use_hmatrix == 1){	/* htools */
    coords = (PetscReal *)realloc(coords,sizeof(PetscReal)*ndim*m);
    for(i=0;i < m;i++)		/* all sources or receiveer coordinates  */
      for(k=0;k < ndim;k++)
	coords[i*ndim+k] = fault[i].x[k];

  }else if(medium->use_hmatrix == 2){	/* H2OPUS */
    fprintf(stderr,"calc_petsc_Isn_matrices: H2OPUS not implemented yet\n");
    exit(-1);
  }else{
#else
    if(medium->use_hmatrix){
      fprintf(stderr,"calc_petsc_Isn_matrices: requesting H matrix but not compiled as such (USE_PETSC_HMAT not set)\n");
      exit(-1);
    }
#endif
    PetscCall(PetscCalloc(m*sizeof(PetscScalar), &avalues));
    PetscCall(PetscCalloc(n*sizeof(PetscInt), &col_idx));
    for (i=0; i < n; i++) 
      col_idx[i] = i;
#ifdef USE_PETSC_HMAT
   }
#endif
  
  if(mode==0){
    ictx->rec_stress_mode = STRIKE;
  }else{
    ictx->rec_stress_mode = NORMAL;
  }

  /* make the matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD, this_mat));
  PetscCall(MatSetSizes(*this_mat, PETSC_DECIDE, PETSC_DECIDE, m, n));
#ifdef USE_PETSC_HMAT
  if(medium->use_hmatrix == 1){
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: using Htool matrix for stress mode %i stress type %i\n",
	      mode,ictx->rec_stress_mode);
    PetscCall(MatSetType(*this_mat, MATHTOOL));
  }else if(medium->use_hmatrix == 2){
    HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: using H2opus matrix for stress mode %i stress type %i\n",
	      mode,ictx->rec_stress_mode);
    PetscCall(MatSetType(*this_mat, MATH2OPUS));
  }else{
#endif
     HEADNODE
      fprintf(stderr,"calc_petsc_Isn_matrices: using dense matrix for stress mode %i stress type %i\n",
	      mode,ictx->rec_stress_mode);
    PetscCall(MatSetType(*this_mat, MATDENSE));
#ifdef USE_PETSC_HMAT
   }
#endif
  PetscCall(MatSetUp(*this_mat));
  PetscCall(MatGetLocalSize(*this_mat, &lm, &ln));
  dn = ln;on = n - ln;
  
  PetscCall(MatSeqAIJSetPreallocation(*this_mat, n, NULL));
  PetscCall(MatMPIAIJSetPreallocation(*this_mat, dn, NULL, on, NULL));
  PetscCall(MatGetOwnershipRange(*this_mat, &medium->rs, &medium->re));
  medium->rn = medium->re  - medium->rs; /* number of local elements */
#ifdef USE_PETSC_HMAT
 if(medium->use_hmatrix==1){
    /* 
       
       H matrix setup HTOOLS
       
    */
    fprintf(stderr,"calc_petsc_Isn_matrices: core %03i/%03i: assigning HTOOL  row %5i to %5i, lm %i ln %i m %i n %i\n",
	    medium->comm_rank,medium->comm_size,medium->rs,medium->re,lm,ln,m,n);
    PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,lm,ln, m, n,ndim,(coords+medium->rs), (coords+medium->re),
				       kernel,(void *)ictx, this_mat));
 }else if(medium->use_hmatrix==2){
    /* 
       
       H matrix setup H2OPUS
       
    */
    fprintf(stderr,"calc_petsc_Isn_matrices: core %03i/%03i: assigning H2OPUS row %5i to %5i, lm %i ln %i m %i n %i\n",
	    medium->comm_rank,medium->comm_size,medium->rs,medium->re,lm,ln,m,n);

    PetscCall(MatCreateH2OpusFromKernel(PETSC_COMM_WORLD, lm, ln, m, n,ndim, (coords+medium->rs),
					PETSC_FALSE, h2opus_kernel,(void *)ictx,
					medium->h2opus_eta, medium->h2opus_leafsize,
					medium->h2opus_basisord, &medium->In));
  
  }else{
#endif
   /* 
       assemble dense matrix 
    */
   fprintf(stderr,"calc_petsc_Isn_matrices: core %03i/%03i: assigning dense row %5i to %5i\n",
	   (int)medium->comm_rank,(int)medium->comm_size,(int)medium->rs,(int)medium->re);
    for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
      GenKEntries_htools(ndim,1,n,&j, col_idx, avalues,ictx);
      PetscCall(MatSetValues(*this_mat, 1, &j, n, col_idx,avalues, INSERT_VALUES));
    }
 #ifdef USE_PETSC_HMAT
 } 
#endif
  PetscCall(MatSetOption(*this_mat, MAT_SYMMETRIC, PETSC_FALSE));
  PetscCall(MatSetFromOptions(*this_mat));
  PetscCall(MatAssemblyBegin(*this_mat, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(*this_mat, MAT_FINAL_ASSEMBLY));
  /* 
     done, now scale (CHECK IF IMPLEMENTED)
  */
  PetscCall(MatScale(*this_mat,scale));
    
 #ifdef USE_PETSC_HMAT
  if(medium->use_hmatrix){
    free(coords);
  }else{
#endif
    PetscCall(PetscFree(avalues));
    PetscCall(PetscFree(col_idx));

 #ifdef USE_PETSC_HMAT
   }
#endif
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}

#endif
