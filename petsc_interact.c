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
#include <petscksp.h>
PetscErrorCode GenEntries(PetscInt , PetscInt , PetscInt ,const PetscInt *, const PetscInt *, PetscScalar *, void *);



void calc_petsc_Isn_matrices(struct med *medium, struct flt *fault,  PetscBool use_h)
{
  struct interact_ctx ictx[1];
  PetscReal   *coords,*avalues=NULL;
  PetscInt    n, m, lm,ln,i,j,k,dn,on, *col_idx=NULL;
  /* defined in interact.c */
  MatHtoolKernelFn *kernel = GenEntries;
  const PetscInt ndim = 3;
  ictx->medium = medium;
  ictx->src_slip_mode = STRIKE;
  ictx->fault = fault;
  medium->use_h = use_h;
  
  m = n = medium->nrflt;

  PetscCall(MatCreate(PETSC_COMM_WORLD, &medium->Is));PetscCall(MatCreate(PETSC_COMM_WORLD, &medium->In));
  PetscCall(MatSetSizes(medium->Is, PETSC_DECIDE, PETSC_DECIDE, m, n));
  PetscCall(MatSetSizes(medium->In, PETSC_DECIDE, PETSC_DECIDE, m, n));  
  if(use_h){
    PetscCall(MatSetType(medium->Is,MATHTOOL));
    PetscCall(MatSetType(medium->In,MATHTOOL));
  }  else { 
    PetscCall(MatSetType(medium->Is, MATDENSE));
    PetscCall(MatSetType(medium->In, MATDENSE));
  }
  PetscCall(MatSetUp(medium->Is));PetscCall(MatSetUp(medium->In));
  PetscCall(MatGetLocalSize(medium->Is, &lm, &ln));PetscCall(MatGetLocalSize(medium->In, &lm, &ln));
  dn = ln;on = n - ln;
  PetscCall(MatSeqAIJSetPreallocation(medium->Is, n, NULL));
  PetscCall(MatMPIAIJSetPreallocation(medium->Is, dn, NULL, on, NULL));
  PetscCall(MatGetOwnershipRange(medium->Is, &medium->rs, &medium->re));
  PetscCall(MatSeqAIJSetPreallocation(medium->In, n, NULL));
  PetscCall(MatMPIAIJSetPreallocation(medium->In, dn, NULL, on, NULL));
  PetscCall(MatGetOwnershipRange(medium->In, &medium->rs, &medium->re));

  medium->rn = medium->re  - medium->rs; /* number of local elements */
  if(use_h){
    coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*n);
    for(i=0;i < m;i++)		/* all sources or receiveer coordinates  */
      for(k=0;k < 3;k++)
	coords[i*ndim+k] = fault[i].x[k];
    fprintf(stderr,"calc_petsc_Isn_matrices: core %03i/%03i: assigning htool row %5i to %5i\n",
	    medium->comm_rank,medium->comm_size,medium->rs,medium->re);
    ictx->rec_stress_mode = STRIKE;
    PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,lm,ln, m, n,
				       ndim,(coords+medium->rs), (coords+medium->re), kernel,ictx, &medium->Is));
    
    ictx->rec_stress_mode = NORMAL;
    PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,lm,ln, m, n,
				       ndim,(coords+medium->rs), (coords+medium->re), kernel,ictx, &medium->In));
    
    free(coords);
  }else{
    /* 
       assemble dense matrix 
    */
    fprintf(stderr,"calc_petsc_Isn_matrices: core %03i/%03i: assigning dense row %5i to %5i\n",
	    medium->comm_rank,medium->comm_size,medium->rs,medium->re);
    PetscCall(PetscCalloc(m*sizeof(PetscScalar), &avalues));
    PetscCall(PetscCalloc(n*sizeof(PetscInt), &col_idx));
    for (i=0; i < n; i++) 
      col_idx[i] = i;
    
    ictx->rec_stress_mode = STRIKE;
    for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
      GenEntries(ndim,1,n,&j, col_idx, avalues,ictx);
      PetscCall(MatSetValues(medium->Is, 1, &j, n, col_idx,avalues, INSERT_VALUES));
    }
    ictx->rec_stress_mode = NORMAL;
    for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
      GenEntries(ndim,1,n,&j, col_idx, avalues,ictx);
      PetscCall(MatSetValues(medium->In, 1, &j, n, col_idx,avalues, INSERT_VALUES));
    }
    
    PetscCall(PetscFree(avalues));
    PetscCall(PetscFree(col_idx));
  }
  PetscCall(MatAssemblyBegin(medium->Is, MAT_FINAL_ASSEMBLY));PetscCall(MatAssemblyEnd(medium->Is, MAT_FINAL_ASSEMBLY));
  PetscCall(MatSetOption(medium->Is, MAT_SYMMETRIC, PETSC_FALSE));PetscCall(MatSetFromOptions(medium->Is));
  PetscCall(MatAssemblyBegin(medium->In, MAT_FINAL_ASSEMBLY));PetscCall(MatAssemblyEnd(medium->In, MAT_FINAL_ASSEMBLY));
  PetscCall(MatSetOption(medium->In, MAT_SYMMETRIC, PETSC_FALSE));PetscCall(MatSetFromOptions(medium->In));

}

#endif
