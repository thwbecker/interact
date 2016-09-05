/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: sparse_solve.c,v 1.4 2002/11/05 17:40:59 tbecker Exp $

*/
#include <math.h>
#include <stdio.h>

#ifdef USE_SUPERLU
// SuperLU includes
 #include "dsp_defs.h"
 #include "util.h"
#endif 

#include "interact.h"
/* 

   solve a linear equation system 

   S x = b
   where S is the CCS sparse representation of A which was n by n, 
   b is n by 1, on return, b is overwritten with x

   uses SuperLU library. part of this routine is copied from a SuperLU
   example

*/

void sparse_driver(unsigned int *si1,unsigned int *si2,
		   A_MATRIX_PREC *a, A_MATRIX_PREC *xsol,
		   A_MATRIX_PREC *b,int n,struct med *medium)
{
#ifndef USE_SUPERLU
  fprintf(stderr,"sparse_driver: the sparse solver is only available when SuperLU support\n");
  fprintf(stderr,"sparse_driver: was compiled in. Uncomment \"include makefile.superlu\" in the makefile,\n");
  fprintf(stderr,"sparse_driver: edit the file accordingly, and recompile if you want to use SuperLU sparse methods.\n");
  exit(-1);
#else
  SuperMatrix slu_a, slu_l, slu_u, slu_b;
  NCformat slu_a_ccs_store[1];
  DNformat slu_b_store[1];
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      nrhs, info, i, permc_spec;
  /* 
     Create matrix A in the format expected by SuperLU. 
  */
  slu_a.Stype = NC;// column-wise, no supernode, NC (CCS)
#ifdef A_MATRIX_SINGLE_PREC
  slu_a.Dtype = _S;// double precision
#else
  slu_a.Dtype = _D;// single precision
#endif
  slu_a.Mtype = GE;// general matrix type
  slu_a.nrow = slu_a.ncol = n;
  slu_a.Store = slu_a_ccs_store;
  // actual CCS type storage, SuperLU NC type
  slu_a_ccs_store[0].nnz = si2[n];
  slu_a_ccs_store[0].nzval = a;// pointer to non-zero values
  slu_a_ccs_store[0].rowind = si1; 
  slu_a_ccs_store[0].colptr = si2; 
  /* 
     Create right-hand side matrix B
  */
  slu_b.Stype = DN;// fortran column format dense
#ifdef A_MATRIX_SINGLE_PREC
  slu_b.Dtype = _S;// single precision
#else
  slu_b.Dtype = _D;// double precision
#endif
  slu_b.Mtype = GE;// general matrix type
  slu_b.nrow = n; slu_b.ncol = 1;
  slu_b.Store = slu_b_store;
  slu_b_store[0].lda = n;
  slu_b_store[0].nzval = b;
  // permutation arrays  
  if ( !(perm_r = intMalloc(n)) ) 
    ABORT("sparse_driver: Malloc fails for perm_r[].");
  if ( !(perm_c = intMalloc(n)) ) 
    ABORT("sparse_driver: Malloc fails for perm_c[].");
  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = 0: use the natural ordering 
   *   permc_spec = 1: use minimum degree ordering on structure of A'*A
   *   permc_spec = 2: use minimum degree ordering on structure of A'+A
   */    	
  permc_spec = 0;
  get_perm_c(permc_spec, &slu_a, perm_c);
  //
  //solve system
  //
#ifdef A_MATRIX_SINGLE_PREC
  sgssv(&slu_a, perm_c, perm_r, &slu_l, &slu_u, &slu_b, &info);
#else
  dgssv(&slu_a, perm_c, perm_r, &slu_l, &slu_u, &slu_b, &info);
#endif
  if(info){
    fprintf(stderr,"sparse_driver: dgssv failed: %i\n",info);
    exit(-1);
  }
  // copy solution from b into x 
  /* 
     xols[n] = slu_b_store[0].nzval[n]
  */
  memcpy(xsol, slu_b_store[0].nzval,sizeof(A_MATRIX_PREC)*n);
  /* De-allocate storage */
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_SuperNode_Matrix(&slu_l);
  Destroy_CompCol_Matrix(&slu_u);
#endif
}



