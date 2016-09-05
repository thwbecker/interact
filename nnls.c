/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: nnls.c,v 2.11 2003/01/13 06:45:08 becker Exp $
*/
#include <math.h>
#include <stdio.h>
#include "interact.h"
/*
  obtain a non-negative solution to the least squares
  problem  A x = b using the Lawson algorithm
  
  A is assumed to be m by n 
  
  input is transposed matrix since we pass to a FORTRAN
  routine. matrix indices are C style 0 ... N-1 

  uses Lawson & Henson's routine as provided in nnls_lawson.F

  WARNING:
  
  solver assumes that matrices are stored in FORTRAN convention,
  not in C style


*/
#define NORM_NNLS_EPS 1.0e-8
#ifdef DEBUG
 #define CHECK_B_NORM 
#endif
void nnls_driver(A_MATRIX_PREC *a,A_MATRIX_PREC *xsol, 
		 A_MATRIX_PREC *b,int m,int n)
{
  int i,mode,*index;
  A_MATRIX_PREC rnorm,*w,*zz;
#ifdef CHECK_B_NORM
  A_MATRIX_PREC xnorm;
#endif
#ifdef DEBUG
  FILE *out1,*out2;
  A_MATRIX_PREC dummy;
  fprintf(stderr,"nnls_driver: writing to %s and %s\n",
	  DEBUG_A_CON_MATRIX_ASCII_OUT,DEBUG_B_CON_VECTOR_ASCII_OUT);
  out1=myopen(DEBUG_A_CON_MATRIX_ASCII_OUT,"w");
  out2=myopen(DEBUG_B_CON_VECTOR_ASCII_OUT,"w");
  print_a_matrix(a,m,n,out1,&dummy,FALSE);
  print_b_vector(b,m,out2,&dummy,FALSE);
  fclose(out1);fclose(out2);
#endif
  if((zz=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*m))==NULL)
    MEMERROR("nnls_driver: 1:");
  if((w=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*n))==NULL)
    MEMERROR("nnls_driver: 2:");
  if((index=(int *)malloc(sizeof(int)*n))==NULL)
  MEMERROR("nnls_driver: 4:");
  /* call lawson non-negative least squares routines */
  law_nnls(a, &m, &n, &n, b, xsol, &rnorm, w, zz, index, &mode);
  if(rnorm > NORM_NNLS_EPS)
    fprintf(stderr,"nnls_driver: ----------> WARNING: res norm  %g > %e <--------- \n",
	    rnorm,NORM_NNLS_EPS);
  free(zz);free(w);free(index);
  /* assign x to b */
  if(mode != 1){
    fprintf(stderr,"nnls_driver: nnls solver problem, mode = %i, no solution obtained\n",mode);
    for(i=0;i<n;i++)
      xsol[i]=0.0;
  }else{
#ifdef CHECK_B_NORM
    xnorm = norm(xsol,n);
#endif
  }
#ifdef CHECK_B_NORM
  if(xnorm < 1e-8)
    fprintf(stderr,"nnls_driver: -----------> WARNING: norm of x is only %g <-----------\n",
	    xnorm);
#endif
#ifdef DEBUG
  fprintf(stderr,"nnls_driver: writing to %s\n",
	  DEBUG_X_CON_VECTOR_ASCII_OUT); 
  out1=myopen(DEBUG_X_CON_VECTOR_ASCII_OUT,"w");
  print_b_vector(xsol,n,out1,&dummy,FALSE);
  fclose(out1);
#endif
}
