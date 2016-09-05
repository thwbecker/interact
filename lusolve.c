/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: lusolve.c,v 1.10 2003/01/13 06:44:58 becker Exp $

*/
#include <math.h>
#include <stdio.h>
#include "interact.h"
/* 

   solve a singular linear equation system by LAPACK LU decomposition

   A x = b
   where A is m by n, b is m by 1
   on return, b is overwritten with x

   WARNING:
  
   solver assumes that matrices are stored in FORTRAN convention,
   not in C style

   A gets overwritten

*/
void lu_driver(A_MATRIX_PREC *a,A_MATRIX_PREC *xsol,
	       A_MATRIX_PREC *b,int m,int n,struct med *medium)
{
  int *ipiv,info;
  static int iunity = 1;
  A_MATRIX_PREC dummy,*bcopy;
  static int print_count=0;
  size_t nsize;
  FILE *out1,*out2;
  char out_string1[STRLEN],out_string2[STRLEN];
  nsize = sizeof(A_MATRIX_PREC)*n;
  /* save original b */
  bcopy = (A_MATRIX_PREC *)malloc(nsize);
  if(!bcopy)MEMERROR("lu_driver");
  memcpy(bcopy,b,nsize);
  /*  */
  if(medium->debug){
    if(n != m){
      fprintf(stderr,"lu_driver: error: maxtrix has to be square: m: %i n: %i\n",
	      m,n);
      exit(-1);
    }
    // debugging output
    sprintf(out_string1,"/tmp/interact/a.%i.dat",print_count);
    sprintf(out_string2,"/tmp/interact/b.%i.dat",print_count);
    fprintf(stderr,"lu_driver: %i times %i system, writing to A to %s and b to %s\n",
	    m,n,out_string1,out_string2);
    out1=myopen(out_string1,"w");
    out2=myopen(out_string2,"w");
    print_a_matrix(a,m,n,out1,&dummy,FALSE);
    print_b_vector(b,m,out2,&dummy,FALSE);
    fclose(out1);fclose(out2);
  }
  //
  // allocate work space
  //
  if((ipiv=(int *)malloc(sizeof(int)*n))==NULL)
    MEMERROR("lu_driver: 1:");
  //
  // solve by LU decomposition 
  //
#ifdef A_MATRIX_SINGLE_PREC
  sgesv_(&n,&iunity,a,&n,ipiv,b,&n,&info);
  if(info != 0){// error in xGESV routine
    fprintf(stderr,"lu_driver: sgesv error code: %i\n",info);
#else
  dgesv_(&n,&iunity,a,&n,ipiv,b,&n,&info);
  if(info != 0){// error in xGESV routine
    fprintf(stderr,"lu_driver: dgesv error code: %i\n",info);
#endif
    fprintf(stderr,"lu_driver: n: %i nrhs: %i lda: %i ldb: %i\n",n,iunity,n,n);
    fprintf(stderr,"lu_driver: is the matrix singular? try the SVD solver\n");
    if(n < 81){// output of small matrices for debugging
      fprintf(stderr,"lu_driver: A after return:\n");
      print_matrix_ftrn(a,m,n,stderr,FALSE);
      fprintf(stderr,"lu_driver: B after return:\n");
      print_vector(b,n,stderr);
    }
    exit(-1);
  }
  //
  // copy b to solution vector xsol
  //
  memcpy(xsol,b,nsize);
  // restore riginal b vector
  memcpy(b,bcopy,nsize);
  free(ipiv);free(bcopy);
  if(medium->debug){
    sprintf(out_string1,"/tmp/interact/x.%i.dat",print_count);
    fprintf(stderr,"lu_driver: writing solution x to %s\n",
	    out_string1);
    out1=myopen(out_string1,"w");
    print_b_vector(xsol,n,out1,&dummy,FALSE);
    fclose(out1);
    print_count++;
  }
}


