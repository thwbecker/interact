/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: sparse.c,v 1.6 2002/11/07 19:45:33 tbecker Exp $
*/

#include "interact.h"
/* 
   
   creation/reading of sparse matrix storage
   in compressed row or compressed column format (the latter is also 
   called  Harwell-Boeing 

   sparse_nr.c holds the numerical recipes type formats

   WE ARE USING C CONVENTION FOR INDICES, IE. ADD ONE FOR ALL INTEGER
   POINTERS, INCLUDING THE LAST ITEM IN IS2 WHICH IS NOW NNZ AND 
   NOT NNZ+1


*/
//
// create a compressed row storage scheme from n by n matrix A
//
//
size_t create_crs_sparse_from_memory(int n, A_MATRIX_PREC *a, 
				     A_MATRIX_PREC coff,
				     unsigned int **si1,
				     unsigned int **si2,
				     A_MATRIX_PREC **val)
{
  unsigned int nnz,nnz1;
  int i,j,os1,rcnt,newrow;
  *si1 = (unsigned int *)malloc(sizeof(unsigned int));
  *si2 = (unsigned int *)calloc(n+1,sizeof(unsigned int));
  *val = (A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC));
  if((! *val) || (!*si1) || (! *si2)){
    fprintf(stderr,"create_crs_from_memory: memerror 1\n");exit(-1);
  }
  nnz = rcnt = 0;
  nnz1 = 1;
  for(i=os1=0;i<n;i++,os1+=n){
    newrow=0;
    for(j=0;j<n;j++)
      if(fabs(a[os1+j]) > coff){
	// store the value of the non-zero Aij entry
	*(*val + nnz) = a[os1+j];
	// store the column in col, fortran convention
	*(*si1 + nnz) = j;// THIS IS J+1 FOR FORTRAN
	if(!newrow){
	  // if we are starting a new A row, store val index in FORTRAN
	  // convention
	  newrow = 1;
	  *(*si2 + rcnt) = nnz;// THIS WOULD BE NNZ1 FOR FORTRAN
	  rcnt++;
	}
	// increment nnz and get more memory
	nnz++;
	nnz1++;
	*si1 = (unsigned int *)realloc(*si1,nnz1*sizeof(unsigned int));
	*val = (A_MATRIX_PREC *)realloc(*val,nnz1*sizeof(A_MATRIX_PREC));
	if((!*val) || (! *si1)){
	  fprintf(stderr,"create_crs_from_memory: memerror 2\n");exit(-1);
	}
      }
  }
  if(rcnt != n){
    fprintf(stderr,"create_crs_from_memory: internal error\n");
    exit(-1);
  }
  // numer of val entries (WOULD BE NNZ+1 FOR FORTRAN)
  *(*si2 + n) = nnz;
  // return memory size
  return (n+1) * sizeof(unsigned int) +
    nnz * ((sizeof(unsigned int) + sizeof(A_MATRIX_PREC)));
}
// like above but read from file
size_t create_crs_sparse_from_file(int n, 
				   A_MATRIX_PREC coff,
				   unsigned int **si1,
				   unsigned int **si2,
				   A_MATRIX_PREC **val,
				   FILE *in)
{
  unsigned int nnz,nnz1;
  int i,j,rcnt,newrow;
  A_MATRIX_PREC ftmp;
  *si1 = (unsigned int *)malloc(sizeof(unsigned int));
  *si2 = (unsigned int *)calloc(n+1,sizeof(unsigned int));
  *val = (A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC));
  if((! *val) || (!*si1) || (! *si2)){
    fprintf(stderr,"create_crs_from_memory: memerror 1\n");exit(-1);
  }
  nnz = rcnt = 0;
  nnz1 = 1;
  for(i=0;i<n;i++){
    newrow=0;
    for(j=0;j<n;j++){
      ftmp = aij_from_file(i,j,n,in);// read i,j element of A 
      if(fabs(ftmp) > coff){
	// store the value of the non-zero Aij entry
	*(*val + nnz) = ftmp;
	// store the column in col, fortran convention
	*(*si1 + nnz) = j;// THIS IS J+1 FOR FORTRAN
	if(!newrow){
	  // if we are starting a new A row, store val index in FORTRAN
	  // convention
	  newrow = 1;
	  *(*si2 + rcnt) = nnz;// THIS WOULD BE NNZ1 FOR FORTRAN
	  rcnt++;
	}
	// increment nnz and get more memory
	nnz++;
	nnz1++;
	*si1 = (unsigned int *)realloc(*si1,nnz1*sizeof(unsigned int));
	*val = (A_MATRIX_PREC *)realloc(*val,nnz1*sizeof(A_MATRIX_PREC));
	if((!*val) || (! *si1)){
	  fprintf(stderr,"create_crs_from_memory: memerror 2\n");exit(-1);
	}
      }
    }
  }
  if(rcnt != n){
    fprintf(stderr,"create_crs_from_memory: internal error\n");
    exit(-1);
  }
  // numer of val entries (WOULD BE NNZ+1 FOR FORTRAN)
  *(*si2 + n) = nnz;
  // return memory size
  return (n+1) * sizeof(unsigned int) +
    nnz * ((sizeof(unsigned int) + sizeof(A_MATRIX_PREC)));
}
//
// create a compressed column storage scheme from n by n matrix A
// also knows as the Harwell-Boeing sparse matrix format
// the CCS format is the CRS format for A^T. 
//
size_t create_ccs_sparse_from_memory(int n, A_MATRIX_PREC *a, 
				     A_MATRIX_PREC coff,
				     unsigned int **si1,
				     unsigned int **si2,
				     A_MATRIX_PREC **val)
{
  unsigned int nnz,nnz1;
  int i,j,os1,rcnt,newrow;
  *si1 = (unsigned int *)malloc(sizeof(unsigned int));
  *si2 = (unsigned int *)calloc(n+1,sizeof(unsigned int));
  *val = (A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC));
  if((! *val) || (!*si1) || (! *si2)){
    fprintf(stderr,"create_crs_from_memory: memerror 1\n");exit(-1);
  }
  nnz = rcnt = 0;
  nnz1 = 1;
  for(i=0;i<n;i++){
    newrow=0;
    for(j=os1=0;j<n;j++,os1+=n)
      if(fabs(a[os1+i]) > coff){
	*(*val + nnz) = a[os1+i];
	*(*si1 + nnz) = j;// J+1 FOR FORTRAN
	if(!newrow){
	  newrow = 1;
	  *(*si2 + rcnt) = nnz;//NNZ1 FOR FORTRAN
	  rcnt++;
	}
	nnz++;
	nnz1++;
	*si1 = (unsigned int *)realloc(*si1,nnz1*sizeof(unsigned int));
	*val = (A_MATRIX_PREC *)realloc(*val,nnz1*sizeof(A_MATRIX_PREC));
	if((!*val) || (! *si1)){
	  fprintf(stderr,"create_crs_from_memory: memerror 2\n");exit(-1);
	}
      }
  }
  if(rcnt != n){
    fprintf(stderr,"create_crs_from_memory: internal error\n");
    exit(-1);
  }
  *(*si2 + n) = nnz;// NNZ1 FOR FORTRAN
  return (n+1) * sizeof(unsigned int) +
    nnz * ((sizeof(unsigned int) + sizeof(A_MATRIX_PREC)));
}
//
// create a compressed column storage scheme from n by n matrix A
// stored on file (see above)
//
size_t create_ccs_sparse_from_file(int n, A_MATRIX_PREC coff,
				   unsigned int **si1,
				   unsigned int **si2,
				   A_MATRIX_PREC **val, FILE *in)
{
  unsigned int nnz,nnz1;
  int i,j,rcnt,newrow;
  A_MATRIX_PREC tflt;
  *si1 = (unsigned int *)malloc(sizeof(unsigned int));
  *si2 = (unsigned int *)calloc(n+1,sizeof(unsigned int));
  *val = (A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC));
  if((! *val) || (!*si1) || (! *si2)){
    fprintf(stderr,"create_crs_from_memory: memerror 1\n");exit(-1);
  }
  nnz = rcnt = 0;
  nnz1 = 1;
  for(i=0;i<n;i++){
    newrow=0;
    for(j=0;j<n;j++){
      tflt = aij_from_file(j,i,n,in);// read j,i element from file
      if(fabs(tflt) > coff){
	*(*val + nnz) = tflt;
	*(*si1 + nnz) = j;// J+1 FOR FORTRAN
	if(!newrow){
	  newrow = 1;
	  *(*si2 + rcnt) = nnz;// NNZ1 FOR FORTRAN
	  rcnt++;
	}
	nnz++;
	nnz1++;
	*si1 = (unsigned int *)realloc(*si1,nnz1*sizeof(unsigned int));
	*val = (A_MATRIX_PREC *)realloc(*val,nnz1*sizeof(A_MATRIX_PREC));
	if((!*val) || (! *si1)){
	  fprintf(stderr,"create_crs_from_memory: memerror 2\n");exit(-1);
	}
      }
    }
  }
  if(rcnt != n){
    fprintf(stderr,"create_crs_from_memory: internal error\n");
    exit(-1);
  }
  *(*si2 + n) = nnz;// NNZ1 FOR FOR FORTRAN
  return (n+1) * sizeof(unsigned int) +
    nnz * ((sizeof(unsigned int) + sizeof(A_MATRIX_PREC)));
}
  



