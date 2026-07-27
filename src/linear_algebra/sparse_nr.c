/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu


*/

#include "interact.h"
/* 
   
   creation/reading of sparse matrix storage
   in numerical recipes scheme

*/

//
// create a sparse matrix in numerical recipes format
size_t create_nrs_sparse(struct med *medium, I_MATRIX_PREC thres,
			 I_MATRIX_PREC *imat,
			 my_boolean matrix_in_memory)
{
  if(matrix_in_memory)// this is for debugging only
    return create_nrs_sparse_from_memory(medium->nmat1, imat,thres,
				     &medium->is1,&medium->val);
  else
    return create_nrs_sparse_from_file(medium->nmat1,thres,
				       &medium->is1,&medium->val,
				       medium->i_mat_in);
}

size_t create_nrs_sparse_from_memory(int n, I_MATRIX_PREC *imat,
				     I_MATRIX_PREC thres,
				     unsigned int **ija,
				     I_MATRIX_PREC **sa)
{
  unsigned int k,k1;
  int i,j;
  size_t ijasize,sasize;
  *sa=(I_MATRIX_PREC *)calloc(n,sizeof(I_MATRIX_PREC));
  *ija=(unsigned int *)calloc(n,sizeof(unsigned int));
  sasize = n * sizeof(I_MATRIX_PREC);
  ijasize = n * sizeof(unsigned int);
  if(! *sa || ! *ija)
    MEMERROR("create_sparse_from_memory: 1:");
  for (j=0;j<n;j++)
    *(*sa+j)=imat[j*n+j];
  k= *(*ija)=(unsigned int)n+1;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      if ((fabs(imat[i*n+j]) > thres) && (i != j)) {
	sasize = (k1=k+1)*sizeof(I_MATRIX_PREC);
	*sa=(I_MATRIX_PREC *)realloc(*sa,sasize);
	ijasize = k1*sizeof(unsigned int);
	*ija=(unsigned int *)realloc(*ija,ijasize);
	if(! *sa || ! *ija)
	  MEMERROR("create_sparse_from_memory: 2:");
	*(*sa+k)  = imat[i*n+j];
	*(*ija+k) = j;
	k++;
      }
    }
    *(*ija+i+1)=k;
  }
  return ijasize + sasize;
}

/*

  create the sparse matrix storage vectors
  ija and sa from an n by n matrix that has been written to a file

  numerical recipes format 

  returns the size of the pointers in bytes
*/
size_t create_nrs_sparse_from_file(int n, I_MATRIX_PREC thres,
				   unsigned int **ija,
				   I_MATRIX_PREC **sa,FILE *in)
{
  unsigned int k,k1;
  int i,j;
  size_t ijasize, sasize;
  I_MATRIX_PREC flttmp;
  *sa=(I_MATRIX_PREC *)calloc(n,sizeof(I_MATRIX_PREC));
  sasize = n * sizeof(I_MATRIX_PREC);
  *ija=(unsigned int *)calloc(n,sizeof(unsigned int));
  ijasize = n * sizeof(unsigned int);
  if(! *sa || ! *ija)
    MEMERROR("create_sparse_from_file: 1:");
  // read in diagonal elements
  for (j=0;j<n;j++)
    *(*sa+j)=aij_from_file(j,j,n,in);
  k= *(*ija)=(unsigned int)n+1;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      flttmp = aij_from_file(i,j,n,in);
      if ((fabs(flttmp) > thres) && (i != j)) {
	sasize = (k1=k+1)*sizeof(I_MATRIX_PREC);
	*sa=(I_MATRIX_PREC *)realloc(*sa,sasize);
	ijasize = k1*sizeof(unsigned int);
	*ija=(unsigned int *)realloc(*ija,ijasize);
	if(! *sa || ! *ija)
	  MEMERROR("create_sparse_from_file: 2:");
	*(*sa+k)  = flttmp;
	*(*ija+k) = j;
	k++;
      }
    }
    *(*ija+i+1)=k;
  }
  return ijasize + sasize;
}

// get element i,j given numerical recipes storaged scheme sparse 
// matrix
I_MATRIX_PREC get_nrs_sparse_el(int get_i,int get_j,
				unsigned int *ija,I_MATRIX_PREC *sa)
{
  unsigned int k;
  if(get_i == get_j)// diagonal element
    return(sa[get_i]);
  else if(ija[get_i] == ija[get_i+1])// no off diagonal elements
    return 0.0;
  else{
    for(k=ija[get_i];(k<ija[get_i+1]) && (ija[k]<get_j);k++)
      ;
    if(k == ija[get_i+1])
      return 0.0;
    else if(ija[k] == get_j)
      return sa[k];
    else
      return 0.0;
  }
}
