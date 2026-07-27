#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "slatec_src.h"

#ifdef SINGLE_PREC
#define COMP_PRECISION float
#define DATA_FORMAT "%f"
#else
#define DATA_FORMAT "%lf"
#define COMP_PRECISION double
#endif
/* 

driver for the nnsl routine from slatec

solve A.x = b with x_i >= 0 

A(m,n), b(m), x(n)

l = n: identical to least squares solution
l = 0: all parameters constrained to be non-negative
 
*/
void nnls_driver(int m,		/* number of data (columns) */
		 int n, 	/* number of parameters (rows) */
		 int l, /* how many of the parameters are free */
		 COMP_PRECISION *a, /* input: design matrix (will not be changed) */
		 COMP_PRECISION *b, /* data */
		 COMP_PRECISION *x, /* output */
		 COMP_PRECISION *rnorm,
		 int *mode)
{
  COMP_PRECISION *atmp,*work;
  COMP_PRECISION prgopt[1]={1.0};
  int me,ma,*iwork,i,nm,itmp1,itmp2;

  // number of exact equations
  me=0;
  // number of least squares equations
  ma=m;
  /* work space sizes */
  itmp1 = me+ma+5*n+20;
  itmp2 = me+ma+n+20;
  work=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*itmp1);  
  iwork=(int *)malloc(sizeof(int)*itmp2);
  iwork[0] = itmp1;
  iwork[1] = itmp2;
  atmp=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*((n+1)*m+20));

  if(!iwork || !work || !atmp){
    fprintf(stderr,"nnls_driver: memory allocation error\n");
    exit(-1);
  }
  /* reassign A and b to atmp */
  nm = n * m;
  /* A matrix */
  memcpy(atmp,a,nm*sizeof(COMP_PRECISION));
  /* b vector */
  memcpy((atmp+nm),b,m*sizeof(COMP_PRECISION));
#ifdef SINGLE_PREC
  wnnls_(atmp, &m, &me, &ma, &n, &l, prgopt, x, rnorm, mode,
	 iwork, work);
#else
  dwnnls_(atmp, &m, &me, &ma, &n, &l, prgopt, x, rnorm, mode,
	  iwork, work);
#endif
  free(atmp);free(iwork);free(work);
}

void nnls_driver_ftn(int *m,int *n, int *l,
		     COMP_PRECISION *a,
		     COMP_PRECISION *b, 
		     COMP_PRECISION *x,
		     COMP_PRECISION *rnorm,
		     int *mode)
{
  nnls_driver(*m,*n,*l,a,b,x,rnorm,mode);
}
void nnls_driver_ftn__(int *m,int *n, int *l,
		      COMP_PRECISION *a,
		      COMP_PRECISION *b, 
		      COMP_PRECISION *x,
		      COMP_PRECISION *rnorm,
		      int *mode)
{

  nnls_driver_ftn(m,n,l,a,b,x,rnorm,mode);

}
