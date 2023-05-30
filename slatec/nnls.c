#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "slatec_src.h"

//
// solve NNLS problem using routine
// provided by Hanson and Haskell
//

#ifdef SINGLE_PREC
#define COMP_PRECISION float
#define DATA_FORMAT "%f"
#else
#define DATA_FORMAT "%lf"
#define COMP_PRECISION double
#endif


int main(int argc,char **argv)
{
  int n,m,ma,me,l,i,j;
  FILE *in;
  COMP_PRECISION *a,*b,*x,rnorm,*work;
  int mode;
  // read in parameters and matrix
  if(argc != 4){
    fprintf(stderr,"NNLS\n");
    fprintf(stderr,"\tusage: \n\n");
    fprintf(stderr,"\t%s m_data n_para l\n\n",argv[0]);
    fprintf(stderr,"\tsolves A x = b in a least squares sense with non-negativity contraints\n");
    fprintf(stderr,"\ton the l+1 ... n_para solution parameters\n\n");
    fprintf(stderr,"\tmatrix A is read from a.dat, b is read from file b.dat\n");
    fprintf(stderr,"\tsolution vector x is written to stdout\n\n");
    fprintf(stderr,"\tif l=n_para: identical to least squares solution\n");
    fprintf(stderr,"\tif l=0:      all parameters constrained to be non-negative\n");
    exit(-1);
  }
  sscanf(argv[1],"%i",&m); 
  sscanf(argv[2],"%i",&n);
  sscanf(argv[3],"%i",&l);
  if(l+1<=n)
    fprintf(stderr,"%s: m(rows): %i n(columns): %i, NN: rows %i to %i\n",
	    argv[0],m,n,l+1,n);
  else
    fprintf(stderr,"%s: m(rows): %i n(columns): %i, no NN rows\n",
	    argv[0],m,n);
    
  if(l<0||l>n){
    fprintf(stderr,"%s: l: %i out of bounds (0<=l<=n)\n",
	    argv[0],l);exit(-1);
  }
  if(l==n)
    fprintf(stderr,"%s: completely unconstrained least squares\n",argv[0]);
  if(l==0)
    fprintf(stderr,"%s: completely non-negative  least squares\n",argv[0]);
  
  a=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n*m);
  x=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n);
  b=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m);
  if(!a || !x || !b){
    fprintf(stderr,"memerr A or x or b\n");exit(-1);
  }
  in=fopen("a.dat","r");
  if(in){
    for(i=0;i<m;i++)
      for(j=0;j<n;j++)
	if(fscanf(in,DATA_FORMAT,(a+j*m+i))!=1){
	  fprintf(stderr,"read error for A\n");
	  exit(-1);
	}
  }else{
    fprintf(stderr,"%s: could not open a.dat\n",
	    argv[0]);
    exit(-1);
  }
  fclose(in);

  in=fopen("b.dat","r");
  if(in){
    j = n * m;
    for(i=0;i<m;i++)
      if(fscanf(in,DATA_FORMAT,(b+i))!=1){
	fprintf(stderr,"read error for b\n");
	exit(-1);
      }
  }else{
    fprintf(stderr,"%s: could not open b.dat\n",
	    argv[0]);
    exit(-1);
  }
  fclose(in);
  /* call the NNLS routine */
  nnls_driver(m,n,l,a,b,x,&rnorm,&mode);
  fprintf(stderr,"%s: dwnnls: rnorm: %11g mode %i\n",
	  argv[0],rnorm,mode);
  for(i=0;i<n;i++)
    printf("%17.10e\n",x[i]);
  free(a);free(b);free(x);
  return 0;
}
