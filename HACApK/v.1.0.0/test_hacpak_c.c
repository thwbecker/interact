#include "HACApK_c_interface.h"
#include <mpi.h>


double vdiff(double *, double *, int );
void axmult_nbyn(double *, double *, int , double *);


int main(int argc, char **argv) {
  
  void *hacapk_struct = NULL;
  int n,i,j,ierr;
  double *x,*y,*z,ztol,*xvh,*bvh,*xvd,*bvd,*Ad;
  ierr = MPI_Init(&argc,&argv);
  /* dimensions */
  n = 20;
  //
  // Allocate  structure
  hacapk_struct = cinit_hacapk_struct(n);
  if(!hacapk_struct){
    fprintf(stderr,"init failed with N = %i\n",n);
    exit(-1);
  }
  printf("init done, N = %i\n",n);
  
  /* coordinates  */
  x = (double *)malloc(sizeof(double)*n);if(!x){fprintf(stderr,"mem error\n");exit(-1);}
  y = (double *)malloc(sizeof(double)*n);if(!y){fprintf(stderr,"mem error\n");exit(-1);}
  z = (double *)malloc(sizeof(double)*n);if(!z){fprintf(stderr,"mem error\n");exit(-1);}
  chacapk_assign_random_coord(x,y,z,&n); /* initialize randomly */
  printf("random init done\n");

  /* 
     init coordinates 
  */
  cset_hacapk_struct_coord(hacapk_struct, x,y,z);
  printf("coord done\n");
  /*  */
  chacapk_set_kernel_par(hacapk_struct,1.0); /* set kernel parameters */

  /* make an H matrix */
  ztol = 1e-6;
  cmake_hacapk_struct_hmat(hacapk_struct, ztol);


  /* generate a dense matrix using same kernel */
  Ad = (double *)malloc(sizeof(double)*n*n);if(!Ad){fprintf(stderr,"mem error\n");exit(-1);}
  chacapk_assemble_dense_mat(hacapk_struct,Ad,n); 
  printf("dense matrix assembled\n");
  //for(i=0;i<n;i++)for(j=0;j<n;j++)printf("%i %i %g\n",i,j,Ad[j*n+i]);

  xvh = (double *)malloc(sizeof(double)*n);if(!xvh){fprintf(stderr,"mem error\n");exit(-1);}
  bvh = (double *)malloc(sizeof(double)*n);if(!bvh){fprintf(stderr,"mem error\n");exit(-1);}
  xvd = (double *)malloc(sizeof(double)*n);if(!xvd){fprintf(stderr,"mem error\n");exit(-1);}
  bvd = (double *)malloc(sizeof(double)*n);if(!bvd){fprintf(stderr,"mem error\n");exit(-1);}
  /* RHS */
  for(i=0;i<n;i++){
    bvh[i] = bvd[i] = 1e3;
  }
  /* inverse solve H matrix version */
  chacapk_solve_Ab_H(hacapk_struct, bvh, xvh, ztol);
  /* inverse solve dense version */
  chacapk_solve_dense(Ad, n, bvd, xvd);
  if(n<30){
    // print solution
    for(i=0;i<n;i++){
      printf("%05i H %12.4e D %12.4e diff %12.4e\n",i+1,xvh[i],xvd[i],fabs(xvh[i]-xvd[i]));
    }
  }
  printf("difference between H and D solve for x: %e\n",vdiff(xvd,xvh,n));
  
  /* vector multiply using H matrix */
  for(i=0;i<n;i++){
    xvh[i] = xvd[i] = 1.;
  }
  chacapk_mult_Ax_H(hacapk_struct, xvh, bvh);
  axmult_nbyn(Ad,xvd,n,bvd);
  if(n<30){
    // print solution
    for(i=0;i<n;i++){
       printf("%05i H %12.4e D %12.4e diff %12.4e\n",i+1,bvh[i],bvd[i],fabs(bvh[i]-bvd[i]));
    }
  }
  printf("difference between H and D multiplication for b: %e\n",vdiff(bvd,bvh,n));


  free(xvd);free(bvd);free(xvh);free(bvh);free(Ad);
  
  free(x);free(y);free(z);
  // Deallocate structure
  cdeallocate_hacapk_struct(hacapk_struct);

  return 0;
}

double vdiff(double *a, double *b, int n) /* vector difference */
{
  int i = 0;
  double sum = 0,tmp;
  while(i<n){
    tmp = a[i]-b[i];
    sum += tmp*tmp;
    i++;
  }
  return sqrt(sum);
}
void axmult_nbyn(double *A, double *x, int n, double *b) /* b = Ax  for n by n */
{
  int i,j;
  for(i=0;i<n;i++)
    b[i] = 0.0;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      b[i] += A[i*n+j]*x[j];
  

}
