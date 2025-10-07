#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
//
// Declare Fortran functions as external C functions
//
extern void init_hacapk_struct(void **obj_ptr, int *n);
extern void deallocate_hacapk_struct(void **obj_ptr);
extern void set_hacapk_struct_coord(void **obj_ptr, double *x, double *y, double *z);
extern void make_hacapk_struct_hmat(void **obj_ptr, double *ztol);
extern void hacapk_mult_Ax(void **obj_ptr, double *x, double *b);
extern void hacapk_solve_Ab(void **obj_ptr, double *b, double *x, double *ztol);
extern void hacapk_assemble_dense_mat(void **obj_ptr, double *Ad, int *n);

int main(int argc, char **argv) {
  
  void *hacapk_struct = NULL;
  int n,i,ierr;
  double *x,*y,*z,ztol,*xv,*bv,*Ad;
  ierr = MPI_Init(&argc,&argv);
  n = 5;
  // Allocate  structure
  init_hacapk_struct(&hacapk_struct, &n);

  /*  */
  
  x = (double *)malloc(sizeof(double)*n);if(!x){fprintf(stderr,"mem error\n");exit(-1);}
  y = (double *)malloc(sizeof(double)*n);if(!y){fprintf(stderr,"mem error\n");exit(-1);}
  z = (double *)malloc(sizeof(double)*n);if(!z){fprintf(stderr,"mem error\n");exit(-1);}

  set_hacapk_struct_coord(&hacapk_struct, x,y,z);


  ztol = 1e-6;
  make_hacapk_struct_hmat(&hacapk_struct, &ztol);
  
  Ad = (double *)malloc(sizeof(double)*n*n);if(!Ad){fprintf(stderr,"mem error\n");exit(-1);}
  hacapk_assemble_dense_mat(&hacapk_struct,Ad,&n);

  
  xv = (double *)malloc(sizeof(double)*n);if(!xv){fprintf(stderr,"mem error\n");exit(-1);}
  bv = (double *)malloc(sizeof(double)*n);if(!bv){fprintf(stderr,"mem error\n");exit(-1);}
  /* RHS */
  for(i=0;i<n;i++)
    bv[i] = 1e3;
  /* inverse solve */
  hacapk_solve_Ab(&hacapk_struct, bv, xv, &ztol);



  /* vector multi */
  hacapk_mult_Ax(&hacapk_struct, xv, bv);
  

  free(xv);free(bv);free(Ad);

  
  free(x);free(y);free(z);
  // Deallocate structure
  deallocate_hacapk_struct(&hacapk_struct);

  return 0;
}
