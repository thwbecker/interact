#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// Declare Fortran functions as external C functions
//

/* 
   
   kernel function, used by both f90 and C code

 */
/* this has to be consistent with the declaration  in m_HACApK_calc_entry_ij.f90 */
double ckernel_func(int , int , void *);

/* 
   
   interface/wrapper functions 

 */
/* initialize hacapk structure with N x N matrix */
extern void *cinit_hacapk_struct(int n, void *ckernel_par);
/* free */
extern void cdeallocate_hacapk_struct(void* c_pointer);

/* set coordinates, should be x[N], y[N], z[N] */
extern void cset_hacapk_struct_coord(void* c_pointer,
				     double *x, double *y, double *z);

/* get coordinates, should be x[N], y[N], z[N] */
extern void *cget_hacapk_struct_coordp(void* c_pointer, int dim);

/* create an H matrix */
extern void cmake_hacapk_struct_hmat(void* c_pointer, double ztol);
/* given an H matrix, multiply A x and return b = A x */
extern void chacapk_mult_Ax_H(void* c_pointer, double *x, double *b);

/* given a dense matrix, solve x = A\b */
extern void chacapk_solve_dense(double *Ad, int n, double *b, double *x);
/* given an H matrix, solve A x = b and return x = A\b */
extern void chacapk_solve_Ab_H(void* c_pointer, double *b, double *x, double ztol);

/* assmeble a dense matrix */
extern void chacapk_assemble_dense_mat(void* c_pointer, double *Ad, int n);




/* from testing routines */
extern void hacapk_assign_random_coord(double *x, double *y, double *z, int *n);

