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
extern void *cinit_hacapk_struct(int , void *);
/* free */
extern void cdeallocate_hacapk_struct(void* );

/* set coordinates, should be x[N], y[N], z[N] */
extern void cset_hacapk_struct_coord(void* ,double *, double *, double *);
extern void cset_hacapk_inorm(void* , int );
extern void cset_hacapk_eta(void* , double );
extern void cset_hacapk_verbosity(void* , int );

/* get coordinates, should be x[N], y[N], z[N] */
extern void *cget_hacapk_struct_coordp(void* , int );

/* create an H matrix */
extern void cmake_hacapk_struct_hmat(void* , double );
/* given an H matrix, multiply A x and return b = A x */
extern void chacapk_mult_Ax_H(void* , double *, double *);

/* given a dense matrix, solve x = A\b */
extern void chacapk_solve_dense(double *, int , double *, double *);
/* given an H matrix, solve A x = b and return x = A\b */
extern void chacapk_solve_Ab_H(void* , double *, double *, double );

/* assmeble a dense matrix */
extern void chacapk_assemble_dense_mat(void* , double *, int );

extern long cget_hacapk_nnz(void *);
extern void chacapk_mult_Ax_H(void *, double *, double *);
PetscErrorCode MatMult_HACApK(Mat , Vec , Vec );


/* from testing routines */
extern void hacapk_assign_random_coord(double *, double *, double *, int *);

