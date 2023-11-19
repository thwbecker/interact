/*
  these are various prototypes for external 
  FORTRAN functions


  $Id: fortran_proto.h,v 1.11 2003/03/02 01:38:27 becker Exp becker $

*/
//
// solver related headers of external functions
//
// $Id: fortran_proto.h,v 1.11 2003/03/02 01:38:27 becker Exp becker $
//
#include "numrec_svd_routines.h"

//
// LAPACK fortran style (not C implementation as in sunperf.h)
// 

extern void ilwvd(int *,int *, int *, int * );

extern void dgelss( int *, int *,int *, double *, int *, 
		    double *, int *,  double *,  double *,
		    int *, double *, int *, int * );

extern void sgelss( int *, int *,int *, float *, int *, 
		    float *, int *,  float *,  float *,
		    int *, float *, int *, int * );

extern void dgelsd( int *, int *,int *, double *, int *, 
		    double *, int *,  double *,  double *,
		    int *, double *, int *, int *,int * );

extern void sgelsd( int *, int *,int *, float *, int *, 
		    float *, int *,  float *,  float *,
		    int *, float *, int *, int *,int * );

extern void dgesvd(char *, char *, int *, int *, 
		    double *, int *, double *, double *, int *, 
		    double *, int *, double *, int *, int *);

extern void sgesvd(char *, char *, int *, int *, 
		   float *, int *, float *, float *, int *, 
		   float *, int *, float *, int *, int *);

extern void dgesv(int *, int *, double *, int *, int *, double *, 
		  int *, int *);

extern void sgesv(int *, int *, float *, int *, int *, float *, 
		  int *, int *);



extern void law_nnls(A_MATRIX_PREC *,int *,int *, int *,
		     A_MATRIX_PREC *,A_MATRIX_PREC *,
		     A_MATRIX_PREC *,A_MATRIX_PREC *,
		     A_MATRIX_PREC *,int *, int *);

/*

  OKADA routines

 */

extern void dc3d(double*,double*,double*,double*,double*,double*,double*,double*,
		 double*,double*,double*,double*,double*,double*,double*,double*,
		 double*,double*,double*,double*,double*,double*,double*,double*,
		 double*,int*);
extern void dc3d0(double*,double*,double*,double*,double*,double*,double*,double*,
		  double*,double*,double*,double*,double*,double*,double*,double*,
		  double*,double*,double*,double*,double*,double*,int*);

extern void palett(int *, float *, float *);

// eispack
extern void  rs(int *, int *, double *,double *, int *, 
		double *, double *,double *, int *);

extern void  s_rs(int *, int *, float *,float *, int *, 
		  float *, float *,float *, int *);

/* crouch and starfield coeff modified */
extern void tdd_coeff(COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,
		      double *,double *,
		      int *,
		      COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,COMP_PRECISION *,
		      COMP_PRECISION *,COMP_PRECISION *,
		      int *);
/*  Nikkhoo, M., Walter, T. R. (2015) half space triangular converted
    to F90 */
extern void tddisphs(COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,
		     COMP_PRECISION *,COMP_PRECISION *,
		     COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *);
extern void tdstresshs(COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,
		       COMP_PRECISION *,COMP_PRECISION *,
		       COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,
		       COMP_PRECISION *);
extern void get_tdcs_base_vectors(COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,
				  COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *);
