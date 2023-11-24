/*

  redefine names for compilers where we need underscores
  when calling fortran routines from C

  $Id: fortran_underscore_defines.h,v 1.4 2004/01/26 01:16:41 becker Exp $


 */
// okada
#define dc3d dc3d_
#define mydc3d mydc3d_
#define dc3d0 dc3d0_
#define dc3d dc3d_
#define dc3d0 dc3d0_
// crouch and starfiled
#define tdd_coeff tdd_coeff_

// pgplot
#define palett palett_
// eispack
#define rs rs_
#define s_rs s_rs_
// LAPACK
#define dgelss dgelss_
#define sgelss sgelss_
#define dgelsd dgelsd_
#define sgelsd sgelsd_
#define dgesv dgesv_
#define sgesv sgesv_
#define dgesvd dgesvd_
#define sgesvd sgesvd_
#define ilwvd ilwvd_
// BLAS
#define ddot ddot_
#define sdot sdot_
#define dgemm dgemm_
#define sgemm sgemm_
#define dgemv dgemv_
#define sgemv sgemv_

#define law_nnls law_nnls_
/*  */
#define tddisphs tddisphs_
#define tddisphs_bird tddisphs_bird_
#define tdstresshs tdstresshs_
#define get_tdcs_base_vectors get_tdcs_base_vectors_
//other
#define my_sincos_deg_ftn my_sincos_deg_ftn_
// weird fix for calling convention
#ifdef LINUX_FORTRAN_CALL_FIX
#define my_sincos_deg_ftn_ my_sincos_deg_ftn__
#define law_nnls_ law_nnls__
#define tdd_coeff_ tdd_coeff__
#endif
