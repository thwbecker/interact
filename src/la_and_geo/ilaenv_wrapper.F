c
c     wrapper for the ilaenv block size determination
c     routine from LAPACK, used for the divide and conquer 
c     method for SVD
c
#ifdef USE_DOUBLE_PRECISION
#include "precision_double.h"
#else
#include "precision_single.h"
#endif
#include "fortran_prec.h"

      subroutine ilwvd(m,n,nrhs,nb)
      integer n,m,nrhs,nb

#ifndef A_MATRIX_PREC_IN_DOUBLE
      nb = ilaenv( 1, 'sgelsd', ' ', m,n,nrhs, -1 )
#else
      nb = ilaenv( 1, 'dgelsd', ' ', m,n,nrhs, -1 )
#endif 
      IF( NB.LE.1 ) NB = MAX( 1, N )
      end
