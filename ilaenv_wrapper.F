c
c     wrapper for the ilaenv block size determination
c     routine from LAPACK, used for the divide and conquer 
c     method for SVD
c
      subroutine ilwvd(m,n,nrhs,nb)
      integer n,m,nrhs,nb

#ifdef A_MATRIX_SINGLE_PREC
      nb = ilaenv( 1, 'sgelsd', ' ', m,n,nrhs, -1 )
#else
      nb = ilaenv( 1, 'dgelsd', ' ', m,n,nrhs, -1 )
#endif 
      IF( NB.LE.1 ) NB = MAX( 1, N )
      end
