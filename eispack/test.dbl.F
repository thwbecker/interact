c
c
c     simple double precision test program for modified EISPACK routines
c
      program main
      parameter(ndim=3)

      real*8 a(ndim*ndim),vec(ndim*ndim),val(ndim),fv1(ndim),fv2(ndim)

      a(1)=1 !a11
      a(2)=2 !a12
      a(3)=3 !a13
      a(5)=4 !a22
      a(6)=5 !a23
      a(9)=6 !a33


      n=ndim
      icvec=1
      call rs(n,n,a,val,icvec,vec,fv1,fv2,ierr)
      if(ierr.ne.0)then
         print *,'error: ',ierr
         stop
      endif
      do i=1,3 
         print *,val(i),vec((i-1)*3+1),vec((i-1)*3+2),vec((i-1)*3+3)
      enddo
c
c     the output should be
c
c -0.515729471589257      -0.736976229099578      -0.327985277605682     
c  0.591009048506104     
c  0.170915188827179       0.591009048506103      -0.736976229099578     
c  0.327985277605682     
c   11.3448142827621       0.327985277605682       0.591009048506104     
c  0.736976229099578     


      end


