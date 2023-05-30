c
c
c     simple single precision test program for modified EISPACK routines
c
      program main
      parameter(ndim=3)

      real*4 a(ndim*ndim),vec(ndim*ndim),val(ndim),fv1(ndim),fv2(ndim)

      a(1)=1 !a11
      a(2)=2 !a12
      a(3)=3 !a13
      a(5)=4 !a22
      a(6)=5 !a23
      a(9)=6 !a33


      n=ndim
      icvec=1
      call s_rs(n,n,a,val,icvec,vec,fv1,fv2,ierr)
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
c     val         v_1            v_2            v_3
c -0.5157297     -0.7369761     -0.3279853      0.5910090    
c  0.1709151      0.5910091     -0.7369761      0.3279852    
c   11.34481      0.3279853      0.5910090      0.7369762    


      end
