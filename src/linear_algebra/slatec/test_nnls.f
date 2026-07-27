      program
      parameter (n=3, m=3)
      real*8 a,prgopt,x,rnorm,work
      integer*4 ma,l,i,mode,iwork
      dimension a(m,n+1),work(m+5*m),iwork(m+n),x(m)
      prgopt=1.0
      for

      call dwnnls(a,m,0,m,n,0,prgopt,x,rnorm,mode,iwork,
     &     work)
      print *,mode
      do i=1,n
         print *,x(i)
      enddo
      
      




      end 
