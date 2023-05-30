      real*4 function s_pythag(a,b)
      real*4 a,b
c
c     finds sqrt(a**2+b**2) without overflow or destructive underflow
c
      real*4 p,r,s,t,u
      p = max(abs(a),abs(b))
      if (p .eq. 0.0) go to 20
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0 + r
         if (t .eq. 4.0) go to 20
         s = r/t
         u = 1.0 + 2.0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 s_pythag = p
      return
      end
