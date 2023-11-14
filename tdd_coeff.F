#ifdef USE_DOUBLE_PRECISION
#include "precision_double.h"
#define EPS_COMP_PREC_FOR_NAN 5.0e-14
#else
#include "precision_single.h"
#define EPS_COMP_PREC_FOR_NAN 5.0e-6
#endif
c
c       this subroutine calculates the stresses and displacemnts
c       in a plane or half-plane given constant slip on a segment
c       it is an implementation of crouch & starfields (5.4 and 7.4)
c       formulae and based on their original fortran code "tdd"
c       the original of this code came from Yuri Fialko in 03/03
c
c       WARNING: normal motion sign convention here is different
c       from the one in interact. we will flip the sign in eval_2dsegment
c
c
c       changes by Thorsten Becker:
c
c       - the atan problem of the f3 function (see C&S, p. 51)
c         is dealt with differently
c       
c       $Id: tdd_coeff.F,v 1.9 2003/03/09 02:13:23 becker Exp $
c
c       input: 
c
c       cx, cy: center of segment 
c       x, y: observational point
c       cosb, sinb: cos and sin of angle beta betwenn x axis and local
c                   x (strike) axis
c       a: half-length
c       ihalf: while-plane: 0 half-plane: 1
c
c       pr: poisson ration nu
c       pr1: (1 - 2 pr) 
c       pr2: 2(1 - pr)
c       con: 1. / (4. * pi * (1. - pr))
c       cons = e/(1. + pr) (this is 2 \mu)
c       chi = 3. - 4. * pr
c
c       the routine returns the stresses sxx, sxy, and syy
c       as well as displacements ux, uy for strike and normal
c       components
c
c
	subroutine tdd_coeff(x,y,cx,cy,a,cosb,sinb,ihalf,
     &                       pr,pr1,pr2,con,cons,chi,
     &                       sxxs,sxxn,syys,syyn,sxys,sxyn,
     &                       uxs,uxn,uys,uyn,iret)
c
	implicit C_PREC(a-h,o-z)
	implicit integer(i-n)
	real*8 cosb,sinb
	
	cos2b = cosb * cosb - sinb * sinb
	sin2b = 2.0d0 * sinb * cosb
	cosb2 = cosb * cosb
	sinb2 = sinb * sinb

c       \bar{x} and \bar{y} (eq. 7.4.2)
	xb =  (x - cx) * cosb + (y - cy) * sinb
	yb = -(x - cx) * sinb + (y - cy) * cosb
c
	xbpa = xb + a
	xbma = xb - a
c
	r1s = (xbma) * (xbma) + yb * yb
	r2s = (xbpa) * (xbpa) + yb * yb
c
c       check for NaN
c
	if((abs(r1s).le.EPS_COMP_PREC_FOR_NAN).or.
     &      (abs(r2s).le.EPS_COMP_PREC_FOR_NAN))then
	   iret=1
	else
	   iret=0
	endif
	if(iret.eq.0)then
c
c       do calculation only for non-Nan quantities 
c       
	   fl1 = 0.5d0 * log(r1s)
	   fl2 = 0.5d0 * log(r2s)
	   fb2 = con * (fl1 - fl2)
	   call det_atan_f3(xbpa,xbma,xb,yb,a,con,fb3)
 20	   fb4 = con * (yb/r1s - yb/r2s)
	   fb5 = con * ((xbma)/r1s - (xbpa)/r2s)
	   fb6 = con * (((xbma)**2-yb*yb)/r1s**2-((xbpa)**2-yb*yb)/r2s**2)
	   fb7 = 2.*con*yb*((xbma)/r1s**2-(xbpa)/r2s**2)
c
c       the `actual' displacements and stresses
c       
	   uxs = -pr1*sinb*fb2+pr2*cosb*fb3+yb*(sinb*fb4-cosb*fb5)
	   uxn = -pr1*cosb*fb2-pr2*sinb*fb3-yb*(cosb*fb4+sinb*fb5)

	   uys =  pr1*cosb*fb2+pr2*sinb*fb3-yb*(cosb*fb4+sinb*fb5)
	   uyn = -pr1*sinb*fb2+pr2*cosb*fb3-yb*(sinb*fb4-cosb*fb5)
c       
	   sxxs = cons*(2.0d0*cosb2*fb4+sin2b*fb5+yb*
     &           (cos2b*fb6-sin2b*fb7))
	   sxxn = cons*(-fb5+yb*(sin2b*fb6+cos2b*fb7))
	   
	   syys = cons*(2.0d0*sinb2*fb4-sin2b*fb5-yb*
     &            (cos2b*fb6-sin2b*fb7))
	   syyn = cons*(-fb5-yb*(sin2b*fb6+cos2b*fb7))
	   
	   sxys = cons*(sin2b*fb4-cos2b*fb5+yb*(sin2b*fb6+cos2b*fb7))
	   sxyn = cons*(-yb*(cos2b*fb6-sin2b*fb7))
	   
c       
c       Now for the image plus surface part in case this is a half-plane
c       
	   if (ihalf .eq. 1) then
	      sin3b = sin2b*cosb + cos2b*sinb
	      cos3b = cos2b*cosb - sin2b*sinb
	      sin4b = 2*sin2b*cos2b
	      cos4b = cos2b**2 - sin2b**2
c
c       new coordinate transformat (eq. 7.4.6)
c      
	      xbp = (x - cx) * cosb - (y + cy) * sinb
	      ybp = (x - cx) * sinb + (y + cy) * cosb
c       
	      xbppa = xbp + a
	      xbpma = xbp - a
c
	      r1sp = (xbpma) * (xbpma) + ybp * ybp
	      r2sp = (xbppa) * (xbppa) + ybp * ybp
	      if((abs(r1sp).le.EPS_COMP_PREC_FOR_NAN).or.
     & 		   (abs(r2sp).le.EPS_COMP_PREC_FOR_NAN))then
		 iret=1
	      endif
	      if(iret.eq.0)then
		 fl1p = 0.5d0 * log(r1sp)
		 fl2p = 0.5d0 * log(r2sp)
		 fbp2 = con * (fl1p - fl2p)
		 call det_atan_f3(xbppa,xbpma,xbp,ybp,a,con,fbp3)
 120		 fbp4 = con * (ybp/r1sp - ybp/r2sp)
		 fbp5 = con * ((xbpma)/r1sp - (xbppa)/r2sp)
		 fbp6 = con * (((xbpma)**2-ybp*ybp)/r1sp**2-
     1                ((xbppa)**2-ybp*ybp)/r2sp**2)
		 fbp7 = 2.0d0*con*ybp*((xbpma)/r1sp**2-(xbppa)/r2sp**2)
		 fbp8 = 2.0d0*con*ybp*((ybp**2-3*(xbpma)**2)/((xbpma)**2+ybp**2)**3-
     1                     (ybp**2-3*(xbppa)**2)/((xbppa)**2+ybp**2)**3)
		 fbp9 = 2.0d0*con*((xbpma)/((xbpma)**2+ybp**2)**2-
     1                 (xbppa)/((xbppa)**2+ybp**2)**2+
     2                 4.0d0*ybp**2*((xbppa)/((xbppa)**2+ybp**2)**3-
     3                           (xbpma)/((xbpma)**2+ybp**2)**3))
c

		 uxpds = pr1*sinb*fbp2-pr2*cosb*fbp3+
     1          (chi*(y*sin2b-yb*sinb)+2.0d0*y*sin2b)*fbp4+
     2          (chi*(y*cos2b-yb*cosb)-y*(1.0d0-2.0d0*cos2b))*fbp5+
     3          2.0d0*y*(y*sin3b-yb*sin2b)*fbp6-
     4          2.0d0*y*(y*cos3b-yb*cos2b)*fbp7

		 uxpdn = pr1*cosb*fbp2+pr2*sinb*fbp3-
     1          (chi*(y*cos2b-yb*cosb)-y)*fbp4+ 
     2          chi*(y*sin2b-yb*sinb)*fbp5-
     3          2.0d0*y*(y*cos3b-yb*cos2b)*fbp6-
     4          2.0d0*y*(y*sin3b-yb*sin2b)*fbp7

		 uypds = -pr1*cosb*fbp2-pr2*sinb*fbp3-
     1          (chi*(y*cos2b-yb*cosb)+y*(1.0d0-2.0d0*cos2b))*fbp4+ 
     2          (chi*(y*sin2b-yb*sinb)-2.0d0*y*sin2b)*fbp5+
     3          2.0d0*y*(y*cos3b-yb*cos2b)*fbp6+
     4          2.0d0*y*(y*sin3b-yb*sin2b)*fbp7
		 uypdn = pr1*sinb*fbp2-pr2*cosb*fbp3-
     1          chi*(y*sin2b-yb*sinb)*fbp4-
     2          (chi*(y*cos2b-yb*cosb)+y)*fbp5+
     3          2.0d0*y*(y*sin3b-yb*sin2b)*fbp6-
     4          2.0d0*y*(y*cos3b-yb*cos2b)*fbp7
c
		 sxxpds = cons*(fbp4-3.0d0*(cos2b*fbp4-sin2b*fbp5)+
     1                 (2.0d0*y*(cosb-3.d0*cos3b)+3.0d0*yb*cos2b)*fbp6+
     2                 (2.0d0*y*(sinb-3.d0*sin3b)+3.0d0*yb*sin2b)*fbp7-
     3          2.0d0*y*(y*cos4b-yb*cos3b)*fbp8-
     4                 2.0d0*y*(y*sin4b-yb*sin3b)*fbp9)
		 sxxpdn = cons*(fbp5+
     1                 (2.0d0*y*(sinb-2.0d0*sin3b)+3.0d0*yb*sin2b)*fbp6-
     2                 (2.0d0*y*(cosb-2.0d0*cos3b)+3.0d0*yb*cos2b)*fbp7-
     3          2.0d0*y*(y*sin4b-yb*sin3b)*fbp8+
     4          2.0d0*y*(y*cos4b-yb*cos3b)*fbp9)

		 syypds = cons*(fbp4-(cos2b*fbp4-sin2b*fbp5)-
     1                 (4.0d0*y*sinb*sin2b-yb*cos2b)*fbp6+
     2                 (4.0d0*y*sinb*cos2b+yb*sin2b)*fbp7+
     3          2.0d0*y*(y*cos4b-yb*cos3b)*fbp8+
     4                 2.0d0*y*(y*sin4b-yb*sin3b)*fbp9)
		 syypdn = cons*(fbp5-
     1                 (2.0d0*y*sinb-yb*sin2b)*fbp6+
     2                 (2.0d0*y*cosb-yb*cos2b)*fbp7+
     3          2.0d0*y*(y*sin4b-yb*sin3b)*fbp8-
     4                  2.0d0*y*(y*cos4b-yb*cos3b)*fbp9)

		 sxypds = cons*(sin2b*fbp4+cos2b*fbp5+
     1                 (2.0d0*y*sinb*(1.0d0+4.0d0*cos2b)-yb*sin2b)*fbp6+
     2                 (2.0d0*y*cosb*(3.0d0-4.0d0*cos2b)+yb*cos2b)*fbp7+
     3          2.0d0*y*(y*sin4b-yb*sin3b)*fbp8-
     4                  2.0d0*y*(y*cos4b-yb*cos3b)*fbp9)
		 sxypdn = cons*(
     1                 (4.0d0*y*sinb*sin2b+yb*cos2b)*fbp6-
     2                 (4.0d0*y*sinb*cos2b-yb*sin2b)*fbp7-
     3          2.0d0*y*(y*cos4b-yb*cos3b)*fbp8-
     4                  2.0d0*y*(y*sin4b-yb*sin3b)*fbp9)
c
c       add to actual displacements
c       
		 uxs = uxs + uxpds
		 uxn = uxn + uxpdn
		 
		 uys = uys + uypds
		 uyn = uyn + uypdn
c       
		 sxxs = sxxs + sxxpds
		 sxxn = sxxn + sxxpdn
		 
		 syys = syys + syypds
		 syyn = syyn + syypdn
		 
		 sxys = sxys + sxypds
		 sxyn = sxyn + sxypdn
	      endif ! end half-plane non-NaN
	   endif ! end half-plane
	endif ! end non-NaN

	return
	end
c
c       deal with special cases for f3 function
c       
c
	subroutine det_atan_f3(xpa,xma,x,y,a,con,f3)
	implicit C_PREC(a-h,o-z)
	implicit integer(i-n)
	data pi/3.14159265358979324d0/


	if (abs(y).le.EPS_COMP_PREC_FOR_NAN) then
	   if (abs(x).lt.a) then
	      if(y.ge.0)then
		 f3 = con * pi
	      else
		 f3 = -con * pi
	      endif
	   else
	      f3 = 0.0d0
	   endif
	else
	   f3 = -con * (atan(xpa/y) - atan(xma/y))
	endif


	return
	end
