/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: mysincos.c,v 1.5 2011/01/09 02:02:43 becker Exp $

  uses SUN's sincos routines that are only available in double precision!?
*/
#include "interact.h"
#include <math.h>
void sincos(double, double *, double *);
void sincosf(float, float *, float *);
void sincosl(long double , long double *, long double *);


void my_sincos_deg(COMP_PRECISION *sin_val,COMP_PRECISION *cos_val,
		   COMP_PRECISION alpha_in_degrees)
{
  COMP_PRECISION tmp;
  tmp = alpha_in_degrees * DEG2RAD;
  my_sincos(sin_val,cos_val,tmp);
}


// fortran version
void my_sincos_deg_ftn(double *sin_val,double *cos_val,
		       double *alpha_in_degrees)
{
  double tmp;
  tmp = *alpha_in_degrees * DEG2RAD;
#if !defined(SUN_TRIG_FUNCTIONS) 
  *sin_val = sin(tmp);
  *cos_val = cos(tmp);
#else
  sincos(tmp,sin_val,cos_val);
#endif
}


void my_sincos(COMP_PRECISION *sin_val,
	       COMP_PRECISION *cos_val,
	       COMP_PRECISION arad)
{
#if !defined(SUN_TRIG_FUNCTIONS) 
  *sin_val = sin(arad);
  *cos_val = cos(arad);
#else
#if defined(USE_DOUBLE_PRECISION)
  sincos(arad,sin_val,cos_val);
#else
  double sin_vald,cos_vald,aradd;
  aradd = (double)arad;
  sincos(aradd,&sin_vald,&cos_vald);
  *sin_val = (COMP_PRECISION)sin_vald;
  *cos_val = (COMP_PRECISION)cos_vald;
#endif
  
#endif
}

