/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: mysincos.c,v 1.5 2011/01/09 02:02:43 becker Exp $

  uses SUN's sincos routines that are only available in double precision!?
*/
#include "interact.h"
#include <math.h>
#ifdef SUN_TRIG_FUNCTIONS
#include <sunmath.h>
#endif

void my_sincos_deg(COMP_PRECISION *sin_val,COMP_PRECISION *cos_val,
		   COMP_PRECISION alpha_in_degrees)
{
#if !defined(SUN_TRIG_FUNCTIONS) || !defined(USE_DOUBLE_PRECISION)
  COMP_PRECISION tmp;
  tmp = alpha_in_degrees * DEG2RAD;
  *sin_val = sin(tmp);
  *cos_val = cos(tmp);
#else
  // using double prec and sun libs
  sincosd(alpha_in_degrees, sin_val, cos_val);
#endif
}


// fortran version
void my_sincos_deg_ftn(double *sin_val,double *cos_val,
		       double *alpha_in_degrees)
{
#if !defined(SUN_TRIG_FUNCTIONS) || !defined(USE_DOUBLE_PRECISION)
  double tmp;
  tmp = *alpha_in_degrees * DEG2RAD;
  *sin_val = sin(tmp);
  *cos_val = cos(tmp);
#else
  // using double prec and sun libs
  sincosd(*alpha_in_degrees, sin_val, cos_val);
#endif
}


void my_sincos(COMP_PRECISION *sin_val,COMP_PRECISION *cos_val,
	       COMP_PRECISION arad)
{
#if !defined(SUN_TRIG_FUNCTIONS) || !defined(USE_DOUBLE_PRECISION)
  *sin_val = sin(arad);
  *cos_val = cos(arad);
#else
  sincos(arad,sin_val,cos_val);
#endif
}

