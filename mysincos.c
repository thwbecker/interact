/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: mysincos.c,v 1.5 2011/01/09 02:02:43 becker Exp $

  uses SUN's sincos routines that are only available in double precision!?
*/
#include "interact.h"

void sincos(double, double *, double *);
void sincosf(float, float *, float *);
void sincosl(long double , long double *, long double *);


/* called with degree, convert */
void my_sincos_deg(COMP_PRECISION *sin_val,COMP_PRECISION *cos_val,
		   COMP_PRECISION alpha_in_degrees)
{
  COMP_PRECISION tmp;
  tmp = alpha_in_degrees * DEG2RAD;
  my_sincos(sin_val,cos_val,tmp);
}


// fortran version
void my_sincos_deg_ftnd(double *sin_val,double *cos_val,double *alpha_in_degrees)
{
  double tmp;
  tmp = *alpha_in_degrees * DEG2RAD;
  sincos(tmp,sin_val,cos_val);
}

void my_sincos_ftn(COMP_PRECISION *sin_val,COMP_PRECISION *cos_val,COMP_PRECISION *alpha)
{
#ifdef USE_DOUBLE_PRECISION
  /* double prec */
 #ifdef USE_MKL_SINCOS		/* this was really slow, it seemed */
  const int unity = 1;
  vdSinCos(unity,alpha,sin_val,cos_val);
 #else
  sincos(*alpha,sin_val,cos_val);
 #endif
#else
  /* single prec */
 #ifdef USE_MKL_SINCOS
  const int unity = 1;
  vsSinCos(unity,alpha,sin_val,cos_val);
 #else
  sincosf(*alpha,sin_val,cos_val);
 #endif
#endif
}

void my_sincos(COMP_PRECISION *sin_val,COMP_PRECISION *cos_val,COMP_PRECISION arad)
{
#ifdef USE_DOUBLE_PRECISION
  /* double prec */
 #ifdef USE_MKL_SINCOS
   const int unity = 1;
   vdSinCos(unity,&arad,sin_val,cos_val);
 #else
   sincos(arad,sin_val,cos_val);
 #endif
#else  /* single prec */
 #ifdef USE_MKL_SINCOS
   const int unity = 1;
   vsSinCos(unity,&arad,sin_val,cos_val);
 #else
   sincosf(arad,sin_val,cos_val);
 #endif
#endif
}


void my_sincosd(double *sin_val,double *cos_val,double arad)
{
#ifdef USE_MKL_SINCOS
  const int unity = 1;
  vdSinCos(unity,&arad,sin_val,cos_val);
#else
  sincos(arad,sin_val,cos_val);
#endif
}
void my_sincos_degd(double *sin_val,double *cos_val,double adeg)
{
  double arad;
  arad = adeg * DEG2RAD;
  my_sincosd(sin_val,cos_val,arad);
}


