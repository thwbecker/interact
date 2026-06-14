/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: randgen.c,v 1.7 2011/01/07 07:19:58 becker Exp $

*/
#include "interact.h"

/*

  assign random values, depending on the mode, dx will be standard deviation
  or range

  output: x
  input:
  x0: mean
  xl: min
  dx: spread

*/
void assign_random(COMP_PRECISION *x,COMP_PRECISION x0, 
		   COMP_PRECISION xl,COMP_PRECISION dx,
		   long *seed,int random_mode)
{
  switch(random_mode){
  case UNIFORM_DISTRIBUTED:
    *x = xl + myrandnr(dx,seed);
    break; 
  case GAUSS_DISTRIBUTED:
    *x = x0 + mygauss_randnr(dx,seed); 
    break;
  default:
    fprintf(stderr,"assign_random: error, mode %i undefined\n",
	    random_mode);
    exit(-1);
  }
}


//
// returns uniformly distributed random number 
// 0 < x < 1
//
COMP_PRECISION myrand(long *seed)
{
  return((COMP_PRECISION)RAND_GEN(seed));
}
//
// this returns a random number between 0 < y < x
// uniformly distributed
//
COMP_PRECISION myrandnr(COMP_PRECISION x,long *seed)
{
  if(x == 0.0)
    return 0.0;
  else
    return((COMP_PRECISION)((double)x*RAND_GEN(seed)));
}
//
// same for an integer between 0 <= j <= i
//
int myrandi(int i,long *seed)
{
  if(i==0)
    return 0;
  else
    return((int)(((double)i)*RAND_GEN(seed)+0.5));
}

//
// Gauss distributed number with zero mean and x 
// standard deviation
//
COMP_PRECISION mygauss_randnr(COMP_PRECISION x,long *seed)
{
  if(x == 0.0)
    return 0.0;
  else
    return((COMP_PRECISION)((double)x*gasdev(seed)));
}


// 
// 
// power law statistics P(x) = C x^n for x [x_0,x_1]
// C = (n+1)/(x_1^(n+1) - x_0^(x+1) from normalization
// the power law distributed quantity is given by 
//
// X = ((x_1^(n+1) - x_0^(n+1)) y + x_0^(n+1))^(1/(n+1))
//
// where y is uniformly distributed
//
// set par_change to unity if the x_0, x_1, and n parameters have changed
// from last call
//
COMP_PRECISION mypower_randnr(COMP_PRECISION x0, COMP_PRECISION x1, 
			      COMP_PRECISION n,int par_change,
			      long *seed)
{
  static COMP_PRECISION a,b,nf1,nf2;
  COMP_PRECISION tmp;
  static int init = 0;
  if((!init) || par_change){
    nf1 = n+1.0;
    if(fabs(nf1) < 1e-7){
      fprintf(stderr,"mypower_randnr: error: exponent needs to be != -1 (%g)\n",
	      n);
      exit(-1);
    }
    nf2=1.0/nf1;
    b = pow(x0,nf1);
    a = pow(x1,nf1) - b;
    init = 1;
  }
  tmp  = a * myrand(seed);
  tmp += b;
  return pow(tmp,nf2);
}


//
// get Gaussian distribution with unity variance (or standard deviation)
//
double gasdev(long *idum)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*RAND_GEN(idum)-1.0;
      v2=2.0*RAND_GEN(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
/*


  ran1 random number generator as of numerical recipes in C
  page 280


*/

#define RAN1_IA 16807
#define RAN1_IM 2147483647
#define RAN1_AM (1.0/RAN1_IM)
#define RAN1_IQ 127773
#define RAN1_IR 2836
#define RAN1_NTAB 32
#define RAN1_NDIV (1+(RAN1_IM-1)/RAN1_NTAB)
#define RAN1_EPS 5.0e-15
#define RAN1_RNMX (1.0-RAN1_EPS)

double ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[RAN1_NTAB];
  double temp;
  
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=RAN1_NTAB+7;j>=0;j--) {
      k=(*idum)/RAN1_IQ;
      *idum=RAN1_IA*(*idum-k*RAN1_IQ)-RAN1_IR*k;
      if (*idum < 0) *idum += RAN1_IM;
      if (j < RAN1_NTAB) iv[j] = *idum;
		}
    iy=iv[0];
  }
  k=(*idum)/RAN1_IQ;
  *idum=RAN1_IA*(*idum-k*RAN1_IQ)-RAN1_IR*k;
  if (*idum < 0) *idum += RAN1_IM;
  j=iy/RAN1_NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=RAN1_AM*iy) > RAN1_RNMX) return RAN1_RNMX;
  else return temp;
}
#undef RAN1_IA
#undef RAN1_IM
#undef RAN1_AM
#undef RAN1_IQ
#undef RAN1_IR
#undef RAN1_NTAB
#undef RAN1_NDIV
#undef RAN1_EPS
#undef RAN1_RNMX

/*


  ran2 number generator from numerical recipes, page 282


 */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 5.0e-15
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
  static int ntabp7 = NTAB + 7;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=ntabp7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



