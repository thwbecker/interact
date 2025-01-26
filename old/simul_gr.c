#include "interact.h"

int main(void)
{
  COMP_PRECISION t,dt,tstop,aspect,tm,m,w,l,rm,
    lmin,lmax,wmin,sd,sdrange,sdmin,mmin,p,mmax,aspecth;
  long seed=-1;
  dt = 0.001;

  
  tstop = 100;
  
  aspect = 3.0;
  aspecth = aspect/2;
  
  /* min and max stress drop (or moment size) */
  sdmin = 1e-5;sdrange=1-sdmin;

  /* min patch size */
  lmin = aspect / 100.0;lmax=aspect-lmin;
  wmin = 1.0 / 100.0;

  /* smallest moment */
  mmin = wmin*wmin*lmin;
  /* largest moment */
  mmax = aspect;

  myrand(&seed);

  tm=rm=0.0;
  t=0.0;
  while(t < tstop){
    /* 
       test quake 
    */
    /* moment pre factor */
    //sd = sdmin + myrandnr(sdrange,&seed);
    sd = 0.5 + mygauss_randnr(1.0,&seed);
    if(sd < sdmin)sd = sdmin;
    if(sd > 1.0)sd = 1.0;
    /* max possible moment */
    m = (tm - rm) * sd;
    /* compute minimum probability */
    p =  m*m*dt*100000/sd;
    /* will earthquake happen? */
    if((m > mmin) && (myrand(&seed) > p)){
      w = 2.0;
      while((w > 1)||(w < wmin)){
	/* pick l */
	//l = lmin + myrandnr(lmax,&seed);
	l = aspecth + mygauss_randnr(aspect,&seed);
	if(l<lmin)l=lmin;
	if(l>aspect)l=aspect;
	w = sqrt(m/l);		/* compute w from m = w^2 l */
      }
      rm += m;
      printf("%12.5e %11g %11g %11g %11g %11g\n",t,l,w,m,rm,tm);
    }
    tm += dt;
    t += dt;
  }
  
  return 0;
}
