/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: optimize.c,v 2.5 2002/10/08 19:24:44 tbecker Exp $
*/
#include "interact.h"

COMP_PRECISION distsq(struct flt *, struct flt *);

void optimize(struct flt *fault,struct med *medium)
{
  int i,j,k,l,nmax,n,start,inc,hit;
  long seed=0;
  COMP_PRECISION oldpen,oldoldpen,newpen,maxdist;
  nmax=100000;

  fprintf(stderr,"optimize: optimization running\n");

  maxdist=max_dist(fault,medium);
  oldpen=penalty_dist(fault,medium,maxdist);
  
  hit=FALSE;
  n=0;
  do{
    i=myrandi(medium->nrflt-1,&seed);
    j=myrandi(medium->nrflt-1,&seed);
    // flip two faults
    fltswap((fault+j),(fault+i));
    // get new penalty
    newpen=penalty_dist(fault,medium,maxdist);
    if(newpen < oldpen)// better?
      oldpen=newpen;// take it
    else
      fltswap((fault+j),(fault+i));// swap back
    if(n % 200==0){
      if(n==0)
	oldoldpen=oldpen;
      else{
	fprintf(stderr,"optimize: n: %6i pen: %11g change: %g\n",
		n,oldpen,(oldoldpen-oldpen)/oldoldpen);
	if((oldoldpen-oldpen)/oldoldpen < 1e-6)
	  hit++;
	oldoldpen=oldpen;
      }
    }
    n++;
  }while(n<nmax && hit<3);
}

COMP_PRECISION distsq(struct flt *a, struct flt *b)
{
  COMP_PRECISION x,tmpd;
  tmpd = a->x[X] - b->x[X];x = SQUARE(tmpd);
  tmpd = a->x[Y] - b->x[Y];x+= SQUARE(tmpd);
  tmpd = a->x[Z] - b->x[Z];x+= SQUARE(tmpd);
  return x;
}
COMP_PRECISION max_dist(struct flt *fault, 
			struct med *medium)
{
  COMP_PRECISION max,d;
  int i,j;
  max=0.0;
  for(i=0;i<medium->nrflt;i++)
    for(j=i+1;j<medium->nrflt;j++){
      d=distsq((fault+i),(fault+j));
      if(d>max)
	max=d;
    }
  return d;
}
COMP_PRECISION penalty_dist(struct flt *fault, 
			    struct med *medium,
			    COMP_PRECISION max_dist)
{
  COMP_PRECISION sum,p1,p2;
  int i,j;
  sum=0.0;
  for(i=0;i<medium->nrflt;i++)
    for(j=i+1;j<medium->nrflt;j++){
      // normalized distance of faults
      p1=distsq((fault+i),(fault+j))/max_dist;
      // normalized distance from diagonal
      
      p2=((COMP_PRECISION)abs(j-i))/((COMP_PRECISION)medium->nrflt);
      sum += p1*p1 * p2;
    }
	
  return sum;
}
