#include <math.h>
#include <stdio.h>
#include "interact.h"
/*


see also the period program in progs/src/misc !


*/

#include "period.h"

//
// analize data sequence with time in first column and y in second
// reads from stdin
//
void find_range(int *,int *,SPECA_PREC ,SPECA_PREC ,SPECA_PREC ,SPECA_PREC ,
		SPECA_PREC *,SPECA_PREC *,SPECA_PREC *,int);

int main(int argc,char **argv)
{
  SPECA_PREC *t,*xm,tmin,tmax,trange,hifac,ofac,f1,f2,*px,*py,prob,t1,t2,dt,
    tmpt,tmpx,frange=0,fstep=0,fmid;
  int i,n,np,it1,it2,nout,jmax,nadd,mode,fslide=0;
  // oversampling factor
  ofac=6.0;
  // max frequency as fraction of nyquist frequency
  hifac=0.75;
  
  f1=0.25;
  f2=0.75;
  dt=0.75;
  mode=1;

  if(argc >= 2)
    sscanf(argv[1],SPPF,&f1);
  if(argc >= 3)
    sscanf(argv[2],SPPF,&f2);
  if(argc >= 4)
    sscanf(argv[3],SPPF,&dt);
  if(argc >= 5)
    sscanf(argv[4],"%i",&mode);

  if((f1<0) && (f2<0)){
    fslide = 1;
    fstep  = -f1;
    frange = -f2/2.0;
  }

  // read in times and moments of events
  t=(SPECA_PREC *)malloc(sizeof(SPECA_PREC));
  xm=(SPECA_PREC *)malloc(sizeof(SPECA_PREC));

  n=0;i=1;nadd=0;
  fprintf(stderr,"%s: reading t f(t) from stdin, min spacing dt: %g mode: %i\n",
	  argv[0],dt,mode);
  while(fscanf(stdin,SPPF2,&tmpt,&tmpx)==2){
    if(n>1){
      while(tmpt > t[n-1]+dt){// fillin zeroes
	t[n] =t[n-1]+dt;
	xm[n]=0.0;
	n++;i++;
	t=(SPECA_PREC *)realloc(t,sizeof(SPECA_PREC)*i);
	xm=(SPECA_PREC *)realloc(xm,sizeof(SPECA_PREC)*i);
	if(!t || !xm)MEMERROR("mspectral, t xm");
	nadd++;
      }
    }
    t[n] =tmpt;
    xm[n]=tmpx;
    n++;i++;
    t=(SPECA_PREC *)realloc(t,sizeof(SPECA_PREC)*i);
    xm=(SPECA_PREC *)realloc(xm,sizeof(SPECA_PREC)*i);
    if(!t || !xm)MEMERROR("mspectral, t xm");
  }
  t=(SPECA_PREC *)realloc(t,sizeof(SPECA_PREC)*n);
  xm=(SPECA_PREC *)realloc(xm,sizeof(SPECA_PREC)*n);
  tmin=t[0];tmax=t[n-1];trange=tmax-tmin;
  fprintf(stderr,"%s: read in %i events from time %g to %g, added %i  zeroes\n",
	  argv[0],n,tmin,tmax,nadd);
  //
  // workspace
  //
  np=(int)((ofac*hifac)/2.0*(SPECA_PREC)n+0.5)*22;
  px=(SPECA_PREC *)malloc(sizeof(SPECA_PREC)*np);
  py=(SPECA_PREC *)malloc(sizeof(SPECA_PREC)*np);
  if(!px || !py)MEMERROR("mspectral, px, py");
  it1=0;it2=n-1;


  if(!fslide){
    // analysis from f1 to f2 as fractions of the whole range, 
    // have to move inward from here
    find_range(&it1,&it2,f1,f2,tmin,trange,t,&t1,&t2,n);
    fprintf(stderr,"%s: calculating from %g (%g%%) to %g (%g%%) with %i samples\n",
	    argv[0],t1,f1*100.,t2,f2*100.,it2-it1+1);
    if(mode == 0)
      period((t+it1-1),(xm+it1-1),(it2-it1+1),ofac,hifac,(px-1),(py-1),np,&nout,
	     &jmax,&prob);
    else
      fasper((t+it1-1),(xm+it1-1),(it2-it1+1),ofac,hifac,(px-1),(py-1),np,&nout,
	     &jmax,&prob);
    for(i=0;i<nout;i++)
      fprintf(stdout,"%g %g\n",px[i],py[i]);
  }else{
    for(fmid=0.0;fmid<=1.0;fmid += fstep){
      f1=fmid-frange;
      f2=fmid+frange;
      find_range(&it1,&it2,f1,f2,tmin,trange,t,&t1,&t2,n);
      fprintf(stderr,"%s: calculating from %g (%g%%) to %g (%g%%) with %i samples\n",
	      argv[0],t1,f1*100.,t2,f2*100.,it2-it1+1);
      if(mode == 0)
	period((t+it1-1),(xm+it1-1),(it2-it1+1),ofac,hifac,(px-1),(py-1),np,&nout,
	       &jmax,&prob);
      else
	fasper((t+it1-1),(xm+it1-1),(it2-it1+1),ofac,hifac,(px-1),(py-1),np,&nout,
	       &jmax,&prob);
      for(i=0;i<nout;i++)
	fprintf(stdout,"%g %g %g\n",px[i],fmid,py[i]);
    }
  }
  exit(0);
}
void find_range(int *it1,int *it2,SPECA_PREC f1,SPECA_PREC f2,
		SPECA_PREC tmin,SPECA_PREC trange,
		SPECA_PREC *t,SPECA_PREC *t1,SPECA_PREC *t2,
		int n)
{
  if(f1<0)f1=0.0;
  if(f2>1)f2=1.0;
  *t1=tmin+trange*f1;
  *t2=tmin+trange*f2;
  *it1=0;
  *it2=n-1;
  while(t[*it1]< *t1)
    *it1 += 1;
  while(t[*it2]> *t2)
    *it2 -= 1;
}

