#include "interact.h"

/*
  
  determine after and foreshock bins
  
 */
#define EVENT_PREC float

int main(int argc,char **argv)
{
  int i,n;
  EVENT_PREC mmin,mmax,tmin,tmax;
  FILE *in;
  struct eve{
    EVENT_PREC time,moment;
  }*event;

  switch(argc){
  case 2:{
    break;
  }
  default:{
    fprintf(stderr,"%s cevents.dat\n",argv[0]);
    fprintf(stderr,"sorts for after and foreshocks\n");
    exit(-1);
  }}
  mmin=tmin=FLT_MAX;mmax=tmax=FLT_MIN;
  n=0;
  event=(struct eve *)calloc(1,sizeof(struct eve)*(n+1));
  in=myopen(argv[1],"r");
  while((i=fscanf(in,"%f %*f %f %*f %*f %*f %*f\n",
		  &(event+n)->time,&(event+n)->moment))==2){
    n++;
    if(event[n].time<tmin)tmin=event[n].time;
    if(event[n].time>tmax)tmax=event[n].time;
    if(event[n].moment>mmax)mmax=event[n].moment;
    if(event[n].moment<mmin)mmin=event[n].moment;
    event=(struct eve *)realloc(event,sizeof(struct eve)*(n+1));
    if(!event)
      MEMERROR("main");
  }
  printf("%i\n",i);
  fprintf(stderr,"%s: read %i events, mom min/max: %12.5e/%12.5e\n",
	  argv[0],n,mmin,mmax);
  fclose(in);
  exit(0);
}
