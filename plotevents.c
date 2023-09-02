/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: plotevents.c,v 2.17 2003/01/07 03:18:21 becker Exp $
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interact.h"

#ifdef USE_PGPLOT

#include "cpgplot.h"

int main(void)
{
  int i,j,n,aflt,*al,naf,oaflt,oc,nriter;
  COMP_PRECISION xmin[3]={1e20,1e20,1e20},
    xmax[3]={-1e20,-1e20,-1e20},dummy,
    otime,mom_thres,optime,tloc[2],alpha,
    sin_dip,cos_dip, corner[4][3];
  float mom,slip[3],time;
  my_boolean del;
  struct flt *fault;
  struct pa *patch;
  FILE *in,*out;
  char tmpstr[100];
  struct med medium[1];
  fault=calloc(1,sizeof(struct flt));
  patch=calloc(1,sizeof(struct pa));
  al=malloc(sizeof(int));

  mom_thres=0.0;
  /* read in geometry */
  in=myopen(GEOMETRY_FILE,"r");
  out=myopen("tmp.dat","w");
  n=0;
  while(fscanf(in,PATCH_CP_FORMAT,&fault[n].x[INT_X],&fault[n].x[INT_Y],&fault[n].x[INT_Z],
	       &fault[n].strike,&fault[n].dip,&fault[n].l,&fault[n].w,
	       &fault[n].group)==8){
    fault=realloc(fault,sizeof(struct flt)*(n+2));
    patch=realloc(patch,sizeof(struct pa)*(n+2));
    alpha= 90.0 - fault[n].strike;
    my_sincos_deg(&fault[n].sin_alpha,&fault[n].cos_alpha,(COMP_PRECISION)alpha);
    my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[n].dip);
    calc_base_vecs(fault[n].t_strike,fault[n].normal,fault[n].t_dip,
		   fault[n].sin_alpha,fault[n].cos_alpha,sin_dip,cos_dip);
    fault[n].area = fault[n].l * fault[n].w * 4.0;
    calculate_corners(corner,(fault+n),&dummy,&dummy);
    for(j=0;j<4;j++){
      patch[n].x[j]=(float)corner[j][INT_X];
      patch[n].y[j]=(float)corner[j][INT_Y];
      patch[n].z[j]=(float)corner[j][INT_Z];
      // extrema
      for(i=0;i<3;i++){
	if(xmax[i]<corner[j][i])
	  xmax[i]=corner[j][i];
	if(xmin[i]>corner[j][i])
	  xmin[i]=corner[j][i];
      }
    }
    patch[n].x[4]=(float)corner[0][INT_X];
    patch[n].y[4]=(float)corner[0][INT_Y];
    patch[n].z[4]=(float)corner[0][INT_Z];
    n++;
  }
  fclose(in);
  medium->nrflt=n;
  fprintf(stderr,"read %i patches from %s\n",n,GEOMETRY_FILE);
  patch=realloc(patch,sizeof(struct pa)*n);
  fault=realloc(fault,sizeof(struct flt)*n);
  /*
    for(j=0;j<3;j++){
    xmin[j]*=1.1;
    xmax[j]*=1.1;
  }
  */
  tloc[INT_X]=(xmax[INT_X]+xmin[INT_X])/2.0;
  tloc[INT_Y]=xmax[INT_Y]*0.9;
  /* create X window */
  if(cpgbeg(i, "/xwindow", 1, 1) != 1)exit(-1);
  // want three pages
  cpgsubp(3,1);
  // device size
  cpgpap(10, .33);
  /* set up axis */
  cpgenv(xmin[INT_X], xmax[INT_X], xmin[INT_Y], xmax[INT_Y], 1, 0);
  cpglab ("x", "y", "x-y");
  cpgenv(xmin[INT_Y], xmax[INT_Y], xmin[INT_Z], xmax[INT_Z], 1, 0);
  cpglab ("y", "z", "y-z");
  cpgenv(xmin[INT_X], xmax[INT_X], xmin[INT_Z], xmax[INT_Z], 1, 0);
  cpglab ("x", "z", "x-z");
  /* 
     draw lines first time around, using the al array 
  */
  for(naf=i=0;i<n;i++) 
    add_to_plotting_list(i,&al,&naf);
  drawset(al,n,1,0,patch);
  /* read in the events */
#ifdef BINARY_PATCH_EVENT_FILE
  in=myopen(EVENT_FILE_BINARY,"r");
#else
  in=myopen(EVENT_FILE_ASCII,"r");
#endif
  naf=0;del=FALSE;
  /* add first event to active event list */
  if(read_patch_event_file(&time,&nriter,&aflt,&mom,slip,in,medium)!=7){
    fprintf(stderr,"could not read events\n");
    exit(-1);
  } 
  printf("%i\n",nriter);
  fprintf(stderr,"%10e %5i %5i %5i %10e %10e %10e %10e\n",
	  time,nriter,aflt,fault[aflt].group,slip[0],slip[1],slip[2],mom);
  fprintf(out,   "%10e %5i %5i %5i %10e %10e %10e %10e\n",
	  time,nriter,aflt,fault[aflt].group,slip[0],slip[1],slip[2],mom);
  
  add_to_plotting_list(aflt,&al,&naf);
  otime=time;
  optime=-1;
  while(read_patch_event_file(&time,&nriter,&aflt,&mom,slip,in,medium)==7){
    fprintf(stderr,"%10e %5i %5i %5i %10e %10e %10e %10e\n",
	    time,nriter,aflt,fault[aflt].group,slip[0],slip[1],slip[2],mom);
    fprintf(out,"%10e %5i %5i %5i %10e %10e %10e %10e\n",
	    time,nriter,aflt,fault[aflt].group,slip[0],slip[1],slip[2],mom);
  
    if(mom>=mom_thres){/* use only if over threshold */
      if(time == otime){
	/* we still add more patches that are active at 
	   the same time */
	if(del){/* still need to delete the old list,
		   ie. only one in new active list*/
	  drawset(al,naf,0,1,patch);
	  drawset(al,naf,1,0,patch);
	  naf=0;
	  add_to_plotting_list(oaflt,&al,&naf);
	  del=FALSE;
	}
	/* add current event to list */
	add_to_plotting_list(aflt,&al,&naf);
      }else{
	/* 
	   now events from a different timestep 
	   update the title of the plot as a preparation
	   for plotting the old active patch list, ie.
	   use old time 
	*/
	if(optime != -1){/* delete old title (time) */
	  cpgqci(&oc);cpgsci(0);
	  sprintf(tmpstr,"time=%10.5f",optime);
	  cpgtext(tloc[INT_X],tloc[INT_Y],tmpstr);
	  cpgsci(oc);
	}
	/* write new time */
	sprintf(tmpstr,"time=%10.5f",otime);
	cpgtext(tloc[INT_X],tloc[INT_Y],tmpstr);
	optime=otime;
	/* draw the old active set  */
	drawset(al,naf,2,1,patch);
	otime=time;
	/* assign the current event to new list after 
	   deleting it */
	oaflt=aflt;
	/* will have to delete it */
	del=TRUE;
	getc(stdin);
      }
    }
  }
  fclose(in);
  cpgend();
  return 0;
}


/* draw a set of faults */
void drawset(int *al, int naf, int color,int fill,
	     struct pa *patch)
{
  int i,oc,j;
  /* get old color */
  cpgqci(&oc);
  /* set newcolor */
  cpgsci(color);
  for(j=0;j<3;j++){
    cpgpanl(1+j, 1);
    /* open buffer */
    cpgbbuf();
    if(color==1){
      if(j==0)
	for(i=0;i<naf;i++)
	  cpgline(5,patch[al[i]].x, patch[al[i]].y);
      else if(j==1)
	for(i=0;i<naf;i++)
	  cpgline(5,patch[al[i]].y, patch[al[i]].z);
      else 
	for(i=0;i<naf;i++)
	  cpgline(5,patch[al[i]].x, patch[al[i]].z);
    }else{
      if(j==0)
	for(i=0;i<naf;i++)
	  cpgpoly(5,patch[al[i]].x, patch[al[i]].y);
      else if(j==1)
	for(i=0;i<naf;i++)
	  cpgpoly(5,patch[al[i]].y, patch[al[i]].z);
      else 
	for(i=0;i<naf;i++)
	  cpgpoly(5,patch[al[i]].x, patch[al[i]].z);
    } 
    cpgebuf();
  }
  /* set color back to old value */
  cpgsci(oc);
}
#endif
