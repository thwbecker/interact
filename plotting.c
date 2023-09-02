/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: plotting.c,v 2.13 2004/03/16 01:11:05 becker Exp $


  several help routines for PGPLOT access


*/
#ifdef USE_PGPLOT
#include "interact.h"
#include "cpgplot.h"
void update_plots(struct med *medium, struct flt *fault)
{
  static float scroll_inc;
  static COMP_PRECISION ticked_time,l_ticked_time;
  static my_boolean init=FALSE;
  int i;
  if(!init){
    scroll_inc=(float)medium->x_scroll_inc;
    ticked_time=0.0;
    l_ticked_time=0.0;
    init=TRUE;
  }
  /* first window */
  if(medium->x_window[0]){
    cpgslct(medium->x_window[0]);
    if(medium->time - medium->x_plot_time >= 
       medium->x_plot_interval){
      /* erase old activated faults, start new */
      for(i=0;i<medium->nr_active_flt;i++){
	plot_patch(medium->active_flt_list[i],medium,fault,1,medium->lw_def);
      }
      /* delete old title */
      plot_time_label(medium, fault,0);
      /* update time and list */
      medium->x_plot_time=medium->time;
      /* write new title */
      plot_time_label(medium, fault,1);
      medium->nr_active_flt=0;
    }
  }
  /* second window */
  if(medium->x_window[1]){
    cpgslct(medium->x_window[1]);
    /* draw tick marks at time axis */
    if(medium->time - medium->time_tic_interval >=
       ticked_time){
      plot_time_tics(medium,fault,1);
      ticked_time= medium->time;
    }
    if(medium->time - medium->l_time_tic_interval >=
       l_ticked_time){
      plot_time_tics(medium,fault,2);
      l_ticked_time= medium->time;
    }
    
    /* scroll screen */
    if(medium->time - medium->x_scroll_time >=
       medium->x_scroll_interval){
      cpgscrl((float)0.0,scroll_inc); 
      medium->x_scroll_time = medium->time;
    }
  }
}


void init_plot_window(struct med *medium, struct flt *fault)
{
  int i,colorscale;
  float bright, contrast, range[2];
  /* create X windows */
  if(medium->nr_xwindows > MAX_NR_X_WINDOWS){
    fprintf(stderr,"too many X windows open attempts, max is %i\n",
		MAX_NR_X_WINDOWS);
    exit(-1);
  }
  /* open all */
  for(i=0;i<medium->nr_xwindows;i++){
    if((i==2 && medium->nrgrp==1)||(i!=2)){
      medium->x_window[i]=cpgopen("/xwindow");
      if(medium->x_window[i]<=0){
	fprintf(stderr,"error opening window %i\n",
		i);
	exit(-1);
      }
    }
  }
  if(medium->x_window[0]){
    /* 
       first window 
    */
    cpgslct(medium->x_window[0]);
    /* obtain base line width for terminal */
    cpgqlw(&medium->lw_def);
    /* set up axes */
    for(i=0;i<2;i++)
      range[i]= (medium->xmax[i] - medium->xmin[i])*0.05;

    cpgenv(medium->xmin[INT_X]-range[INT_X], 
	   medium->xmax[INT_X]+range[INT_X], 
	   medium->xmin[INT_Y]-range[INT_X], 
	   medium->xmax[INT_Y]+range[INT_Y], 1, 0);
    cpglab("x", "y", "map view");
    /* determine position of time label */
    medium->tloc[INT_X]=((medium->xmax[INT_X]+medium->xmin[INT_X])/2.0);
    medium->tloc[INT_Y]=(medium->xmax[INT_Y]+range[INT_Y]*2);
    for(i=0;i<medium->nrflt;i++)
      plot_patch(i,medium, fault,1,medium->lw_def);
    /* write title */
    plot_time_label(medium, fault,1);
    /* initialize active fault list */
    medium->active_flt_list=malloc(sizeof(int));
    medium->nr_active_flt=0;
  }
  if(medium->x_window[1]){
    /* 
       second window 
    */
    cpgslct(medium->x_window[1]);
    cpgenv(medium->xmin[INT_X], medium->xmax[INT_X], 
	   0.0, medium->stop_time, 0, 2);
    
    cpglab("x", "time", "x-t view");
    
  }
  if(medium->x_window[2] && medium->nrgrp==1){
    /* 
       third window 
    */
    cpgslct(medium->x_window[2]);
    fprintf(stderr,"init_plot_window: initialize moment window with bounds x(%g, %g) y(%g, %g) for n: %i m: %i\n",
	    0.0,((float)medium->grp0_m)*medium->lmean,
	    0.0,((float)medium->grp0_n)*medium->wmean,medium->grp0_n,medium->grp0_m);
    cpgenv(0.0,((float)medium->grp0_m+1)*2.0*medium->lmean,
	   0.0,((float)medium->grp0_n+1)*2.0*medium->wmean, 1, 0);
    cpglab("x", "z", "total moment release");
    if((medium->momrel=(float *)calloc(medium->grp0_n*medium->grp0_m,
				       sizeof(float)))== NULL)
      MEMERROR("init_plot_window");
    medium->moment_array_init=TRUE;

    bright = 0.5;
    contrast  = 1.0;
    colorscale = 3;
    palett(&colorscale, &contrast, &bright);

  }
  if(medium->x_window[0])
    cpgslct(medium->x_window[0]);
}


void close_plot_window(struct med *medium, 
		       struct flt *fault)
{
  int i;
  for(i=0;i<medium->nr_xwindows;i++){
    cpgslct(medium->x_window[i]);
    cpgclos();
  }
}

void plot_time_label(struct med *medium, 
		     struct flt *fault,int color)
{
  char tmpstr[STRLEN];
  int oc;
  cpgqci(&oc);
  cpgsci(color);
  sprintf(tmpstr,"time=%10.5f",medium->x_plot_time);
  cpgtext(medium->tloc[INT_X],medium->tloc[INT_Y],tmpstr);
  cpgsci(oc);
}
void plot_patch(int flt,struct med *medium, struct flt *fault,int color, int width)
{
  float x[2],y[2],cos_alpha,sin_alpha;
  /* determine geometry */
  cos_alpha=fault[flt].cos_alpha*fault[flt].l;
  sin_alpha=fault[flt].sin_alpha*fault[flt].l;
  x[0]=fault[flt].x[INT_X]+cos_alpha;
  y[0]=fault[flt].x[INT_Y]+sin_alpha;
  x[1]=fault[flt].x[INT_X]-cos_alpha;
  y[1]=fault[flt].x[INT_Y]-sin_alpha;
  /* draw patch */
  cpgsci(color);
  cpgslw(width);
  cpgline(2,x,y);
}
void plot_projected_patch(int flt,struct med *medium, struct flt *fault,int color,
			  int width)
{
  float x[2],y[2],cos_alpha;
  /* determine geometry of window */
  cpgqwin(x,(x+1),y,y+1);
  /* y coordinates are at the minimum */
  y[1]=y[0]=y[0]+(float)medium->x_scroll_inc/2.0;
  /* x coordinates depend on fault */
  cos_alpha=fault[flt].cos_alpha * fault[flt].l;
  x[0]=fault[flt].x[INT_X]-cos_alpha;
  x[1]=fault[flt].x[INT_X]+cos_alpha;
  /* draw patch */
  cpgsci(color);
  cpgslw(width);
  cpgline(2,x,y);
}
void plot_time_tics(struct med *medium, struct flt *fault,float fac)
{
  float x1,x2,y1,y2,xr,x[2],y[2];
  cpgqwin(&x1,&x2,&y1,&y2);
  xr=(x2-x1)/100.0*fac;
  x[0]=x1;
  x[1]=x1+xr;
  y[0]=y[1]=y1+(float)medium->x_scroll_inc/2.0;;
  cpgsci(1);
  cpgline(2,x,y);
  x[0]=x2-xr;
  x[1]=x2;
  cpgline(2,x,y);
}
/* add an active patch to list */
void add_to_plotting_list(int aflt,int **al,int *naf)
{
  *al=realloc(*al,sizeof(int)*(*naf+2));
  (*al)[*naf]=aflt;
  *naf += 1;
}

void plot_quake(int r_flt,struct med *medium, 
		struct flt *fault)
{
  if(medium->x_window[0]){
    /* add to list */
    
    add_to_plotting_list(r_flt,
			 &medium->active_flt_list,
			 &medium->nr_active_flt);
    
    /* plot on first window */
    cpgslct(medium->x_window[0]);
    plot_patch(r_flt,medium,fault,2,medium->lw_def);
  }
  if(medium->x_window[1]){
    /* and on second */
    cpgslct(medium->x_window[1]);
    plot_projected_patch(r_flt,medium,fault,2,
			 medium->lw_def);
  }
}



void plot_moment_array(struct med *medium,float *farray,int window)
{
  static float tr[6];
  static my_boolean initialized=FALSE;
  static int asize;
  float mmax,mmin;
  int i;
  if(!initialized){// boundaries for array
    //
    //The transformation matrix TR is used to calculate the world
    //coordinates of the center of the "cell" that represents each
    //array element. The world coordinates of the center of the cell
    //corresponding to array element A(I,J) are given by:
    //   X = TR[0] + TR[1]*I + TR[2]*J
    //   Y = TR[3] + TR[4]*I + TR[5]*J
    tr[0]= -medium->lmean/2.0;
    tr[1]=  medium->lmean;
    tr[2]=  0.0;
    tr[3]= -medium->wmean/2.0;
    tr[4]=  0.0;
    tr[5]=  medium->wmean;
    asize=medium->grp0_n*medium->grp0_m;
    initialized=TRUE;
  }
  // find extrema
  mmin= FLT_MAX;mmax=-FLT_MAX;
  for(i=0;i<asize ;i++){
    if(farray[i] > mmax)mmax=farray[i];
    if(farray[i] < mmin)mmin=farray[i];
  }
  if(medium->x_window[window]){
    cpgslct(medium->x_window[window]);
    cpgimag(farray,medium->grp0_m,medium->grp0_n,
	    1,medium->grp0_m,1,medium->grp0_n,mmin,mmax,tr);
  }
}
/*
  
  plot a stick plot given x, y vectors with n elements

*/
void psticks(float *x, float *y, int n)
{
  int i;
  float xmin,xmax,ymin,ymax;
  cpgqwin(&xmin,&xmax,&ymin,&ymax);
  for(i=0;i<n;i++){
    cpgmove(x[i],ymin);
    cpgdraw(x[i],y[i]);
  }
}

#endif
