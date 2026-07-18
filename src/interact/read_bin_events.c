#include "interact.h"
/*

  $Id: read_bin_events.c,v 1.11 2003/02/13 22:45:12 becker Exp $

  read the binary event file written by interact

*/
#include "cpgplot.h"



int main(int argc, char **argv)
{
  FILE *in;
  struct med *medium;
  struct flt *fault;
  int aflt,nriter,i,j,offset,nrxwindows=3,colorscale,asize,gn;
  float slip[3],time,mom,*slipa,bright,contrast,told,xmin,xmax,ymin,ymax,
    aspect,y[2],z[2],*group_moment_m,*group_moment_time,gmmin,gmmax,gmtick;
  my_boolean plot=TRUE,init=FALSE;
  char tmpchar;
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(argc!=1){
    fprintf(stderr,"%s\nreads single patch activations from binary file %s and patch geometry from %s\n",
	    argv[0],EVENT_FILE_BINARY,GEOMETRY_FILE);
    fprintf(stderr,"output is in the following format to stdout:\n");
    fprintf(stderr,"time nr_iteration nr_patch group(patch) loc_x loc_y slip_strike slip_dip slip_normal moment\n\n");
    exit(-1);
  }
  read_geometry(GEOMETRY_FILE,&medium,&fault,FALSE,FALSE,FALSE,FALSE);
  // read group events file
  read_moment_file(&group_moment_time,&group_moment_m,&gmmin,&gmmax,&gn,TRUE);
  gmtick=(gmmax-gmmin)*0.05+gmmin;
#ifdef USE_PGPLOT
  if(plot){// init all windows
     if(nrxwindows > MAX_NR_X_WINDOWS){
	fprintf(stderr,"too many X windows open attempts, max is %i\n",
		MAX_NR_X_WINDOWS);
	exit(-1);
     }
     for(i=0;i < nrxwindows;i++){
       // create
       medium->x_window[i]=cpgopen("/xwindow");
       if(medium->x_window[i] <= 0){
	 fprintf(stderr,"error opening window %i\n",i);
	 exit(-1);
       }
       // select
       cpgslct(medium->x_window[i]);
       if(i<2){
	 // size of fault
	 xmin = 0.0;xmax = ((float)medium->grp0_m)*medium->lmean;
	 ymin = 0.0;ymax = ((float)medium->grp0_n)*medium->wmean;
	 fprintf(stderr,"init_plot_window: initializing window %i with bounds x(%g, %g) y(%g, %g) for m: %i n: %i\n",
		 i,xmin,xmax,ymin,ymax,medium->grp0_m,medium->grp0_n);
	 aspect=xmax/ymax;
	 // scaling factors for arrays
	 // scaling for arrays such that pos variables get mapped into 0..n-1
	 y[STRIKE]=((float)medium->grp0_m-1.0)/(2.0*(float)aspect-medium->lmean);
	 z[STRIKE]=((float)medium->grp0_m-1.0)/2.0;
	 y[DIP]=(1.0-(float)medium->grp0_n)/(2.0*((float)medium->wmean-1.0));
	 z[DIP]=((float)medium->grp0_n-1.0)/2.0;
	 // environment
	 cpgenv(xmin,xmax,ymin,ymax, 1, 0);
	 // label
	 if(i==0){
	   asize = medium->grp0_n*medium->grp0_m;
	   cpglab("x/L", "y/W", "total moment release");
	   if((medium->momrel=(float *)calloc(asize,sizeof(float)))== NULL)
	     MEMERROR("init_plot_window");
	   medium->moment_array_init=TRUE;
	 }else{
	   // slip at certain time
	   cpglab("x/L","y/W","event slip");
	   if((slipa=(float *)malloc(asize*sizeof(float)))== NULL)
	     MEMERROR("init_plot_window");
	 }
	 bright = 0.5;
	 contrast  = 1.0;
	 colorscale = 3;
	 palett(&colorscale, &contrast, &bright);
       }else{// windows higher than 2
	 // y logscale plot for group moments
	 cpgenv(0.0,group_moment_time[gn-1]*1.05,gmmin,gmmax, 0, 0);
	 cpglab("time","ln(moment)","moment release");
	 psticks(group_moment_time,group_moment_m,gn);
       }
     }
  }
#endif
  told = -FLT_MAX;
  
  // read patch events
  in=myopen(EVENT_FILE_BINARY,"r");
  while(read_patch_event_file(&time,&nriter,&aflt,&mom,slip,in,medium)==7){
#ifdef USE_PGPLOT
    if(plot){
      if(time > told + .05){// coarse graining
	if(init){
	  // update plots if time has progressed, ie. new activation
	  for(i=0;i < nrxwindows;i++){
	    if(i == 0)
	      plot_moment_array(medium,medium->momrel,i);
	    else if(i==1)
	      plot_moment_array(medium,slipa,i);
	    else if(i==2){
	      cpgslct(medium->x_window[i]);
	      // erase old and add new tick
	      cpgsci(0);cpgmove(told,gmmin);cpgdraw(told,gmtick);
	      cpgsci(2);cpgmove(time,gmmin);cpgdraw(time,gmtick);
	      cpgsci(1);
	    }else
	      ;
	  }
	}else
	  init = TRUE;
	if(nrxwindows > 1){
	  // set slip arr to zero
	  for(i=0;i<asize;)
	    slipa[i++]=0.0;
	}
	told = time;
	// wait for input
	scanf("%c",&tmpchar);
      }
      i=(int)(fault[aflt].pos[DIP]*y[DIP] + z[DIP]+0.5);
      j=(int)(fault[aflt].pos[STRIKE]*y[STRIKE] + z[STRIKE]+0.5);
      offset  = i * medium->grp0_m + j;
      medium->momrel[offset] += mom;
      if(nrxwindows > 1)
	slipa[offset] += fabs(slip[STRIKE]);
    }
    /*
      printf("%10e %5i %5i %5i %12.5e %12.5e  %12.5e %12.5e %12.5e %12.5e\n",
      time,nriter,aflt,	    fault[aflt].group,
      fault[aflt].pos[0],fault[aflt].pos[1],
      slip[0],slip[1],slip[2],mom);
    */
#else
    printf("%10e %5i %5i %5i %12.5e %12.5e  %12.5e %12.5e %12.5e %12.5e\n",
	    time,nriter,aflt,	    fault[aflt].group,
	    fault[aflt].pos[0],fault[aflt].pos[1],
	    slip[0],slip[1],slip[2],mom);

#endif
  }
  fclose(in);
#ifdef USE_PGPLOT
  for(i=0;i < nrxwindows;i++){
    cpgslct(medium->x_window[i]);
    cpgclos();
  }
#endif
}
