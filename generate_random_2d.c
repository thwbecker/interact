#include "interact.h"
#include "properties.h"
/* 


generate a random stress field by superposition of 2-D dislocations



(might also want to take a look at randomflt.c, which is more sophisticated in terms
of checking for overlaps and such

$Id: generate_random_2d.c,v 1.1 2005/07/27 00:19:05 becker Exp $



*/

int main(int argc, char **argv)
{
  int nx,ny,i,j,k,n,random_mode,err;
  long int seed=-1;
  double alpha;
  COMP_PRECISION x[3],xmin,xrange,ymin,yrange,dx,dy;
  COMP_PRECISION fxmin,fxrange,fxmean,lmin,lrange,lmean,tmpdbl,s[3][3],sl[3][3],ul[3];
  struct flt *fault;
  /* 

  parameters

  
  */
  if(argc < 2){
    n = 10000;			/* number of faults */
  }else{
    sscanf(argv[1],"%i",&n);
  }
  fprintf(stderr,"%s: generating stress field for %i random faults\n",
	  argv[0],n);

  nx = ny = 101;		/* number of samples in x and y direction */
  
  xmin = -10;xrange=20;/* geographic bounds */
  yrange=xrange;
  ymin=xmin;
  /* derived spacings */
  dx = xrange/(nx-1);
  dy = yrange/(ny-1);

  /* 
     
     fault properties 

  */
  random_mode = UNIFORM_DISTRIBUTED;
  lmin = 0.001;			/* half length */
  lrange = xrange/2;		/* range */
  lmean = xrange/4;		/* for Gauss distributed */
  /* location */
  fxmin = xmin;			/* mid point location */
  fxrange = xrange;
  fxmean = 0;

  // intialize random number generator
  RAND_GEN(&seed);

  
  /* make room for faults */
  if((fault=malloc(sizeof(struct flt)*n))==NULL)
    MEMERROR("main: 2");
  
  if(n==1){
    /* single fault for debugging */
    fault[0].x[INT_X] = fault[0].x[INT_Y] = fault[0].x[INT_Z] = 0.0;
    fault[0].l = 2;
    fault[0].w = 0;
    fault[0].strike = 0;
    alpha=90.0-(double)fault[0].strike;
    my_sincos_degd(&fault[0].sin_alpha,&fault[0].cos_alpha,alpha);
    fault[0].u[STRIKE] = 2e-4 * fault[0].l;
    fault[0].u[NORMAL] = fault[0].u[DIP] =0.0;
  }else{
    /* 
       generate random faults 
    */
    for(i=0;i < n;i++){
      /* 
	 length 
      */
      assign_random(&fault[i].l,lmin,lmean,lrange,&seed,random_mode);
      fault[i].w = 0.0;		/* width is null */
      /* 
	 location 
      */
      assign_random(&fault[i].x[INT_X],fxmin,fxmean,fxrange,&seed,random_mode);
      assign_random(&fault[i].x[INT_Y],fxmin,fxmean,fxrange,&seed,random_mode);
      fault[i].x[INT_Z] = 0.0;
      /* 
       strike 
      */
      assign_random(&tmpdbl,0,180,360,&seed,random_mode);
      fault[i].strike = tmpdbl;	/* strike is a float variable in degree */
      /* derived quantities, need them for base vectors */
      alpha=90.0-(COMP_PRECISION)fault[i].strike;
      // get angle cosines
      my_sincos_degd(&fault[i].sin_alpha,&fault[i].cos_alpha,alpha);
      /* 
	 displacement
      */
      fault[i].u[STRIKE] = 2e-4 * fault[i].l;
      fault[i].u[NORMAL] = fault[i].u[DIP] = 0.0;
    } /* end fault loop */
  }

  /* 

  evaluate stress

  */
  x[INT_Z] = 0.0;
  for(i=0,x[INT_Y]=ymin;i<ny;i++,x[INT_Y]+=dy){	/* y loop */
    for(j=0,x[INT_X]=xmin;j<nx;j++,x[INT_X]+=dx){ /* x loop */
      /* init the important stress componenents */
      s[INT_X][INT_X] = s[INT_X][INT_Y] = s[INT_Y][INT_Y] = 0.0;
      for(k=0;k < n;k++){		
	/* 
	   sum over faults 
	*/
	/* get contribution from fault k */
	eval_2dsegment_plane_strain(x,(fault+k),fault[k].u,ul,sl,&err);
	if(!err){
	  /* non-infinite, add only three components */
	  s[INT_X][INT_X] += sl[INT_X][INT_X];
	  s[INT_X][INT_Y] += sl[INT_X][INT_Y];
	  s[INT_Y][INT_Y] += sl[INT_Y][INT_Y];
	}else{
	  fprintf(stderr,"%s: solution infinite at location x,y: %g, %g, fault %i\n",
		  argv[0],x[INT_X],x[INT_Y],k+1);
	}
      }	/* end fault loop */
      /* output */
      fprintf(stdout,"%11g %11g\t%11g %11g %11g\n",
	      x[INT_X],x[INT_Y],s[INT_X][INT_X],s[INT_X][INT_Y],s[INT_Y][INT_Y]);
    }
  }
  free(fault);

  return 0;
}
