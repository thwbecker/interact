#include "interact.h"
/*

  calculate spatial correlation coefficients

  input:  

  x, n, dim, y, bins

  x[n*dim] vector with the spatial coordinates of n 
  data points with scalar value y[n]

          
  output: r[bins], the correlation at a mid 
  distance range of xr[bins]
  and the number of pairs in this bin, nr[bins]

  all vectors have to be initialized outside this routine

  range_frac is the part of the total range that is used for bins

  cr is the characteristic correlation (say 0.1), 
  and xcr the value for which f(xcr) < cr 


  WARNING: the routine assumes that the distance between points (i.e.l
  the geometry) doesn't change between subsequent calls, only the
  scalar properties to be analyzed

*/
#define FAULT_BASED

#ifdef FAULT_BASED
//
// based on faults, ie. pos[2] will replace the coordinate vector x
// and fault[].s[component] the y vector
//
void calc_spatial_correlation(struct flt *fault,int n, int dim, int component,
			      int bins, COMP_PRECISION **r,COMP_PRECISION **xr,
			      int **nr,COMP_PRECISION range_frac,
			      COMP_PRECISION cr, COMP_PRECISION *xcr,
			      float **dist,float **xrl,float **xrr)
#else
// general version
void calc_spatial_correlation(float *x,int n, int dim, COMP_PRECISION *y,
			      int bins, COMP_PRECISION **r,COMP_PRECISION **xr,
			      int **nr,COMP_PRECISION range_frac,
			      COMP_PRECISION cr, COMP_PRECISION *xcr,
			      float **dist,float **xrl,float **xrr)
#endif
{
  COMP_PRECISION *xc=NULL,*yc=NULL;
  int i,j,k,l,m,o1,o2,o3,o4,dsize;
  static COMP_PRECISION maxdist,mindist,dx,tx1,tx2,tx3,range;
  static my_boolean init = FALSE;
  if(!init){
    /*

      init bins and distances only for first call

    */
    if(bins < 1 || n < 2 || dim < 1){
      fprintf(stderr,"calc_spatial_correlation: range error bins/n/dim: %i/%i/%i\n",
	    bins,n,dim);
      exit(-1);
    }
    // allocate memory
    dsize=(int)((COMP_PRECISION)SQUARE(n));
    *dist = (float *)malloc(sizeof(float)*dsize);
    *r=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*bins);
    *nr=(int *)malloc(sizeof(int)*bins);
    *xrl=(float *)malloc(sizeof(float)*bins);
    *xrr=(float *)malloc(sizeof(float)*bins);
    *xr=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*bins);
    if(! *dist|| ! *r || ! *nr || ! *xrl || ! *xrr || ! *xr)
      MEMERROR("calc_spatial_correlation");
    //
    // find min and max distance between points
    //
    mindist=  FLT_MAX;maxdist= -FLT_MIN;
#ifdef FAULT_BASED
    for(m=i=0;i < n;i++)
      for(j=i+1;j<n;j++){
	*(*dist+m) = 
	  distance_float(fault[i].pos,fault[j].pos,(int)2);
#else
    for(m=o1=i=0;i < n;i++,o1 += dim)
      for(j=i+1,o2=j*dim;j<n;j++,o2 += dim){
	*(*dist+m) = distance_float((x+o1),(x+o2),dim);
#endif
	if(*(*dist+m) < mindist)mindist = (COMP_PRECISION)*(*dist+m);
	if(*(*dist+m) > maxdist)maxdist = (COMP_PRECISION)*(*dist+m);
	m++;
      }
    if(range_frac < 0 || range_frac > 1){
      fprintf(stderr,"calc_spatial_correlation: range_frac (%g) error\n",
	      range_frac);
      exit(-1);
    }
    // distance range to explore for correlation
    range = (maxdist - mindist)*range_frac;
    maxdist = mindist + range + EPS_COMP_PREC;
    mindist -= EPS_COMP_PREC;
    range = maxdist - mindist;
    //
    // create bins
    dx = range / (COMP_PRECISION)bins;
    tx1 = mindist;// left bound of interval
    tx2 = tx1 + dx/2.0;// center
    tx3 = tx1 + dx;// right bound of interval
    for(i=0;i<bins;i++){
      *(*xrl+i) = (float)tx1;
      *(*xr+i)  = tx2;
      *(*xrr+i) = (float)tx3;
      tx1 += dx;tx2 += dx; tx3 += dx;
    }
    init = TRUE;
  }
    //
  // loop through distance bins and calculate correlation
  for(i=0;i<bins;i++){
    //
    // find pairs with distances in range
    //
    // reset this bin
    *(*nr+i) = o1 = 0;// counter of points in range
    o4 = dim;
    // correlation vectors
    xc=(COMP_PRECISION *)realloc(xc,sizeof(COMP_PRECISION)*o4);
    yc=(COMP_PRECISION *)realloc(yc,sizeof(COMP_PRECISION)*o4);
    //
    // loop through all possible combinations of points
    //
    for(o2=m=j=0;j < n;j++,o2+=dim)
      for(k=j+1,o3=k*dim;k < n;k++,o3+=dim){
	if((*(*dist+m) >= *(*xrl+i))&&// check if distance is within range
	   (*(*dist+m) < *(*xrr+i))){
	  //
	  // add function value to correlation vector
	  //
	  for(l=0;l<dim;l++){
#ifdef FAULT_BASED
	    if(component < 3){// fault stresses
	      xc[o1+l] = fault[j].s[component];
	      yc[o1+l] = fault[k].s[component];
	    }else if(component == 3){// stat fric
	      xc[o1+l] = (COMP_PRECISION)fault[j].mu_s;
	      yc[o1+l] = (COMP_PRECISION)fault[k].mu_s;
	    }else if(component == 4){// dyn fric
	      xc[o1+l] = (COMP_PRECISION)fault[j].mu_d;
	      yc[o1+l] = (COMP_PRECISION)fault[k].mu_d;
	    }else{
	      fprintf(stderr,"calc_spatial_corr: component error, %i undefined\n",
		      component);
	      exit(-1);
	    }
#else
	    xc[o1+l]=y[j];yc[o1+l]=y[k];
#endif
	  }
	  // increment counters
	  *(*nr+i) += 1;
	  o1 += dim;o4 += dim;
	  xc=(COMP_PRECISION *)realloc(xc,sizeof(COMP_PRECISION)*o4);
	  yc=(COMP_PRECISION *)realloc(yc,sizeof(COMP_PRECISION)*o4);
	  if(!xc || !yc)MEMERROR("calc_spatial_correlation: xc/yc");
	}// end add-to-bin branch if in range
	m++;
      }
    // assign the correlation to this distance range
    *(*r+i) = correlation_coefficient(xc,yc,*(*nr+i));
  }
  //
  // find the cr cutoff, ie. the xcr value where the correlation is below
  // cr
  //
  for(i=0;i<bins;i++)
    if(*(*r+i) < cr){
      if(i==0)
	*xcr = *(*xr+i);
      else{// linearly interpolate
	*xcr = *(*xr+i-1) + ((*(*r+i-1)-cr)/(*(*r+i-1)- *(*r+i)))*dx;
      }
      break;
    }
  free(xc);free(yc);
}
  /*

    calculate the correlation coefficient for two x[n] y[n] vectors
    if n<2, returns NaN
    
   */
COMP_PRECISION correlation_coefficient(COMP_PRECISION *x,
				       COMP_PRECISION *y,
				       int n)
{
  COMP_PRECISION mx,my,s1,s2,s3,dx,dy;
  int i;
  if(n >= 2){
    mx=0.0;my=0.0;
    for(i=0;i<n;i++){
      mx += x[i];my += y[i];
    }
    mx /= (COMP_PRECISION)n;my /= (COMP_PRECISION)n;
    s1=0.0;s2=0.0;s3=0.0;
    for(i=0;i<n;i++){
      dx=x[i]-mx;dy=y[i]-my;
      s1 += (dx*dy);s2 += (dx*dx);s3 += (dy*dy);
    }
    return s1/(sqrt(s2)*sqrt(s3));
  }else{
    //
    // not enough entries
    //
    return NaN;
  }
}
