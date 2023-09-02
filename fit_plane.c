#include "interact.h"



/*


subroutine to fit a rectangular plane to n (input) points in 3-D,
stored in xq[n*3] vector (input)

uses LAPACK SVD routine and finds the normal vector to the plane that
contains the rectangle from the eigenvector of the smallest eigenvalue

t_strike, t_dip, and normal (output) are [3] vectors
  
l and w (output) are the half length and width of the best fit
rectangle respectively

mx[3] (output) are the coordinates of the center of the plane
  
sin_alpha, cos_alpha, fdip, and fstrike (output) refer to the fault
angles (sin(alpha) and cos(alpha) as opposed to alpha directly) as
defined in read_geometry

if n = 4, will expect a quad input and adjust area such that original
and new area are the same, based on a triangular approximation of the,
possibly, irregular initial rectangle
  

$Id: fit_plane.c,v 1.13 2003/06/27 22:14:28 becker Exp $

*/
void fit_plane(int n,
	       COMP_PRECISION *xq,COMP_PRECISION *t_strike, 
	       COMP_PRECISION *t_dip,
	       COMP_PRECISION *normal, 
	       COMP_PRECISION *l, COMP_PRECISION *w,
	       COMP_PRECISION *sin_alpha, 
	       COMP_PRECISION *cos_alpha,
	       float *fstrike, float *fdip, 
	       COMP_PRECISION *mx,my_boolean verbose,
	       my_boolean adjust_area)
{
  char c1='O',c2='N';		/* for SVD routine */
  COMP_PRECISION *p,*work,s[3],*dummy=NULL,scale;
  COMP_PRECISION orig_area=0,cx[3],aspect,strike,dip,alpha,
    sin_dip,cos_dip,*xqc,dist_to_plane,*ll,*ww,llm,wwm;
  int mn,m=3,lwork,i,j,info;
  //
  // array sizes
  mn = 3*n;
  lwork=((MAX(3*MIN(3,n)+MAX(3,n),5*MIN(3,n)-4))*4);// factor is for improved speed
  //
  // allocate memory
  ll=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n);
  ww=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n);
  p=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*mn);
  xqc=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*mn);
  work=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*lwork);
  if(!ll || !ww || !p || !xqc || !work)MEMERROR("fit_plane");
  //
  // orginal area of deformed quad
  //
  if(adjust_area){
    if(n == 4){
      orig_area = quad_area(xq);
      if(verbose)
	fprintf(stderr,"fit_plane: adjusting area for quad input\n");
    }else{
      fprintf(stderr,"fit_plane: adjusting area only works for quads (n=4) input, n: %i\n",n);
      exit(-1);
    }
  }
  //
  // calculate mean location (mx) from coordinates and 
  // remove from coordinates. the new coordinates are xqc
  //
  calc_mean_quad_coord(xq,mx);
  for(scale=0.0,i=0;i<mn;i += 3) /* shift to new origin */
    c_eq_a_minus_b_3d((xqc+i),(xq+i),mx);
  //
  // P matrix has points in columns, copy of xqc 
  // (since P gets overwritten with the vectors of the SVD)
  //
  a_equals_b_vector(p,xqc,mn);
  //
  // compute SVD to find normal vector associated with 
  // smallest singular value
  // first vector is g, second h, last normal. 
  // g, h, and normal are normalized
  //
  // ratio of s[1] to s[0] is like mean w/l
  //
  //fprintf(stderr,"ok1\n");
#ifdef USE_DOUBLE_PRECISION
  dgesvd_(&c1,&c2,&m,&n,p,&m,s,dummy,&m,dummy,&n,
	  work,&lwork,&info);
#else
  sgesvd_(&c1,&c2,&m,&n,p,&m,s,dummy,&m,dummy,&n,
	  work,&lwork,&info);
#endif
  //fprintf(stderr,"ok2\n");

  /*
    
    on output:
    p[0,1,2] is the strike (g) vector
    p[3,4,5] is the dip (h) vector
    p[6,7,8] is the normal vector

  */
  if(info){
    fprintf(stderr,"fit_plane: svd: error in dgesvd, code: %i\n",info);
    exit(-1);
  }
  //
  // aspect ratio = W/L ~ s[1]/s[0]
  if((s[0] == 0.0) || (s[1] == 0.0)){
    fprintf(stderr,"fit_plane: svd: error: sigular values 0,1,2: %g %g %g\n",
	    s[0],s[1],s[2]);
    exit(-1);
  }
  for(i=0;i<mn;i+=3){// project all points onto the best-fit plane
    dist_to_plane = dotp_3d((p+6),(xqc+i));// (p+6) is the normal vector
    for(j=0;j<3;j++)// reassign xqc to the projected locations
      xqc[i+j] -= dist_to_plane * p[6+j];
  }
  // determine centroid within newly projected locations on plane
  calc_centroid_quad(xqc,cx);
  // add cx correction to mean fault coordinate
  add_b_to_a_vector_3d(mx,cx);
  /*
    
    determine fault angles from normal vector as given by SVD
    
  */
#ifdef DEBUG
  if(fabs(norm_3d((p+6))-1.0) > EPS_COMP_PREC){
    fprintf(stderr,"fit_plane: error: normal vector not normalized\n");
    exit(-1);
  }
#endif
  dip   = RAD2DEGF(acos(p[6+INT_Z]));
  alpha = RAD2DEGF(atan2(-p[6+INT_X],p[6+INT_Y]));
  my_sincos_deg(sin_alpha,cos_alpha,(COMP_PRECISION)alpha);
  my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)dip);
  if(dip != 0.0){// check which sign of dip we need
    if(*sin_alpha != 0.0){ // can we use sin_alha to check sign?
      if((p[6+INT_X]/(*sin_alpha))/sin_dip < 0){// should be positive
	dip = -dip;
	sin_dip= -sin_dip;
      }
    }else{// have to use n_y and cos_alpga, should be negative
      if((p[6+INT_Y]/(*cos_alpha))/sin_dip > 0){
	dip = -dip;
	sin_dip= -sin_dip;
      }
    }
  }
  calc_base_vecs(t_strike,normal,t_dip,
		 *sin_alpha,*cos_alpha,sin_dip,cos_dip);
  if(verbose)
    fprintf(stderr,"fit_plane: S1: %10.3e S2: %10.3e S3: %10.3e g: %6.3f %6.3f %6.3f h: %6.3f %6.3f %6.3f n: %6.3f %6.3f %6.3f\n",
	    s[0],s[1],s[2],p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8]);
  //
  // get patch type angles
  //
  strike=90.0-alpha;
  check_angles(&dip,&strike);
  *fdip=(float)dip;
  *fstrike=(float)strike;
  //
  // get best fit projections
  //
  for(llm=wwm=0.0,
	i=j=0;i<mn;i+=3,j++){
    // correct xqc by centroid
    sub_b_from_a_vector_3d((xqc+i),cx);
    // get projection of distance to mid point
    ll[j]=fabs(dotp_3d(t_strike,(xqc+i)));
    ww[j]=fabs(dotp_3d(t_dip,(xqc+i)));
    llm += ll[j];
    wwm += ww[j];
  }
  llm /= (COMP_PRECISION)n;
  wwm /= (COMP_PRECISION)n;
  //
  // sort the points (inverted)
  qsort(ll,n,sizeof(COMP_PRECISION),inv_compare_flt);
  qsort(ww,n,sizeof(COMP_PRECISION),inv_compare_flt);
  //
  // get aspect ratio corrected for area
  //

  if(adjust_area){// this only works for four nodes and quads!

    aspect=(ww[0]+ww[1])/(ll[0]+ll[1]);
    *l = 0.5 * sqrt(orig_area/aspect);
    *w = aspect * *l;
  }else{
    //
    //*l = llm; *w = wwm;
    *l = s[0]/2.0; *w = s[1]/2.0;
  }
  free(ll);free(ww);free(p);free(xqc);free(work);
}




/*
  
  convert quad forming points in xq[12] into fault patch
  
  calculates angles and base vectors using the fit_plane subroutine
  from above, additionally sets the correct fault code 

  if adjust_area is set, will modify length and width such that the
  original (approximated by two triangles and final area are the same)
  
*/
void points2patch(struct flt *fault,COMP_PRECISION *xq,
		  my_boolean adjust_area)
{
#ifdef ALLOW_NON_3DQUAD_GEOM
  fault->type = RECTANGULAR_PATCH;
#endif
  /*
    fprintf(stderr,"before: x: %g %g %g; %g %g %g; %g %g %g; %g %g %g\n",
    xq[INT_X],xq[INT_Y],xq[INT_Z], xq[INT_X+3],xq[INT_Y+3],xq[INT_Z+3], xq[INT_X+6],xq[INT_Y+6],xq[INT_Z+6], 
    xq[INT_X+9],xq[INT_Y+9],xq[INT_Z+9]);
  */
  // fit to plane and calculate g and h vectors
  fit_plane(4,xq,fault->t_strike,fault->t_dip,fault->normal,
	    &fault->l,&fault->w,&fault->sin_alpha,&fault->cos_alpha,
	    &fault->strike,&fault->dip,fault->x,FALSE,adjust_area);
}
