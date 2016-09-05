/*
  interact: model fault interactions using dislocations in a 
  halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu
  
  $Id: geometry.c,v 2.44 2003/12/20 21:53:59 becker Exp $
*/
#include "interact.h"
#include <math.h>

/*

  calculate the lower hemisphere projection of a fault geometry
  given the dip (in degrees from horizontal downward) and the strike
  (in degree CW from North) 

  output is projection in x[3] cartesian coordinates
*/
void calc_lhemi_proj(COMP_PRECISION dip, COMP_PRECISION strike, 
		     COMP_PRECISION *x)
{
  COMP_PRECISION plunge,sinpl,cospl;
  plunge= (90.0-dip) * DEG2RAD;
  x[Z]  = sin((PIHALF-plunge)/2.0);
  my_sincos_deg(&sinpl,&cospl,90.0+strike);
  x[X] = x[Z] * sinpl;
  x[Y] = x[Z] * cospl;
}

/* 
   obtain a traction vector given matrix sm and
   normal vector on the plane 
*/

void resolve_force(COMP_PRECISION *norm,COMP_PRECISION sm[3][3],
		   COMP_PRECISION *trac)
{
  int i,j;
  for(i=0;i<3;i++)
    for(trac[i]=0.0,j=0;j<3;j++)
      trac[i] += norm[j] * sm[i][j];
}
/* 
   given a fault patch with specified 
   angles strike and dip (resp. their cosines and sines),
   calculate the normal, dip, and strike vectors of that
   patch
   
   INPUT:
   
   the sines and cosines of the following angles:
   
   alpha: counterclockwise angle from east
   dip:   positive upward from horizontal, ie. a vertical
   fault has dip 90 degrees
   
   OUTPUT:
   strike[3], normal[3], and dip[3]
   
*/

void calc_base_vecs(COMP_PRECISION *strike,COMP_PRECISION *normal,
		    COMP_PRECISION *dip,
		    COMP_PRECISION sin_alpha,
		    COMP_PRECISION cos_alpha,
		    COMP_PRECISION sin_dip, 
		    COMP_PRECISION cos_dip)
{
  /* tangential vector in strike direction, by definition 
     this vector is in the x-y plane */
  strike[X]=  cos_alpha;
  strike[Y]=  sin_alpha;
  strike[Z]=  0.0;
  /* normal vector on fault plane */
  normal[X]=   sin_alpha * sin_dip;
  normal[Y]=  -cos_alpha * sin_dip;
  normal[Z]=   cos_dip;
  /* tangential vector in dip direction */
  dip[X]=  -sin_alpha * cos_dip;
  dip[Y]=   cos_alpha * cos_dip;
  dip[Z]=   sin_dip;
}

/*

  wrapper for the routine below with f factors being calculated locally

*/
void get_maxsdir_stress_drops2(COMP_PRECISION *stress,
			       COMP_PRECISION stress_drop,
			       COMP_PRECISION *tstress_drop)
{
  COMP_PRECISION f,f1,f2;
  f = (fabs(stress[DIP]) < EPS_COMP_PREC)?(-1.0):
    (fabs(stress[STRIKE]/stress[DIP]));
  get_maxsdir_stress_drops(stress,stress_drop,tstress_drop,f,&f1,&f2);
}
/*

  determines the scalar stress drop values for strike and dip
  direction given the stress state as specified in
  stress[STIKE,DIP,NORMAL], the scalar stress drop as given by
  stress_drop, and the f factor as defined below

  tstress_drop[STRIKE,DIP] is the output, it's normal component
  remains unspecified


  f= \sigma_s / \sigma_d for sigma_d != 0 and f=-1 else
  
  IS INPUT based on the initial stress configuration like so:
  
  f = (fabs(bstress[DIP]) < EPS_COMP_PREC)?(-1.0):
           (fabs(bstress[STRIKE]/bstress[DIP]));
	     
  output is 

  f1=f/sqrt(1+f^2)

  f2=1/sqrt(1+f^2)

  then

  \Delta \sigma_s = f1 * \Delta s
  \Delta \sigma_d = f2 * \Delta s

*/
void get_maxsdir_stress_drops(COMP_PRECISION *stress,
			      COMP_PRECISION stress_drop,
			      COMP_PRECISION *tstress_drop, 
			      COMP_PRECISION f,
			      COMP_PRECISION *f1, COMP_PRECISION *f2)
{
  COMP_PRECISION local_f[2];// leave coding like that to allow passing
#ifdef DEBUG
  if(stress_drop < 0){
    fprintf(stderr,"get_maxsdir_stress_drops: error, stress drop negative on input (%g)\n",
	    stress_drop);
    exit(-1);
  }
#endif
  // of dummy arguments
  if(f == -1.0){// zero dip stress 
    tstress_drop[STRIKE] = (stress[STRIKE]>0)?(-stress_drop):(stress_drop);
    tstress_drop[DIP] = 0.0;
    local_f[0] = 1.0; local_f[1] = 0.0;
  }else{
    //
    // calculate how the stress drop should be distributed given the
    // vector direction of the shear stress ie., drop in fractions as
    // given by the shear stress vector in the strike/dip plane
    //
    local_f[1] = 1.0/hypot((double)1.0,(double)f);
    local_f[0] = f * local_f[1];
    tstress_drop[STRIKE] =  stress_drop * local_f[0];
    tstress_drop[DIP]    =  stress_drop * local_f[1];
    // take care of sign
    if(stress[STRIKE] > 0)
      tstress_drop[STRIKE] = -tstress_drop[STRIKE];
    if(stress[DIP] > 0)
      tstress_drop[DIP]    = -tstress_drop[DIP];
  }
  *f1 = local_f[0];
  *f2 = local_f[1];
}

/* 
   calculate the positions of the corners of a fault patch if
   fault is a 3-D rectangle

   if fault is a segment, only calculate the first two corners, 
   which will be the left and right endpoint of the segment 

*/
void calculate_corners(COMP_PRECISION corner[4][3],
		       struct flt *fault,
		       COMP_PRECISION *l, COMP_PRECISION *w)
{
 
#ifdef ALLOW_NON_3DQUAD_GEOM
  if(patch_is_2d(fault->type)){
    calculate_seg_corners(corner,fault,1.0);
    *l = fault->l;*w = fault->w;
  }else if(fault->type == RECTANGULAR_PATCH){
    calculate_quad_corners(corner,fault,1.0);
    *l = fault->l;*w = fault->w;
  }else if(fault->type == POINT_SOURCE)
    calculate_point_source_corners(corner,fault,1.0,l,w);
  else{
    fprintf(stderr,"calculate_corners: error, type %i is not implemented\n",
	    fault->type);
    exit(-1);
  }
#else
  calculate_quad_corners(corner,fault,1.0);
  *l = fault->l;*w = fault->w;
#endif
}
/*
  blows up the faults length and width by a factor 
  leeway and calculates the corners
  if leeway = 1 same as above
  
*/
void calculate_bloated_corners(COMP_PRECISION corner[4][3],
			       struct flt *fault,
			       COMP_PRECISION leeway)
{
#ifdef ALLOW_NON_3DQUAD_GEOM
  COMP_PRECISION lloc,wloc;
  if(patch_is_2d(fault->type))
    calculate_seg_corners(corner,fault,leeway);
  else if(fault->type == RECTANGULAR_PATCH)
    calculate_quad_corners(corner,fault,leeway);
  else if(fault->type == POINT_SOURCE)
    calculate_point_source_corners(corner,fault,leeway,&lloc,
				   &wloc);
  else{
    fprintf(stderr,"calculate_bloated_corners: error, fault type %i not implemented\n",
	    fault->type);
    exit(-1);
  }
#else
  calculate_quad_corners(corner,fault,leeway);
#endif
}

/*
  
  calculate the corners of a rectangular patch

 */
void calculate_quad_corners(COMP_PRECISION corner[4][3],struct flt *fault,
			    COMP_PRECISION leeway)
{
  int i;
  COMP_PRECISION sx,dx;
  for(i=0;i<3;i++){
    sx = fault->t_strike[i] * (COMP_PRECISION)fault->l * leeway;
    dx = fault->t_dip[i]    * (COMP_PRECISION)fault->w * leeway;
    // lower left
    corner[0][i]=fault->x[i]-sx-dx;
    // lower right
    corner[1][i]=fault->x[i]+sx-dx;
    // upper right
    corner[2][i]=fault->x[i]+sx+dx;
    // upper left
    corner[3][i]=fault->x[i]-sx+dx;
  }
}
/*
  
  calculate the "corners" of a point source patch

*/
void calculate_point_source_corners(COMP_PRECISION corner[4][3],
				    struct flt *fault,
				    COMP_PRECISION leeway,
				    COMP_PRECISION *l,
				    COMP_PRECISION *w)
{
  int i;
  COMP_PRECISION sx,dx;
  if(fault->l > 0){
    fprintf(stderr,"calculate_point_source_corners: fault->l should be < 0 and \"aspect ratio\"\n");
    exit(-1);
  }
  // W = sqrt(area/(4 aspect ratio)) = sqrt(W'/(4 -L'))
  *w = sqrt(((COMP_PRECISION)fault->w)/
	    (-4.0*(COMP_PRECISION)fault->l));
  // L = aspect ratio * W
  *l = (COMP_PRECISION)-fault->l * (*w);
  for(i=0;i<3;i++){
    sx = fault->t_strike[i] * (*l) * leeway;
    dx = fault->t_dip[i]    * (*w) * leeway;
    // lower left
    corner[0][i]=fault->x[i]-sx-dx;
    // lower right
    corner[1][i]=fault->x[i]+sx-dx;
    // upper right
    corner[2][i]=fault->x[i]+sx+dx;
    // upper left
    corner[3][i]=fault->x[i]-sx+dx;
  }
}
/*
  
  calculate the two corners of a 2-D segment

  corners 3 and four remain undefined
*/
void calculate_seg_corners(COMP_PRECISION corner[4][3],struct flt *fault,
			   COMP_PRECISION leeway)
{
  int i;
  COMP_PRECISION sx;
  for(i=0;i<2;i++){
    sx = fault->t_strike[i] * (COMP_PRECISION)fault->l * leeway;
    // lower left
    corner[0][i]=fault->x[i]-sx;
    // lower right
    corner[1][i]=fault->x[i]+sx;
  }
  corner[0][Z] = corner[1][Z] = 0.0;
}

/*
  
  calculate the area of an (irregular) quad using two triangles
  we take the two possible divisions and average ...
*/
COMP_PRECISION quad_area(COMP_PRECISION *xq)
{
  COMP_PRECISION area,area2,gvec[3],hvec[3];
  // first possibility
  c_eq_a_minus_b_3d(gvec,(xq+3),xq);
  c_eq_a_minus_b_3d(hvec,(xq+9),xq);
  area = triangle_area_gh(gvec,hvec);
  c_eq_a_minus_b_3d(gvec,(xq+9),(xq+6));
  c_eq_a_minus_b_3d(hvec,(xq+3),(xq+6));
  area += triangle_area_gh(gvec,hvec);
  // second possibility
  c_eq_a_minus_b_3d(gvec,xq,(xq+3));
  c_eq_a_minus_b_3d(hvec,(xq+6),(xq+3));
  area2  = triangle_area_gh(gvec,hvec);
  c_eq_a_minus_b_3d(gvec,xq,(xq+9));
  c_eq_a_minus_b_3d(hvec,(xq+6),(xq+9));
  area2 += triangle_area_gh(gvec,hvec);
  
  return ((area+area2)/2.0);
}

/*
  calculates the area of a triangle
  given xt which has xt[i*3 + j] i=0,1,2 nodes in 
  j=0,1,2 dimensions
*/
COMP_PRECISION triangle_area(COMP_PRECISION *xt)
{
  COMP_PRECISION gvec[3],hvec[3];
  get_gh_tri_vec(xt,gvec,hvec);
  return(triangle_area_gh(gvec,hvec));
}

/*
  calculates the area of a triangle
  given the loval g (x_2-x-1) and h (x_3-x_1) vectors
*/
COMP_PRECISION triangle_area_gh(COMP_PRECISION *gvec,COMP_PRECISION *hvec)
{
  COMP_PRECISION cross[3];
  cross_product(gvec,hvec,cross);
  return norm_3d(cross)/2.0;
}
/*
  
  determine g and h vectors, see above
  points are in FE ordering, counterclockwise
  
  3
  \
  ^  \
  |   \
  |    \
  |     \
  h      \
  
  1 g---> 2
  
  
*/

void get_gh_tri_vec(COMP_PRECISION *xt,// 3 points in FE ordering
		    COMP_PRECISION *gvec,// output
		    COMP_PRECISION *hvec)// output
{
  gvec[X]=xt[3+X]-xt[ +X];
  gvec[Y]=xt[3+Y]-xt[ +Y];
  gvec[Z]=xt[3+Z]-xt[ +Z];
  hvec[X]=xt[6+X]-xt[ +X];
  hvec[Y]=xt[6+Y]-xt[ +Y];
  hvec[Z]=xt[6+Z]-xt[ +Z];
}

/*
  
  determine g and h for rectangular element
  
  
  3 -----------b----------- 2
  |            ^            |
  |            |            |
  |            |h           |
  |            x----------->a
  |                  g      |
  |                         |
  0 ----------------------- 1
  
  such that |g| is L (half lenght) and |h| is W (half-width)
  x is the centroid
  
*/
void get_gh_quad_vec(COMP_PRECISION *xq,// four points in FE ordering
		     COMP_PRECISION *xc,// centroid point
		     COMP_PRECISION *gvec,// output vectors
		     COMP_PRECISION *hvec)
{
  COMP_PRECISION a[3],b[3];
  c_eq_a_plus_b_3d(a,(xq+6),(xq+3));scale_vector_3d(a,0.5);
  c_eq_a_plus_b_3d(b,(xq+9),(xq+6));scale_vector_3d(b,0.5);
  c_eq_a_minus_b_3d(gvec,a,xc);
  c_eq_a_minus_b_3d(hvec,b,xc);
  if(fabs(gvec[Z])>fabs(hvec[Z])) // dip should indicate depth
    swap_ab_vector_3d(gvec,hvec);
  else if(gvec[Z] == hvec[Z]){// in plane
    if(fabs(vec_to_strike(hvec)) > fabs(vec_to_strike(gvec)))
      swap_ab_vector_3d(gvec,hvec);
  }
}

/*
  
  check if four points in FE ordering are on a plane
  
*/
my_boolean check_planar(COMP_PRECISION *x)
{
  COMP_PRECISION vec[4][3],normal[3];
  static COMP_PRECISION eps=EPS_COMP_PREC * 1000.0;
  c_eq_a_minus_b_3d(&vec[0][X],(x+3+X),(x+X));// get vectors along the edges
  c_eq_a_minus_b_3d(&vec[1][X],(x+6+X),(x+3+X));
  c_eq_a_minus_b_3d(&vec[2][X],(x+9+X),(x+6+X));
  c_eq_a_minus_b_3d(&vec[3][X],(x+X),(x+9+X));
  cross_product(&vec[0][X],&vec[1][X],normal);
  if((project_vector(&vec[2][X],normal) > eps)||
     (project_vector(&vec[2][X],normal) > eps)){
    /* 
       fprintf(stderr,"check_planar: rectangular element points possibly not on plane\n");
       fprintf(stderr,"check_planar: points: %g %g %g; %g %g %g; %g %g %g; %g %g %g\n",
       x[X],x[Y],x[Z], x[X+3],x[Y+3],x[Z+3], x[X+6],x[Y+6],x[Z+6], 
       x[X+9],x[Y+9],x[Z+9]);
    */
    return(FALSE);
  }
  // check right angles
  if((fabs(project_vector(&vec[0][X],&vec[1][X])) > eps)||
     (fabs(project_vector(&vec[1][X],&vec[2][X])) > eps)||
     (fabs(project_vector(&vec[2][X],&vec[3][X])) > eps)){
    /*
      fprintf(stderr,"check_planar: rectangular element points possibly not all 90 degree angles\n");
      fprintf(stderr,"check_planar: points: %g %g %g; %g %g %g; %g %g %g; %g %g %g\n",
      x[X],x[Y],x[Z], x[X+3],x[Y+3],x[Z+3], x[X+6],x[Y+6],x[Z+6], 
      x[X+9],x[Y+9],x[Z+9]);
      fprintf(stderr,"check_planar: dotps: 01: %g 12: %g 23: %g\n",
      project_vector(&vec[0][X],&vec[1][X]),
      project_vector(&vec[1][X],&vec[2][X]),
      project_vector(&vec[2][X],&vec[3][X]));
    */
    return(FALSE);
  }
  return TRUE;
}

/*
  
  determine alpha and dip vectors given non-normalized
  g and h vectors, also determines area
  
*/
void get_alpha_dip_tri_gh(COMP_PRECISION *xt,COMP_PRECISION *sin_alpha,
			  COMP_PRECISION *cos_alpha,COMP_PRECISION *dip,
			  COMP_PRECISION *area)
{
  COMP_PRECISION normal[3],nl,alpha,gvec[3],hvec[3];
  int i;
  get_gh_tri_vec(xt,gvec,hvec);
  //fprintf(stderr,"g: %g %g %g\n",gvec[X],gvec[Y],gvec[Z]);
  //fprintf(stderr,"h: %g %g %g\n",hvec[X],hvec[Y],hvec[Z]);
  cross_product(gvec,hvec,normal);
  //fprintf(stderr,"n: %g %g %g\n",normal[X],normal[Y],normal[Z]);
  nl=norm_3d(normal);
  if(nl == 0.0){
    fprintf(stderr,"get_alpha_dip_tri_gh: triangle degenerate:\n");
    for(i=0;i<3;i++)
      fprintf(stderr,"%i: x: %g y: %g z: %g\n",i+1,
	      xt[i*3+X],xt[i*3+Y],xt[i*3+Z]);
    exit(-1);
  }
  *area= nl/2.0;
  //normal[X]/=nl; not needed
  normal[Y]/=nl;
  normal[Z]/=nl;
  *dip= acos(normal[Z]);
  alpha=acos(normal[Y]);
  // might have to check for bounds
  check_angles(dip,&alpha);
  my_sincos(sin_alpha,cos_alpha,(COMP_PRECISION)alpha);
}



/*
  
  determine several geometrical quantities of 
  a fault group that consists of several patches
  (assumes planar faults)
  
*/

void calc_group_geometry(struct med *medium,struct flt *fault,
			 struct geog *grp)
{
  int i,j,k,clim;
  COMP_PRECISION dx[3],corner[4][3],fac,pos[2],
    l,w;
  /* 
     determine center of mass and average strike 
     and dip vectors for each patch group 
  */
  for(i=0;i<medium->nrflt;i++){
    grp[fault[i].group].nrflt++;
    add_b_to_a_vector_3d(grp[fault[i].group].center,fault[i].x);
    add_b_to_a_vector_3d(grp[fault[i].group].strike_vec,fault[i].t_strike);
    add_b_to_a_vector_3d(grp[fault[i].group].dip_vec,fault[i].t_dip);
  }
  for(i=0;i<medium->nrgrp;i++){
    if(grp[i].nrflt){// get average
      fac=1.0/((COMP_PRECISION)grp[i].nrflt);
      scale_vector_3d(grp[i].center,fac);
      scale_vector_3d(grp[i].strike_vec,fac);
      normalize_3d(grp[i].strike_vec);
      scale_vector_3d(grp[i].dip_vec,fac);
      normalize_3d(grp[i].dip_vec);
    }
    for(j=0;j<2;j++){
      grp[i].pmin[j]= FLT_MAX;
      grp[i].pmax[j]=-FLT_MAX;
    }
  }
  /* 

     assign local position to patches in terms of distance in strike
     and dip direction

  */
  for(i=0;i<medium->nrflt;i++){
    // determine distance of center of patch from 
    // center of group of patches
    for(j=0;j<3;j++)
      dx[j] = fault[i].x[j] - grp[fault[i].group].center[j];
    fault[i].pos[STRIKE]=(float)
      project_vector(dx,grp[fault[i].group].strike_vec);
    fault[i].pos[DIP]=(float)
      project_vector(dx,grp[fault[i].group].dip_vec);

#ifdef ATZ_NATZ_VOODOO_GNATZ    
    fprintf(stderr,"mx: %g %g %g ms: %g %g %g md: %g %g %g P: %g %g\n",
	    grp[fault[i].group].center[X],grp[fault[i].group].center[Y],
	    grp[fault[i].group].center[Z],
	    grp[fault[i].group].strike_vec[X],
	    grp[fault[i].group].strike_vec[Y],
	    grp[fault[i].group].strike_vec[Z],
	    grp[fault[i].group].dip_vec[X],
	    grp[fault[i].group].dip_vec[Y],
	    grp[fault[i].group].dip_vec[Z],
	    fault[i].pos[STRIKE],fault[i].pos[DIP]);
#endif
    // determine actual extent of corners of patch
    calculate_corners(corner,(fault+i),&l,&w);
#ifdef ALLOW_NON_3DQUAD_GEOM
    if(patch_is_2d(fault->type))
      clim = 2;
    else if((fault->type == RECTANGULAR_PATCH)||
	    (fault->type == POINT_SOURCE))
      clim = 4;
    else{
      fprintf(stderr,"calc_group_geometry: geometries other than 2-D, point source, and quad not implemented yet\n");
      exit(-1);
    }
#else
    clim=4;
#endif
    for(j=0;j<clim;j++){
      for(k=0;k<3;k++)
	dx[k]=corner[j][k] - grp[fault[i].group].center[k];
      pos[STRIKE] = project_vector(dx,grp[fault[i].group].strike_vec);
      pos[DIP] = project_vector(dx,grp[fault[i].group].dip_vec);
      // determine min/max
      for(k=0;k<2;k++){
	if(grp[fault[i].group].pmin[k] > pos[k])
	  grp[fault[i].group].pmin[k] =  pos[k];
	if(grp[fault[i].group].pmax[k] < pos[k])
	  grp[fault[i].group].pmax[k] =  pos[k];
      }
    }
  }
  // determine the half ranges in strike and dip direction
  for(i=0;i<medium->nrgrp;i++){
    if(grp[i].nrflt){
      grp[i].prange[DIP]=   (grp[i].pmax[DIP]-   grp[i].pmin[DIP])/2.0;
      grp[i].prange[STRIKE]=(grp[i].pmax[STRIKE]-grp[i].pmin[STRIKE])/2.0;
    }else{
      grp[i].prange[DIP]=grp[i].prange[STRIKE]=medium->nan;
    }
  }
}


/*

  determine the strike and dip angles given a vector, dip and strike
  are output in degrees

*/
void vec_to_angles(COMP_PRECISION *x,COMP_PRECISION *dip, 
		   COMP_PRECISION *strike)
{
  *dip =    vec_to_dip(x); 
  *strike = vec_to_strike(x);
  check_angles(dip,strike);
}
/*
  go the other way, from dip and strike angles (in deg) to a normalized vector
*/
void angles_to_vec(COMP_PRECISION dip, COMP_PRECISION strike, COMP_PRECISION *x)
{
  COMP_PRECISION cos_dip,sin_dip,sin_strike,cos_strike;
  my_sincos_deg(&sin_dip,&cos_dip,dip);
  my_sincos_deg(&sin_strike,&cos_strike,strike);

  x[X] = cos_dip * sin_strike;
  x[Y] = cos_dip * cos_strike;
  x[Z] = sin_dip;
}
/*
  determine the strike angle in degree from any given 3D vector
*/
COMP_PRECISION vec_to_strike(COMP_PRECISION *x)
{
  return(atan2(x[X],x[Y])*RAD2DEG);
}
/*
  determine the dip angle in degree from any given 3D vector
*/
COMP_PRECISION vec_to_dip(COMP_PRECISION *x)
{
  return(atan2(x[Z],hypot(x[X],x[Y]))*RAD2DEG);
}
/* 
   check for the dip being between 0 and 90 degrees
*/
void check_fault_angles(struct flt *fault)
{
  COMP_PRECISION dip,strike;
  dip=(COMP_PRECISION)fault->dip;strike=(COMP_PRECISION)fault->strike;
  check_angles(&dip,&strike);
  fault->dip = (float)dip; fault->strike = (float)strike;
}
/*
  make sure angles are in the right range
*/
void check_angles(COMP_PRECISION *dip,COMP_PRECISION *strike)
{
  if(*dip > 90.0){
    *dip = 180.0 - *dip;
    *strike += 180.0;
  }
  if(*dip < 0.0){
    *dip = -(*dip);
    *strike += 180.0;
  }
  fix_azimuth(strike);
}
// same for azimuth, should be between 0 and 360 degrees
void fix_azimuth(COMP_PRECISION *azi)
{
  if(*azi >= 360.0)
    *azi -= 360.0;
  if(*azi < 0.0)
    *azi = 360.0 + *azi;
}

/*
  
  returns the global coordinate x[3] given three nodes
  on the triangle (FE convention for numbering) in xt
  and local coordinates g and h
  
*/
void globalx(COMP_PRECISION *xt, COMP_PRECISION g,COMP_PRECISION h,
	     COMP_PRECISION *x)
{
  COMP_PRECISION tmp;
  tmp=1.0-g-h;
  x[X] = tmp * xt[X] + g * xt[3+X] + h * xt[6+X];
  x[Y] = tmp * xt[Y] + g * xt[3+Y] + h * xt[6+Y];
  x[Z] = tmp * xt[Z] + g * xt[3+Z] + h * xt[6+Z];
}
//
// determines the centroid of a triangular element
// given a set of points in xt
// output is xc
void calc_centroid_tri(COMP_PRECISION *xt,COMP_PRECISION *xc)
{
  xc[X] = (xt[  +X] + xt[3+X] + xt[6+X])/3.0;
  xc[Y] = (xt[  +Y] + xt[3+Y] + xt[6+Y])/3.0;
  xc[Z] = (xt[  +Z] + xt[3+Z] + xt[6+Z])/3.0;
}
//
// determines the mean coordinates of a rectangular element
// given a set of points xt  output is xc
//
void calc_mean_quad_coord(COMP_PRECISION *xq,COMP_PRECISION *xc)
{
  xc[X] = (xq[  +X] + xq[3+X] + xq[6+X] + xq[9+X])/4.0;
  xc[Y] = (xq[  +Y] + xq[3+Y] + xq[6+Y] + xq[9+Y])/4.0;
  xc[Z] = (xq[  +Z] + xq[3+Z] + xq[6+Z] + xq[9+Z])/4.0;
}
/*
  
  calculate the centroid of a quad assuming that it is in a plane
  by means of adding the area weighted centroids of two triangles
  
*/
void calc_centroid_quad(COMP_PRECISION *xq, COMP_PRECISION *xc)
{
  COMP_PRECISION area[2],xt[9],c[3];
  int i,j;
  xc[X]=xc[Y]=xc[Z]=0.0;
  for(j=0;j<2;j++){// two triangles
    for(i=0;i<3;i++){
      xt[  i]=xq[(j==1?6:0)+i];
      xt[3+i]=xq[(j==1?9:3)+i];
      xt[6+i]=xq[(j==1?0:6)+i];
    }
    area[j] = triangle_area(xt);
    calc_centroid_tri(xt,c);
    for(i=0;i<3;i++)
      xc[i] += area[j] * c[i];
  }
  area[0] += area[1];
  for(i=0;i<3;i++)
    xc[i] /= area[0];
}





#ifdef ALLOW_NON_3DQUAD_GEOM
// decide if patch is 2-D in geometry
my_boolean patch_is_2d(MODE_TYPE type)
{
  if((type == TWO_DIM_SEGMENT_PLANE_STRESS)||
     (type == TWO_DIM_SEGMENT_PLANE_STRAIN)||
     (type == TWO_DIM_HALFPLANE_PLANE_STRAIN))
    return TRUE;
  else
    return FALSE;
}
#endif


/* 

   determine position of patch in fault group in along strike and
   along dip coordinates assuming a rectangular patch

   the output is the local coordinate of the patch center, 
   normalized by the half width and half length 


   WARNING: this also initializes the group structure 
*/

void calculate_position_of_patch(struct med *medium, struct flt *fault)
{
  int i,j,nlength;
  struct geog *grp;
#ifdef USE_PGPLOT
  COMP_PRECISION group_aspect;
#endif
#ifdef DEBUG
  COMP_PRECISION pmin[2]={FLT_MAX,FLT_MAX},pmax[2]={-FLT_MAX,-FLT_MAX};
#endif
  if(!(grp=(struct geog *)calloc(medium->nrgrp,sizeof(struct geog))))
    MEMERROR("calculate_position_of_patch");
  //
  // determine average geometrical quantities of the groups
  // and calculate the patches pos[STRIKE] and pos[DIP] relative
  // to the average strike and dip vectors of the group
  //
  calc_group_geometry(medium,fault,grp);
  /*
    
     normalize the position of a patch in both strike and dip direction 
     by maximum half width W of the fault group

     if a is the aspect ratio, ie. total length/width, and n(m) the number 
     of patches in strike(dip) direction, then pos will be within the range
     -a + L/n ... a - L/n and 
     -1 + W/m ... 1 - W/m

     FOR 2D:

     normalize by L

  */
  nlength = DIP; // usually, normalize by width
#ifdef ALLOW_NON_3DQUAD_GEOM
  for(i=0;i<medium->nrflt;i++){
    if(patch_is_2d(fault[i].type)){
      nlength = STRIKE;
      break;
    }
  }
#endif
  for(i=0;i<medium->nrflt;i++){
    for(j=0;j<2;j++){
      if(fabs(grp[fault[i].group].prange[nlength])>EPS_COMP_PREC){
	fault[i].pos[j] /= grp[fault[i].group].prange[nlength];
      }
      if(fabs(fault[i].pos[j]) < EPS_COMP_PREC)
	fault[i].pos[j] = 0.0;
#ifdef DEBUG
      if(fault[i].pos[j] < pmin[j])
	pmin[j] = fault[i].pos[j];
      if(fault[i].pos[j] > pmax[j])
	pmax[j] = fault[i].pos[j];
#endif
    }
  }
#ifdef DEBUG
  for(i=0;i<2;i++)
    fprintf(stderr,"calculate_position_of_patch: dir: %i pmin: %g pmax: %g\n",
	    i,pmin[i],pmax[i]);
#endif
#ifdef USE_PGPLOT
  //
  // determine fault group array sizes
  // for plotting output, ie. reconstruct the number of patches in group
  if(medium->nrgrp == 1){
    if(grp[0].prange[DIP] && grp[0].prange[STRIKE]){
      group_aspect=(float)((int)(grp[0].prange[STRIKE]/grp[0].prange[nlength]));
      medium->grp0_n=(int)(sqrt(medium->nrflt/
				(COMP_PRECISION)group_aspect)+0.5);
      if(medium->grp0_n)
	medium->grp0_m=medium->nrflt/medium->grp0_n;
    }else if(grp[0].prange[DIP]==0 && grp[0].prange[STRIKE]!=0){
      medium->grp0_m=medium->nrflt;medium->grp0_n=1;
    }else{
      medium->grp0_n=medium->nrflt;medium->grp0_m=1;
    }
  }
#endif
  medium->pos_init=TRUE;
  free(grp);
}

/* 

compute the cartesian representation of a slip vector

 */
void compute_cartesian_slip(COMP_PRECISION *ux,
			    COMP_PRECISION *us,
			    struct flt *flt)
{
  int i;
  for(i=0;i<3;i++){
    ux[i]  = flt->t_strike[i] * us[STRIKE];
    ux[i] += flt->t_dip[i]    * us[DIP];
    ux[i] += flt->normal[i]   * us[NORMAL];
  }
}
