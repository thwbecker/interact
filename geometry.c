/*
  interact: model fault interactions using dislocations in a 
  halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu
  
  $Id: geometry.c,v 2.44 2003/12/20 21:53:59 becker Exp $
*/
#include "interact.h"
#include "properties.h"


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
  x[INT_Z]  = sin((PIHALF-plunge)/2.0);
  my_sincos_deg(&sinpl,&cospl,90.0+strike);
  x[INT_X] = x[INT_Z] * sinpl;
  x[INT_Y] = x[INT_Z] * cospl;
}

/* 
   obtain a traction vector given matrix sm and normal vector on the
   plane, left multiply
*/

void resolve_force(COMP_PRECISION *norm,COMP_PRECISION sm[3][3],
		   COMP_PRECISION *trac)
{
  int i,j;
  for(i=0;i<3;i++){
    trac[i] = 0.0;
    for(j=0;j < 3;j++)
      trac[i] += norm[j] * sm[i][j];
  }
}
/* 
   given a quad fault patch with specified angles strike and dip
   (resp. their cosines and sines), calculate the normal, dip, and
   strike vectors of that patch
   
   INPUT:
   
   the sines and cosines of the following angles:
   
   alpha: counterclockwise angle from east
   dip:   positive upward from horizontal, ie. a vertical
   fault has dip 90 degrees
   
   OUTPUT:
   strike[3], normal[3], and dip[3]
   
*/

void calc_quad_base_vecs(COMP_PRECISION *strike,COMP_PRECISION *normal,
			 COMP_PRECISION *dip,
			 COMP_PRECISION sin_alpha,
			 COMP_PRECISION cos_alpha,
			 COMP_PRECISION sin_dip, 
			 COMP_PRECISION cos_dip)
{
  /* tangential vector in strike direction, by definition 
     this vector is in the x-y plane */
  strike[INT_X]=  cos_alpha;
  strike[INT_Y]=  sin_alpha;
  strike[INT_Z]=  0.0;
  /* normal vector on fault plane */
  normal[INT_X]=   sin_alpha * sin_dip;
  normal[INT_Y]=  -cos_alpha * sin_dip;
  normal[INT_Z]=   cos_dip;
  /* tangential vector in dip direction */
  dip[INT_X]=  -sin_alpha * cos_dip;
  dip[INT_Y]=   cos_alpha * cos_dip;
  dip[INT_Z]=   sin_dip;
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

   this will also return total (not half) fault length and width

*/
void calculate_corners(COMP_PRECISION corner[4][3],
		       struct flt *fault,
		       COMP_PRECISION *l, COMP_PRECISION *w)
{
#ifdef ALLOW_NON_3DQUAD_GEOM
  if(patch_is_2d(fault->type)){	/* 2D */
    calculate_seg_corners(corner,fault,1.0);
    *l = fault->l*2;*w = NAN;
  }else if(fault->type == RECTANGULAR_PATCH){ /* quad */
    calculate_quad_corners(corner,fault,1.0);
    *l = fault->l*2;*w = fault->w*2;
  }else if(fault->type == POINT_SOURCE){ /* point source */
    calculate_point_source_corners(corner,fault,1.0,l,w);
    *l = *w = NAN;
  }else if(fault->type == TRIANGULAR){
    calculate_tri_corners(corner,fault,1.0);
    *l = *w = NAN;
  }else{
    fprintf(stderr,"calculate_corners: error, element type %i is not implemented\n",
	    fault->type);
    exit(-1);
  }
#else
  calculate_quad_corners(corner,fault,1.0);
  *l = fault->l*2;*w = fault->w*2;
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
  else if(fault->type == TRIANGULAR)
    calculate_tri_corners(corner,fault,leeway);
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
/* determine the number of vertices in each patch */
int ncon_of_patch(struct flt *fault)
{
#ifdef ALLOW_NON_3DQUAD_GEOM
  if(patch_is_2d(fault->type))
    return 2;
  else if(fault->type == RECTANGULAR_PATCH)
    return 4;
  else if(fault->type == TRIANGULAR)
    return 3;
  else if(fault->type == POINT_SOURCE)
    return 1;
  else{
    fprintf(stderr,"ncon_of_patch: mode %i undefined\n",fault->type);
    exit(-1);
  }
#else
  return 4;
#endif
}
/* determine the VTK code of the patch type */
int vtk_type_of_patch(struct flt *fault)
{
#ifdef ALLOW_NON_3DQUAD_GEOM
  if(patch_is_2d(fault->type))
    return 3;			/* vtk line */
  else if(fault->type == RECTANGULAR_PATCH)
    return 9;			     /* vtk quad */
  else if(fault->type == TRIANGULAR) /*  */
    return 5;			     /* vtk tri */
  else if(fault->type == POINT_SOURCE)
    return 1;			/* vertex */
  else{
    fprintf(stderr,"ncon_of_patch: mode %i undefined\n",fault->type);
    exit(-1);
  }
#else
  return 9;
#endif
}
/*
  
  calculate the corners of a rectangular patch, sorted FE CCW style
  from lower left

  3 --- 2
  |     |
  0 --- 1 

 */
void calculate_quad_corners(COMP_PRECISION corner[4][3],struct flt *fault,
			    COMP_PRECISION leeway)
{
  int i;
  COMP_PRECISION sx,dx;
  for(i=0;i < 3;i++){
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
#ifdef ALLOW_NON_3DQUAD_GEOM

/* 

   calculate the three vertices for a triangle those are identical to
   the ->xt array, but for visualization purposes, we also allow for
   shrinking a patch, in which case they are computed by reduced
   vectors from centroid
*/
void calculate_tri_corners(COMP_PRECISION corner[4][3],struct flt *fault,
			   COMP_PRECISION leeway)
{
  static my_boolean init = FALSE;
  COMP_PRECISION vec[3];
  int i,j;
  if(leeway == 1.0){
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	corner[i][j] = fault->xt[i*3+j];
   
  }else{
    if(!init)
      fprintf(stderr,"calculate_tri_corners: WARNING: leeway != 1 only approximate\n");
    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
	vec[j] = fault->xt[i*3+j] - fault->x[j]; /* diff from centroid */
	corner[i][j] = fault->x[j] + leeway * vec[j];
      }
    }
  }
  for(j=0;j<3;j++)
    corner[3][j] = NAN;
  init = TRUE;
}
#endif
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
  int i,j;
  COMP_PRECISION sx;
  for(i=0;i<2;i++){
    sx = fault->t_strike[i] * (COMP_PRECISION)fault->l * leeway;
    // lower left
    corner[0][i]=fault->x[i]-sx;
    // lower right
    corner[1][i]=fault->x[i]+sx;
  }
  corner[0][INT_Z] = corner[1][INT_Z] = 0.0;
  for(i=3;i<4;i++)
    for(j=0;j<3;j++)
      corner[i][j] = NAN;
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
  gvec[INT_X]=xt[3+INT_X]-xt[INT_X];
  gvec[INT_Y]=xt[3+INT_Y]-xt[INT_Y];
  gvec[INT_Z]=xt[3+INT_Z]-xt[INT_Z];
  
  hvec[INT_X]=xt[6+INT_X]-xt[INT_X];
  hvec[INT_Y]=xt[6+INT_Y]-xt[INT_Y];
  hvec[INT_Z]=xt[6+INT_Z]-xt[INT_Z];
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
  if(fabs(gvec[INT_Z])>fabs(hvec[INT_Z])) // dip should indicate depth
    swap_ab_vector_3d(gvec,hvec);
  else if(gvec[INT_Z] == hvec[INT_Z]){// in plane
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
  c_eq_a_minus_b_3d(&vec[0][INT_X],(x+3+INT_X),(x+INT_X));// get vectors along the edges
  c_eq_a_minus_b_3d(&vec[1][INT_X],(x+6+INT_X),(x+3+INT_X));
  c_eq_a_minus_b_3d(&vec[2][INT_X],(x+9+INT_X),(x+6+INT_X));
  c_eq_a_minus_b_3d(&vec[3][INT_X],(x+INT_X),(x+9+INT_X));
  cross_product(&vec[0][INT_X],&vec[1][INT_X],normal);
  if((dotp_3d(&vec[2][INT_X],normal) > eps)||
     (dotp_3d(&vec[2][INT_X],normal) > eps)){
    /* 
       fprintf(stderr,"check_planar: rectangular element points possibly not on plane\n");
       fprintf(stderr,"check_planar: points: %g %g %g; %g %g %g; %g %g %g; %g %g %g\n",
       x[INT_X],x[INT_Y],x[INT_Z], x[INT_X+3],x[INT_Y+3],x[INT_Z+3], x[INT_X+6],x[INT_Y+6],x[INT_Z+6], 
       x[INT_X+9],x[INT_Y+9],x[INT_Z+9]);
    */
    return(FALSE);
  }
  // check right angles
  if((fabs(dotp_3d(&vec[0][INT_X],&vec[1][INT_X])) > eps)||
     (fabs(dotp_3d(&vec[1][INT_X],&vec[2][INT_X])) > eps)||
     (fabs(dotp_3d(&vec[2][INT_X],&vec[3][INT_X])) > eps)){
    /*
      fprintf(stderr,"check_planar: rectangular element points possibly not all 90 degree angles\n");
      fprintf(stderr,"check_planar: points: %g %g %g; %g %g %g; %g %g %g; %g %g %g\n",
      x[INT_X],x[INT_Y],x[INT_Z], x[INT_X+3],x[INT_Y+3],x[INT_Z+3], x[INT_X+6],x[INT_Y+6],x[INT_Z+6], 
      x[INT_X+9],x[INT_Y+9],x[INT_Z+9]);
      fprintf(stderr,"check_planar: dotps: 01: %g 12: %g 23: %g\n",
      dotp_3d(&vec[0][INT_X],&vec[1][INT_X]),
      dotp_3d(&vec[1][INT_X],&vec[2][INT_X]),
      dotp_3d(&vec[2][INT_X],&vec[3][INT_X]));
    */
    return(FALSE);
  }
  return TRUE;
}



/*
  
  determine several geometrical quantities of a fault group that
  consists of several patches (assumes planar faults)
  
*/

void calc_group_geometry(struct med *medium,struct flt *fault,
			 struct geog *grp)
{
  int i,j,k,clim,igrp;
  COMP_PRECISION dx1[3],dx2[3],corner[4][3],fac,pos[2],l,w,dist_max,dist;
  static int group_geom_mode = 1;
#ifdef ALLOW_NON_3DQUAD_GEOM
  COMP_PRECISION global_dip_rad,sin_global_dip_rad,cos_global_dip_rad,
    gnormal[3],gstrike[3],gdip[3];  
#endif
  /*
     determine center of mass and average strike 
     and dip vectors for each patch group 
     
  */
  for(i=0;i<medium->nrflt;i++){
    igrp = fault[i].group;
    grp[igrp].nrflt++;
    
    add_b_to_a_vector_3d(grp[igrp].center,fault[i].x);
#ifdef ALLOW_NON_3DQUAD_GEOM
    if(fault[i].type == TRIANGULAR){
      /* 
	 compute appropriate projection vectors 

	 for triangle, global
      */
      global_dip_rad   = DEG2RADF((COMP_PRECISION)fault[i].dip);
      my_sincos(&sin_global_dip_rad,&cos_global_dip_rad,global_dip_rad);

      calc_quad_base_vecs(gstrike, gnormal, gdip,
			  fault[i].sin_alpha, fault[i].cos_alpha,
			  sin_global_dip_rad, cos_global_dip_rad);
      //fprintf(stderr,"%g %g %g %g %g\n",fault[i].strike,fault[i].dip,gstrike[0],gstrike[1],gstrike[2]);
      a_equals_b_vector_3d(dx1,gstrike);
      a_equals_b_vector_3d(dx2,gdip);
    }else{
      /* for quad, local */
      a_equals_b_vector_3d(dx1,fault[i].t_strike);
      a_equals_b_vector_3d(dx2,fault[i].t_dip);
    }
#else
    a_equals_b_vector_3d(dx1,fault[i].t_strike);
    a_equals_b_vector_3d(dx2,fault[i].t_dip);
#endif
    if(group_geom_mode == 1){
      /* new mode, weighted by patch area */
      scale_vector_3d(dx1,fault[i].area);
      scale_vector_3d(dx2,fault[i].area);
    }
    add_b_to_a_vector_3d(grp[igrp].strike_vec,dx1);
    add_b_to_a_vector_3d(grp[igrp].dip_vec,dx2);
  }
  for(i=0;i < medium->nrgrp;i++){
    if(grp[i].nrflt){// get average
      fac=1.0/((COMP_PRECISION)grp[i].nrflt);
      scale_vector_3d(grp[i].center,fac);
      /* don't need to normalize those by number or weight, will be unity
	 normalized */
      normalize_3d(grp[i].strike_vec);
      normalize_3d(grp[i].dip_vec);
    }
    for(j=0;j<2;j++){
      grp[i].pmin[j] = FLT_MAX;
      grp[i].pmax[j] =-FLT_MAX;
    }
  }
  /* 

     assign local position to patches in terms of distance in strike
     and dip direction

  */
  dist_max = 0;
  for(i=0;i < medium->nrflt;i++){
    // determine distance of center of patch from 
    // center of group of patches
    c_eq_a_minus_b_3d(dx1,fault[i].x,grp[fault[i].group].center);
    dist = norm_3d(dx1);
    if(dist > dist_max)
      dist_max = dist;
    fault[i].pos[STRIKE]=(float)
      dotp_3d(dx1,grp[fault[i].group].strike_vec);
    fault[i].pos[DIP]=(float)
      dotp_3d(dx1,grp[fault[i].group].dip_vec);
  }

  for(i=0;i < medium->nrflt;i++){
#ifdef ATZ_NATZ_VOODOO_GNATZ
    fprintf(stderr,"mx: %g %g %g ms: %g %g %g md: %g %g %g P: %6.3f %6.3f\n",
	    grp[fault[i].group].center[INT_X],grp[fault[i].group].center[INT_Y],
	    grp[fault[i].group].center[INT_Z],
	    grp[fault[i].group].strike_vec[INT_X],
	    grp[fault[i].group].strike_vec[INT_Y],
	    grp[fault[i].group].strike_vec[INT_Z],
	    grp[fault[i].group].dip_vec[INT_X],
	    grp[fault[i].group].dip_vec[INT_Y],
	    grp[fault[i].group].dip_vec[INT_Z],
	    fault[i].pos[STRIKE]/dist_max,fault[i].pos[DIP]/dist_max);
#endif
    // determine actual extent of corners of patch
    calculate_corners(corner,(fault+i),&l,&w);
    clim = ncon_of_patch((fault+i));
    /* here, we treate point source as having quasi corners (?!) */
    if(clim == 1)
      clim = 4;
    for(j=0;j < clim;j++){
      for(k=0;k<3;k++)
	dx1[k] = corner[j][k] - grp[fault[i].group].center[k];
      pos[STRIKE] = dotp_3d(dx1,grp[fault[i].group].strike_vec);
      pos[DIP] = dotp_3d(dx1,grp[fault[i].group].dip_vec);
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
  double alpha;
  *dip =    vec_to_dip(x); 
  alpha = vec_to_strike(x);
  check_angles(dip,&alpha);
  *strike = (COMP_PRECISION)alpha;
}
/*
  go the other way, from dip and strike angles (in deg) to a normalized vector
*/
void angles_to_vec(COMP_PRECISION dip, COMP_PRECISION strike, COMP_PRECISION *x)
{
  COMP_PRECISION cos_dip,sin_dip,sin_strike,cos_strike;
  my_sincos_deg(&sin_dip,&cos_dip,dip);
  my_sincos_deg(&sin_strike,&cos_strike,strike);

  x[INT_X] = cos_dip * sin_strike;
  x[INT_Y] = cos_dip * cos_strike;
  x[INT_Z] = sin_dip;
}
/*
  determine the strike angle in degree from any given 3D vector
*/
COMP_PRECISION vec_to_strike(COMP_PRECISION *x)
{
  return(atan2(x[INT_X],x[INT_Y])*RAD2DEG);
}
/*
  determine the dip angle in degree from any given 3D vector
*/
COMP_PRECISION vec_to_dip(COMP_PRECISION *x)
{
  return(atan2(x[INT_Z],hypot(x[INT_X],x[INT_Y]))*RAD2DEG);
}
/* 
   check for the dip being between 0 and 90 degrees
*/
void check_fault_angles(struct flt *fault)
{
  COMP_PRECISION dip;
  double strike;
  dip=(COMP_PRECISION)fault->dip;strike=(COMP_PRECISION)fault->strike;
  check_angles(&dip,&strike);
  fault->dip = (float)dip; fault->strike = (COMP_PRECISION)strike;
}
/*
  make sure angles are in the right range
*/
void check_angles(COMP_PRECISION *dip,double *strike)
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
void fix_azimuth(double *azi)
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
  x[INT_X] = tmp * xt[INT_X] + g * xt[3+INT_X] + h * xt[6+INT_X];
  x[INT_Y] = tmp * xt[INT_Y] + g * xt[3+INT_Y] + h * xt[6+INT_Y];
  x[INT_Z] = tmp * xt[INT_Z] + g * xt[3+INT_Z] + h * xt[6+INT_Z];
}
//
// determines the centroid of a triangular element
// given a set of points in xt
// output is xc
void calc_centroid_tri(COMP_PRECISION *xt,COMP_PRECISION *xc)
{
  xc[INT_X] = (xt[  +INT_X] + xt[3+INT_X] + xt[6+INT_X])/3.0;
  xc[INT_Y] = (xt[  +INT_Y] + xt[3+INT_Y] + xt[6+INT_Y])/3.0;
  xc[INT_Z] = (xt[  +INT_Z] + xt[3+INT_Z] + xt[6+INT_Z])/3.0;
}
//
// determines the mean coordinates of a rectangular element
// given a set of points xt  output is xc
//
void calc_mean_quad_coord(COMP_PRECISION *xq,COMP_PRECISION *xc)
{
  xc[INT_X] = (xq[  +INT_X] + xq[3+INT_X] + xq[6+INT_X] + xq[9+INT_X])/4.0;
  xc[INT_Y] = (xq[  +INT_Y] + xq[3+INT_Y] + xq[6+INT_Y] + xq[9+INT_Y])/4.0;
  xc[INT_Z] = (xq[  +INT_Z] + xq[3+INT_Z] + xq[6+INT_Z] + xq[9+INT_Z])/4.0;
}
/*
  
  calculate the centroid of a quad assuming that it is in a plane
  by means of adding the area weighted centroids of two triangles
  
*/
void calc_centroid_quad(COMP_PRECISION *xq, COMP_PRECISION *xc)
{
  COMP_PRECISION area[2],xt[9],c[3];
  int i,j;
  xc[INT_X]=xc[INT_Y]=xc[INT_Z]=0.0;
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
    PMEMERROR("calculate_position_of_patch");
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
      if(fabs(grp[fault[i].group].prange[nlength]) > EPS_COMP_PREC){
	fault[i].pos[j] /= grp[fault[i].group].prange[nlength];
      }
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
    HEADNODE
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

/*

  calculate the background stress (sm[3][3] matrix) at location x[3]
  and time "time" given the constant stress matrix factors a[6] and the 
  loading rates b[6] as well as the pressure "pressure"

  

*/

void background_stress(COMP_PRECISION sm[3][3], COMP_PRECISION *x, 
		       COMP_PRECISION time,COMP_PRECISION *a,
		       COMP_PRECISION *b,COMP_PRECISION pressure)
{
  COMP_PRECISION locp;
#ifdef HYDROSTATIC_PRESSURE
  locp = -(x[INT_Z]/HYDROSTATIC_PRESSURE) * pressure; 
#else
  locp = pressure; 
#endif
  /* isotropic elements, compression negative */
  sm[INT_X][INT_X] = a[0] + time * b[0] - locp;
  sm[INT_Y][INT_Y] = a[3] + time * b[3] - locp;
  sm[INT_Z][INT_Z] = a[5] + time * b[5] - locp;
  /* off diagonal elements  */  
  sm[INT_X][INT_Y]=sm[INT_Y][INT_X] = a[1] + time * b[1];
  sm[INT_X][INT_Z]=sm[INT_Z][INT_X] = a[2] + time * b[2];
  sm[INT_Y][INT_Z]=sm[INT_Z][INT_Y] = a[4] + time * b[4];
}
void background_disp(COMP_PRECISION *u, COMP_PRECISION *x, 
		     struct med *medium,COMP_PRECISION *a,
		     COMP_PRECISION *b)
{
  /* the characteristic strain rate for simple shear 
     is given by 

     characteristic stressing rate
     --------------------------
     2 mu

     integration gives the characteristic strain

     characteristic stressing rate
     ----------------------------- y_location
     mu
     
  */
  int i;
  my_boolean hit=FALSE;
  for(i=0;i<6;i++)
    if(a[i]!=0.0||((i!=1)&&(b[i]!=0.0))){hit=TRUE;break;}
  if(hit){
    fprintf(stderr,"background_disp: EXITING: background displacement is inaccurate since no simple shear stressing\n");
    exit(-1);
  }
  u[INT_X]=medium->time * (b[1]/SHEAR_MODULUS)*u[INT_Y];
  u[INT_Y]=u[INT_Z]=0.0;
}
/*

  obtain the local coordinates given the base vectors vec_1 and vec_2


*/
void get_local_x_on_plane(COMP_PRECISION *xl,COMP_PRECISION *x,
			  COMP_PRECISION *flt_mean_x,COMP_PRECISION *vec_1,
			  COMP_PRECISION *vec_2)
{
  int i;
  for(i=0;i<3;i++){
    xl[i]  = flt_mean_x[i];
    xl[i] += vec_1[i] * x[INT_X];
    xl[i] += vec_2[i] * x[INT_Y];
  }
}
// this geometry routine moved to eval_triangle
// void get_tri_prop_based_on_gh(struct flt *fault)
/*



  obtain average faulk plane vectors and location
  on return flt_mean_x will hold the mean location and vec_1 and vec_2
  the mean strike and dip or the mean strike and normal vectors, depending
  on the n[INT_Z] flag, -1 or -2 
  
  doesn't make sense for triangular

 */
void get_fault_plane_basevec(COMP_PRECISION *flt_mean_x,
			     COMP_PRECISION *vec_1,COMP_PRECISION *vec_2,
			     struct flt *fault,struct med *medium)
{
  int n,i,j;
  // get average fault plane vectors
  // and mean location of patches
  for(i=0;i<3;i++)
    vec_1[i]=vec_2[i]=flt_mean_x[i]=0.0;
  if(medium->n[INT_Z] == -1){
    fprintf(stderr,"get_fault_plane_basevec: base vectors are average strike and dip of fault group 0\n");
    for(n=i=0;i<medium->nrflt;i++)
      if(fault[i].group == 0){
	n++;
	for(j=0;j<3;j++){
	  flt_mean_x[j]   += fault[i].x[j];
	  vec_1[j]        += fault[i].t_strike[j];
	  vec_2[j]        += fault[i].t_dip[j];
	}
      }
  }else if(medium->n[INT_Z] == -2){
    fprintf(stderr,"get_fault_plane_basevec: base vectors are average strike and normal of fault group 0\n");
    for(n=i=0;i<medium->nrflt;i++)
      if(fault[i].group == 0){
	n++;
	for(j=0;j<3;j++){
	  flt_mean_x[j]   += fault[i].x[j];
	  vec_1[j]        += fault[i].t_strike[j];
	  vec_2[j]        += fault[i].normal[j];
	}
      }
  }else{
    fprintf(stderr,"get_fault_plane_basevec: medium->n[Z] has to be -1 or -2 but is %i\n",
	    medium->n[INT_Z]);
    exit(-1);
  }
  if(n)
    for(i=0;i<3;i++){
      flt_mean_x[i]  /=(COMP_PRECISION)n;
      if(fabs(flt_mean_x[i])<EPS_COMP_PREC)
	flt_mean_x[i]=0.0;
      vec_1[i]  /=(COMP_PRECISION)n;
      if(fabs(vec_1[i])<EPS_COMP_PREC)
	vec_1[i]=0.0;
      vec_2[i]     /=(COMP_PRECISION)n;
      if(fabs(vec_2[i])<EPS_COMP_PREC)
	vec_2[i]=0.0;
    }
  normalize_3d(vec_1);normalize_3d(vec_2);
}
/*

  calculate the deviator stress and pressure
  given a stress matrix sm
  
  output is dm and pressure (of original tensor), and second invariant
  of deviatoric tensor

*/

void calc_deviatoric_stress(COMP_PRECISION sm[3][3],COMP_PRECISION dm[3][3],
			    COMP_PRECISION *pressure, COMP_PRECISION *s2)
{
  int i,j;
  *pressure= -(sm[INT_X][INT_X] + sm[INT_Y][INT_Y] + sm[INT_Z][INT_Z])/3.0;
 
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      dm[i][j] = sm[i][j] + ((i==j)?(*pressure):(0.0));
  /* second invariant of deviatoric tensor */
  *s2 = sqrt(0.5* (dm[INT_X][INT_X] * dm[INT_X][INT_X] + 
		   dm[INT_X][INT_Y] * dm[INT_X][INT_Y] * 2.0 + 
		   dm[INT_Y][INT_Y] * dm[INT_Y][INT_Y] + 
		   dm[INT_Y][INT_Z] * dm[INT_Y][INT_Z] * 2.0 + 
		   dm[INT_Z][INT_Z] * dm[INT_Z][INT_Z] + 
		   dm[INT_X][INT_Z] * dm[INT_X][INT_Z] * 2.0));
}


#ifdef ALLOW_NON_3DQUAD_GEOM

/* 
   compute the global basis vectors for a triangular patch as
   specified by its strike and dip using the rectangular patch
   convention

*/
void calc_global_strike_dip_from_local(struct flt *fault,
				       COMP_PRECISION *gstrike, COMP_PRECISION *gnormal, COMP_PRECISION *gdip)
{
  COMP_PRECISION global_dip_rad,sin_global_dip_rad,cos_global_dip_rad;
#ifdef DEBUG
  if(fault->type != TRIANGULAR){
    fprintf(stderr,"calc_global_strike_dip_from_local: called for a non triangular patch?!\n");
    exit(-1);
  }
#endif
  /* 
     compute appropriate projection vectors, triangular fault holds
     global dip and strike as sin/cos alpha
  */
  global_dip_rad   = DEG2RADF((COMP_PRECISION)fault->dip);
  my_sincos(&sin_global_dip_rad,&cos_global_dip_rad,global_dip_rad);
  calc_quad_base_vecs(gstrike, gnormal, gdip,
		      fault->sin_alpha, fault->cos_alpha,
		      sin_global_dip_rad,   cos_global_dip_rad);
}
/* 
   compute the projection of fault local properties to global assuming
   projection vectors already computed
   input slip and tractions are fault local

   vectors are sorted STRIKE, NORMAL, DIP
*/
void calc_global_slip_and_traction_from_local(struct flt *fault,COMP_PRECISION *slip,COMP_PRECISION *trac,
					      COMP_PRECISION *gstrike, COMP_PRECISION *gnormal, COMP_PRECISION *gdip,
					      COMP_PRECISION *gslip,COMP_PRECISION *gtrac,my_boolean only_traction)
{
  int j;
  COMP_PRECISION lslip[3],ltrac[3];
#ifdef DEBUG
  if(fault->type != TRIANGULAR){
    fprintf(stderr,"calc_global_strike_dip_from_local: called for a non triangular patch?!\n");
    exit(-1);
  }
#endif
  if(only_traction){
    for(j=0;j<3;j++){		/* global slip and strike vectors */
      ltrac[j] = fault->t_strike[j] * trac[STRIKE] + fault->t_dip[j] * trac[DIP];// + fault->normal[j] * trac[NORMAL];
    }
  }else{
    for(j=0;j<3;j++){		/* global slip and strike vectors */
      lslip[j] = fault->t_strike[j] * slip[STRIKE] + fault->t_dip[j] * slip[DIP];//+ fault->normal[j] * slip[NORMAL];
      ltrac[j] = fault->t_strike[j] * trac[STRIKE] + fault->t_dip[j] * trac[DIP];// + fault->normal[j] * trac[NORMAL];
    }
  }
  if(!only_traction){
    /* project into global coordinate system */
    gslip[0] = dotp_3d(lslip,gstrike);
    //gslip[1] = dotp_3d(lslip,gnormal);
    gslip[1] = slip[NORMAL];
    gslip[2] = dotp_3d(lslip,gdip);
  }
  gtrac[0] = dotp_3d(ltrac,gstrike);
  //gtrac[1] = dotp_3d(ltrac,gnormal);
  gtrac[1] = trac[NORMAL];
  gtrac[2] = dotp_3d(ltrac,gdip);
}
#endif
