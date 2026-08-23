/*
  interact: model fault interactions using dislocations in a 
            halfspace

	    (C) Thorsten Becker, thwbecker@post.harvard.edu

  part of interact, (C) Thorsten Becker; see README.md and COPYRIGHT

  evaluate 2-D segment given the solution from Crouch and Starfield
  (1974), p. 91

  @Book{crouch83,
  author =	 "Crouch, S. L. and Starfield, A. M.",
  title =	 "Boundary Element Methods in Solid Mechanics. With
                  Applications in Rock Mechanics",
  publisher =	 "Allen and Unwin",
  year =	 1983,
  address =	 "London"
  }

  with respect to Dy of Crouch & Starfield, we flip the sign for
  normal direction slip


  in this case, our coordinate system is x-y where y <= 0 for a half
  space setup, and strike is counted in degrees clockwise from up, +y

  slip is typically in plane shear (mode II) or in plane normal mode
  (I) but can be out of plane, screw disclocation type (mode III)


*/
#include "interact.h"
#include "properties.h"


void eval_2dsegment_plane_strain(COMP_PRECISION *x,
				 struct flt *fault,
				 COMP_PRECISION *disp,
				 COMP_PRECISION *u_global, 
				 COMP_PRECISION sm_global[3][3],
				 int *iret, MODE_TYPE mode,
				 struct el_par elastic)
{
  COMP_PRECISION u[3],x_local[3],dx[3],sm_local[3][3],f2,f3,f4,f5,f6,f7,
    c1,c12,c2,c22,c3,c32,c4,c42,y2;
  COMP_PRECISION 
    pfac2 = (1.0 -       elastic.poisson), 
    pfac1 = 4.0 * PI * pfac2,
    pfac3 = (1.0 - 2.0 * elastic.poisson);

#ifdef DEBUG
  COMP_PRECISION l,w,corners[MAX_NR_EL_VERTICES*3];
  if(x[INT_Z] != 0.0){
    fprintf(stderr,"eval_2dsegment_plane_strain: z coordinate has to be zero: z: %g\n",
	    x[INT_Z]);
    exit(-1);
  }
#endif
  /* 
     shift and rotate observational point into local reference frame. 
     first, move fault to origin
  */
  dx[INT_X] = x[INT_X] - fault->x[INT_X];
  dx[INT_Y] = x[INT_Y] - fault->x[INT_Y];
  //
  // then line up fault strike with the x-axis, ie. rotate dx into the
  // fault local system x_local
  //
  rotate_vec2d(dx,x_local,fault->cos_alpha,fault->sin_alpha);
  //
  // geometrical factors
  //
  get_2dseg_geo(x_local,fault->l,&c1, &c2,&c12, &c22,&c3, 
		&c4,&c32, &c42,&y2,iret);
  if(*iret){
    //
    // solution is infinite
    //
#ifdef DEBUG
    calculate_vertices(corners,fault,&l,&w);
    fprintf(stderr,
	    "eval_2dsegment_plane_strain: SINGULAR: x: (%g, %g) segment: (%g, %g) to (%g, %g)\n",
	    x_local[INT_X],x_local[INT_Y],corners[0*3+INT_X],corners[0*3+INT_Y],
	    corners[1*3+INT_X],corners[1*3+INT_Y]);
#endif
    set_stress_and_disp_nan(sm_global,u_global,mode);
  }else{
    // get f factors
    get_2dseg_ffac(&f2,&f3,&f4,&f5,&f6,&f7,c1,c2,c12,c22,c3,c4,
		   c32,c42,pfac1,y2,x_local);
    // displacements
    get_2dseg_disp(u,disp,x_local,f2,f3,f4,f5,pfac2,pfac3);
    //
    // plane strain solution
    //
    // stresses
    get_2dseg_stress(sm_local,disp,x_local,f4,f5,f6,f7,elastic);
    // no shear stresses with z for plain strain since exz and eyz = 0
    sm_local[INT_X][INT_Z] = sm_local[INT_Z][INT_X] = sm_local[INT_Z][INT_Y] = sm_local[INT_Y][INT_Z] = 0.0;
    // from plane strain: szz = \nu/ (sxx + syy)
    sm_local[INT_Z][INT_Z] = elastic.poisson * (sm_local[INT_X][INT_X] + sm_local[INT_Y][INT_Y]);
    /* 

       rotate displacements back into global frame this subroutine
       takes the first three components of u, ie. u_x u_y u_z

    */
    if(mode != GC_STRESS_ONLY){
      rotate_vec2d(u,u_global,fault->cos_alpha,-fault->sin_alpha);
      u_global[INT_Z] = 0.0;
    }
    if(mode != GC_DISP_ONLY){
      /* 
	 rotate stress matrix back into global field 
      */
      rotate_mat_z(sm_local,sm_global,fault->cos_alpha,-fault->sin_alpha);
    }
  }
}
//
// repeat exercise for plane stress which we can obtain by
// letting  \nu --> \nu/(1+\nu) in the above formulation for plane strain
//
void eval_2dsegment_plane_stress(COMP_PRECISION *x,struct flt *fault,
				 COMP_PRECISION *disp,
				 COMP_PRECISION *u_global, 
				 COMP_PRECISION sm_global[3][3],int *iret,
				 MODE_TYPE mode,struct el_par elastic)
{  
  COMP_PRECISION u[3],x_local[3],dx[3],sm_local[3][3],f2,f3,f4,f5,f6,f7,
    c1,c12,c2,c22,c3,c32,c4,c42,y2;
  COMP_PRECISION 
    // the last factror, pfac5, uses the real poisson ratio for e_zz
    pfac0 = elastic.poisson/(1.0+elastic.poisson),
    pfac1 = 4.0 * PI * (1.0 - (elastic.poisson/(1.0+elastic.poisson))),
    pfac2 = 1.0 -       pfac0, 
    pfac3 = 1.0 - 2.0 * pfac0,
    pfac5 = elastic.poisson/(2.0*elastic.shear*(1.0+elastic.poisson));
#ifdef DEBUG
  COMP_PRECISION l,w,corners[MAX_NR_EL_VERTICES*3];
  if(x[INT_Z] != 0.0){
    fprintf(stderr,"eval_2dsegment_plane_stress: z coordinate has to be zero: z: %g\n",
	    x[INT_Z]);
    exit(-1);
  }
#endif
  dx[INT_X] = x[INT_X] - fault->x[INT_X];
  dx[INT_Y] = x[INT_Y] - fault->x[INT_Y];
  rotate_vec2d(dx,x_local,fault->cos_alpha,fault->sin_alpha);
  get_2dseg_geo(x_local,fault->l,&c1, &c2,&c12, &c22,&c3, 
		&c4,&c32, &c42,&y2,iret);
  if(*iret){
#ifdef DEBUG
    calculate_vertices(corners,fault,&l,&w);
    fprintf(stderr,
	    "eval_2dsegment_plane_stress: SINGULAR: x: (%g, %g) segment: (%g, %g) to (%g, %g)\n",
	    x_local[INT_X],x_local[INT_Y],corners[0*3+INT_X],corners[0*3+INT_Y],
	    corners[1*3+INT_X],corners[1*3+INT_Y]);
#endif
    set_stress_and_disp_nan(sm_global,u_global,mode);
  }else{
    get_2dseg_ffac(&f2,&f3,&f4,&f5,&f6,&f7,c1,c2,c12,c22,
		   c3,c4,c32,c42,pfac1,y2,x_local);
    get_2dseg_disp(u,disp,x_local,f2,f3,f4,f5,pfac2,pfac3);
    get_2dseg_stress(sm_local,disp,x_local,f4,f5,f6,f7,elastic);
    // no shear stresses with z 
    sm_local[INT_X][INT_Z] = sm_local[INT_Z][INT_X] =
      sm_local[INT_Z][INT_Y] = sm_local[INT_Y][INT_Z] = 0.0;
    sm_local[INT_Z][INT_Z] = 0.0;
    if(mode != GC_STRESS_ONLY)
      rotate_vec2d(u,u_global,fault->cos_alpha,-fault->sin_alpha);
    if(mode != GC_DISP_ONLY){
      rotate_mat_z(sm_local,sm_global,fault->cos_alpha,-fault->sin_alpha);
      // u_z is undefined without plate thickness, but we shall volunteer
      // e_zz = - pfac5 (s_xx + s_yy) instead
      if(mode != GC_STRESS_ONLY)
	u_global[INT_Z] = -pfac5 * (sm_global[INT_X][INT_X] + sm_global[INT_Y][INT_Y]);
    }
  }

}
/*

here follow copies for faults at origin and strike angle 90

*/

void eval_2dsegment_plane_strain_basic(COMP_PRECISION *x,
				       struct flt *fault,
				       COMP_PRECISION *disp,
				       COMP_PRECISION *u_global, 
				       COMP_PRECISION sm_global[3][3],
				       int *iret,struct el_par elastic)
{
  COMP_PRECISION f2,f3,f4,f5,f6,f7,c1,c12,c2,c22,c3,c32,c4,
    c42,y2;
  COMP_PRECISION
    pfac0 = 1.0 - elastic.poisson,
    pfac1 = 4.0 * PI * pfac0,
    pfac3 = (1.0 - 2.0 * elastic.poisson);
#ifdef DEBUG
  COMP_PRECISION l,w,corners[MAX_NR_EL_VERTICES*3];
  if(x[INT_Z] != 0.0){
    fprintf(stderr,"eval_2dsegment_plane_strain_basic: z coordinate has to be zero: z: %g\n",
	    x[INT_Z]);exit(-1);
  }
  if(norm(fault->x,2) > EPS_COMP_PREC){
    fprintf(stderr,"eval_2dsegment_plane_strain_basic: error, fault origin should be 0,0 but is %g, %g\n",
	    fault->x[INT_X],fault->x[INT_Y]);exit(-1);
  }
#endif
  if(fault->strike != 90.0){
    fprintf(stderr,"eval_2dsegment_plane_strain_basic: error, fault strike should be 90 but is %g\n",
	    fault->strike);exit(-1);
  }
  get_2dseg_geo(x,fault->l,&c1,&c2,&c12,&c22,&c3,&c4,&c32,
		&c42,&y2,iret);
  if(*iret){
#ifdef DEBUG
    calculate_vertices(corners,fault,&l,&w);
    fprintf(stderr,
	    "eval_2dsegment_plane_strain_basic: SINGULAR: x: (%g, %g) segment: (%g, %g) to (%g, %g)\n",
	    x[INT_X],x[INT_Y],corners[0*3+INT_X],corners[0*3+INT_Y],
	    corners[1*3+INT_X],corners[1*3+INT_Y]);
#endif
    set_stress_and_disp_nan(sm_global,u_global,GC_DISP_AND_STRESS);
  }else{
    get_2dseg_ffac(&f2,&f3,&f4,&f5,&f6,&f7,c1,c2,c12,c22,c3,c4,
		   c32,c42,pfac1,y2,x);
    get_2dseg_disp(u_global,disp,x,f2,f3,f4,f5,pfac0,pfac3);
    get_2dseg_stress(sm_global,disp,x,f4,f5,f6,f7,elastic);
    sm_global[INT_X][INT_Z] = sm_global[INT_Z][INT_X] = sm_global[INT_Z][INT_Y] =
      sm_global[INT_Y][INT_Z] = 0.0;
    sm_global[INT_Z][INT_Z] = elastic.poisson * (sm_global[INT_X][INT_X] + sm_global[INT_Y][INT_Y]);
    u_global[INT_Z] = 0.0;
  }
}
void eval_2dsegment_plane_stress_basic(COMP_PRECISION *x,
				       struct flt *fault,
				       COMP_PRECISION *disp,
				       COMP_PRECISION *u_global, 
				       COMP_PRECISION sm_global[3][3],
				       int *iret,struct el_par elastic)
{  
  COMP_PRECISION f2,f3,f4,f5,f6,f7,c1,c12,c2,c22,c3,c32,c4,c42,
    y2;
  COMP_PRECISION
    pfac0 = elastic.poisson/(1.0+elastic.poisson),
    pfac2 = (1.0 - pfac0),
    pfac1 = 4.0 * PI * pfac2,
    pfac3 = (1.0 - 2.0 * pfac0),
    pfac5 = elastic.poisson/(elastic.mu2*(1.0+elastic.poisson));
#ifdef DEBUG
  COMP_PRECISION l,w,corners[MAX_NR_EL_VERTICES*3];
  if(x[INT_Z] != 0.0){
    fprintf(stderr,"eval_2dsegment: z coordinate has to be zero: z: %g\n",
	    x[INT_Z]);exit(-1);
  }  
  if(norm(fault->x,2) > EPS_COMP_PREC){
    fprintf(stderr,"eval_2dsegment_plane_stress_basic: error, fault origin should be 0,0 but is %g, %g\n",
	    fault->x[INT_X],fault->x[INT_Y]);exit(-1);
  }
#endif
  if(fault->strike != 90.0){
    fprintf(stderr,"eval_2dsegment_plane_stress_basic: error, fault strike should be 90 but is %g\n",
	    fault->strike);exit(-1);
  }
  get_2dseg_geo(x,fault->l,&c1,&c2,&c12,&c22,&c3,&c4,&c32,
		&c42,&y2,iret);
  if(*iret){
#ifdef DEBUG
    calculate_vertices(corners,fault,&l,&w);
    fprintf(stderr,
	    "eval_2dsegment_plane_stress_basic: SINGULAR: x: (%g, %g) segment: (%g, %g) to (%g, %g)\n",
	    x[INT_X],x[INT_Y],corners[0*3+INT_X],corners[0*3+INT_Y],corners[1*3+INT_X],
	    corners[1*3+INT_Y]);
#endif
    set_stress_and_disp_nan(sm_global,u_global,GC_DISP_AND_STRESS);
  }else{
    get_2dseg_ffac(&f2,&f3,&f4,&f5,&f6,&f7,c1,c2,c12,c22,
		   c3,c4,c32,c42,pfac1,y2,x);
    get_2dseg_disp(u_global,disp,x,f2,f3,f4,f5,pfac2,pfac3);
    get_2dseg_stress(sm_global,disp,x,f4,f5,f6,f7,elastic);
    sm_global[INT_X][INT_Z] = sm_global[INT_Z][INT_X] =
      sm_global[INT_Z][INT_Y] = sm_global[INT_Y][INT_Z] = 0.0;
    sm_global[INT_Z][INT_Z] = 0.0;
    u_global[INT_Z] = -pfac5 * (sm_global[INT_X][INT_X] + sm_global[INT_Y][INT_Y]);
  }

}


/*

get fault geometry dependent c factors, and y2

*/
void get_2dseg_geo(COMP_PRECISION *x,COMP_PRECISION l,
		   COMP_PRECISION *c1,COMP_PRECISION *c2,
		   COMP_PRECISION *c12,COMP_PRECISION *c22,
		   COMP_PRECISION *c3,COMP_PRECISION *c4,
		   COMP_PRECISION *c32,COMP_PRECISION *c42,
		   COMP_PRECISION *y2, int *iret)
{
  *c1  = x[INT_X] - l;// x - a
  *c2  = x[INT_X] + l;// x + a
  *y2  = x[INT_Y] * x[INT_Y];
  *c12 = (*c1) * (*c1);
  *c22 = (*c2) * (*c2);
  *c3 = (*c12) + (*y2);
  *c4 = (*c22) + (*y2);
  *iret = 0;
  if(*c3 <= EPS_COMP_PREC_FOR_NAN){
    *iret = 1;
    *c32 = 0.0;
  }else{
    *c32 = (*c3)*(*c3);
  }
  if(*c4 <= EPS_COMP_PREC_FOR_NAN){
    *iret = 1;
    *c42 = 0.0;
  }else{
    *c42 = (*c4)*(*c4);
  }
}
/*
   calculate F factors, see Crouch & starfield, p.58 and 91


   TO DO: make sure that the atan2 works the same on all
   platforms

*/
void get_2dseg_ffac(COMP_PRECISION *f2,  COMP_PRECISION *f3,
		    COMP_PRECISION *f4,  COMP_PRECISION *f5,
		    COMP_PRECISION *f6,  COMP_PRECISION *f7,
		    COMP_PRECISION c1,   COMP_PRECISION c2,
		    COMP_PRECISION c12,  COMP_PRECISION c22,
		    COMP_PRECISION c3,   COMP_PRECISION c4,
		    COMP_PRECISION c32,  COMP_PRECISION c42,
		    COMP_PRECISION pfac1,COMP_PRECISION y2,
		    COMP_PRECISION *x)
{
  *f2 =  0.5*(log(c3)              -              log(c4))/pfac1;
  // this is problematic for y=0 and x=a, see c & S, p. 51
  //
  // lim(atan(y/(x-a)) - atan(y/(x+a)) =  0  for |x|>a, y=0
  //                                   =  pi for |x|<a y=0+
  //                                   = -pi for |x|<a,y=0-
  //
  *f3 =     -(atan2(x[INT_Y],c1)       -       atan2(x[INT_Y],c2))/pfac1;
  *f4 =      (x[INT_Y]/c3              -              x[INT_Y]/c4)/pfac1;
  *f5 =      (        c1/c3        -                c2/c4)/pfac1;
  *f6 =      ((c12 - y2)/c32       -         (c22-y2)/c42)/pfac1;
  *f7 =      2.0*x[INT_Y]*(c1/c32-                     c2/c42)/pfac1;
}
		   
/*
  
assign part of the displacements, u[INT_Z] depends on approximation

*/
void get_2dseg_disp(COMP_PRECISION *u,COMP_PRECISION *disp,
		    COMP_PRECISION *x,
		    COMP_PRECISION f2, COMP_PRECISION f3,
		    COMP_PRECISION f4, COMP_PRECISION f5,
		    COMP_PRECISION pfac2, COMP_PRECISION pfac3)
{
  // displacements, disp[NORMAL] = -Dy
  COMP_PRECISION tmp[4];
  tmp[0] = 2.0*f3*pfac2;
  tmp[1] = f2*pfac3;
  tmp[2] = f5*x[INT_Y];
  tmp[3] = f4*x[INT_Y];
  u[INT_X]  = disp[STRIKE] * (     tmp[0] - tmp[2]);
  u[INT_X] -= disp[NORMAL] * (-    tmp[1] - tmp[3]);
  u[INT_Y]  = disp[STRIKE] * (     tmp[1] - tmp[3]);
  u[INT_Y] -= disp[NORMAL] * (     tmp[0] + tmp[2]);
}
/*
  
assign part of the stress matrix, rest depends on plane strain/stress
approximation

*/
void get_2dseg_stress(COMP_PRECISION sm[3][3], 
		      COMP_PRECISION *disp,COMP_PRECISION *x,
		      COMP_PRECISION f4,COMP_PRECISION f5,
		      COMP_PRECISION f6,COMP_PRECISION f7,
		      struct el_par  elastic)
{
  COMP_PRECISION tmp[4];
  tmp[0] = x[INT_Y]*f6;
  tmp[1] = x[INT_Y]*f7;
  tmp[2] = elastic.mu2*disp[STRIKE];
  tmp[3] = elastic.mu2*disp[NORMAL];
  sm[INT_X][INT_X] =tmp[2]*(2.0*f4 + tmp[0]);
  sm[INT_X][INT_X]-=tmp[3]*(   -f5 + tmp[1]);
  sm[INT_Y][INT_Y] =tmp[2]*(       - tmp[0]);
  sm[INT_Y][INT_Y]-=tmp[3]*(   -f5 - tmp[1]);
  sm[INT_X][INT_Y] =tmp[2]*(   -f5 + tmp[1]);
  sm[INT_X][INT_Y]-=tmp[3]*(       - tmp[0]);
  sm[INT_Y][INT_X] = sm[INT_X][INT_Y];
}
/*  

driver for the original Crouch and Starfield method
of calculating the stresses and strains, tdd_coeff
ihalf = 0 : full plane
ihalf = 1 : half plane


WE FLIP THE SIGN OF THE NORMAL MOTION TO MAKE THINGS
LIKE IN INTERACT

*/
void eval_2dsegment_plane_strain_tdd(COMP_PRECISION *x,
				     struct flt *fault,
				     COMP_PRECISION *disp,
				     COMP_PRECISION *u_global, 
				     COMP_PRECISION sm_global[3][3],
				     int ihalf,int *iret, MODE_TYPE mode,
				     struct el_par elastic)
{
  COMP_PRECISION pr0,pr,pr1,pr2,con,cons,chi;
  COMP_PRECISION sxx[3],sxy[3],syy[3],ux[3],uy[3];
  pr = (COMP_PRECISION)elastic.poisson;
  pr0 =  1.0 - pr;
  pr1 =  (1.0 - 2.0 * pr);
  pr2=   2.0*pr0;
  con=   1.0 / (4.0 * PI * pr0);
  cons = (COMP_PRECISION)elastic.mu2;
  chi =  3.0 - 4.0 * pr;

  if((ihalf)&&(x[INT_Y] > 0)){	/* half-plane in air */
    set_stress_and_disp_nan(sm_global,u_global,mode);
  }else{
    tdd_coeff((x+INT_X),(x+INT_Y),&fault->x[INT_X],&fault->x[INT_Y],
	      &fault->l,&fault->cos_alpha,&fault->sin_alpha,
	      &ihalf,&pr,&pr1,&pr2,&con,&cons,&chi,
	      (sxx+STRIKE),(sxx+NORMAL),(syy+STRIKE),(syy+NORMAL),
	      (sxy+STRIKE),(sxy+NORMAL),(ux+STRIKE), (ux+NORMAL),
	      (uy+STRIKE),(uy+NORMAL),iret);
    if(!(*iret)){
      if(mode != GC_STRESS_ONLY){
	/* add up contribution from strike and normal motion */
	u_global[INT_X] = ux[STRIKE] * disp[STRIKE] - 
	  ux[NORMAL] * disp[NORMAL];
	u_global[INT_Y] = uy[STRIKE] * disp[STRIKE] -
	  uy[NORMAL] * disp[NORMAL];
	u_global[INT_Z] = 0.0;
      }
      if(mode != GC_DISP_ONLY){
	/* stresses */
	sm_global[INT_X][INT_X] = sxx[STRIKE] * disp[STRIKE] -
	  sxx[NORMAL] * disp[NORMAL];
	sm_global[INT_Y][INT_X] = sm_global[INT_X][INT_Y] = 
	  sxy[STRIKE] * disp[STRIKE] -
	  sxy[NORMAL] * disp[NORMAL];
	sm_global[INT_Y][INT_Y] = syy[STRIKE] * disp[STRIKE] -
	  syy[NORMAL] * disp[NORMAL];
	sm_global[INT_X][INT_Z] = sm_global[INT_Z][INT_X] = sm_global[INT_Z][INT_Y] = 
	  sm_global[INT_Y][INT_Z] = 0.0;
	sm_global[INT_Z][INT_Z] = pr * (sm_global[INT_X][INT_X] + sm_global[INT_Y][INT_Y]);
      }
    }else{
      set_stress_and_disp_nan(sm_global,u_global,mode);
    }
  }
}
/*
  
  anti plane routines, from here on down based on Claude code, above TWB

  
  geometry and sign conventions are as for the in-plane routines: the
  segment is centered on fault->x, has half length fault->l and is
  rotated by alpha (fault->cos_alpha, fault->sin_alpha) with respect
  to the global x axis. all observational points have to lie in the
  z = 0 plane. z is the out of plane direction, y is up, and, for the
  half plane case, the free surface is at y = 0 with the medium at
  y <= 0, i.e. the same convention as for the ihalf switch of
  eval_2dsegment_plane_strain_tdd.

  the anti-plane problem is governed by

    \nabla^2 u_z = 0,  s_xz = \mu u_z,x,  s_yz = \mu u_z,y

  and all other displacement and stress components are zero. the
  solution is therefore independent of the Poisson ratio and there is
  no distinction between plane strain and plane stress, only
  elastic.shear enters.

  for a constant displacement discontinuity D over the element, the
  full plane solution is that of two screw dislocations of strength
  -D and +D at the two element tips,

    u_z = -D (\theta_1 - \theta_2)/(2 \pi)

  with \theta_1 = atan2(y_l, x_l - a), \theta_2 = atan2(y_l, x_l + a)
  in the element local frame and a = fault->l. as for the in-plane
  disp[STRIKE] component of Crouch and Starfield, the sign convention
  is

    D = u_z(y_l = 0^-) - u_z(y_l = 0^+)

  i.e. the same sense of motion as a positive disp[STRIKE], but out of
  the plane. the half plane solution follows from the image method,

    u_z^{hs}(x,y) = u_z^{fs}(x,y) + u_z^{fs}(x,-y),

  which is regular in y < 0 and gives s_yz = 0 on y = 0 exactly. note
  that the image term is evaluated by mirroring the observational
  point about the free surface, so the local gradient of the image
  contribution has to be rotated back into the global frame before the
  sign of its y component is flipped.

  the singularity handling, the iret flag, and the mode switches are
  as for the in-plane routines.

*/

/*

  out of plane shear (mode III) segment at arbitrary location and
  orientation

  ihalf = 0 : full plane
  ihalf = 1 : half plane, free surface at y = 0, medium at y <= 0

*/
void eval_2dsegment_antiplane(COMP_PRECISION *x,
			      struct flt *fault,
			      COMP_PRECISION *disp,
			      COMP_PRECISION *u_global, 
			      COMP_PRECISION sm_global[3][3],
			      int ihalf,int *iret, MODE_TYPE mode,
			      struct el_par elastic)
{
  COMP_PRECISION dx[3],grad[3],grad_img[3],uz,uz_img,mud;
#ifdef DEBUG
  if(x[INT_Z] != 0.0){
    fprintf(stderr,"eval_2dsegment_antiplane: z coordinate has to be zero: z: %g\n",
	    x[INT_Z]);
    exit(-1);
  }
  if(ihalf){
    if(fault->x[INT_Y] + fabs(fault->l * fault->sin_alpha) > EPS_COMP_PREC)
      fprintf(stderr,"eval_2dsegment_antiplane: WARNING: segment reaches above the free surface, y_c: %g l: %g\n",
	      fault->x[INT_Y],fault->l);
  }
#endif
  *iret = 0;
  if((ihalf)&&(x[INT_Y] > 0)){	/* half-plane in air */
    set_stress_and_disp_nan(sm_global,u_global,mode);
    return;
  }
  /* 
     source contribution, shift the observational point into the
     segment local frame within get_2dseg_antiplane_terms
  */
  dx[INT_X] = x[INT_X] - fault->x[INT_X];
  dx[INT_Y] = x[INT_Y] - fault->x[INT_Y];
  get_2dseg_antiplane_terms(dx,fault,&uz,grad,iret);
  if(*iret){
    /* 
       solution is infinite 
    */
    set_stress_and_disp_nan(sm_global,u_global,mode);
    return;
  }
  if(ihalf){
    /* 
       image contribution, mirror the observational point about the
       free surface 
    */
    dx[INT_X] =   x[INT_X] - fault->x[INT_X];
    dx[INT_Y] = - x[INT_Y] - fault->x[INT_Y];
    get_2dseg_antiplane_terms(dx,fault,&uz_img,grad_img,iret);
    if(*iret){
      set_stress_and_disp_nan(sm_global,u_global,mode);
      return;
    }
    /* 
       u^{hs}(x,y) = u(x,y) + u(x,-y), the chain rule flips the sign
       of the y derivative of the image term
    */
    uz  += uz_img;
    grad[INT_X] += grad_img[INT_X];
    grad[INT_Y] -= grad_img[INT_Y];
  }
  if(mode != GC_STRESS_ONLY){
    /* 
       only the out of plane displacement is non-zero, and it is a
       scalar, no rotation needed
    */
    u_global[INT_X] = u_global[INT_Y] = 0.0;
    u_global[INT_Z] = disp[STRIKE] * uz;
  }
  if(mode != GC_DISP_ONLY){
    /* 
       s_xz = \mu u_z,x and s_yz = \mu u_z,y, all other components
       vanish for anti-plane strain
    */
    mud = elastic.shear * disp[STRIKE];
    sm_global[INT_X][INT_X] = sm_global[INT_Y][INT_Y] = 
      sm_global[INT_Z][INT_Z] = 0.0;
    sm_global[INT_X][INT_Y] = sm_global[INT_Y][INT_X] = 0.0;
    sm_global[INT_X][INT_Z] = sm_global[INT_Z][INT_X] = mud * grad[INT_X];
    sm_global[INT_Y][INT_Z] = sm_global[INT_Z][INT_Y] = mud * grad[INT_Y];
  }
}
/*

  copy for a segment at the origin with strike angle 90, i.e. no shift
  and no rotation, following eval_2dsegment_plane_strain_basic

*/
void eval_2dsegment_antiplane_basic(COMP_PRECISION *x,
				    struct flt *fault,
				    COMP_PRECISION *disp,
				    COMP_PRECISION *u_global, 
				    COMP_PRECISION sm_global[3][3],
				    int ihalf,int *iret,
				    struct el_par elastic)
{
  COMP_PRECISION xm[3],f3,f4,f5,f3i,f4i,f5i,c1,c12,c2,c22,c3,c32,c4,c42,y2,mud,
    pfac1 = 2.0 * PI;
#ifdef DEBUG
  if(x[INT_Z] != 0.0){
    fprintf(stderr,"eval_2dsegment_antiplane_basic: z coordinate has to be zero: z: %g\n",
	    x[INT_Z]);exit(-1);
  }
  if(norm(fault->x,2) > EPS_COMP_PREC){
    fprintf(stderr,"eval_2dsegment_antiplane_basic: error, fault origin should be 0,0 but is %g, %g\n",
	    fault->x[INT_X],fault->x[INT_Y]);exit(-1);
  }
#endif
  if(fault->strike != 90.0){
    fprintf(stderr,"eval_2dsegment_antiplane_basic: error, fault strike should be 90 but is %g\n",
	    fault->strike);exit(-1);
  }
  *iret = 0;
  if((ihalf)&&(x[INT_Y] > 0)){	/* half-plane in air */
    set_stress_and_disp_nan(sm_global,u_global,GC_DISP_AND_STRESS);
    return;
  }
  get_2dseg_geo(x,fault->l,&c1,&c2,&c12,&c22,&c3,&c4,&c32,&c42,&y2,iret);
  if(*iret){
    set_stress_and_disp_nan(sm_global,u_global,GC_DISP_AND_STRESS);
    return;
  }
  get_2dseg_antiplane_ffac(&f3,&f4,&f5,c1,c2,c3,c4,pfac1,x);
  if(ihalf){
    /* 
       for a segment at the origin with alpha = 0, the local and
       global frames coincide and mirroring only affects y
    */
    xm[INT_X] =   x[INT_X];
    xm[INT_Y] = - x[INT_Y];
    xm[INT_Z] = 0.0;
    get_2dseg_geo(xm,fault->l,&c1,&c2,&c12,&c22,&c3,&c4,&c32,&c42,&y2,iret);
    if(*iret){
      set_stress_and_disp_nan(sm_global,u_global,GC_DISP_AND_STRESS);
      return;
    }
    get_2dseg_antiplane_ffac(&f3i,&f4i,&f5i,c1,c2,c3,c4,pfac1,xm);
    f3 += f3i;
    f4 += f4i;
    f5 -= f5i;			/* sign flip for the y derivative */
  }
  u_global[INT_X] = u_global[INT_Y] = 0.0;
  u_global[INT_Z] = disp[STRIKE] * f3;
  mud = elastic.shear * disp[STRIKE];
  sm_global[INT_X][INT_X] = sm_global[INT_Y][INT_Y] = 
    sm_global[INT_Z][INT_Z] = 0.0;
  sm_global[INT_X][INT_Y] = sm_global[INT_Y][INT_X] = 0.0;
  sm_global[INT_X][INT_Z] = sm_global[INT_Z][INT_X] =   mud * f4;
  sm_global[INT_Y][INT_Z] = sm_global[INT_Z][INT_Y] = - mud * f5;
}
/*

  displacement and displacement gradient of a single, full plane
  anti-plane segment, per unit displacement discontinuity

  input:  dx     offset of the observational point from the segment
                 center, in the global frame
          fault  segment
  output: uz     u_z / D
          grad   (u_z,x , u_z,y) / D, rotated back into the global frame
          iret   set if the point coincides with one of the tips

*/
void get_2dseg_antiplane_terms(COMP_PRECISION *dx,
			       struct flt *fault,
			       COMP_PRECISION *uz,
			       COMP_PRECISION *grad,
			       int *iret)
{
  COMP_PRECISION x_local[3],grad_local[3],f3,f4,f5,
    c1,c12,c2,c22,c3,c32,c4,c42,y2;
#ifdef DEBUG
  COMP_PRECISION l,w,corners[MAX_NR_EL_VERTICES*3];
#endif
  /*
    line up the segment with the local x-axis
  */
  rotate_vec2d(dx,x_local,fault->cos_alpha,fault->sin_alpha);
  /* 
     geometrical factors, same as for the in-plane solution 
  */
  get_2dseg_geo(x_local,fault->l,&c1,&c2,&c12,&c22,&c3,
		&c4,&c32,&c42,&y2,iret);
  if(*iret){
#ifdef DEBUG
    calculate_vertices(corners,fault,&l,&w);
    fprintf(stderr,
	    "get_2dseg_antiplane_terms: SINGULAR: x: (%g, %g) segment: (%g, %g) to (%g, %g)\n",
	    x_local[INT_X],x_local[INT_Y],corners[0*3+INT_X],corners[0*3+INT_Y],
	    corners[1*3+INT_X],corners[1*3+INT_Y]);
#endif
    return;
  }
  get_2dseg_antiplane_ffac(&f3,&f4,&f5,c1,c2,c3,c4,2.0*PI,x_local);
  *uz = f3;
  grad_local[INT_X] =   f4;
  grad_local[INT_Y] = - f5;
  /* 
     the gradient is a vector in the x-y plane and rotates back like
     the displacements do for the in-plane solution 
  */
  rotate_vec2d(grad_local,grad,fault->cos_alpha,-fault->sin_alpha);
}
/*

  the anti-plane F factors. these are the f3, f4, and f5 of
  get_2dseg_ffac, but with pfac1 = 2 pi instead of 4 pi (1-nu) since
  the Poisson ratio does not enter the anti-plane problem. with

    f3 = -(\theta_1 - \theta_2)/(2 \pi)

  we get u_z = D f3, u_z,x = D f4, and u_z,y = -D f5 in the segment
  local frame.

  the atan2 comment of get_2dseg_ffac applies here as well: for y = 0
  the difference of the two angles is 0 outside and +/- pi inside the
  segment, which is the intended displacement discontinuity.

*/
void get_2dseg_antiplane_ffac(COMP_PRECISION *f3, COMP_PRECISION *f4,
			      COMP_PRECISION *f5,
			      COMP_PRECISION c1, COMP_PRECISION c2,
			      COMP_PRECISION c3, COMP_PRECISION c4,
			      COMP_PRECISION pfac1,COMP_PRECISION *x)
{
  *f3 =    -(atan2(x[INT_Y],c1)  -  atan2(x[INT_Y],c2))/pfac1;
  *f4 =     (x[INT_Y]/c3         -         x[INT_Y]/c4)/pfac1;
  *f5 =     (        c1/c3       -               c2/c4)/pfac1;
}
