/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: eval_2dsegment.c,v 1.13 2003-03-01 23:33:40-08 becker Exp tbecker $

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


  WARNING: FOR ALL CHANGES, REMEMBER TO FIX THE BASIC VERSIONS!

*/
#include "interact.h"
#include "properties.h"

//
// be generous with those
//
#ifdef USE_DOUBLE_PRECISION
#define EPS_COMP_PREC_FOR_NAN 5.0e-14
#else
#define EPS_COMP_PREC_FOR_NAN 5.0e-6
#endif

void eval_2dsegment_plane_strain(COMP_PRECISION *x,
				 struct flt *fault,
				 COMP_PRECISION *disp,
				 COMP_PRECISION *u_global, 
				 COMP_PRECISION sm_global[3][3],
				 int *iret)
{
  COMP_PRECISION u[3],x_local[3],dx[3],sm_local[3][3],f2,f3,f4,f5,f6,f7,
    c1,c12,c2,c22,c3,c32,c4,c42,y2;
  static COMP_PRECISION 
    pfac1 = 4.0 * 3.141592653589793 * (1.0 - POISSON_NU),
    pfac2 = (1.0 -       POISSON_NU), 
    pfac3 = (1.0 - 2.0 * POISSON_NU),
    pfac4 = POISSON_NU;
#ifdef DEBUG
  COMP_PRECISION l,w,corners[4][3];
  if(x[Z] != 0.0){
    fprintf(stderr,"eval_2dsegment_plane_strain: z coordinate has to be zero: z: %g\n",
	    x[Z]);
    exit(-1);
  }
#endif
  /* 
     shift and rotate observational point into local reference frame. 
     first, move fault to origin
  */
  dx[X] = x[X] - fault->x[X];
  dx[Y] = x[Y] - fault->x[Y];
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
    calculate_corners(corners,fault,&l,&w);
    fprintf(stderr,
	    "eval_2dsegment_plane_strain: SINGULAR: x: (%g, %g) segment: (%g, %g) to (%g, %g)\n",
	    x_local[X],x_local[Y],corners[0][X],corners[0][Y],
	    corners[1][X],corners[1][Y]);
#endif
    set_stress_and_disp_nan(sm_global,u_global);
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
    get_2dseg_stress(sm_local,disp,x_local,f4,f5,f6,f7);
    // no shear stresses with z for plain strain since exz and eyz = 0
    sm_local[X][Z] = sm_local[Z][X] = sm_local[Z][Y] = sm_local[Y][Z] = 0.0;
    // from plane strain: szz = \nu/ (sxx + syy)
    sm_local[Z][Z] = pfac4 * (sm_local[X][X] + sm_local[Y][Y]);
    /* 

       rotate displacements back into global frame this subroutine
       takes the first three components of u, ie. u_x u_y u_z

    */
    rotate_vec2d(u,u_global,fault->cos_alpha,-fault->sin_alpha);
    u_global[Z] = 0.0;
    /* 
       rotate stress matrix back into global field 
    */
    rotate_mat_z(sm_local,sm_global,fault->cos_alpha,-fault->sin_alpha);
  }
}
//
// repeat exercise for plane stress which we can obtain by
// letting  \nu --> \nu/(1+\nu) in the above formulation for plane strain
//
void eval_2dsegment_plane_stress(COMP_PRECISION *x,struct flt *fault,
				 COMP_PRECISION *disp,
				 COMP_PRECISION *u_global, 
				 COMP_PRECISION sm_global[3][3],int *iret)
{  
  COMP_PRECISION u[3],x_local[3],dx[3],sm_local[3][3],f2,f3,f4,f5,f6,f7,
    c1,c12,c2,c22,c3,c32,c4,c42,y2;
  static COMP_PRECISION // in pfac1 ... pfac4, \nu has been replaced by  \nu/(1+\nu)
    // the last factror, pfac5, uses the real poisson ratio for e_zz
    pfac1 = 4.0 * 3.141592653589793 * (1.0 - (POISSON_NU/(1.0+POISSON_NU))),
		       pfac2 = (1.0 -       (POISSON_NU/(1.0+POISSON_NU))), 
    pfac3 = (1.0 - 2.0 * (POISSON_NU/(1.0+POISSON_NU))),
    pfac5 = POISSON_NU/(TWO_TIMES_SHEAR_MODULUS*(1.0+POISSON_NU));
#ifdef DEBUG
  COMP_PRECISION l,w,corners[4][3];
  if(x[Z] != 0.0){
    fprintf(stderr,"eval_2dsegment_plane_stress: z coordinate has to be zero: z: %g\n",
	    x[Z]);
    exit(-1);
  }
#endif
  dx[X] = x[X] - fault->x[X];
  dx[Y] = x[Y] - fault->x[Y];
  rotate_vec2d(dx,x_local,fault->cos_alpha,fault->sin_alpha);
  get_2dseg_geo(x_local,fault->l,&c1, &c2,&c12, &c22,&c3, 
		&c4,&c32, &c42,&y2,iret);
  if(*iret){
#ifdef DEBUG
    calculate_corners(corners,fault,&l,&w);
    fprintf(stderr,
	    "eval_2dsegment_plane_stress: SINGULAR: x: (%g, %g) segment: (%g, %g) to (%g, %g)\n",
	    x_local[X],x_local[Y],corners[0][X],corners[0][Y],
	    corners[1][X],corners[1][Y]);
#endif
    set_stress_and_disp_nan(sm_global,u_global);
  }else{
    get_2dseg_ffac(&f2,&f3,&f4,&f5,&f6,&f7,c1,c2,c12,c22,
		   c3,c4,c32,c42,pfac1,y2,x_local);
    get_2dseg_disp(u,disp,x_local,f2,f3,f4,f5,pfac2,pfac3);
    get_2dseg_stress(sm_local,disp,x_local,f4,f5,f6,f7);
    // no shear stresses with z 
    sm_local[X][Z] = sm_local[Z][X] = sm_local[Z][Y] = sm_local[Y][Z] = 0.0;
    sm_local[Z][Z] = 0.0;
    rotate_vec2d(u,u_global,fault->cos_alpha,-fault->sin_alpha);
    rotate_mat_z(sm_local,sm_global,fault->cos_alpha,-fault->sin_alpha);
    // u_z is undefined without plate thickness, but we shall volunteer
    // e_zz = - pfac5 (s_xx + s_yy) instead
    u_global[Z] = -pfac5 * (sm_global[X][X] + sm_global[Y][Y]);
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
				       int *iret)
{
  COMP_PRECISION f2,f3,f4,f5,f6,f7,c1,c12,c2,c22,c3,c32,c4,
    c42,y2;
  static COMP_PRECISION 
    pfac1 = 4.0 * 3.141592653589793 * (1.0 - POISSON_NU),
    pfac2 = (1.0 -       POISSON_NU), 
    pfac3 = (1.0 - 2.0 * POISSON_NU),
    pfac4 = POISSON_NU;
#ifdef DEBUG
  COMP_PRECISION l,w,corners[4][3];
  if(x[Z] != 0.0){
    fprintf(stderr,"eval_2dsegment_plane_strain_basic: z coordinate has to be zero: z: %g\n",
	    x[Z]);exit(-1);
  }
  if(norm(fault->x,2) > EPS_COMP_PREC){
    fprintf(stderr,"eval_2dsegment_plane_strain_basic: error, fault origin should be 0,0 but is %g, %g\n",
	    fault->x[X],fault->x[Y]);exit(-1);
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
    calculate_corners(corners,fault,&l,&w);
    fprintf(stderr,
	    "eval_2dsegment_plane_strain_basic: SINGULAR: x: (%g, %g) segment: (%g, %g) to (%g, %g)\n",
	    x[X],x[Y],corners[0][X],corners[0][Y],
	    corners[1][X],corners[1][Y]);
#endif
    set_stress_and_disp_nan(sm_global,u_global);
  }else{
    get_2dseg_ffac(&f2,&f3,&f4,&f5,&f6,&f7,c1,c2,c12,c22,c3,c4,
		   c32,c42,pfac1,y2,x);
    get_2dseg_disp(u_global,disp,x,f2,f3,f4,f5,pfac2,pfac3);
    get_2dseg_stress(sm_global,disp,x,f4,f5,f6,f7);
    sm_global[X][Z] = sm_global[Z][X] = sm_global[Z][Y] = sm_global[Y][Z] = 0.0;
    sm_global[Z][Z] = pfac4 * (sm_global[X][X] + sm_global[Y][Y]);
    u_global[Z] = 0.0;
  }
}
void eval_2dsegment_plane_stress_basic(COMP_PRECISION *x,
				       struct flt *fault,
				       COMP_PRECISION *disp,
				       COMP_PRECISION *u_global, 
				       COMP_PRECISION sm_global[3][3],
				       int *iret)
{  
  COMP_PRECISION f2,f3,f4,f5,f6,f7,c1,c12,c2,c22,c3,c32,c4,c42,
    y2;
  static COMP_PRECISION // in pfac1 ... pfac4, \nu has been replaced by  \nu/(1+\nu)
    // the last factror, pfac5, uses the real poisson ratio for e_zz
    pfac1 = 4.0 * 3.141592653589793 * (1.0 - (POISSON_NU/(1.0+POISSON_NU))),
		       pfac2 = (1.0 -       (POISSON_NU/(1.0+POISSON_NU))), 
    pfac3 = (1.0 - 2.0 * (POISSON_NU/(1.0+POISSON_NU))),
    pfac5 = POISSON_NU/(TWO_TIMES_SHEAR_MODULUS*(1.0+POISSON_NU));
#ifdef DEBUG
  COMP_PRECISION l,w,corners[4][3];
  if(x[Z] != 0.0){
    fprintf(stderr,"eval_2dsegment: z coordinate has to be zero: z: %g\n",
	    x[Z]);exit(-1);
  }  
  if(norm(fault->x,2) > EPS_COMP_PREC){
    fprintf(stderr,"eval_2dsegment_plane_stress_basic: error, fault origin should be 0,0 but is %g, %g\n",
	    fault->x[X],fault->x[Y]);exit(-1);
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
    calculate_corners(corners,fault,&l,&w);
    fprintf(stderr,
	    "eval_2dsegment_plane_stress_basic: SINGULAR: x: (%g, %g) segment: (%g, %g) to (%g, %g)\n",
	    x[X],x[Y],corners[0][X],corners[0][Y],corners[1][X],
	    corners[1][Y]);
#endif
   set_stress_and_disp_nan(sm_global,u_global);
  }else{
    get_2dseg_ffac(&f2,&f3,&f4,&f5,&f6,&f7,c1,c2,c12,c22,
		   c3,c4,c32,c42,pfac1,y2,x);
    get_2dseg_disp(u_global,disp,x,f2,f3,f4,f5,pfac2,pfac3);
    get_2dseg_stress(sm_global,disp,x,f4,f5,f6,f7);
    sm_global[X][Z] = sm_global[Z][X] = sm_global[Z][Y] = sm_global[Y][Z] = 0.0;
    sm_global[Z][Z] = 0.0;
    u_global[Z] = -pfac5 * (sm_global[X][X] + sm_global[Y][Y]);
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
  *c1  = x[X] - l;// x - a
  *c2  = x[X] + l;// x + a
  *y2  = x[Y] * x[Y];
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
  *f3 =     -(atan2(x[Y],c1)       -       atan2(x[Y],c2))/pfac1;
  *f4 =      (x[Y]/c3              -              x[Y]/c4)/pfac1;
  *f5 =      (        c1/c3        -                c2/c4)/pfac1;
  *f6 =      ((c12 - y2)/c32       -         (c22-y2)/c42)/pfac1;
  *f7 =      2.0*x[Y]*(c1/c32-                     c2/c42)/pfac1;
}
		   
/*
  
assign part of the displacements, u[Z] depends on approximation

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
  tmp[2] = f5*x[Y];
  tmp[3] = f4*x[Y];
  u[X]  = disp[STRIKE] * (     tmp[0] - tmp[2]);
  u[X] -= disp[NORMAL] * (-    tmp[1] - tmp[3]);
  u[Y]  = disp[STRIKE] * (     tmp[1] - tmp[3]);
  u[Y] -= disp[NORMAL] * (     tmp[0] + tmp[2]);
}
/*
  
assign part of the stress matrix, rest depends on plane strain/stress
approximation

*/
void get_2dseg_stress(COMP_PRECISION sm[3][3], 
		      COMP_PRECISION *disp,COMP_PRECISION *x,
		      COMP_PRECISION f4,COMP_PRECISION f5,
		      COMP_PRECISION f6,COMP_PRECISION f7)
{
  COMP_PRECISION tmp[4];
  tmp[0] = x[Y]*f6;
  tmp[1] = x[Y]*f7;
  tmp[2] = TWO_TIMES_SHEAR_MODULUS*disp[STRIKE];
  tmp[3] = TWO_TIMES_SHEAR_MODULUS*disp[NORMAL];
  sm[X][X] =tmp[2]*(2.0*f4 + tmp[0]);
  sm[X][X]-=tmp[3]*(   -f5 + tmp[1]);
  sm[Y][Y] =tmp[2]*(       - tmp[0]);
  sm[Y][Y]-=tmp[3]*(   -f5 - tmp[1]);
  sm[X][Y] =tmp[2]*(   -f5 + tmp[1]);
  sm[X][Y]-=tmp[3]*(       - tmp[0]);
  sm[Y][X] = sm[X][Y];
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
				     int ihalf,int *iret)
{
  static my_boolean init=FALSE;
  static COMP_PRECISION pr,pr1,pr2,con,cons,chi;
  COMP_PRECISION sxx[3],sxy[3],syy[3],ux[3],uy[3];
  if(!init){		/* initialize parameters */
    pr = (COMP_PRECISION)POISSON_NU;
    pr1 = (1.0 - 2.0 * pr);
    pr2=  2.0*(1.0 - pr);
    con= 1.0 / (4.0 * PI * (1.0 - pr));
    cons = (COMP_PRECISION)YOUNG_MODULUS/(1.0 + pr);
    chi = 3.0 - 4.0 * pr;
    init = TRUE;
  }
  if((ihalf)&&(x[Y] > 0)){	/* half-plane in air */
    set_stress_and_disp_nan(sm_global,u_global);
  }else{
    tdd_coeff((x+X),(x+Y),&fault->x[X],&fault->x[Y],
	      &fault->l,&fault->cos_alpha,&fault->sin_alpha,
	      &ihalf,&pr,&pr1,&pr2,&con,&cons,&chi,
	      (sxx+STRIKE),(sxx+NORMAL),(syy+STRIKE),(syy+NORMAL),
	      (sxy+STRIKE),(sxy+NORMAL),(ux+STRIKE), (ux+NORMAL),
	      (uy+STRIKE),(uy+NORMAL),iret);
    if(!(*iret)){
      /* add up contribution from strike and normal motion */
      u_global[X] = ux[STRIKE] * disp[STRIKE] - 
	ux[NORMAL] * disp[NORMAL];
      u_global[Y] = uy[STRIKE] * disp[STRIKE] -
	uy[NORMAL] * disp[NORMAL];
      u_global[Z] = 0.0;
      /* stresses */
      sm_global[X][X] = sxx[STRIKE] * disp[STRIKE] -
	sxx[NORMAL] * disp[NORMAL];
      sm_global[Y][X] = sm_global[X][Y] = 
	sxy[STRIKE] * disp[STRIKE] -
	sxy[NORMAL] * disp[NORMAL];
      sm_global[Y][Y] = syy[STRIKE] * disp[STRIKE] -
	syy[NORMAL] * disp[NORMAL];
      sm_global[X][Z] = sm_global[Z][X] = sm_global[Z][Y] = 
	sm_global[Y][Z] = 0.0;
      sm_global[Z][Z] = pr * (sm_global[X][X] + sm_global[Y][Y]);
    }else{
      set_stress_and_disp_nan(sm_global,u_global);
    }
  }
}
