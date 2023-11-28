/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu


  wraps C calls to the Okada (1992) halfspace routines for Greens
  functions given constant slip on a rectangle or for a point source


WARNING: FOR ALL CHANGES: REMEMBER TO ALSO CHANGE THE BASIC VERSION
BELOW!


*/
#include "interact.h"
#include "properties.h"
// offsets for the indices in the Okada routines
#define OKUX 0
#define OKUY 1
#define OKUZ 2
#define OKUXX 3
#define OKUYX 4
#define OKUZX 5
#define OKUXY 6
#define OKUYY 7
#define OKUZY 8
#define OKUXZ 9
#define OKUYZ 10
#define OKUZZ 11
extern void dc3d(double*,double*,double*,double*,double*,double*,double*,double*,
		 double*,double*,double*,double*,double*,double*,double*,double*,
		 double*,double*,double*,double*,double*,double*,double*,double*,
		 double*,int*);
extern void dc3d0(double*,double*,double*,double*,double*,double*,double*,double*,
		  double*,double*,double*,double*,double*,double*,double*,double*,
		  double*,double*,double*,double*,double*,double*,int*);

/*
  
  evaluate stresses and displacements due to slip on a 

  RECTANGULAR FAULT PATCH 

  using the Okada code 

  input is the fault geometry (mid-point, half-length and width) 
  in the fault structure, the observational point x, and the 
  applied slip[3] in the disp vector [STRIKE, DIP, NORMAL]
  
  output is the displacement in u_global and the stress 
  matrix [6] in sm_global 


  (also see eval_rectangle_basic below, which only works for
  strike = 90)

*/
void eval_rectangle(COMP_PRECISION *x,struct flt *fault,
		    COMP_PRECISION *disp,
		    COMP_PRECISION *u_global, 
		    COMP_PRECISION sm_global[3][3],int *iret,
		    MODE_TYPE mode)
{

  COMP_PRECISION iso,dx[3],sm_local[3][3];
#ifndef USE_DOUBLE_PRECISION
  int i;
  double x_local_d[3],u_d[12],disp_d[3];
#endif
  COMP_PRECISION x_local[3],u[12];
  static double medium_alpha = ALPHA_CONST;
  double al1,al2,aw1,aw2,depth,cpdip;
#ifdef DEBUG
  COMP_PRECISION corners[4][3],l,w;
#endif
  cpdip=(double)fault->dip;
#ifdef DEBUG
  if(cpdip < 0 || cpdip > 90){
    fprintf(stderr,"eval_rectangle: dip (%g) should be between 0 and 90 (maybe)\n",
	    cpdip);
    exit(-1);
  }
  if(fault->x[INT_Z] >0){
    fprintf(stderr,"eval_rectangle: error: fault z: %g (needs to be <=0)\n",
	    fault->x[INT_Z]);
    exit(-1);
  }
#endif
  /* 
     shift and rotate observational point 
     into local reference frame. first, move fault to origin
  */
  dx[INT_X]=x[INT_X] - fault->x[INT_X];
  dx[INT_Y]=x[INT_Y] - fault->x[INT_Y];
  dx[INT_Z]=x[INT_Z];
  // then line up fault strike with the x-axis
  rotate_vec(dx,x_local,fault->cos_alpha,fault->sin_alpha);
  /* depth has to be positive */
  depth= (double)-fault->x[INT_Z];
  // half length
  al1 = (double)-fault->l;al2 = (double)fault->l;
  // half width
  aw1 = (double)-fault->w;aw2 = (double)fault->w;
#ifdef USE_DOUBLE_PRECISION
  // call to Okada routine
  dc3d(&medium_alpha,(x_local+INT_X),(x_local+INT_Y),(x_local+INT_Z),
       &depth,&cpdip,&al1,&al2,&aw1,&aw2,
       (disp+STRIKE),(disp+DIP),(disp+NORMAL),
       (u+OKUX),(u+OKUY),(u+OKUZ),(u+OKUXX),(u+OKUYX),(u+OKUZX),
       (u+OKUXY),(u+OKUYY),(u+OKUZY),(u+OKUXZ),(u+OKUYZ),(u+OKUZZ),iret);
#else  /* single prec */
  for(i=0;i<3;i++){
    disp_d[i] = (double)disp[i];
    x_local_d[i] = (double)x_local[i];
  }

  dc3d(&medium_alpha,(x_local_d+INT_X),(x_local_d+INT_Y),(x_local_d+INT_Z),
       &depth,&cpdip,&al1,&al2,&aw1,&aw2,
       (disp_d+STRIKE),(disp_d+DIP),(disp_d+NORMAL),
       (u_d+OKUX),(u_d+OKUY),(u_d+OKUZ),(u_d+OKUXX),(u_d+OKUYX),(u_d+OKUZX),
       (u_d+OKUXY),(u_d+OKUYY),(u_d+OKUZY),(u_d+OKUXZ),(u_d+OKUYZ),(u_d+OKUZZ),
       iret);
  for(i=0;i < 12;i++){
    u[i] = (COMP_PRECISION)u_d[i];
  }
#endif
  if(*iret){
#ifdef DEBUG
    calculate_corners(corners,fault,&l,&w);
    fprintf(stderr,
	    "eval_rectangle: SINGULAR: x: (%g, %g, %g) fault corners: (%g, %g, %g) (%g, %g, %g) (%g, %g, %g) (%g, %g, %g)\n",
	    x_local[INT_X],x_local[INT_Y],x_local[INT_Z],
	    corners[0][INT_X],corners[0][INT_Y],corners[0][INT_Z],
	    corners[1][INT_X],corners[1][INT_Y],corners[1][INT_Z],
	    corners[2][INT_X],corners[2][INT_Y],corners[2][INT_Z],
	    corners[3][INT_X],corners[3][INT_Y],corners[3][INT_Z]);
#endif
    set_stress_and_disp_nan(sm_global,u_global,mode);
  }else{
    if(mode != GC_STRESS_ONLY){	/* may not be allocated */
      /* 
	 rotate displacements back into global frame 
	 this subroutine takes the first three components 
	 of u, ie. u_x u_y u_z
      */
      rotate_vec(u,u_global,fault->cos_alpha,-fault->sin_alpha);
    }
    if(mode != GC_DISP_ONLY){
      /* 
	 convert displacement derivatives into 
	 strains and stresses 
      */
      /* this is the isotropic component */
      iso= LAMBDA_CONST*(u[OKUXX]+u[OKUYY]+u[OKUZZ]);
      /* 
	 assign sxx,sxy,sxz,syy,syz,szz 
	 according to 
	 s_ij = \lambda \sum_i e_ii \delta_ij + 2 \mu e_ij
      */
      sm_local[INT_X][INT_X]=               iso + TWO_TIMES_SHEAR_MODULUS*u[OKUXX];
      sm_local[INT_X][INT_Y]=sm_local[INT_Y][INT_X]=SHEAR_MODULUS*(u[OKUXY]+u[OKUYX]);
      sm_local[INT_X][INT_Z]=sm_local[INT_Z][INT_X]=SHEAR_MODULUS*(u[OKUXZ]+u[OKUZX]);
      sm_local[INT_Y][INT_Y]=               iso + TWO_TIMES_SHEAR_MODULUS*u[OKUYY];
      sm_local[INT_Y][INT_Z]=sm_local[INT_Z][INT_Y]=SHEAR_MODULUS*(u[OKUYZ]+u[OKUZY]);
      sm_local[INT_Z][INT_Z]=               iso+TWO_TIMES_SHEAR_MODULUS*u[OKUZZ];
      /* 
	 rotate stress matrix back into global field 
      */
      rotate_mat_z(sm_local,sm_global,fault->cos_alpha,
		   -fault->sin_alpha);
    }
  }
}
/*

this is the same as above, but for fault coordinates 
x: x,y = 0,0 and strike = 90

input is L, W, dip, and depth (>0) as opposed to fault structure

will compute both displacement and stress

*/
void eval_rectangle_basic(COMP_PRECISION *x,
			  COMP_PRECISION l, COMP_PRECISION w,
			  COMP_PRECISION dip,
			  COMP_PRECISION depth,
			  COMP_PRECISION *disp,
			  COMP_PRECISION *u_global, 
			  COMP_PRECISION sm_global[3][3],
			  int *iret)
{
#ifndef USE_DOUBLE_PRECISION
  double depth_d,x_d[3],disp_d[3],u_d[12],dip_d,u_global_d[3];
  int i;
#endif
  double medium_alpha = ALPHA_CONST;
  COMP_PRECISION u[12],iso;
  double al1,al2,aw1,aw2;
  al1 = (double)-l;al2 = (double)l;
  aw1 = (double)-w;aw2 = (double)w;
  //#ifdef DEBUG
  if((depth < 0)||(x[INT_Z]>0)){
    fprintf(stderr,"eval_rectangle_basic: error: depth: %g (has to be >0) z: %g (has to be <0)\n",
	    depth,x[INT_Z]);
    exit(-1);
  }
  //#endif
#ifdef USE_DOUBLE_PRECISION
  dc3d(&medium_alpha,(x+INT_X),(x+INT_Y),(x+INT_Z),&depth,&dip,
       &al1,&al2,&aw1,&aw2,(disp+STRIKE),(disp+DIP),
       (disp+NORMAL),(u_global+INT_X),(u_global+INT_Y),(u_global+INT_Z),
       (u+OKUXX),(u+OKUYX),(u+OKUZX),(u+OKUXY),(u+OKUYY),
       (u+OKUZY),(u+OKUXZ),(u+OKUYZ),(u+OKUZZ),iret);
#else
  for(i=0;i<3;i++){
    x_d[i] = (double)x[i];
    disp_d[i] = (double)disp[i];
  }
  depth_d = (double)depth;
  dip_d = (double)dip;
  dc3d(&medium_alpha,(x_d+INT_X),(x_d+INT_Y),(x_d+INT_Z),&depth_d,&dip_d,
       &al1,&al2,&aw1,&aw2,(disp_d+STRIKE),(disp_d+DIP),
       (disp_d+NORMAL),(u_global_d+INT_X),(u_global_d+INT_Y),(u_global_d+INT_Z),
       (u_d+OKUXX),(u_d+OKUYX),(u_d+OKUZX),(u_d+OKUXY),(u_d+OKUYY),
       (u_d+OKUZY),(u_d+OKUXZ),(u_d+OKUYZ),(u_d+OKUZZ),iret);
  for(i=0;i<3;i++)
    u_global[i] = (COMP_PRECISION)u_global_d[i];
  for(i=0;i<12;i++)
    u[i] = (COMP_PRECISION)u_d[i];
#endif
  if(*iret){
    set_stress_and_disp_nan(sm_global,u_global,GC_DISP_AND_STRESS);
  }else{
    iso= LAMBDA_CONST*(u[OKUXX]+u[OKUYY]+u[OKUZZ]);
    sm_global[INT_X][INT_X]=               iso + TWO_TIMES_SHEAR_MODULUS*u[OKUXX];
    sm_global[INT_X][INT_Y]=sm_global[INT_Y][INT_X]=SHEAR_MODULUS*(u[OKUXY]+u[OKUYX]);
    sm_global[INT_X][INT_Z]=sm_global[INT_Z][INT_X]=SHEAR_MODULUS*(u[OKUXZ]+u[OKUZX]);
    sm_global[INT_Y][INT_Y]=               iso + TWO_TIMES_SHEAR_MODULUS*u[OKUYY];
    sm_global[INT_Y][INT_Z]=sm_global[INT_Z][INT_Y]=SHEAR_MODULUS*(u[OKUYZ]+u[OKUZY]);
    sm_global[INT_Z][INT_Z]=               iso + TWO_TIMES_SHEAR_MODULUS*u[OKUZZ];
  }
}


/* 
   evaluate the displacements and stresses due to a 

   POINT SOURCE

   at the location as specified in the fault structure 
   using the Okada code 

   input is the fault geometry (mid-point, half-length and half-width)
   in the fault structure, the observational point x

   output is the displacement in u_global and the stress 
   matrix [6] in sm_global */


void eval_point(COMP_PRECISION *x,struct flt *fault,
		COMP_PRECISION *disp,COMP_PRECISION *u_global, 
		COMP_PRECISION sm_global[3][3],int *iret,
		MODE_TYPE mode)
{
  eval_point_short(x,fault->x,fault->area,fault->sin_alpha,
		   fault->cos_alpha,(COMP_PRECISION)fault->dip,
		   disp,u_global,sm_global,iret,mode);
}
/*

  slightly modified from above, will not take fault array
  but xf[3] and area for point source location and fault `area'
  sin(alpha), cos(alpha), dip
  
*/
void eval_point_short(COMP_PRECISION *x,COMP_PRECISION *xf,COMP_PRECISION area,
		      COMP_PRECISION sin_alpha,COMP_PRECISION cos_alpha,
		      COMP_PRECISION dip,COMP_PRECISION *disp,
		      COMP_PRECISION *u_global, 
		      COMP_PRECISION sm_global[3][3],
		      int *iret,MODE_TYPE mode)
{
  double medium_alpha = ALPHA_CONST;
  COMP_PRECISION mu_alpha = SHEAR_MODULUS/LAMBDA_CONST;
  COMP_PRECISION u[12],
    iso,x_local[3],dx[3],sm_local[3][3];
  double depth,potency[4];
#ifndef USE_DOUBLE_PRECISION
  double u_d[12],x_local_d[3],dip_d;
  int i;
#endif
  dx[INT_X]=x[INT_X] - xf[INT_X];
  dx[INT_Y]=x[INT_Y] - xf[INT_Y];
  dx[INT_Z]=x[INT_Z];
  rotate_vec(dx,x_local,cos_alpha,sin_alpha);
  depth= (double)-xf[INT_Z];
  potency[0]=potency[1]=potency[2]=potency[3]=0.0;
  // potency is normally moment/myu, ie. area
  if(disp[STRIKE] != 0.)
    potency[0]=(double)area * disp[STRIKE];
  if(disp[DIP] != 0.)
    potency[1]=(double)area * disp[DIP];
  if(disp[NORMAL] != 0.)
    potency[2]=(double)area * disp[NORMAL]* mu_alpha;
#ifdef USE_DOUBLE_PRECISION
  dc3d0(&medium_alpha,(x_local+INT_X),(x_local+INT_Y),(x_local+INT_Z),
	&depth,&dip,
	(potency+0),(potency+1),(potency+2),(potency+3),
	(u+OKUX),(u+OKUY),(u+OKUZ),
	(u+OKUXX),(u+OKUYX),(u+OKUZX),
	(u+OKUXY),(u+OKUYY),(u+OKUZY),
	(u+OKUXZ),(u+OKUYZ),(u+OKUZZ),iret);
#else
  dip_d = (double)dip;
  for(i=0;i<3;i++)
    x_local_d[i] = (double)x_local[i];
  dc3d0(&medium_alpha,(x_local_d+INT_X),(x_local_d+INT_Y),(x_local_d+INT_Z),
	&depth,&dip_d,
	(potency+0),(potency+1),(potency+2),(potency+3),
	(u_d+OKUX),(u_d+OKUY),(u_d+OKUZ),
	(u_d+OKUXX),(u_d+OKUYX),(u_d+OKUZX),
	(u_d+OKUXY),(u_d+OKUYY),(u_d+OKUZY),
	(u_d+OKUXZ),(u_d+OKUYZ),(u_d+OKUZZ),iret);
  for(i=0;i<12;i++)
    u[i] = (COMP_PRECISION)u_d[i];
#endif
  if(*iret){
    set_stress_and_disp_nan(sm_global,u_global,mode);
  }else{
    if(mode !=  GC_STRESS_ONLY)
      rotate_vec(u,u_global,cos_alpha,-sin_alpha);
    if(mode != GC_DISP_ONLY){
      iso= LAMBDA_CONST*(u[OKUXX]+u[OKUYY]+u[OKUZZ]);

      sm_local[INT_X][INT_X]=iso+TWO_TIMES_SHEAR_MODULUS*u[OKUXX];
      sm_local[INT_X][INT_Y]=sm_local[INT_Y][INT_X]=
	SHEAR_MODULUS*(u[OKUXY]+u[OKUYX]);
      sm_local[INT_X][INT_Z]=sm_local[INT_Z][INT_X]=
	SHEAR_MODULUS*(u[OKUXZ]+u[OKUZX]);
      sm_local[INT_Y][INT_Y]=iso+TWO_TIMES_SHEAR_MODULUS*u[OKUYY];
      sm_local[INT_Y][INT_Z]=sm_local[INT_Z][INT_Y]=
	SHEAR_MODULUS*(u[OKUYZ]+u[OKUZY]);
      sm_local[INT_Z][INT_Z]=iso+TWO_TIMES_SHEAR_MODULUS*u[OKUZZ];
      rotate_mat_z(sm_local,sm_global,cos_alpha,-sin_alpha);
    }
  }
}
#undef OKUX 
#undef OKUY 
#undef OKUZ 

#undef OKUXX 
#undef OKUYX 
#undef OKUZX 

#undef OKUXY 
#undef OKUYY 
#undef OKUZY 

#undef OKUXZ 
#undef OKUYZ 
#undef OKUZZ 



void set_stress_and_disp_nan(COMP_PRECISION s[3][3],COMP_PRECISION *u,
			     MODE_TYPE mode)
{
  if(mode != GC_STRESS_ONLY){
    u[INT_X]=u[INT_Y]=u[INT_Z]=NAN;
  }
  if(mode != GC_DISP_ONLY){
    s[INT_X][INT_X]=s[INT_X][INT_Y]=s[INT_X][INT_Z]=NAN;
    s[INT_Y][INT_X]=s[INT_Y][INT_Y]=s[INT_Y][INT_Z]=NAN;
    s[INT_Z][INT_X]=s[INT_Z][INT_Y]=s[INT_Z][INT_Z]=NAN;
  }
}
