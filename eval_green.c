/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: eval_green.c,v 2.17 2003-03-01 23:33:59-08 becker Exp tbecker $


  routine for selecting appropriate Greens functions routines given
  the patch geometry

*/
#include "interact.h"
#include "properties.h"

/* 

   evaluate the displacements and stresses based on various half-space
   and 2-D Greens functions

   will call point, rectangle or triangle routines from here

*/


void eval_green(COMP_PRECISION *x,struct flt *fault,
		COMP_PRECISION *disp,COMP_PRECISION *u_global, 
		COMP_PRECISION sm_global[3][3],int *iret)
{
#ifdef ALLOW_NON_3DQUAD_GEOM
  switch(fault->type){
  case POINT_SOURCE:{
    eval_point(x,fault,disp,u_global,sm_global,iret);
    break;
  }
  case TRIANGULAR:{
    eval_triangle(x,fault,disp,u_global,sm_global,iret);
    break;
  }
  case RECTANGULAR_PATCH:{
    eval_rectangle(x,fault,disp,u_global,sm_global,iret);
    break;
  }
  case TWO_DIM_SEGMENT_PLANE_STRAIN:{
    eval_2dsegment_plane_strain(x,fault,disp,u_global,sm_global,iret);
    //eval_2dsegment_plane_strain_tdd(x,fault,disp,u_global,sm_global, 
    //0,iret); 
    break;
  }
  case TWO_DIM_SEGMENT_PLANE_STRESS:{
    eval_2dsegment_plane_stress(x,fault,disp,u_global,sm_global,iret);
    break;
  }
  case TWO_DIM_HALFPLANE_PLANE_STRAIN:{
    eval_2dsegment_plane_strain_tdd(x,fault,disp,u_global,sm_global, 
     				    1,iret); 
    break;
  }
  default:{
    fprintf(stderr,"eval_green: fault patch type %i undefined\n",
	    fault->type);
    exit(-1);
  }}
#else
  eval_rectangle(x,fault,disp,u_global,sm_global,iret);
#endif
}
/* 
   same as above but assume that the faults coordinates
   are in the origin (x: 0, y: 0) and the strike is 0 

*/
void eval_green_basic(COMP_PRECISION *x,struct flt *fault,
		      COMP_PRECISION *disp,
		      COMP_PRECISION *u_global, 
		      COMP_PRECISION sm_global[3][3],int *iret)
{
#ifdef DEBUG
  if((fault->strike != 90)||(norm(fault->x,2) > EPS_COMP_PREC)){
    fprintf(stderr,"eval_green_basic: fault should have strike=90 and be at origin\n");
    fprintf(stderr,"eval_green_basic: strike: %g x,y: %g,%g\n",
	    fault->strike,fault->x[INT_X],fault->x[INT_Y]);
    exit(-1);
  }
#endif
#ifdef ALLOW_NON_3DQUAD_GEOM
  switch(fault->type){
  case RECTANGULAR_PATCH:{
    eval_rectangle_basic(x,
			 (COMP_PRECISION)fault->l,
			 (COMP_PRECISION)fault->w,
			 (COMP_PRECISION)fault->dip,
			 (COMP_PRECISION) -fault->x[INT_Z],
			 disp,u_global,sm_global,iret);
    break;
  }
  case TWO_DIM_SEGMENT_PLANE_STRAIN:{
    eval_2dsegment_plane_strain_basic(x,fault,disp,
				      u_global,sm_global,iret);
    break;
  }
  case TWO_DIM_SEGMENT_PLANE_STRESS:{
    eval_2dsegment_plane_stress_basic(x,fault,disp,u_global,
				      sm_global,iret);
    break;
  }
  default:{
    fprintf(stderr,"eval_green_basic: fault patch type %i here undefined\n",
	    fault->type);
    exit(-1);
  }}
#else
  eval_rectangle(x,fault,disp,u_global,sm_global,iret);
#endif
}
