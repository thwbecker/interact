/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu


  routine for selecting appropriate Greens functions routines given
  the patch geometry

*/
#include "interact.h"
#include "properties.h"

/* 
   evaluate the greens function and project to fault stress 
   ireceive is the receiving fault
   slip the slip on islip fault
   s[] is normally fault[ireceive].s[]
   
*/

void eval_green_and_project_stress_to_fault(struct flt *fault, int ireceive,
					    int islip, COMP_PRECISION *slip,
					    COMP_PRECISION *s)
{
  //#define SUPER_DUPER_DEBUG
  int iret;
  COMP_PRECISION trac[3],sm[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},u[3]={0.,0.,0.};
  /* evaluate the greens function for slipping fault islip at receiver
     location centroids of fault ireceive */
#ifdef SUPER_DUPER_DEBUG
  fprintf(stderr,"add_quake_stress_3: slip: %10.3e %10.3e %10.3e evaluated at %10.3e %10.3e %10.3e \n",
	  slip[0],slip[1],slip[2],fault[ireceive].x[0],fault[ireceive].x[1],fault[ireceive].x[2]);
#endif
  eval_green(fault[ireceive].x,(fault+islip),slip,u,sm,&iret,GC_STRESS_ONLY,(ireceive==islip)?(TRUE):(FALSE));
  if(!iret){
    /* project the stresses */
    resolve_force(fault[ireceive].normal,sm,trac); /* convert to local
						      traction
						      vector */
    s[STRIKE]  += dotp_3d(trac,fault[ireceive].t_strike);
    s[DIP]     += dotp_3d(trac,fault[ireceive].t_dip);
    s[NORMAL]  += dotp_3d(trac,fault[ireceive].normal);
#ifdef SUPER_DUPER_DEBUG
    fprintf(stderr,"add_quake_stress_3: rec %03i rup %03i t %10.3e %10.3e %10.3e with sv (%9.2e,%9.2e,%9.2e) %9.2e (t:%9.2e) dv (%9.2e,%9.2e,%9.2e) %9.2e (t:%9.2e) nv (%9.2e,%9.2e,%9.2e) %9.2e (t:%9.2e)\n",ireceive,islip,
	    trac[0],trac[1],trac[2],
	    fault[ireceive].t_strike[0],fault[ireceive].t_strike[1],fault[ireceive].t_strike[2],
	    dotp_3d(trac,fault[ireceive].t_strike),s[STRIKE],
	    fault[ireceive].t_dip[0],fault[ireceive].t_dip[1],fault[ireceive].t_dip[2],
	    dotp_3d(trac,fault[ireceive].t_dip),s[DIP],
	    fault[ireceive].normal[0],fault[ireceive].normal[1],fault[ireceive].normal[2],
	    dotp_3d(trac,fault[ireceive].normal),s[NORMAL]);
#endif
  }
}

/* 

   evaluate the displacements and stresses based on various
   half-space, triangular, and 2-D Greens functions

   will call point, rectangle or triangle routines from here

   mode type can be GC_DISP_AND_STRES, GC_DISP_ONLY, GC_STRESS_ONLY
   depending on asking for displacement and stress, displacement, or
   stress only evaluation
   
   on_element is true when called from 
   eval_green_and_project_stress_to_fault to evaluate on
   centroid, for example
*/


void eval_green(COMP_PRECISION *x,struct flt *fault,
		COMP_PRECISION *disp,COMP_PRECISION *u_global, 
		COMP_PRECISION sm_global[3][3],int *iret,
		MODE_TYPE mode,my_boolean eval_on_itself)
/* mode type only for triangular */
{
#ifdef ALLOW_NON_3DQUAD_GEOM
  switch(fault->type){
  case POINT_SOURCE:{
    eval_point(x,fault,disp,u_global,sm_global,iret,mode);
    break;
  }
  case TRIANGULAR:{
#ifdef USE_TGF_TRIANGLE
    eval_triangle_tgf(x,fault,disp,u_global,sm_global,iret,mode,eval_on_itself);
#else
    eval_triangle_nw(x,fault,disp,u_global,sm_global,iret,mode);
#endif
    break;
  }
  case IQUAD:{
    eval_iquad(x,fault,disp,u_global,sm_global,iret,mode);
    break;
  }
  case OKADA_PATCH:{
    eval_okada(x,fault,disp,u_global,sm_global,iret,mode);
    break;
  }
  case TWO_DIM_SEGMENT_PLANE_STRAIN:{
    eval_2dsegment_plane_strain(x,fault,disp,u_global,sm_global,iret,mode);
    break;
  }
  case TWO_DIM_SEGMENT_PLANE_STRESS:{
    eval_2dsegment_plane_stress(x,fault,disp,u_global,sm_global,iret,mode);
    break;
  }
  case TWO_DIM_HALFPLANE_PLANE_STRAIN:{
    eval_2dsegment_plane_strain_tdd(x,fault,disp,u_global,sm_global, 1,iret,
				    mode); 
    break;
  }
  default:{
    fprintf(stderr,"eval_green: fault patch type %i undefined\n",
	    fault->type);
    exit(-1);
  }}
#else
  eval_okada(x,fault,disp,u_global,sm_global,iret,mode);
#endif
}
/* 
   same as above but assume that the faults coordinates
   are in the origin (x: 0, y: 0) and the strike is 0 

   will compute both stress and diaplacement

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
  case OKADA_PATCH:{
    eval_okada_basic(x,
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
  eval_okada(x,fault,disp,u_global,sm_global,iret,GC_DISP_AND_STRESS);
#endif
}
