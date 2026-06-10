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
  eval_green_at_receiver(fault,ireceive,islip,slip,u,sm,&iret,GC_STRESS_ONLY);
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


   the M236, M244, and HYBR evaluations are taken from 

   Noda, H., Inconsistency of a single-point evaluation of traction on
   a fault discretized with 2 triangular elements and several improved
   methods, GJI, 2025.
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
    /* 
       all triangular modes: single (collocation) point evaluation at
       the requested location x

       NOTE: the Noda (2025) multi-point evaluations (M236, M244,
       HYBR) are a property of the RECEIVER element, not of the
       source, and are therefore handled in eval_green_at_receiver()
       which is used for all fault-on-fault evaluations; the mixture
       has to be applied to the contributions of ALL sources (the
       artifact it cancels stems from the neighboring elements, cf.
       Noda's eqs. 9-11 where the focal element does not contribute),
       so it cannot be implemented here where only the source element
       is known
    */
  case TRIANGULAR:
  case TRIANGULAR_M236:
  case TRIANGULAR_M244:
  case TRIANGULAR_HYBR:{
    eval_triangle_general(x,fault,disp,u_global,sm_global,iret,mode,eval_on_itself);
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

#ifdef ALLOW_NON_3DQUAD_GEOM
void eval_triangle_general(COMP_PRECISION *x,struct flt *fault,
			   COMP_PRECISION *slip,COMP_PRECISION *u, 
			   COMP_PRECISION sm[3][3],int *giret,MODE_TYPE mode,
			   my_boolean is_self)
{
#ifdef INT_USE_TGF_TRIANGLE
  fprintf(stderr,"eval_triangle_tgf: WARNING: not tested well\n");
  eval_triangle_tgf(x,fault,slip,u,sm,giret,mode,is_self);
#else
  eval_triangle_nw(x,fault,slip,u,sm,giret,mode);
#endif
}
/* 
   Noda (2025, EPS) multi-point receiver evaluation machinery.

   get_noda_points() returns the number of evaluation points (4 for
   M236/M244, 7 for HYBR), their normalized locations eta (such that
   the evaluation point is sum_i xn_i/eta_i, cf. calc_tri_bary_coord),
   and the mixing weights w; returns 0 if the receiver type does not
   use multi-point evaluation
*/
static int get_noda_points(int rtype,COMP_PRECISION eta[7][3],COMP_PRECISION w[7])
{
  const COMP_PRECISION alpha = -4.1;	/* HYB mixing parameter, Noda eq. 38 */
  const COMP_PRECISION two_sqrt3 = 2.0*sqrt(3.0),
    fac = 1.0/(3.0 - 2.0*sqrt(3.0));
  static const COMP_PRECISION eta236[4][3]={{3.,3.,3.},{2.,3.,6.},{3.,6.,2.},{6.,2.,3.}};
  static const COMP_PRECISION eta244[4][3]={{3.,3.,3.},{2.,4.,4.},{4.,2.,4.},{4.,4.,2.}};
  int i,j,n;
  switch(rtype){
  case TRIANGULAR_M236:{	/* Noda eq. 31 */
    for(i=0;i<4;i++)
      for(j=0;j<3;j++)
	eta[i][j] = eta236[i][j];
    w[0] = 4.0;w[1] = w[2] = w[3] = -1.0;
    n = 4;
    break;
  }
  case TRIANGULAR_M244:{	/* Noda eq. 28 */
    for(i=0;i<4;i++)
      for(j=0;j<3;j++)
	eta[i][j] = eta244[i][j];
    w[0] = -two_sqrt3 * fac;w[1] = w[2] = w[3] = fac;
    n = 4;
    break;
  }
  case TRIANGULAR_HYBR:{	/* Noda eq. 38, sharing the centroid */
    for(j=0;j<3;j++){
      eta[0][j] = eta244[0][j];	          /* centroid */
      for(i=1;i<4;i++){
	eta[i  ][j] = eta244[i][j];
	eta[3+i][j] = eta236[i][j];
      }
    }
    w[0] = alpha * (-two_sqrt3 * fac) + (1.0-alpha) * 4.0;
    w[1] = w[2] = w[3] = alpha * fac;
    w[4] = w[5] = w[6] = (1.0-alpha) * (-1.0);
    n = 7;
    break;
  }
  default:
    n = 0;
    break;
  }
  return n;
}
#endif

/* 
   evaluate the Greens function (stress and, depending on mode,
   displacement) of source fault isrc at the receiver fault irec

   for regular patch types, this evaluates at the receiver centroid,
   as before; for the Noda (2025) triangular multi-point types (M236,
   M244, HYBR) of the RECEIVER, the source contribution is evaluated
   at the receiver's multiple interior points and mixed with the Noda
   weights - this has to be done for ALL sources (not only the
   self-interaction) for the artifact cancellation to work, since the
   element-scale oscillation stems from the neighboring elements

   all fault-on-fault evaluations should go through this function
*/
void eval_green_at_receiver(struct flt *fault,int irec,int isrc,
			    COMP_PRECISION *slip,COMP_PRECISION *u,
			    COMP_PRECISION sm_global[3][3],int *iret,
			    MODE_TYPE mode)
{
#ifdef ALLOW_NON_3DQUAD_GEOM
  COMP_PRECISION eta[7][3],w[7],xl[3],u_loc[3],sm_loc[3][3];
  int np,p,i;
  np = get_noda_points((int)fault[irec].type,eta,w);
  if(np){
    zero_3x3_matrix(sm_global);
    u[INT_X] = u[INT_Y] = u[INT_Z] = 0.0;
    for(p=0;p < np;p++){
      calc_tri_bary_coord(fault[irec].xn,xl,eta[p][0],eta[p][1],eta[p][2]);
      eval_green(xl,(fault+isrc),slip,u_loc,sm_loc,iret,mode,
		 (irec==isrc)?(TRUE):(FALSE));
      if(*iret)
	return;			/* singular evaluation */
      add_ay_to_3x3_matrix(sm_global,sm_loc,w[p]);
      for(i=0;i<3;i++)
	u[i] += w[p] * u_loc[i];
    }
  }else{
#endif
    eval_green(fault[irec].x,(fault+isrc),slip,u,sm_global,iret,mode,
	       (irec==isrc)?(TRUE):(FALSE));
#ifdef ALLOW_NON_3DQUAD_GEOM
  }
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
    fprintf(stderr,"eval_green_basic: fault patch type %i  undefined here (might work for eval_green)\n",
	    fault->type);
    exit(-1);
  }}
#else
  eval_okada(x,fault,disp,u_global,sm_global,iret,GC_DISP_AND_STRESS);
#endif
}
