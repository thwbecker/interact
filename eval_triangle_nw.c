/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu

  
  program to evaluate the displacement and stresses due to a
  triangular element using the F90 conversion of 

  Nikkhoo, M., Walter, T. R. (2015): Triangular dislocation: an
  analytical, artefact-free solution. - Geophysical Journal
  International, 201, 1117-1139. doi: 10.1093/gji/ggv035

  input is the observational point x, the triangular
  coordinates xt, and the displacement in slip

  slip is in strike, dip, and normal format

  output is u_global (displacements) and sm_global
  (stress matrix)

*/
#include "interact.h"
#include "properties.h"

#ifdef ALLOW_NON_3DQUAD_GEOM
void eval_triangle_nw(COMP_PRECISION *x,struct flt *fault,
		      COMP_PRECISION *slip,COMP_PRECISION *u, 
		      COMP_PRECISION sm[3][3],int *giret,
		      MODE_TYPE mode)
{
 
  COMP_PRECISION stress[6];
#ifdef USE_HBI_TDDEF
  double mu = SHEAR_MODULUS,lambda = LAMBDA_CONST;
#else
  COMP_PRECISION strain[6];
#endif
    
  if(mode != GC_STRESS_ONLY){
    /* input is x y z, vertices, and slip as strike, dip, normal
       displacements are output as east, north, up
    */
    tddisphs(x,  &(fault->xn[0]),&(fault->xn[3]),&(fault->xn[6]),(slip),(slip+1),(slip+2),u);
    //tddisphs_bird((x+0),(x+1),(x+2),&(fault->xn[0]),&(fault->xn[3]),&(fault->xn[6]),(slip),(slip+1),(slip+2),&nu,(u+0),(u+1),(u+2));
#ifdef CRAZY_DEBUG
    fprintf(stderr,"eval_triangle_nw: xt %g %g %g\t%g %g %g\t%g %g %g\tslip %g %g %g\tx %g %g %g\tu %g %g %g\n", 
	    fault->xn[0],fault->xn[1],fault->xn[2], 
	    fault->xn[3],fault->xn[4],fault->xn[5], 
	    fault->xn[6],fault->xn[7],fault->xn[8], 
	    slip[0],slip[1],slip[2],x[0],x[1],x[2],
	    u[0],u[1],u[2]);
    fprintf(stderr,"eval_triangle_nw: u: %g %g %g\n",u[0],u[1],u[2]);
#endif
    if(!finite(u[0])||(!finite(u[1]))||(!finite(u[2]))){
      set_stress_and_disp_nan(sm,u,mode);
      *giret = 1;
      return;
    }
    /* don't even access u if not desired, might have gotten called as
     * NULL */
  }
  if(mode !=  GC_DISP_ONLY){
    /* 
       stress and strain are given as: Sxx, Syy, Szz, Sxy, Sxz and Syz
    */
#ifdef USE_HBI_TDDEF
    hbi_tdstresshs_(x,(x+1),(x+2),&(fault->xn[0]),&(fault->xn[3]),&(fault->xn[6]),slip,(slip+1),(slip+2),&mu,&lambda,
		    (stress),(stress+1),(stress+2),(stress+3),(stress+4),(stress+5));
#else
    tdstresshs(x,&(fault->xn[0]),&(fault->xn[3]),&(fault->xn[6]),
	       (slip),(slip+1),(slip+2),stress,strain);
#ifdef CRAZY_DEBUG
    fprintf(stderr,"eval_triangle_nw: stress: %g %g %g %g %g %g\n",strain[0],strain[1],strain[2],strain[3],strain[4],strain[5]);
#endif
#endif
    if((!finite(stress[0]))||(!finite(stress[1]))||(!finite(stress[2]))||
       (!finite(stress[3]))||(!finite(stress[4]))||(!finite(stress[5]))){
      set_stress_and_disp_nan(sm,u,mode);
      *giret = 1;
      return;
    }
    sm[INT_X][INT_X] = stress[0];
    sm[INT_Y][INT_Y] = stress[1];
    sm[INT_Z][INT_Z] = stress[2];
    sm[INT_X][INT_Y] = stress[3];
    sm[INT_X][INT_Z] = stress[4];
    sm[INT_Y][INT_Z] = stress[5];
    
    sm[INT_Y][INT_X]=sm[INT_X][INT_Y];
    sm[INT_Z][INT_X]=sm[INT_X][INT_Z];
    sm[INT_Z][INT_Y]=sm[INT_Y][INT_Z];
  }else{
    set_stress_and_disp_nan(sm,u,mode);
  }
  *giret = 0;
}


/*
  
  compute local coordinate system, strike and dip given initialized
  triangular points xt 

  moved here such that geometry.c can compile without the f90 routines

*/
void get_tri_prop_based_on_gh(struct flt *fault)
{
  COMP_PRECISION dip;
  double strike,alpha;
  calc_centroid_tri(fault->xn,fault->x); /* assign centroid to x */
  /* F90 routine */
  get_tdcs_base_vectors((fault->xn+0),(fault->xn+3),(fault->xn+6),
			fault->t_strike,fault->t_dip,fault->normal,
			&fault->area);
  fault->l = fault->w = sqrt(fault->area);

  strike = atan2(fault->t_strike[INT_X],fault->t_strike[INT_Y])*RAD2DEG;  
  /*  */
  dip = (COMP_PRECISION)asin((double)fault->t_dip[INT_Z])*RAD2DEG;
  
  check_angles(&dip,&strike);
  /*  */
  fault->strike = (COMP_PRECISION)strike;
  fault->dip    = dip;
  /*  */
  alpha = 90 - fault->strike;
  my_sincos_degd(&(fault->sin_alpha),&(fault->cos_alpha),alpha);

}
/* end of triangle geometry allowed branch */
#endif
