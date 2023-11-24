/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@post.harvard.edu


  program to evaluate the displacement and stresses due to a
  triangular element using the F90 conversion of 

  Nikkhoo, M., Walter, T. R. (2015): Triangular dislocation: an
  analytical, artefact-free solution. - Geophysical Journal
  International, 201, 1117-1139. doi: 10.1093/gji/ggv035

  input is the observational point x, the triangular
  coordinates xt, and the displacement in disp

  disp is in strike, dip, and normal format

  output is u_global (displacements) and sm_global
  (stress matrix)

*/
#include "interact.h"
#include "properties.h"

#ifdef ALLOW_NON_3DQUAD_GEOM
void eval_triangle(COMP_PRECISION *x,struct flt *fault,
		   COMP_PRECISION *slip,COMP_PRECISION *u, 
		   COMP_PRECISION sm[3][3],int *giret,
		   MODE_TYPE mode)
{
 
  COMP_PRECISION stress[6],strain[6];
  const COMP_PRECISION nu =  POISSON_NU;

  if(mode != GC_STRESS_ONLY){
    /* input is x y z, vertices, and slip as strike, dip, normal
       displacements are output as east, north, up
    */
    tddisphs(x,  &(fault->xt[0]),&(fault->xt[3]),&(fault->xt[6]),(slip),(slip+1),(slip+2),&nu,u);
    //tddisphs_bird((x+0),(x+1),(x+2),&(fault->xt[0]),&(fault->xt[3]),&(fault->xt[6]),(slip),(slip+1),(slip+2),&nu,(u+0),(u+1),(u+2));
#ifdef CRAZY_DEBUG
    fprintf(stderr,"eval_triangle: xt %g %g %g\t%g %g %g\t%g %g %g\tslip %g %g %g\tx %g %g %g\tu %g %g %g\n", 
	    fault->xt[0],fault->xt[1],fault->xt[2], 
	    fault->xt[3],fault->xt[4],fault->xt[5], 
	    fault->xt[6],fault->xt[7],fault->xt[8], 
	    slip[0],slip[1],slip[2],x[0],x[1],x[2],
	    u[0],u[1],u[2]);
    fprintf(stderr,"eval_triangle: u: %g %g %g\n",u[0],u[1],u[2]);
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
    tdstresshs(x,&(fault->xt[0]),&(fault->xt[3]),&(fault->xt[6]),
	       (slip),(slip+1),(slip+2),stress,strain);
#ifdef CRAZY_DEBUG
    fprintf(stderr,"eval_triangle: stress: %g %g %g %g %g %g\n",strain[0],strain[1],strain[2],strain[3],strain[4],strain[5]);
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
#endif
