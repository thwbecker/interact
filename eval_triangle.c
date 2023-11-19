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
		   COMP_PRECISION sm[3][3],int *giret)
{
 
  COMP_PRECISION stress[6],strain[6];
  const COMP_PRECISION nu =  POISSON_NU, lambda = LAMBDA,mu = SHEAR_MODULUS;

  /* input is x y z, vertices, and slip as strike, dip, normal
     displacements are output as east, north, up
  */
  tddisphs(x,  &(fault->xt[0]),&(fault->xt[3]),&(fault->xt[6]),(slip),(slip+1),(slip+2),&nu,u);
  /* 
     stress and strain are given as: Sxx, Syy, Szz, Sxy, Sxz and Syz
  */
  tdstresshs(x,&(fault->xt[0]),&(fault->xt[3]),&(fault->xt[6]),
	     (slip),(slip+1),(slip+2),&mu,&lambda,
	     stress,strain);

  sm[INT_X][INT_X] = stress[0];
  sm[INT_Y][INT_Y] = stress[1];
  sm[INT_Z][INT_Z] = stress[2];
  sm[INT_X][INT_Y] = stress[3];
  sm[INT_X][INT_Z] = stress[4];
  sm[INT_Y][INT_Z] = stress[5];

  sm[INT_Y][INT_X]=sm[INT_X][INT_Y];
  sm[INT_Z][INT_X]=sm[INT_X][INT_Z];
  sm[INT_Z][INT_Y]=sm[INT_Y][INT_Z];
  *giret = 0;
}
#endif
