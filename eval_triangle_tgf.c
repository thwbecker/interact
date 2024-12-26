/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu

  
  program to evaluate the displacement and stresses due to a
  triangular element using the f90 code of

  Gimbutas, Zydrunas, et al. "On the Calculation of Displacement,
  Stress, and Strain Induced by Triangular Dislocations." Bulletin of
  the Seismological Society of America 102.6 (2012): 2776-2780.
  (see tgf/COPYRIGHT and readme)

  input is the observational point x, the triangular
  coordinates xt, and the displacement in slip

  slip is in strike, dip, and normal format

  output is u_global (displacements) and sm_global
  (stress matrix)

  ________________________________________________________________________________
  THIS HAS NOT BEEN TESTED YET
  ________________________________________________________________________________

*/
#include "interact.h"
#include "properties.h"

#ifdef ALLOW_NON_3DQUAD_GEOM
#ifdef USE_TGF_TRIANGLE
void eval_triangle_tgf(COMP_PRECISION *x,struct flt *fault,
		       COMP_PRECISION *slip,COMP_PRECISION *u, 
		       COMP_PRECISION sm[3][3],int *giret,int mode, my_boolean is_self)
{
 
  COMP_PRECISION etracel;
  const double mu = SHEAR_MODULUS,lambda = LAMBDA_CONST,two_mu = 2.0*(SHEAR_MODULUS);
  double fstrain[9],strain0[9],u0[3],triangle[9],*dnormal,*dslip,*du,*dx;
  int numfunev,unity=1;
  int i,j,i3,j3;
  int calc_strain;
#ifdef USE_DOUBLE_PRECISION
  dslip = slip;
  dnormal = fault->normal;
  dx=x;
  du=u;
#else
  dnormal=(double *)malloc(sizeof(double)*3);
  dslip=(double *)malloc(sizeof(double)*3);
  du=(double *)malloc(sizeof(double)*3);
  dx=(double *)malloc(sizeof(double)*3);
  for(i=0;i<3;i++){
    dx[i] = (double)x[i];
    du[i] = (double)u[i];
    dslip[i] = (double)slip[i];
    dnormal[i] = (double)fault->normal[i];
  }
#endif

  
  if(mode== GC_DISP_ONLY)
    calc_strain = 0;
  else
    calc_strain = 1;
  
  for(i=i3=0;i<3;i++,i3+=3){
    for(j=j3=0;j<3;j++,j3+=3)
      triangle[j3+i] = fault->xn[i3+j]; /* convert to fortran convention */
  }
  
  if(is_self)
    eltst3triadirectself(&lambda, &mu, &unity, &unity, triangle, dslip, dnormal, dx, du, &calc_strain, fstrain);
  else
    eltst3triadirecttarg(&lambda, &mu, &unity, triangle, dslip, dnormal, dx, du, &calc_strain, fstrain);
  elth3triaadap(         &lambda, &mu, &unity, triangle, dslip, dnormal, dx, &unity, u0, &calc_strain, strain0, &numfunev);

  for(i=i3=0;i < 3;i++,i3+=3){
    u[i] += u0[i];
    if(calc_strain)
      for(j=0;j < 3;j++)
	fstrain[i3+j] += strain0[i3+j];
  }
  if(calc_strain){
    etracel = (fstrain[0*3+0]+fstrain[1*3+1]+fstrain[2*3+2]) * lambda;
    sm[0][0] = two_mu * fstrain[0*3+0] + etracel;
    sm[1][1] = two_mu * fstrain[1*3+1] + etracel;
    sm[2][2] = two_mu * fstrain[2*3+2] + etracel;
    /* flip the ordering here from fortran */
    sm[0][1]=sm[1][0] = two_mu * fstrain[1*3+0];
    sm[0][2]=sm[2][0] = two_mu * fstrain[2*3+0];
    sm[1][2]=sm[2][1] = two_mu * fstrain[2*3+1];
  }
  
  *giret = 0;
#ifndef USE_DOUBLE_PRECISION
  free(dnormal);free(dslip);free(dx);free(du);
#endif
}


#endif
#endif
