/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@post.harvard.edu


  this computes the displacements and stresses due an irregular quad,
  realized with three triangles

 2--------------1
 |\            /|
 | \          / |
 |  \   X    /  |
 |   \      /   |
 |    \    /    |
 | N1  \  / N2  |
 |      \/      |
 3 ---- 0 ----- 4
 

 
 input is the observational point x, the coordinates in xn, and the
 displacement in disp
 
 disp is in strike, dip, and normal format
 
 output is u_global (displacements) and sm_global
 (stress matrix)



*/
#include "interact.h"
#include "properties.h"

#ifdef ALLOW_NON_3DQUAD_GEOM
void eval_iquad(COMP_PRECISION *x,struct flt *fault,
		COMP_PRECISION *gslip,
		COMP_PRECISION *ug, 
		COMP_PRECISION smg[3][3],int *giret,
		MODE_TYPE mode)
{
  COMP_PRECISION sm[3][3],u[3],slip[3];
  struct flt afault;
  static my_boolean afault_init=FALSE;
  int j;
  if(!afault_init){
    afault.xn = (COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*9);
    afault_init = TRUE;
  }
  /* main triangle */
  eval_triangle(x,fault,gslip,u,sm,giret,mode);
  //#define SUPER_DUPER_DEBUG

#ifdef SUPER_DUPER_DEBUG
  fprintf(stderr,"eval_iquad: MT: xt %g %g %g\t%g %g %g\t%g %g %g\tslip %g %g %g\tx %g %g %g\tu %g %g %g\n", 
	  fault->xn[0],fault->xn[1],fault->xn[2], 
	  fault->xn[3],fault->xn[4],fault->xn[5], 
	  fault->xn[6],fault->xn[7],fault->xn[8], 
	  gslip[0],gslip[1],gslip[2],x[0],x[1],x[2],
	  u[0],u[1],u[2]);
  fprintf(stderr,"eval_triangle: u: %g %g %g\n",u[0],u[1],u[2]);
#endif

  
  ug[0]=u[0];  ug[1]=u[1];  ug[2]=u[2];
  smg[0][0]=sm[0][0];  smg[0][1]=sm[0][1];  smg[0][2]=sm[0][2];
  smg[1][0]=sm[1][0];  smg[1][1]=sm[1][1];  smg[1][2]=sm[1][2];
  smg[2][0]=sm[2][0];  smg[2][1]=sm[2][1];  smg[2][2]=sm[2][2];

  /* first aux triangle */
  slip[0]=slip[1]=slip[2]=0.0;
  for(j=0;j < 3;j++){
    afault.xn[0*3+j] = fault->xn[3*3+j]; /* nodes */
    afault.xn[1*3+j] = fault->xn[0*3+j];
    afault.xn[2*3+j] = fault->xn[2*3+j];
    slip[0] += gslip[j] * fault->xn[5*3+j]; /* projection */
    slip[1] += gslip[j] * fault->xn[6*3+j];
    slip[2] += gslip[j] * fault->xn[7*3+j];
  }
  eval_triangle(x,&afault,slip,u,sm,giret,mode);
  ug[0] += u[0];  ug[1] += u[1];  ug[2] += u[2];
  smg[0][0] += sm[0][0];  smg[0][1] += sm[0][1];  smg[0][2] += sm[0][2];
  smg[1][0] += sm[1][0];  smg[1][1] += sm[1][1];  smg[1][2] += sm[1][2];
  smg[2][0] += sm[2][0];  smg[2][1] += sm[2][1];  smg[2][2] += sm[2][2];
 
  /* second aux triangle */
  slip[0]=slip[1]=slip[2]=0.0;
  for(j=0;j < 3;j++){
    afault.xn[0*3+j] = fault->xn[0*3+j]; /* nodes */
    afault.xn[1*3+j] = fault->xn[4*3+j];
    afault.xn[2*3+j] = fault->xn[1*3+j];
    slip[0] += gslip[j] * fault->xn[ 8*3+j]; /* projection */
    slip[1] += gslip[j] * fault->xn[ 9*3+j];
    slip[2] += gslip[j] * fault->xn[10*3+j];
  }
  eval_triangle(x,&afault,slip,u,sm,giret,mode);
  ug[0] += u[0];  ug[1] += u[1];  ug[2] += u[2];
  smg[0][0] += sm[0][0];  smg[0][1] += sm[0][1];  smg[0][2] += sm[0][2];
  smg[1][0] += sm[1][0];  smg[1][1] += sm[1][1];  smg[1][2] += sm[1][2];
  smg[2][0] += sm[2][0];  smg[2][1] += sm[2][1];  smg[2][2] += sm[2][2];
}


/* end of triangle geometry allowed branch */
#endif
