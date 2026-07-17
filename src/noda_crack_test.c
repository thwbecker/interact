#include "interact.h"
#include "properties.h"
/* 
   Noda (2025) style convergence test of the triangular multi-point
   receiver evaluations through interact's own Greens function
   machinery: Eshelby circular shear crack of unit radius and unit
   stress drop, discretized with a uniform equilateral mesh in a
   vertical plane, placed deep in the half space so that the analytic
   full space solution applies

   prescribes the analytic slip at the element centroids and compares
   the evaluated strike traction with the analytic value of -1 (in
   units of the stress drop)
*/
int main(int argc, char **argv)
{
  struct flt *fault;
  COMP_PRECISION dx,h,z0,x0,zp0,zp1,xoff,r,cen[2],disp[3],s[3],
    smax,dt,e[4],e09[4],u[3],sm[3][3],dummy[2];
  COMP_PRECISION *sval;
  my_boolean full_space = FALSE;
  int q = 3,nx,nz,i,j,k,n,nf,n09,iret;
  char *tname[4]={"-tv 0 (CTR)","-tv 1 (M244)","-tv 2 (M236)","-tv 3 (HYB)"};
  if(argc > 1)
    sscanf(argv[1],"%i",&q);
  dx = pow(2.0,-(COMP_PRECISION)q);
  h = dx*sqrt(3.0)/2.0;
  z0 = -200.0;			/* crack center depth, >> radius */
  nx = (int)ceil(1.2/dx);nz = (int)ceil(1.2/h);
  /* count and build mesh: vertical plane y=0, local coords (x, zp) */
  fault = (struct flt *)calloc(2*(2*nx+2)*(2*nz+2),sizeof(struct flt));
  nf = 0;
  for(j=-nz;j < nz;j++){
    zp0 = j*h;zp1 = (j+1)*h;
    xoff = 0.5*dx*(COMP_PRECISION)(abs(j)%2);
    for(i=-nx;i < nx;i++){
      x0 = i*dx - xoff;
      for(k=0;k < 2;k++){
	if(k==0){		/* upward triangle, ordered for +y normal */
	  fault[nf].xn = (COMP_PRECISION *)malloc(9*sizeof(COMP_PRECISION));
	  fault[nf].xn[0]=x0;       fault[nf].xn[1]=0;fault[nf].xn[2]=z0+zp0;
	  fault[nf].xn[3]=x0+0.5*dx;fault[nf].xn[4]=0;fault[nf].xn[5]=z0+zp1;
	  fault[nf].xn[6]=x0+dx;    fault[nf].xn[7]=0;fault[nf].xn[8]=z0+zp0;
	}else{			/* downward triangle */
	  fault[nf].xn = (COMP_PRECISION *)malloc(9*sizeof(COMP_PRECISION));
	  fault[nf].xn[0]=x0+dx;    fault[nf].xn[1]=0;fault[nf].xn[2]=z0+zp0;
	  fault[nf].xn[3]=x0+0.5*dx;fault[nf].xn[4]=0;fault[nf].xn[5]=z0+zp1;
	  fault[nf].xn[6]=x0+1.5*dx;fault[nf].xn[7]=0;fault[nf].xn[8]=z0+zp1;
	}
	cen[0] = (fault[nf].xn[0]+fault[nf].xn[3]+fault[nf].xn[6])/3.0;
	cen[1] = (fault[nf].xn[2]+fault[nf].xn[5]+fault[nf].xn[8])/3.0 - z0;
	r = sqrt(cen[0]*cen[0]+cen[1]*cen[1]);
	if(r < 1.0){
	  get_tri_prop_based_on_gh(fault+nf);
	  fault[nf].group = 0;
	  nf++;
	}else{
	  free(fault[nf].xn);
	}
      }
    }
  }
  /* analytic slip at centroids, scaled for unit stress drop */
  smax = 8.0/PI*(1.0-POISSON_NU)/(2.0-POISSON_NU)/SHEAR_MODULUS;
  sval = (COMP_PRECISION *)malloc(nf*sizeof(COMP_PRECISION));
  for(i=0;i < nf;i++){
    r = sqrt(fault[i].x[INT_X]*fault[i].x[INT_X] +
	     (fault[i].x[INT_Z]-z0)*(fault[i].x[INT_Z]-z0));
    sval[i] = smax*sqrt(1.0-r*r);
  }
  fprintf(stderr,"noda_crack_test: q=%i dx=1/%i nf=%i normal[y]=%g/%g t_strike[x]=%g/%g\n",
	  q,(int)(1.0/dx),nf,fault[0].normal[INT_Y],fault[1].normal[INT_Y],
	  fault[0].t_strike[INT_X],fault[1].t_strike[INT_X]);
  /* loop over evaluation schemes via the run-level flag */
  for(k=0;k < 4;k++){
    for(i=0;i < nf;i++)
      fault[i].type = TRIANGULAR + k;
    e[k] = e09[k] = 0.0;n09 = 0;
    for(i=0;i < nf;i++){	/* receivers */
      s[0]=s[1]=s[2]=0.0;
      for(j=0;j < nf;j++){	/* sources */
	disp[STRIKE] = sval[j];disp[DIP] = disp[NORMAL] = 0.0;
	eval_green_and_project_stress_to_fault(fault,i,j,disp,s,TRUE,full_space);
      }
      dt = fabs(s[STRIKE] + 1.0); /* analytic: -1 */
      e[k] += dt;
      r = sqrt(fault[i].x[INT_X]*fault[i].x[INT_X] +
	       (fault[i].x[INT_Z]-z0)*(fault[i].x[INT_Z]-z0));
      if(r < 0.9){
	e09[k] += dt;n09++;
      }
    }
    e[k] /= (COMP_PRECISION)nf;e09[k] /= (COMP_PRECISION)n09;
    fprintf(stdout,"q=%i N=%5i %16s  E = %10.5f  E0.9 = %10.5f\n",
	    q,nf,tname[k],e[k],e09[k]);
  }
  /* 
     verify -tv 4 (split scheme): operator assembly calls
     (multi_point_eval FALSE) must reproduce CTR, evaluation calls
     (multi_point_eval TRUE) must reproduce HYB
  */
  for(i=0;i < nf;i++)
    fault[i].type = TRIANGULAR + 4;
  for(k=0;k < 2;k++){		/* k=0: FALSE (operator), k=1: TRUE (evaluation) */
    dummy[k] = 0.0;
    for(i=0;i < nf;i++){
      s[0]=s[1]=s[2]=0.0;
      for(j=0;j < nf;j++){
	disp[STRIKE] = sval[j];disp[DIP] = disp[NORMAL] = 0.0;
	eval_green_and_project_stress_to_fault(fault,i,j,disp,s,(k)?(TRUE):(FALSE),full_space);
      }
      dummy[k] += fabs(s[STRIKE] + 1.0);
    }
    dummy[k] /= (COMP_PRECISION)nf;
  }
  fprintf(stdout,"-tv 4 operator path:    E = %10.5f  (must equal CTR E = %10.5f): %s\n",
	  dummy[0],e[0],(fabs(dummy[0]-e[0]) < 1e-12)?("PASS"):("FAIL"));
  fprintf(stdout,"-tv 4 evaluation path:  E = %10.5f  (must equal HYB E = %10.5f): %s\n",
	  dummy[1],e[3],(fabs(dummy[1]-e[3]) < 1e-12)?("PASS"):("FAIL"));
  return 0;
}
