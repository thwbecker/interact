/*

  unit test for the visco-elastic Prony kernel machinery
  (src/ve/prony_kernel.c); see rsf_ve_design.md

  reads a geometry from geom.in and, for a set of source-receiver
  pairs, checks that

  (1) the sum of all amplitude terms (constant plus exponentials)
      reproduces the elastic kernel at rounding level (the t = 0,
      instantaneous limit)

  (2) the held-out Laplace sample residual is at rounding level,
      confirming the pole set is complete for the kernels in use

  (3) for the stress kernel, the constant (relaxed) amplitude is at
      rounding level relative to the elastic kernel (a homogeneous
      Maxwell medium relaxes fully)

  (4) for the displacement kernel at surface observation points, the
      tau_m amplitude is at rounding level (G cancels in
      displacements) and the constant term matches a direct
      evaluation at a very small Laplace variable (the permanent,
      relaxed deformation)

  exit 0 with a PASS line if all checks hold, exit -1 otherwise

  usage: ve_check [tau_m, 1] [tolerance, 1e-9]

  part of interact

*/
#include "interact.h"
#include "properties.h"
#include "prony_kernel.h"

int main(int argc, char **argv)
{
  struct flt *fault;
  struct med *medium;
  struct prony_spec spec;
  struct el_par ep_small;
  COMP_PRECISION C[VE_MAX_NP+1][3][3],D[VE_MAX_NP+1][3];
  COMP_PRECISION slip[3],u[3],sm[3][3],x[3],up_small[3],up2[2][3];
  COMP_PRECISION tau_m,tol,res,scl,dev,worst[4];
  COMP_PRECISION s_small;
  int i,j,k,it,ipair,iret,isrc,irec,nrflt;
  int src_list[3],rec_list[3],nsrc,mode;
  my_boolean ok;
  tau_m = 1.0;
  tol = 1e-9;
  if(argc > 1)sscanf(argv[1],ONE_CP_FORMAT,&tau_m);
  if(argc > 2)sscanf(argv[2],ONE_CP_FORMAT,&tol);
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("ve_check");
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25);
  nrflt = medium->nrflt;
  ve_spec_homogeneous(&spec,medium->elastic.shear,
		      medium->elastic.poisson,tau_m);
  fprintf(stderr,"ve_check: %i patches, np %i (tau %g %g %g) + const, ns %i, tol %g\n",
	  nrflt,spec.np,spec.tau[0],spec.tau[1],spec.tau[2],spec.ns,tol);
  nsrc = (nrflt >= 3)?(3):(1);
  src_list[0] = 0;src_list[1] = nrflt/2;src_list[2] = nrflt-1;
  rec_list[0] = nrflt-1;rec_list[1] = nrflt/3;rec_list[2] = 0;
  for(i=0;i < 4;i++)worst[i] = 0.0;
  ok = TRUE;
  for(mode=STRIKE;mode <= NORMAL;mode++){
    slip[STRIKE] = slip[DIP] = slip[NORMAL] = 0.0;
    slip[mode] = 1.0;
    for(ipair=0;ipair < nsrc;ipair++){
      isrc = src_list[ipair];irec = rec_list[ipair];
      /* stress kernel amplitudes and held-out residual (check 2) */
      res = ve_prony_amplitudes_stress(&spec,medium,fault,irec,isrc,
				       slip,C);
      if(res > worst[1])worst[1] = res;
      /* elastic reference (check 1 and scale for check 3) */
      eval_green_at_receiver(fault,irec,isrc,slip,u,sm,&iret,
			     GC_STRESS_ONLY,TRUE,medium->full_space,
			     medium->elastic);
      scl = 0.0;
      for(i=0;i < 3;i++)
	for(j=0;j < 3;j++)
	  if(fabs(sm[i][j]) > scl)scl = fabs(sm[i][j]);
      if(scl < EPS_COMP_PREC)continue;
      for(i=0;i < 3;i++)
	for(j=0;j < 3;j++){
	  dev = 0.0;
	  for(it=0;it < spec.nterm;it++)
	    dev += C[it][i][j];
	  dev = fabs(dev - sm[i][j])/scl;
	  if(dev > worst[0])worst[0] = dev;
	  /* relaxed stress amplitude should vanish (check 3) */
	  dev = fabs(C[spec.np][i][j])/scl;
	  if(dev > worst[2])worst[2] = dev;
	}
    }
  }
  /* displacement checks at a surface observation point (check 4);
     strike slip on the middle patch */
  slip[STRIKE] = 1.0;slip[DIP] = slip[NORMAL] = 0.0;
  isrc = nrflt/2;
  x[INT_X] = fault[isrc].x[INT_X] + 2.0*fault[isrc].l;
  x[INT_Y] = 0.0;
  x[INT_Z] = 0.0;
  res = ve_prony_amplitudes_disp(&spec,medium,fault,x,isrc,slip,D);
  if(res > worst[1])worst[1] = res;
  scl = 0.0;
  for(i=0;i < 3;i++)
    for(it=0;it < spec.nterm;it++)
      if(fabs(D[it][i]) > scl)scl = fabs(D[it][i]);
  if(scl > EPS_COMP_PREC){
    /* tau_m amplitude of displacements should vanish */
    for(i=0;i < 3;i++){
      dev = fabs(D[0][i])/scl;
      if(dev > worst[3])worst[3] = dev;
    }
    /* constant term vs direct evaluation extrapolated to s -> 0
       (two small-s evaluations, linear Richardson, since the
       small-s value itself is O(s) off the limit) */
    for(k=0;k < 2;k++){
      s_small = ((k == 0)?(1e-6):(1e-7))/tau_m;
      ve_effective_elpar(spec.g0,spec.nu0,tau_m,s_small,&ep_small);
      eval_green(x,(fault+isrc),slip,up2[k],sm,&iret,
		 GC_DISP_ONLY,FALSE,medium->full_space,ep_small);
    }
    for(i=0;i < 3;i++){
      /* extrapolate: U0 = (s1 U(s2) - s2 U(s1))/(s1 - s2) with
	 s1 = 1e-6, s2 = 1e-7 */
      up_small[i] = (1e-6*up2[1][i] - 1e-7*up2[0][i])/(1e-6 - 1e-7);
      dev = fabs(D[spec.np][i] - up_small[i])/scl;
      if(dev > worst[3])worst[3] = dev;
    }
  }
  fprintf(stderr,"ve_check: worst elastic-sum dev %8.1e, held-out resid %8.1e, relaxed-stress amp %8.1e, disp checks %8.1e\n",
	  worst[0],worst[1],worst[2],worst[3]);
  for(i=0;i < 4;i++)
    if(!(worst[i] < tol))
      ok = FALSE;
  if(ok){
    fprintf(stderr,"ve_check: PASS (tolerance %g)\n",tol);
    return 0;
  }else{
    fprintf(stderr,"ve_check: FAIL (tolerance %g)\n",tol);
    return -1;
  }
}
