/*

  step 7a validation driver: antiplane elastic plate over a Maxwell
  substrate, general rigidity contrast, exercised through the
  dip = -90 antiplane elements (eval_2dsegment_antiplane) via the
  image-series machinery of prony_kernel.c

  the geometry read from geom.in must be antiplane elements forming a
  vertical fault strictly inside the plate; the driver uses the
  middle element as the source and surface observation points

  checks:

  (1) elastic anchor: the s -> inf sample (Gamma -> Gamma_0) against
      the Rybicki (1971) elastic layered arctan series
  (2) relaxed anchor: the s -> 0 sample (Gamma -> 1) against the
      free-bottom (fully relaxed) Rybicki series
  (3) elastic layered cross-check at intermediate Gamma values: the
      numeric image machinery against the analytic series, two
      independent implementations of the same solution
  (4) THE step 7a measurement: Prony term count. for np = 2..6, fit
      the sampled s-domain kernel and compare the reconstructed time
      history against the independent Nur and Mavko (1974)
      time-domain series, reporting both the held-out s-domain
      residual and the worst time-domain deviation relative to the
      coseismic amplitude

  usage: ve_layered_check [plate_h, 2] [t_M, 1] [g2_over_g1, 1] [n_img, 60]

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
  COMP_PRECISION slip[3],u[3],sm[3][3],x[3];
  COMP_PRECISION plate_h,t_M,g2fac,gamma0,rate_b,xs,c1,c2,uana,gam,s;
  COMP_PRECISION smp[VE_MAX_NS],amp[VE_MAX_NP+1],tv[8];
  COMP_PRECISION dev,worst,res,scl,model,t,coseis;
  int n_img,np,i,k,it,isrc,nt,iobs;
  plate_h = 2.0;t_M = 1.0;g2fac = 1.0;n_img = 60;
  if(argc > 1)sscanf(argv[1],ONE_CP_FORMAT,&plate_h);
  if(argc > 2)sscanf(argv[2],ONE_CP_FORMAT,&t_M);
  if(argc > 3)sscanf(argv[3],ONE_CP_FORMAT,&g2fac);
  if(argc > 4)sscanf(argv[4],"%i",&n_img);
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("ve_layered_check");
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25);
  ve_layered_gamma_pars(medium->elastic.shear,
			g2fac*medium->elastic.shear,t_M,
			&gamma0,&rate_b);
  slip[STRIKE] = 1.0;slip[DIP] = slip[NORMAL] = 0.0;
  isrc = medium->nrflt/2;
  /* segment depth interval and horizontal position, from the element
     (vertical: strike 0, so the along-segment direction is y) */
  xs = fault[isrc].x[INT_X];
  c1 = -(fault[isrc].x[INT_Y] + fault[isrc].l);
  c2 = -(fault[isrc].x[INT_Y] - fault[isrc].l);
  if((c1 < 0)||(c2 > plate_h)){
    fprintf(stderr,"ve_layered_check: source element [%g, %g] outside plate [0, %g]\n",
	    c1,c2,plate_h);
    exit(-1);
  }
  fprintf(stderr,"ve_layered_check: H %g t_M %g g2/g1 %g (Gamma0 %.4f, b %.4f), source depths [%g, %g], n_img %i\n",
	  plate_h,t_M,g2fac,gamma0,rate_b,c1,c2,n_img);
  x[INT_Y] = 0.0;x[INT_Z] = 0.0;
  worst = 0.0;
  for(iobs=0;iobs < 3;iobs++){
    x[INT_X] = xs + 0.75 + 1.25*(COMP_PRECISION)iobs;
    /* (1) elastic anchor, Gamma(s -> inf) = Gamma_0 */
    ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
				gamma0,n_img,plate_h,GC_DISP_ONLY,u,sm);
    uana = ve_rybicki_antiplane_uz(x[INT_X],xs,c1,c2,plate_h,gamma0,n_img);
    dev = fabs(u[INT_Z]-uana)/fabs(uana);
    if(dev > worst)worst = dev;
    if(iobs == 0)
      fprintf(stderr,"elastic anchor (Gamma = Gamma0):        rel dev %8.1e\n",dev);
    /* (2) relaxed anchor, Gamma = 1 */
    ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
				1.0,n_img,plate_h,GC_DISP_ONLY,u,sm);
    uana = ve_rybicki_antiplane_uz(x[INT_X],xs,c1,c2,plate_h,1.0,n_img);
    dev = fabs(u[INT_Z]-uana)/fabs(uana);
    if(dev > worst)worst = dev;
    if(iobs == 0)
      fprintf(stderr,"relaxed anchor (Gamma = 1):             rel dev %8.1e\n",dev);
    /* (3) intermediate Gamma cross-checks */
    for(k=1;k <= 3;k++){
      gam = 0.25*(COMP_PRECISION)k;
      ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
				  gam,n_img,plate_h,GC_DISP_ONLY,u,sm);
      uana = ve_rybicki_antiplane_uz(x[INT_X],xs,c1,c2,plate_h,gam,n_img);
      dev = fabs(u[INT_Z]-uana)/fabs(uana);
      if(dev > worst)worst = dev;
    }
  }
  fprintf(stderr,"machinery vs Rybicki, all anchors + Gamma sweep, 3 obs: worst rel dev %8.1e\n",
	  worst);
  if(worst > 1e-8){
    fprintf(stderr,"ve_layered_check: FAIL (elastic image machinery)\n");
    return -1;
  }
  /* 
     (4) Prony term count: fit the s-domain samples per obs point,
     reconstruct the time history, compare against Nur-Mavko; times
     span the layered relaxation clock 1/b 
  */
  tv[0] = 0.2/rate_b;tv[1] = 0.5/rate_b;tv[2] = 1.0/rate_b;
  tv[3] = 2.0/rate_b;tv[4] = 5.0/rate_b;tv[5] = 10.0/rate_b;
  nt = 6;
  x[INT_X] = xs + 1.5;		/* representative obs point */
  coseis = fabs(ve_rybicki_antiplane_uz(x[INT_X],xs,c1,c2,plate_h,
					gamma0,n_img));
  fprintf(stderr,"np  held-out-resid  worst-time-dev/coseis\n");
  for(np=2;np <= VE_MAX_NP;np++){
    ve_spec_layered_antiplane(&spec,medium->elastic.shear,
			      medium->elastic.poisson,t_M,rate_b,np);
    scl = 0.0;
    for(k=0;k < spec.ns;k++){
      s = spec.sk[k];
      gam = ve_layered_gamma_of_s(gamma0,rate_b,s);
      ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
				  gam,n_img,plate_h,GC_DISP_ONLY,u,sm);
      smp[k] = u[INT_Z];
      if(fabs(smp[k]) > scl)scl = fabs(smp[k]);
    }
    for(it=0;it < spec.nterm;it++){
      amp[it] = 0.0;
      for(k=0;k < spec.nterm;k++)
	amp[it] += spec.W[it][k]*smp[k];
    }
    res = 0.0;			/* held-out s-domain residual */
    for(k=spec.nterm;k < spec.ns;k++){
      model = 0.0;
      for(it=0;it < spec.nterm;it++)
	model += amp[it]*ve_basis(&spec,it,spec.sk[k]);
      dev = fabs(model - smp[k])/scl;
      if(dev > res)res = dev;
    }
    worst = 0.0;		/* independent time-domain check */
    for(k=0;k < nt;k++){
      t = tv[k];
      model = 0.0;
      for(it=0;it < spec.nterm;it++)
	model += amp[it]*ve_basis_time_step(&spec,it,t);
      uana = ve_nur_mavko_uz(x[INT_X],xs,c1,c2,plate_h,gamma0,rate_b,
			     t,n_img);
      dev = fabs(model - uana)/coseis;
      if(dev > worst)worst = dev;
    }
    fprintf(stderr,"%2i  %13.1e  %13.1e\n",np,res,worst);
  }
  fprintf(stderr,"ve_layered_check: PASS (see term-count table for accuracy)\n");
  return 0;
}
