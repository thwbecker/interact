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

  usage: ve_layered_check [plate_h, 2] [t_M, 1] [g2_over_g1, 1] [n_img, 60] [profile_prefix]

  with a fifth argument, additionally writes surface-profile data
  files for plotting, with ALL elements slipping by one unit (the
  full fault): {prefix}_elastic.dat (x/H, then per Gamma in {1, .75,
  .5, 0, -.5, -1}: u_analytic u_numeric), {prefix}_ve_disp.dat and
  {prefix}_ve_vel.dat (x/H, then per t' = t/tau_r in {0, .5, 1, 5,
  10}: analytic numeric), tau_r = 1/b; velocities are scaled by
  tau_r. analytic = Rybicki / Nur-Mavko series, numeric = the image
  machinery through the antiplane elements (elastic) and the np = 6
  Prony reconstruction fitted from its s-domain samples (VE)

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
  COMP_PRECISION gamv[6],tpv[5],cmin,cmax,xnorm,uana2,vana,vnum;
  COMP_PRECISION vamp[VE_MAX_NP],vmat[VE_MAX_NP][VE_MAX_NP+1],uel,piv;
  int ip,jp,kp;
  FILE *fpe,*fpd,*fpv;
  char fname[512];
  int ig,ix,nx;
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
  if(argc > 5){
    /* 
       surface-profile output for plotting: all elements slip by one
       unit; the analytic series uses the union depth interval of the
       elements (assumed contiguous and vertical)
    */
    cmin = 1e30;cmax = -1e30;
    for(i=0;i < medium->nrflt;i++){
      res = -(fault[i].x[INT_Y] + fault[i].l);
      if(res < cmin)cmin = res;
      res = -(fault[i].x[INT_Y] - fault[i].l);
      if(res > cmax)cmax = res;
    }
    fprintf(stderr,"profiles: full fault [%g, %g], files %s_{elastic,ve_disp,ve_vel}.dat\n",
	    cmin,cmax,argv[5]);
    gamv[0] = 1.0;gamv[1] = 0.75;gamv[2] = 0.5;
    gamv[3] = 0.0;gamv[4] = -0.5;gamv[5] = -1.0;
    tpv[0] = 0.0;tpv[1] = 0.5;tpv[2] = 1.0;tpv[3] = 5.0;tpv[4] = 10.0;
    snprintf(fname,512,"%s_elastic.dat",argv[5]);fpe = myopen(fname,"w");
    snprintf(fname,512,"%s_ve_disp.dat",argv[5]);fpd = myopen(fname,"w");
    snprintf(fname,512,"%s_ve_vel.dat",argv[5]);fpv = myopen(fname,"w");
    fprintf(fpe,"# x/H then per Gamma in {1, .75, .5, 0, -.5, -1}: u_ana u_num\n");
    fprintf(fpd,"# x/H then per t' in {0, .5, 1, 5, 10}: u_ana u_num (t' = t/tau_r, tau_r = 1/b)\n");
    fprintf(fpv,"# x/H then per t' in {0, .5, 1, 5, 10}: v_ana v_num (times tau_r)\n");
    nx = 50;
    ve_spec_layered_antiplane(&spec,medium->elastic.shear,
			      medium->elastic.poisson,t_M,rate_b,VE_MAX_NP);
    for(ix=0;ix < nx;ix++){
      xnorm = 0.1 + 9.9*(COMP_PRECISION)ix/((COMP_PRECISION)(nx-1));
      x[INT_X] = xs + xnorm * plate_h;
      /* elastic profiles over the Gamma family */
      fprintf(fpe,"%10.4f",xnorm);
      for(ig=0;ig < 6;ig++){
	uana = ve_rybicki_antiplane_uz(x[INT_X],xs,cmin,cmax,plate_h,
				       gamv[ig],n_img);
	model = 0.0;
	for(isrc=0;isrc < medium->nrflt;isrc++){
	  ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
				      gamv[ig],n_img,plate_h,
				      GC_DISP_ONLY,u,sm);
	  model += u[INT_Z];
	}
	fprintf(fpe," %12.5e %12.5e",uana,model);
      }
      fprintf(fpe,"\n");
      /* visco-elastic: fit the summed-source samples once per x */
      for(k=0;k < spec.ns;k++){
	s = spec.sk[k];
	gam = ve_layered_gamma_of_s(gamma0,rate_b,s);
	smp[k] = 0.0;
	for(isrc=0;isrc < medium->nrflt;isrc++){
	  ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
				      gam,n_img,plate_h,
				      GC_DISP_ONLY,u,sm);
	  smp[k] += u[INT_Z];
	}
      }
      for(it=0;it < spec.nterm;it++){
	amp[it] = 0.0;
	for(k=0;k < spec.nterm;k++)
	  amp[it] += spec.W[it][k]*smp[k];
      }
      /* 
	 velocity-consistent fit: the samples are the elastic solution
	 at the effective moduli, i.e. s ubar(s), so the velocity
	 transform is vbar(s) = s ubar(s) - u(0+) = sample - elastic
	 value, with basis transforms tau_p/(1 + s tau_p); fitting
	 these directly avoids differentiating the displacement fit,
	 whose amplitude oscillation the derivative amplifies
      */
      ve_layered_antiplane_sample(medium,fault,0,slip,x,gamma0,n_img,
				  plate_h,GC_DISP_ONLY,u,sm);
      uel = 0.0;
      for(isrc=0;isrc < medium->nrflt;isrc++){
	ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
				    gamma0,n_img,plate_h,
				    GC_DISP_ONLY,u,sm);
	uel += u[INT_Z];
      }
      for(ip=0;ip < spec.np;ip++){	/* np x np system on sk[1..np] */
	s = spec.sk[ip+1];
	for(jp=0;jp < spec.np;jp++)
	  vmat[ip][jp] = spec.tau[jp]/(1.0 + s*spec.tau[jp]);
	vmat[ip][spec.np] = smp[ip+1] - uel;
      }
      for(ip=0;ip < spec.np;ip++){	/* small Gauss, partial pivot */
	piv = fabs(vmat[ip][ip]);kp = ip;
	for(jp=ip+1;jp < spec.np;jp++)
	  if(fabs(vmat[jp][ip]) > piv){piv = fabs(vmat[jp][ip]);kp = jp;}
	if(kp != ip)
	  for(jp=ip;jp <= spec.np;jp++){
	    piv = vmat[ip][jp];vmat[ip][jp] = vmat[kp][jp];vmat[kp][jp] = piv;
	  }
	for(jp=ip+1;jp < spec.np;jp++){
	  piv = vmat[jp][ip]/vmat[ip][ip];
	  for(kp=ip;kp <= spec.np;kp++)
	    vmat[jp][kp] -= piv*vmat[ip][kp];
	}
      }
      for(ip=spec.np-1;ip >= 0;ip--){
	vamp[ip] = vmat[ip][spec.np];
	for(jp=ip+1;jp < spec.np;jp++)
	  vamp[ip] -= vmat[ip][jp]*vamp[jp];
	vamp[ip] /= vmat[ip][ip];
      }
      fprintf(fpd,"%10.4f",xnorm);
      fprintf(fpv,"%10.4f",xnorm);
      for(k=0;k < 5;k++){
	t = tpv[k]/rate_b;
	uana2 = ve_nur_mavko_uz(x[INT_X],xs,cmin,cmax,plate_h,gamma0,
				rate_b,t,n_img);
	vana = ve_nur_mavko_vz(x[INT_X],xs,cmin,cmax,plate_h,gamma0,
			       rate_b,t,n_img)/rate_b;
	model = 0.0;vnum = 0.0;
	for(it=0;it < spec.nterm;it++)
	  model += amp[it]*ve_basis_time_step(&spec,it,t);
	for(it=0;it < spec.np;it++)	/* velocity from its own fit */
	  vnum += vamp[it]*exp(-t/spec.tau[it]);
	vnum /= rate_b;
	fprintf(fpd," %12.5e %12.5e",uana2,model);
	fprintf(fpv," %12.5e %12.5e",vana,vnum);
      }
      fprintf(fpd,"\n");
      fprintf(fpv,"\n");
    }
    fclose(fpe);fclose(fpd);fclose(fpv);
  }
  return 0;
}
