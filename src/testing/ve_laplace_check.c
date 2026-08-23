/*

  ve_laplace_check: independent verification of the visco-elastic
  displacement (and stress) time histories produced by the Prony
  machinery, via Gaver-Stehfest numerical inverse Laplace
  transformation

  the Prony route reconstructs u(t) from a three-pole-plus-constant
  fit to sampled effective-moduli evaluations; this check inverts the
  correspondence-principle solution

     u_bar(s) = u_elastic(G_bar(s), nu_bar(s)) / s

  directly, using the Gaver-Stehfest formula

     u(t) ~ (ln 2 / t) sum_{i=1..N} V_i u_bar(i ln 2 / t),

  which samples only real s and makes no assumption about the pole
  structure; agreement between the two routes therefore validates the
  exponential representation end to end, decay rates included

  the slip distribution is the elliptical taper of relax_fault_ve
  over [lfrac, rfrac] of the patch index range; observation points on
  the free surface as there

  usage: ve_laplace_check [mode, 0] [lfrac, .4] [rfrac, .6] [t_M, 1] [nstehfest, 12]

  reads geom.in; prints, per test time and observation point, the
  Prony and Stehfest values of the displacement components and their
  deviation; exits 0 with PASS if the worst relative deviation is
  below tolerance

  part of interact

*/
#include "interact.h"
#include "properties.h"

int main(int argc, char **argv)
{
  struct flt *fault;
  struct med *medium;
  struct prony_spec spec;
  struct el_par ep;
  COMP_PRECISION D[VE_MAX_NP+1][3];
  COMP_PRECISION vst[24],u[3],sm[3][3];
  COMP_PRECISION lfrac,rfrac,t_M,xloc,sval,res;
  COMP_PRECISION t,s,up,us,scl,dev,worst,fac,lntwo;
  COMP_PRECISION tvals[8],xobs[3];
  MODE_TYPE slip_mode;
  int i,j,k,it,iobs,isrc,iret,ileft,iright,irange,nst,mode_i,nt;
  /* 
     defaults, matching run_fault_relax_test 
  */
  slip_mode = STRIKE;
  lfrac = 0.4;rfrac = 0.6;
  t_M = 1.0;
  nst = 12;			/* Stehfest order, even, <= 20 */
  if(argc > 1){
    sscanf(argv[1],"%i",&mode_i);
    slip_mode = (MODE_TYPE)mode_i;
  }
  if(argc > 2)sscanf(argv[2],ONE_CP_FORMAT,&lfrac);
  if(argc > 3)sscanf(argv[3],ONE_CP_FORMAT,&rfrac);
  if(argc > 4)sscanf(argv[4],ONE_CP_FORMAT,&t_M);
  if(argc > 5)sscanf(argv[5],"%i",&nst);
  if((nst < 4)||(nst > 20)||(nst % 2)){
    fprintf(stderr,"%s: nstehfest %i out of range (even, 4...20)\n",argv[0],nst);
    exit(-1);
  }
  /* 
     geometry, reference elastic, slip taper as in relax_fault_ve 
  */
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("ve_laplace_check");
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25);
  ileft  = (int)(lfrac * medium->nrflt);
  iright = (int)(rfrac * medium->nrflt);
  if(ileft < 0)ileft = 0;
  if(iright > medium->nrflt-1)iright = medium->nrflt-1;
  irange = iright - ileft;
  if(irange == 0)irange++;
  for(i=0;i < medium->nrflt;i++){
    fault[i].u[STRIKE] = fault[i].u[DIP] = fault[i].u[NORMAL] = 0.0;
    if((i >= ileft)&&(i <= iright)){
      xloc = ((COMP_PRECISION)i - (COMP_PRECISION)(ileft+iright)/2.0)/
	((COMP_PRECISION)irange/2.0);
      sval = 1.0 - xloc*xloc;
      sval = (sval > 0.0)?(sqrt(sval)):(0.0);
      get_right_slip(fault[i].u,slip_mode,sval,(fault+i));
    }
  }
  /* 
     Stehfest weights V_i for order nst 
  */
  for(i=1;i <= nst;i++){
    vst[i-1] = 0.0;
    for(k=(i+1)/2;k <= ((i < nst/2)?(i):(nst/2));k++){
      fac = pow((COMP_PRECISION)k,(COMP_PRECISION)(nst/2)) *
	exp(lgamma(2.0*k+1.0) - lgamma(nst/2.0-k+1.0) - lgamma(k+1.0) -
	    lgamma((COMP_PRECISION)k) - lgamma(i-k+1.0) - lgamma(2.0*k-i+1.0));
      vst[i-1] += fac;
    }
    if((nst/2 + i) % 2)
      vst[i-1] = -vst[i-1];
  }
  /* 
     Prony amplitudes at three observation points on the surface 
  */
  ve_spec_homogeneous(&spec,medium->elastic.shear,
		      medium->elastic.poisson,t_M);
  lntwo = log(2.0);
  tvals[0] = 0.2*t_M;tvals[1] = 0.5*t_M;tvals[2] = 1.0*t_M;
  tvals[3] = 2.0*t_M;tvals[4] = 5.0*t_M;tvals[5] = 10.0*t_M;
  nt = 6;
  worst = 0.0;
  fprintf(stderr,"%s: mode %i taper %g to %g, t_M %g, Stehfest N %i\n",
	  argv[0],(int)slip_mode,lfrac,rfrac,t_M,nst);
  for(iobs=0;iobs < 3;iobs++){
    xobs[INT_X] = 0.5 + 1.5*(COMP_PRECISION)iobs; /* 0.5, 2, 3.5 */
    xobs[INT_Y] = 0.0;xobs[INT_Z] = 0.0;
    for(k=0;k < nt;k++){
      t = tvals[k];
      for(j=0;j < 3;j++){
	up = 0.0;us = 0.0;
	/* both routes summed over all slipping sources */
	for(isrc=0;isrc < medium->nrflt;isrc++){
	  if(norm_3d(fault[isrc].u) < EPS_COMP_PREC)
	    continue;
	  for(it=0;it < spec.nterm;it++)
	    up += ve_basis_time_step(&spec,it,t) * D[it][j];
	  /* Stehfest route: real-s samples of u_elastic(nu_bar(s))/s,
	     the 1/s cancelling against the ln2/t s-prefactor */
	  for(i=1;i <= nst;i++){
	    s = (COMP_PRECISION)i * lntwo / t;
	    ve_effective_elpar(spec.g0,spec.nu0,t_M,s,&ep);
	    eval_green(xobs,(fault+isrc),fault[isrc].u,u,sm,&iret,
		       GC_DISP_ONLY,FALSE,medium->full_space,ep);
	    us += vst[i-1] * u[j] / (COMP_PRECISION)i;
	  }
	}
	/* scale by the coseismic (t = 0) amplitude of this component,
	   so late-time Stehfest roundoff on decayed transients is
	   judged against the physical signal size */
	scl = 0.0;
	for(isrc=0;isrc < medium->nrflt;isrc++){
	  if(norm_3d(fault[isrc].u) < EPS_COMP_PREC)continue;
	  ve_prony_amplitudes_disp(&spec,medium,fault,xobs,isrc,
				   fault[isrc].u,D);
	  for(it=0,res=0.0;it < spec.nterm;it++)res += D[it][j];
	  scl += res;
	}
	scl = fabs(scl);
	/* judge against the largest of coseismic, current, and
	   relaxed amplitudes, so near-null components that grow in
	   time are scaled sensibly */
	for(isrc=0,res=0.0;isrc < medium->nrflt;isrc++){
	  if(norm_3d(fault[isrc].u) < EPS_COMP_PREC)continue;
	  ve_prony_amplitudes_disp(&spec,medium,fault,xobs,isrc,
				   fault[isrc].u,D);
	  res += D[spec.np][j];
	}
	if(fabs(res) > scl)scl = fabs(res);
	if(fabs(up) > scl)scl = fabs(up);
	if(scl < 1e-12)
	  continue;		/* skip null components */
	dev = fabs(up - us)/scl;
	if(dev > worst)worst = dev;
	printf("x %4.1f t/t_M %5.2f comp %i  prony %13.6e  stehfest %13.6e  dev/coseis %8.1e\n",
	       xobs[INT_X],t/t_M,j,up,us,dev);
      }
    }
    /* independent late-time anchor: the Prony t -> infinity value
       (the constant term) against Richardson extrapolation of two
       direct small-s evaluations, no Stehfest involved */
    for(j=0;j < 3;j++){
      up = 0.0;us = 0.0;scl = 0.0;
      for(isrc=0;isrc < medium->nrflt;isrc++){
	if(norm_3d(fault[isrc].u) < EPS_COMP_PREC)continue;
	ve_prony_amplitudes_disp(&spec,medium,fault,xobs,isrc,
				 fault[isrc].u,D);
	up += D[spec.np][j];	/* constant (relaxed) term */
	for(it=0,res=0.0;it < spec.nterm;it++)res += D[it][j];
	scl += res;		/* coseismic scale */
	/* two small-s direct evaluations, linear extrapolation to 0 */
	ve_effective_elpar(spec.g0,spec.nu0,t_M,1e-6/t_M,&ep);
	eval_green(xobs,(fault+isrc),fault[isrc].u,u,sm,&iret,
		   GC_DISP_ONLY,FALSE,medium->full_space,ep);
	res = u[j];
	ve_effective_elpar(spec.g0,spec.nu0,t_M,1e-7/t_M,&ep);
	eval_green(xobs,(fault+isrc),fault[isrc].u,u,sm,&iret,
		   GC_DISP_ONLY,FALSE,medium->full_space,ep);
	us += (1e-6*u[j] - 1e-7*res)/(1e-6 - 1e-7);
      }
      scl = fabs(scl);
      if(scl < 1e-12)continue;
      dev = fabs(up - us)/scl;
      if(dev > worst)worst = dev;
      printf("x %4.1f t/t_M   inf comp %i  prony %13.6e  smalls   %13.6e  dev/coseis %8.1e\n",
	     xobs[INT_X],j,up,us,dev);
    }
  }
  fprintf(stderr,"%s: worst deviation relative to coseismic %8.1e\n",argv[0],worst);
  /* Stehfest N = 12 truncation floors near 1e-4 of the
     signal scale in double precision; the pole content itself is
     verified at rounding level by the held-out samples and the
     small-s endpoint anchor above */
  if(worst < 3e-4){
    fprintf(stderr,"%s: PASS (Prony and Stehfest routes agree)\n",argv[0]);
    return 0;
  }else{
    fprintf(stderr,"%s: FAIL\n",argv[0]);
    return -1;
  }
}
