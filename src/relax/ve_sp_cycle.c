/*

  ve_sp_cycle: Savage-Prescott style earthquake-cycle testbed for the
  2-D antiplane elastic plate over a Maxwell visco-elastic half space
  (step 7a of rsf_ve_implementation_plan.md, cycle stage)

  geometry: a vertical antiplane fault (w = 0 elements, dip -90,
  strike 0, x = 0) cutting the elastic plate from the free surface
  (y = 0) to the interface (y = -H), read from geom.in, e.g. from

    makefault -z 0 -n <n> -y -<H/2> -l <H/2> | \
        gawk '{print($1,$2,$3,$4,-$5,$6,0,$8)}' > geom.in

  cycle: uniform slip du = V_pl * T on all patches every T, no other
  loading; the relaxed (constant) part of the layered kernel carries
  the secular plate motion, so the spun-up far field moves at
  +/- V_pl/2 and the substrate relaxation re-loads the fault
  interseismically, which is the Savage and Prescott (1978)
  construction

  three routes are computed and compared:

  route A, the MACHINERY UNDER TEST: per-receiver Prony amplitudes
  (np-term ladder plus relaxed constant, ve_spec_layered_antiplane)
  fitted from s-domain samples of the layered elastic image kernel
  through the actual antiplane elements, advanced through the cycle
  by the exact-exponential per-pole state update (the rsf_solve
  step 5 h mode: h_p += a_p du at an event, exponential decay
  between events); surface velocity is the h-state decay sum
  -sum_p h_p/tau_p and fault stress the resolved h-state sum

  route B, the exact-in-time reference: the same element-based image
  families X_n(receiver) but with the EXACT Nur-Mavko time weights
  W_n(t) (displacement/stress) and Wdot_n(t) (velocity) summed
  explicitly over all past events; differences A - B isolate the
  Prony fit error with the spatial part held identical

  route C, fully analytic (arctan image series, independent
  implementation), lives in the companion script sp_analytic.py and
  is compared against the files this program writes

  internal checks (program exits nonzero on violation): the
  sequential h recurrence against its geometric closed form
  (rounding); spin-up convergence; the spun-up cycle-mean surface
  velocity (interseismic drift plus coseismic jump) against the
  plate velocity V_pl/2 at the far observation point; route A vs
  route B deviation gates at fit-accuracy consistent levels

  usage:
    ve_sp_cycle [H, 2] [t_M, 1] [g2_over_g1, 1] [T_over_taur, 1]
                [np, 6] [n_img_fit, 60] [prefix, sp]

  taur = 1/b is the layered relaxation time (2 t_M at equal
  rigidity); output files {prefix}_vel_profiles.dat,
  {prefix}_disp_profiles.dat, {prefix}_stress_hist.dat,
  {prefix}_surf_hist.dat, {prefix}_checks.dat

  units: velocities in V_pl, displacements in du = V_pl T, times in
  T, distances in H, stress raw (G = 1, event slip du with V_pl = 1)

  part of interact - this program only ADDS a binary; no shared code
  is modified

*/
#include "interact.h"
#include "properties.h"
#include "prony_kernel.h"

#define VESP_NIMG_B 600		/* max image families (sampling and exact route B) */
#define VESP_NT 9		/* output cycle fractions */
#define VESP_NHIST 240		/* history samples over the final two cycles */
#define VESP_NOBS 40		/* surface profile points */
#define VESP_NDEPTH 3		/* stress history receiver depths */

/* 

   element-based spatial value of image family n for one receiver:
   sum over all source patches (unit strike slip) of the
   reference-moduli elastic kernel of the source copy shifted down by
   2 n H; mode GC_DISP_ONLY (u_z at point x) or GC_STRESS_ONLY
   (stress resolved in STRIKE mode on receiver patch irec, proper
   self-interaction for the unshifted family)

*/
static COMP_PRECISION family_value(struct med *medium, struct flt *fault,
				   int n, COMP_PRECISION plate_h,
				   MODE_TYPE mode,
				   COMP_PRECISION *x, int irec)
{
  struct flt img[1];
  COMP_PRECISION slip[3],u[3],sm[3][3],val;
  int isrc,iret;
  slip[STRIKE] = 1.0;slip[DIP] = slip[NORMAL] = 0.0;
  val = 0.0;
  for(isrc=0;isrc < medium->nrflt;isrc++){
    if((n == 0) && (mode == GC_STRESS_ONLY)){
      eval_green_at_receiver(fault,irec,isrc,slip,u,sm,&iret,
			     GC_STRESS_ONLY,TRUE,FALSE,
			     medium->elastic);
    }else{
      *img = fault[isrc];
      img->x[INT_Y] = fault[isrc].x[INT_Y] - 2.0*(COMP_PRECISION)n*plate_h;
      eval_green((mode == GC_STRESS_ONLY)?(fault[irec].x):(x),
		 img,slip,u,sm,&iret,mode,FALSE,FALSE,
		 medium->elastic);
    }
    if(iret)
      fprintf(stderr,"family_value: WARNING: singular, family %i, src %i\n",
	      n,isrc);
    if(mode == GC_STRESS_ONLY)
      val += resolve_stress_on_fault(sm,(fault+irec),DIP);
    else
      val += u[INT_Z];
    if(n > 0){			/* second (surface-mirrored) image
				   family, see ve_rybicki_antiplane_uz */
      *img = fault[isrc];
      img->x[INT_Y] = -fault[isrc].x[INT_Y] - 2.0*(COMP_PRECISION)n*plate_h;
      eval_green((mode == GC_STRESS_ONLY)?(fault[irec].x):(x),
		 img,slip,u,sm,&iret,mode,FALSE,FALSE,
		 medium->elastic);
      if(iret)
	fprintf(stderr,"family_value: WARNING: singular, mirror family %i, src %i\n",
		n,isrc);
      if(mode == GC_STRESS_ONLY)
	val += resolve_stress_on_fault(sm,(fault+irec),DIP);
      else
	val += u[INT_Z];
    }
  }
  return val;
}
/* 

   exact Nur-Mavko weights W_n (step response) and Wdot_n (rate) for
   all n = 0..nmax at one age: upward recurrence for Gamma_0 = 0
   (W_n = P(n, b age) Erlang ramp, Wdot_n Erlang density), binomial
   routines otherwise

*/
static void nm_weights_all(int nmax, COMP_PRECISION gamma0,
			   COMP_PRECISION b, COMP_PRECISION age,
			   COMP_PRECISION *w, COMP_PRECISION *wdot)
{
  COMP_PRECISION xx,eterm,psum;
  int n;
  if(fabs(gamma0) < 1e-14){
    xx = b * age;
    w[0] = 1.0;wdot[0] = 0.0;
    eterm = exp(-xx);		/* exp(-x) x^k/k! at k = 0 */
    psum = 0.0;
    for(n=1;n <= nmax;n++){
      psum += eterm;		/* now sum_{k<n} */
      w[n] = 1.0 - psum;
      if(w[n] < 0.0)w[n] = 0.0;	/* rounding guard */
      wdot[n] = b * eterm;	/* b exp(-x) x^{n-1}/(n-1)! */
      eterm *= xx/((COMP_PRECISION)n);
    }
  }else{
    for(n=0;n <= nmax;n++){
      w[n]    = ve_nur_mavko_weight(n,gamma0,b,age);
      wdot[n] = ve_nur_mavko_weight_dot(n,gamma0,b,age);
    }
  }
}
/* 

   route B: event-summed weights at age t after the most recent
   event, K events T apart; skip_last skips the most recent event
   (for evaluating the previous cycle); sw and swd have nmax+1 entries

*/
static void route_b_weights(int nmax, COMP_PRECISION gamma0,
			    COMP_PRECISION b, COMP_PRECISION t,
			    COMP_PRECISION T, int K, int skip_last,
			    COMP_PRECISION *sw, COMP_PRECISION *swd)
{
  COMP_PRECISION wv[VESP_NIMG_B+1],wdv[VESP_NIMG_B+1],age;
  int n,j;
  for(n=0;n <= nmax;n++){sw[n] = 0.0;swd[n] = 0.0;}
  for(j=(skip_last?1:0);j < K;j++){
    age = t + (COMP_PRECISION)j*T;
    nm_weights_all(nmax,gamma0,b,age,wv,wdv);
    for(n=0;n <= nmax;n++){
      sw[n] += wv[n];swd[n] += wdv[n];
    }
  }
}
/* 

   route A: value and rate at age t after the most recent event from
   the h states (which hold the state at age 0+, i.e. just after that
   event); skip_last removes the most recent event's contribution
   first (previous cycle); a = the receiver's amplitude row, du the
   event slip

*/
static void route_a_val(struct prony_spec *spec, COMP_PRECISION *hrow,
			COMP_PRECISION cacc_r, COMP_PRECISION *a,
			COMP_PRECISION du, COMP_PRECISION t,
			int skip_last,
			COMP_PRECISION *val, COMP_PRECISION *rate)
{
  COMP_PRECISION h;
  int ip;
  *val = cacc_r - (skip_last?(a[spec->np]*du):(0.0));
  *rate = 0.0;
  for(ip=0;ip < spec->np;ip++){
    h = hrow[ip] - (skip_last?(a[ip]*du):(0.0));
    h *= exp(-t/spec->tau[ip]);
    *val += h;
    *rate -= h/spec->tau[ip];
  }
}

int main(int argc, char **argv)
{
  struct flt *fault;
  struct med *medium;
  struct prony_spec spec,spec_s;
  COMP_PRECISION plate_h,t_M,g2fac,gamma0,rate_b,taur,T,du,Vpl,T_over_taur;
  COMP_PRECISION slip[3],u[3],sm[3][3],x[3];
  COMP_PRECISION cmin,cmax,res,scl,s,gam,dev,worst_res,tslow;
  COMP_PRECISION t,geo,efac,tt;
  COMP_PRECISION *amp_d,*amp_s,*xn_d,*xn_s,*hp,*cacc,*arow;
  COMP_PRECISION swv[VESP_NIMG_B+1],swdv[VESP_NIMG_B+1];
  COMP_PRECISION tfrac[VESP_NT] = {0.02,0.05,0.1,0.2,0.35,0.5,0.7,0.85,0.98};
  COMP_PRECISION xobs[VESP_NOBS],zdep[VESP_NDEPTH] = {0.25,0.5,0.9};
  COMP_PRECISION worst_ab_v,worst_ab_u,worst_ab_s,mean_vA,mean_vB;
  COMP_PRECISION sscl,vA,uA,uA0,vB,uB,uB0,sA,sB;
  int np,n_img,nib,i,j,k,n,it,ip,irec,nrflt,nobs,K,idep[VESP_NDEPTH];
  int nrec_d,nrec_s,worst_rec,iret,isrc,io[2],jp,kp;
  COMP_PRECISION *smp_d,*vamp_d,vmat[VE_MAX_NP][VE_MAX_NP+1],piv;
  COMP_PRECISION vA2,qp,worst_ab_v2;
  FILE *fv,*fd,*fs,*fh,*fc;
  char fname[512],prefix[256];
  struct flt img[1];
  /* defaults */
  plate_h = 2.0;t_M = 1.0;g2fac = 1.0;T_over_taur = 1.0;
  np = VE_MAX_NP;n_img = VESP_NIMG_B;
  snprintf(prefix,256,"sp");
  Vpl = 1.0;
  if(argc > 1)sscanf(argv[1],ONE_CP_FORMAT,&plate_h);
  if(argc > 2)sscanf(argv[2],ONE_CP_FORMAT,&t_M);
  if(argc > 3)sscanf(argv[3],ONE_CP_FORMAT,&g2fac);
  if(argc > 4)sscanf(argv[4],ONE_CP_FORMAT,&T_over_taur);
  if(argc > 5)sscanf(argv[5],"%i",&np);
  if(argc > 6)sscanf(argv[6],"%i",&n_img);
  if(argc > 7)snprintf(prefix,256,"%s",argv[7]);
  if(n_img > VESP_NIMG_B)n_img = VESP_NIMG_B;
  nib = n_img;			/* image depth, sampling AND route B */
  /* geometry, reference moduli, layered rates */
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("ve_sp_cycle");
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25);
  ve_layered_gamma_pars(medium->elastic.shear,
			g2fac*medium->elastic.shear,t_M,
			&gamma0,&rate_b);
  taur = 1.0/rate_b;
  T = T_over_taur * taur;
  du = Vpl * T;
  nrflt = medium->nrflt;
  /* through-plate geometry check */
  cmin = 1e30;cmax = -1e30;
  for(i=0;i < nrflt;i++){
    res = -(fault[i].x[INT_Y] + fault[i].l);
    if(res < cmin)cmin = res;
    res = -(fault[i].x[INT_Y] - fault[i].l);
    if(res > cmax)cmax = res;
  }
  if((fabs(cmin) > 1e-8) || (fabs(cmax - plate_h) > 1e-8)){
    fprintf(stderr,"ve_sp_cycle: fault depth range [%g, %g] does not span the plate [0, %g]\n",
	    cmin,cmax,plate_h);
    fprintf(stderr,"ve_sp_cycle: the Savage-Prescott cycle requires a through-plate fault\n");
    exit(-1);
  }
  fprintf(stderr,"ve_sp_cycle: H %g t_M %g g2/g1 %g (Gamma0 %.4f b %.4f taur %.4f) T/taur %.4f np %i n_img %i/%i, %i patches\n",
	  plate_h,t_M,g2fac,gamma0,rate_b,taur,T/taur,np,n_img,
	  nib,nrflt);
  /* observation profile x/H log-spaced in [0.1, 10] */
  nobs = VESP_NOBS;
  for(i=0;i < nobs;i++)
    xobs[i] = plate_h * 0.1 *
      pow(100.0,(COMP_PRECISION)i/((COMP_PRECISION)(nobs-1)));
  /* nearest patches to the stress history depths */
  for(k=0;k < VESP_NDEPTH;k++){
    idep[k] = 0;res = 1e30;
    for(i=0;i < nrflt;i++){
      dev = fabs(-fault[i].x[INT_Y] - zdep[k]*plate_h);
      if(dev < res){res = dev;idep[k] = i;}
    }
  }
  nrec_d = nobs;nrec_s = nrflt;
  /* 
     route A assembly: s-domain samples -> per-receiver amplitudes;
     the spec uses a wider ladder and lower minimum sample than
     ve_spec_layered_antiplane because the spun-up cycle accumulates
     slow image content: rates on [b/45, 1.5 b], samples log-spaced
     on [b/100, 5 b]; the sampling image depth n_img must satisfy
     Gamma(s_min)^n_img << 1, i.e. n_img >> b/s_min, or the relaxed
     term is biased low (through-plate far-field tail goes like
     x/(2 pi H n_img))
  */
  spec.g0 = medium->elastic.shear;spec.nu0 = medium->elastic.poisson;
  spec.t_M = t_M;spec.bulk0 = bulk_mod_from_G_nu(spec.g0,spec.nu0);
  spec.np = np;
  for(i=0;i < np;i++)		/* geometric ladder [b/45, 1.5 b] */
    spec.tau[i] = 1.0/(rate_b*(1.0/45.0)*
		       pow(67.5,(COMP_PRECISION)i/
			   ((np > 1)?((COMP_PRECISION)(np-1)):(1.0))));
  spec.has_const = TRUE;
  spec.nterm = np + 1;
  spec.ns = spec.nterm + VE_NHELD;
  for(i=0;i < spec.nterm;i++)	/* fit samples [b/100, 5 b] */
    spec.sk[i] = rate_b*(1.0/100.0)*
      pow(500.0,(COMP_PRECISION)i/((COMP_PRECISION)spec.nterm - 1.0));
  spec.sk[spec.nterm]   = rate_b/30.0;	/* held out */
  spec.sk[spec.nterm+1] = rate_b*0.6;
  ve_solve_weights(&spec);
  /* stress spec: same ladder, NO constant basis (the through-plate
     relaxed stress is exactly zero; fitting a constant would hand
     its truncation noise a factor-of-K secular amplification in the
     spun-up event sum) */
  spec_s = spec;
  spec_s.has_const = FALSE;
  spec_s.nterm = np;
  spec_s.ns = spec_s.nterm + VE_NHELD;
  for(i=0;i < spec_s.nterm;i++)	/* same span, one fewer point */
    spec_s.sk[i] = rate_b*(1.0/100.0)*
      pow(500.0,(COMP_PRECISION)i/((COMP_PRECISION)spec_s.nterm - 1.0));
  spec_s.sk[spec_s.nterm]   = rate_b/30.0;
  spec_s.sk[spec_s.nterm+1] = rate_b*0.6;
  ve_solve_weights(&spec_s);
  amp_d = (COMP_PRECISION *)calloc((size_t)spec.nterm*nrec_d,sizeof(COMP_PRECISION));
  amp_s = (COMP_PRECISION *)calloc((size_t)spec.nterm*nrec_s,sizeof(COMP_PRECISION));
  if((!amp_d)||(!amp_s))MEMERROR("ve_sp_cycle: amps");
  slip[STRIKE] = 1.0;slip[DIP] = slip[NORMAL] = 0.0;
  worst_res = 0.0;worst_rec = 0;
  x[INT_Y] = x[INT_Z] = 0.0;
  smp_d  = (COMP_PRECISION *)calloc((size_t)spec.ns*nrec_d,sizeof(COMP_PRECISION));
  vamp_d = (COMP_PRECISION *)calloc((size_t)spec.np*nrec_d,sizeof(COMP_PRECISION));
  if((!smp_d)||(!vamp_d))MEMERROR("ve_sp_cycle: vfit");
  {
    COMP_PRECISION smp[VE_MAX_NS],wfac;
    for(irec=0;irec < nrec_d + nrec_s;irec++){
      struct prony_spec *sp = (irec < nrec_d)?(&spec):(&spec_s);
      for(k=0;k < sp->ns;k++){
	s = sp->sk[k];
	gam = ve_layered_gamma_of_s(gamma0,rate_b,s);
	smp[k] = 0.0;
	if(irec < nrec_d){	/* displacement receiver */
	  x[INT_X] = xobs[irec];
	  for(isrc=0;isrc < nrflt;isrc++){
	    ve_layered_antiplane_sample(medium,fault,isrc,slip,x,gam,
					nib,plate_h,GC_DISP_ONLY,u,sm);
	    smp[k] += u[INT_Z];
	  }
	}else{			/* stress receiver */
	  int jr = irec - nrec_d;
	  for(isrc=0;isrc < nrflt;isrc++){
	    wfac = 1.0;
	    for(n=0;n <= nib;n++){
	      if(fabs(wfac) < 1e-16)break;
	      if(n == 0){
		eval_green_at_receiver(fault,jr,isrc,slip,u,sm,&iret,
				       GC_STRESS_ONLY,TRUE,FALSE,
				       medium->elastic);
	      }else{
		*img = fault[isrc];
		img->x[INT_Y] = fault[isrc].x[INT_Y] -
		  2.0*(COMP_PRECISION)n*plate_h;
		eval_green(fault[jr].x,img,slip,u,sm,&iret,
			   GC_STRESS_ONLY,FALSE,FALSE,medium->elastic);
	      }
	      smp[k] += wfac * resolve_stress_on_fault(sm,(fault+jr),DIP);
	      if(n > 0){	/* second image family */
		*img = fault[isrc];
		img->x[INT_Y] = -fault[isrc].x[INT_Y] -
		  2.0*(COMP_PRECISION)n*plate_h;
		eval_green(fault[jr].x,img,slip,u,sm,&iret,
			   GC_STRESS_ONLY,FALSE,FALSE,medium->elastic);
		smp[k] += wfac * resolve_stress_on_fault(sm,(fault+jr),DIP);
	      }
	      wfac *= gam;
	    }
	  }
	}
      }
      if(irec < nrec_d)		/* keep samples for the velocity fit */
	for(k=0;k < spec.ns;k++)
	  smp_d[(size_t)spec.ns*irec+k] = smp[k];
      /* amplitude rows are (np + 1) wide for BOTH kinds; the stress
	 rows carry an explicit zero constant */
      arow = (irec < nrec_d)?(amp_d + (size_t)spec.nterm*irec):
	(amp_s + (size_t)spec.nterm*(irec - nrec_d));
      for(it=0;it < sp->nterm;it++){
	arow[it] = 0.0;
	for(k=0;k < sp->nterm;k++)
	  arow[it] += sp->W[it][k]*smp[k];
      }
      if(irec >= nrec_d)
	arow[spec.np] = 0.0;	/* exact zero relaxed stress */
      scl = 0.0;
      for(k=0;k < sp->ns;k++)
	if(fabs(smp[k]) > scl)scl = fabs(smp[k]);
      if(scl > 0.0){
	for(k=sp->nterm;k < sp->ns;k++){
	  res = 0.0;
	  for(it=0;it < sp->nterm;it++)
	    res += arow[it]*ve_basis(sp,it,sp->sk[k]);
	  dev = fabs(res - smp[k])/scl;
	  if(dev > worst_res){worst_res = dev;worst_rec = irec;}
	}
      }
    }
  }
  fprintf(stderr,"ve_sp_cycle: assembly done, worst held-out residual %.1e (receiver %i of %i)\n",
	  worst_res,worst_rec,nrec_d+nrec_s);
  /* 
     route B assembly: exact family values per receiver 
  */
  xn_d = (COMP_PRECISION *)calloc((size_t)(VESP_NIMG_B+1)*nrec_d,sizeof(COMP_PRECISION));
  xn_s = (COMP_PRECISION *)calloc((size_t)(VESP_NIMG_B+1)*nrec_s,sizeof(COMP_PRECISION));
  if((!xn_d)||(!xn_s))MEMERROR("ve_sp_cycle: families");
  x[INT_Y] = x[INT_Z] = 0.0;
  for(n=0;n <= nib;n++){
    for(i=0;i < nrec_d;i++){
      x[INT_X] = xobs[i];
      xn_d[(size_t)n*nrec_d+i] = family_value(medium,fault,n,plate_h,
					      GC_DISP_ONLY,x,0);
    }
    for(i=0;i < nrec_s;i++)
      xn_s[(size_t)n*nrec_s+i] = family_value(medium,fault,n,plate_h,
					      GC_STRESS_ONLY,x,i);
  }
  fprintf(stderr,"ve_sp_cycle: exact families assembled (%i)\n",nib+1);
  /* 
     velocity-consistent s-domain fit per displacement receiver:
     vbar(s) = sample(s) - elastic, basis tau_p/(1 + s tau_p), np x np
     system on the interior fit samples (see ve_layered_check)
  */
  for(i=0;i < nrec_d;i++){
    COMP_PRECISION uel;
    COMP_PRECISION wf = 1.0;
    uel = 0.0;			/* exact elastic: Gamma0^n weighted families */
    for(n=0;n <= nib;n++){
      uel += wf*xn_d[(size_t)n*nrec_d+i];
      wf *= gamma0;
      if(fabs(wf) < 1e-16)break;
    }
    for(ip=0;ip < spec.np;ip++){
      s = spec.sk[ip+1];
      for(jp=0;jp < spec.np;jp++)
	vmat[ip][jp] = spec.tau[jp]/(1.0 + s*spec.tau[jp]);
      vmat[ip][spec.np] = smp_d[(size_t)spec.ns*i+ip+1] - uel;
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
      vamp_d[(size_t)spec.np*i+ip] = vmat[ip][spec.np];
      for(jp=ip+1;jp < spec.np;jp++)
	vamp_d[(size_t)spec.np*i+ip] -= vmat[ip][jp]*vamp_d[(size_t)spec.np*i+jp];
      vamp_d[(size_t)spec.np*i+ip] /= vmat[ip][ip];
    }
  }
  /* spin-up length from the slowest content */
  tslow = spec.tau[0];
  for(ip=1;ip < spec.np;ip++)
    if(spec.tau[ip] > tslow)tslow = spec.tau[ip];
  /* route B completeness: the order-n Erlang front arrives at age
     n taur; cover nib fronts plus margin, and 30 e-folds of the
     slowest ladder pole */
  K = (int)(((COMP_PRECISION)nib +
	     20.0*sqrt((COMP_PRECISION)nib))*taur/T + 30.0*tslow/T) + 3;
  fprintf(stderr,"ve_sp_cycle: spinning up over %i cycles\n",K);
  /* 
     route A: sequential exact-exponential h updates over K events 
  */
  hp   = (COMP_PRECISION *)calloc((size_t)(nrec_d+nrec_s)*spec.np,sizeof(COMP_PRECISION));
  cacc = (COMP_PRECISION *)calloc((size_t)(nrec_d+nrec_s),sizeof(COMP_PRECISION));
  if((!hp)||(!cacc))MEMERROR("ve_sp_cycle: states");
  for(j=0;j < K;j++){
    for(irec=0;irec < nrec_d+nrec_s;irec++){
      arow = (irec < nrec_d)?(amp_d + (size_t)spec.nterm*irec):
	(amp_s + (size_t)spec.nterm*(irec-nrec_d));
      for(ip=0;ip < spec.np;ip++)
	hp[(size_t)irec*spec.np+ip] += arow[ip]*du;
      cacc[irec] += arow[spec.np]*du;	/* zero for stress rows */
    }
    if(j < K-1){
      for(ip=0;ip < spec.np;ip++){
	efac = exp(-T/spec.tau[ip]);
	for(irec=0;irec < nrec_d+nrec_s;irec++)
	  hp[(size_t)irec*spec.np+ip] *= efac;
      }
    }
  }
  /* recurrence vs geometric closed form */
  dev = 0.0;
  for(irec=0;irec < nrec_d+nrec_s;irec++){
    arow = (irec < nrec_d)?(amp_d + (size_t)spec.nterm*irec):
      (amp_s + (size_t)spec.nterm*(irec-nrec_d));
    for(ip=0;ip < spec.np;ip++){
      efac = exp(-T/spec.tau[ip]);
      geo = arow[ip]*du*(1.0 - pow(efac,(COMP_PRECISION)K))/(1.0 - efac);
      res = fabs(hp[(size_t)irec*spec.np+ip] - geo)/(fabs(geo) + 1e-30);
      if(res > dev)dev = res;
    }
  }
  fprintf(stderr,"ve_sp_cycle: h recurrence vs geometric closed form: %.1e %s\n",
	  dev,(dev < 1e-10)?("(ok)"):("(FAIL)"));
  if(dev > 1e-10)exit(-1);
  /* spin-up remainder */
  res = exp(-((COMP_PRECISION)K)*T/tslow);
  fprintf(stderr,"ve_sp_cycle: spin-up remainder (slowest content) %.1e %s\n",
	  res,(res < 1e-9)?("(ok)"):("(FAIL)"));
  if(res > 1e-9)exit(-1);
  /* 
     outputs; state is just AFTER the K-th event, t = age since it 
  */
  snprintf(fname,512,"%s_vel_profiles.dat",prefix);fv = myopen(fname,"w");
  snprintf(fname,512,"%s_disp_profiles.dat",prefix);fd = myopen(fname,"w");
  snprintf(fname,512,"%s_stress_hist.dat",prefix);fs = myopen(fname,"w");
  snprintf(fname,512,"%s_surf_hist.dat",prefix);fh = myopen(fname,"w");
  snprintf(fname,512,"%s_checks.dat",prefix);fc = myopen(fname,"w");
  fprintf(fv,"# spun-up interseismic surface velocity / V_pl; H %g Gamma0 %g b %g T/taur %g np %i n_img %i K %i\n",
	  plate_h,gamma0,rate_b,T/taur,np,n_img,K);
  fprintf(fv,"# x/H then per t/T in {");
  for(k=0;k < VESP_NT;k++)fprintf(fv,"%g ",tfrac[k]);
  fprintf(fv,"}: v_A(h-state) v_A2(velocity-fit) v_B(exact)\n");
  fprintf(fd,"# spun-up postseismic surface displacement (u(t)-u(0+))/du; columns as the velocity file\n");
  fprintf(fs,"# resolved fault stress over the final two cycles (G = 1, event slip du = V_pl T)\n");
  fprintf(fs,"# t/T (negative = previous cycle) then per depth z/H in {%g %g %g}: s_A s_B\n",
	  zdep[0],zdep[1],zdep[2]);
  fprintf(fh,"# postseismic surface displacement (u(t)-u(0+)) at x/H = 1 and 3, units du\n");
  fprintf(fh,"# t/T u_A(1) u_B(1) u_A(3) u_B(3)\n");
  worst_ab_v = worst_ab_v2 = worst_ab_u = worst_ab_s = 0.0;
  /* stress deviations are collected raw and scaled at the end by the
     maximum |route B| stress of the history (the cycle stress swing;
     the coseismic stress of UNIFORM through-plate slip nearly
     vanishes at interior patches and is no scale) */
  sscl = 0.0;
  /* (a) velocity and displacement profiles at the cycle fractions */
  for(i=0;i < nobs;i++){
    fprintf(fv,"%12.6f",xobs[i]/plate_h);
    fprintf(fd,"%12.6f",xobs[i]/plate_h);
    arow = amp_d + (size_t)spec.nterm*i;
    route_a_val(&spec,hp+(size_t)i*spec.np,cacc[i],arow,du,1e-12*T,0,
		&uA0,&vA);
    route_b_weights(nib,gamma0,rate_b,1e-12*T,T,K,0,swv,swdv);
    uB0 = 0.0;
    for(n=0;n <= nib;n++)
      uB0 += du*swv[n]*xn_d[(size_t)n*nrec_d+i];
    for(k=0;k < VESP_NT;k++){
      t = tfrac[k]*T;
      route_a_val(&spec,hp+(size_t)i*spec.np,cacc[i],arow,du,t,0,&uA,&vA);
      /* velocity-consistent fit, spun up by the geometric event sum */
      vA2 = 0.0;
      for(ip=0;ip < spec.np;ip++){
	qp = exp(-T/spec.tau[ip]);
	vA2 += vamp_d[(size_t)spec.np*i+ip]*exp(-t/spec.tau[ip])*
	  (1.0 - pow(qp,(COMP_PRECISION)K))/(1.0 - qp);
      }
      vA2 *= du;
      route_b_weights(nib,gamma0,rate_b,t,T,K,0,swv,swdv);
      uB = vB = 0.0;
      for(n=0;n <= nib;n++){
	uB += du*swv[n]*xn_d[(size_t)n*nrec_d+i];
	vB += du*swdv[n]*xn_d[(size_t)n*nrec_d+i];
      }
      fprintf(fv," %14.6e %14.6e %14.6e",vA/Vpl,vA2/Vpl,vB/Vpl);
      fprintf(fd," %14.6e %14.6e",(uA-uA0)/du,(uB-uB0)/du);
      dev = fabs(vA-vB)/Vpl;if(dev > worst_ab_v)worst_ab_v = dev;
      dev = fabs(vA2-vB)/Vpl;if(dev > worst_ab_v2)worst_ab_v2 = dev;
      dev = fabs((uA-uA0)-(uB-uB0))/du;if(dev > worst_ab_u)worst_ab_u = dev;
    }
    fprintf(fv,"\n");fprintf(fd,"\n");
  }
  /* (b) cycle-mean velocity check: interseismic drift plus the
     coseismic jump must equal the plate velocity V_pl/2; evaluated
     at x/H ~ 4 because the image-truncation tail of both routes
     grows like x/(2 pi H n_img) */
  i = 0;
  {
    COMP_PRECISION xtarget = 4.0*plate_h,best = 1e30;
    for(k=0;k < nobs;k++)
      if(fabs(xobs[k] - xtarget) < best){best = fabs(xobs[k] - xtarget);i = k;}
  }
  arow = amp_d + (size_t)spec.nterm*i;
  route_a_val(&spec,hp+(size_t)i*spec.np,cacc[i],arow,du,1e-12*T,0,&uA0,&vA);
  route_a_val(&spec,hp+(size_t)i*spec.np,cacc[i],arow,du,T,0,&uA,&vA);
  res = 0.0;			/* elastic (coseismic) jump from the fit */
  for(it=0;it < spec.nterm;it++)
    res += arow[it];
  mean_vA = ((uA - uA0) + res*du)/T;
  route_b_weights(nib,gamma0,rate_b,1e-12*T,T,K,0,swv,swdv);
  uB0 = 0.0;
  for(n=0;n <= nib;n++)uB0 += du*swv[n]*xn_d[(size_t)n*nrec_d+i];
  route_b_weights(nib,gamma0,rate_b,T,T,K,0,swv,swdv);
  uB = 0.0;
  for(n=0;n <= nib;n++)uB += du*swv[n]*xn_d[(size_t)n*nrec_d+i];
  {				/* exact elastic jump: Gamma0^n weights */
    COMP_PRECISION wf = 1.0,el = 0.0;
    for(n=0;n <= nib;n++){
      el += wf*xn_d[(size_t)n*nrec_d+i];
      wf *= gamma0;
      if(fabs(wf) < 1e-16)break;
    }
    mean_vB = ((uB - uB0) + el*du)/T;
  }
  fprintf(stderr,"ve_sp_cycle: cycle-mean velocity at x/H %g: A %.6f B %.6f (target %.6f)\n",
	  xobs[i]/plate_h,mean_vA/Vpl,mean_vB/Vpl,0.5);
  /* (c) stress history over the final two cycles */
  for(k=0;k < VESP_NHIST;k++){
    tt = -1.0 + 2.0*((COMP_PRECISION)k + 0.5)/((COMP_PRECISION)VESP_NHIST);
    t = (tt < 0.0)?((tt+1.0)*T):(tt*T);
    fprintf(fs,"%10.5f",tt);
    for(j=0;j < VESP_NDEPTH;j++){
      irec = nrec_d + idep[j];
      arow = amp_s + (size_t)spec.nterm*idep[j];
      if(tt >= 0.0){
	route_a_val(&spec,hp+(size_t)irec*spec.np,cacc[irec],arow,du,t,0,
		    &sA,&vA);
	route_b_weights(nib,gamma0,rate_b,t,T,K,0,swv,swdv);
      }else{			/* previous cycle: drop the last event,
				   ages then count from the (K-1)-th */
	route_a_val(&spec,hp+(size_t)irec*spec.np,cacc[irec],arow,du,t,1,
		    &sA,&vA);
	route_b_weights(nib,gamma0,rate_b,t,T,K,1,swv,swdv);
      }
      sB = 0.0;			/* transient weights W_n - 1; the
					   event count cancels the exact
					   zero relaxed stress */
      {
	COMP_PRECISION nev = (COMP_PRECISION)((tt >= 0.0)?K:(K-1));
	for(n=0;n <= nib;n++)
	  sB += du*(swv[n] - nev)*xn_s[(size_t)n*nrec_s+idep[j]];
      }
      fprintf(fs," %14.6e %14.6e",sA,sB);
      if(fabs(sB) > sscl)sscl = fabs(sB);
      dev = fabs(sA-sB);if(dev > worst_ab_s)worst_ab_s = dev;
    }
    fprintf(fs,"\n");
  }
  /* (d) surface displacement history at x/H = 1 and 3 */
  for(j=0;j < 2;j++){
    COMP_PRECISION xt = (j == 0)?(1.0):(3.0);
    io[j] = 0;res = 1e30;
    for(i=0;i < nobs;i++){
      dev = fabs(xobs[i]/plate_h - xt);
      if(dev < res){res = dev;io[j] = i;}
    }
  }
  {
    COMP_PRECISION uref_a[2],uref_b[2];
    for(j=0;j < 2;j++){		/* u(0+) references */
      arow = amp_d + (size_t)spec.nterm*io[j];
      route_a_val(&spec,hp+(size_t)io[j]*spec.np,cacc[io[j]],arow,du,
		  1e-12*T,0,&uref_a[j],&vA);
      route_b_weights(nib,gamma0,rate_b,1e-12*T,T,K,0,swv,swdv);
      uref_b[j] = 0.0;
      for(n=0;n <= nib;n++)
	uref_b[j] += du*swv[n]*xn_d[(size_t)n*nrec_d+io[j]];
    }
    for(k=0;k < VESP_NHIST/2;k++){
      t = T*((COMP_PRECISION)k + 0.5)/((COMP_PRECISION)(VESP_NHIST/2));
      fprintf(fh,"%10.5f",t/T);
      for(j=0;j < 2;j++){
	arow = amp_d + (size_t)spec.nterm*io[j];
	route_a_val(&spec,hp+(size_t)io[j]*spec.np,cacc[io[j]],arow,du,t,0,
		    &uA,&vA);
	route_b_weights(nib,gamma0,rate_b,t,T,K,0,swv,swdv);
	uB = 0.0;
	for(n=0;n <= nib;n++)
	  uB += du*swv[n]*xn_d[(size_t)n*nrec_d+io[j]];
	fprintf(fh," %14.6e %14.6e",(uA-uref_a[j])/du,(uB-uref_b[j])/du);
      }
      fprintf(fh,"\n");
    }
  }
  fclose(fv);fclose(fd);fclose(fs);fclose(fh);
  /* summary and gates */
  fprintf(fc,"# ve_sp_cycle: H %g Gamma0 %g b %g T/taur %g np %i n_img %i K %i patches %i\n",
	  plate_h,gamma0,rate_b,T/taur,np,n_img,K,nrflt);
  fprintf(fc,"spinup_cycles %i\n",K);
  fprintf(fc,"depth_0 %g\n",-fault[idep[0]].x[INT_Y]);
  fprintf(fc,"depth_1 %g\n",-fault[idep[1]].x[INT_Y]);
  fprintf(fc,"depth_2 %g\n",-fault[idep[2]].x[INT_Y]);
  fprintf(fc,"xobs_1 %g\n",xobs[io[0]]);
  fprintf(fc,"xobs_3 %g\n",xobs[io[1]]);
  fprintf(fc,"n_img %i\n",nib);
  fprintf(fc,"held_out_residual %e\n",worst_res);
  fprintf(fc,"mean_velocity_A %e\n",mean_vA/Vpl);
  fprintf(fc,"mean_velocity_B %e\n",mean_vB/Vpl);
  fprintf(fc,"worst_AB_velocity_dev_Vpl %e\n",worst_ab_v);
  fprintf(fc,"worst_AB_velocityfit_dev_Vpl %e\n",worst_ab_v2);
  fprintf(fc,"worst_AB_displacement_dev_du %e\n",worst_ab_u);
  fprintf(fc,"worst_AB_stress_dev_rel %e\n",worst_ab_s/((sscl > 0.0)?sscl:1.0));
  fclose(fc);
  if(sscl > 0.0)worst_ab_s /= sscl;
  fprintf(stderr,"ve_sp_cycle: A vs B worst dev: velocity(h-state) %.2e / velocity(fit, unused for gate) %.2e V_pl, displacement %.2e du, stress %.2e of cycle swing\n",
	  worst_ab_v,worst_ab_v2,worst_ab_u,worst_ab_s);
  if((worst_res > 1e-3)||(worst_ab_u > 5e-3)||(worst_ab_v > 1e-2)||
     (worst_ab_s > 2e-2)||
     (fabs(mean_vB/Vpl - 0.5) > 3e-3)||(fabs(mean_vA/Vpl - 0.5) > 5e-3)){
    fprintf(stderr,"ve_sp_cycle: FAIL (see gates at the end of main)\n");
    return -1;
  }
  fprintf(stderr,"ve_sp_cycle: PASS\n");
  return 0;
}
