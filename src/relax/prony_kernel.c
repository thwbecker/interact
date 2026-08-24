/*

  visco-elastic Prony kernel machinery, core routines

  see src/includes/prony_kernel.h and rsf_ve_design.md /
  rsf_ve_implementation_plan.md for background

  the central object is a prony_spec holding relaxation times tau_p,
  Laplace sample points s_k, effective elastic parameters at each
  sample point, and the weight matrix W that maps kernel samples to
  exponential amplitudes; ve_prony_amplitudes_* evaluate the existing
  elastic Green's functions at the sampled effective moduli and
  return per-term amplitude sets plus a held-out verification
  residual

  conventions: t_M and the tau_p share one time unit (the caller's
  choice); Laplace sample points are in the inverse of that unit; the
  basis function of an exponential term at sample s is s/(s + 1/tau_p)
  and that of the constant (relaxed) term is 1

  part of interact

*/
#include "interact.h"
#include "properties.h"

/* 

   effective elastic parameters at Laplace variable s for a shear
   Maxwell / elastic bulk material with reference (g0, nu0) and
   Maxwell time t_M

*/
void ve_effective_elpar(COMP_PRECISION g0, COMP_PRECISION nu0,
			COMP_PRECISION t_M, COMP_PRECISION s,
			struct el_par *ep)
{
  COMP_PRECISION bulk,gb,nub;
  bulk = bulk_mod_from_G_nu(g0,nu0);
  /* scaled shear moulus */
  gb  = g0 * s/(s + 1.0/t_M);
  /* poisson from scaled G and bulk */
  nub =  nu_from_G_bulk(gb,bulk);
  calc_medium_elastic_parameters(ep,gb,nub);
}

/* 

   fill a prony_spec for the homogeneous shear Maxwell / elastic bulk
   material: three simple poles at the material rates (see
   rsf_ve_design.md), plus the constant term for generality, sample
   points spread geometrically around the rates, VE_NHELD extra
   points held out for verification

*/
void ve_spec_homogeneous(struct prony_spec *spec,
			 COMP_PRECISION g0, COMP_PRECISION nu0,
			 COMP_PRECISION t_M)
{
  int i;
  spec->g0 = g0;
  spec->nu0 = nu0;
  spec->t_M = t_M;
  /*  */
  spec->bulk0 =  bulk_mod_from_G_nu(g0,nu0);
  spec->np = 3;
  spec->tau[0] = t_M;
  spec->tau[1] = t_M * 3.0*(1.0 - nu0)/(1.0 + nu0);
  spec->tau[2] = t_M * 3.0/(2.0*(1.0 + nu0));
  spec->has_const = TRUE;
  spec->nterm = spec->np + 1;
  spec->ns = spec->nterm + VE_NHELD;
  /* fit samples bracketing the rates, held-out samples interleaved */
  spec->sk[0] = 0.05/t_M;
  spec->sk[1] = 0.40/t_M;
  spec->sk[2] = 1.30/t_M;
  spec->sk[3] = 5.00/t_M;
  spec->sk[4] = 0.15/t_M;	/* held out */
  spec->sk[5] = 2.50/t_M;	/* held out */
  for(i=0;i < spec->ns;i++)
    ve_effective_elpar(g0,nu0,t_M,spec->sk[i],&spec->ep[i]);
  ve_solve_weights(spec);
}
/* 

   basis function of term iterm (exponential or constant) at Laplace
   sample s

*/
COMP_PRECISION ve_basis(struct prony_spec *spec, int iterm,
			COMP_PRECISION s)
{
  if(iterm < spec->np)
    return s/(s + 1.0/spec->tau[iterm]);
  else				/* constant (relaxed) term */
    return 1.0;
}
/* 

   solve the small nterm x nterm system for the weight matrix W such
   that amplitude[i] = sum_k W[i][k] sample[k] reproduces the basis
   at the first nterm sample points; plain Gauss elimination with
   partial pivoting on the transposed basis matrix

*/
void ve_solve_weights(struct prony_spec *spec)
{
  COMP_PRECISION a[VE_MAX_NP+1][2*(VE_MAX_NP+1)],fac,tmp;
  int n,i,j,k,ip;
  n = spec->nterm;
  /* augmented [M | I], M[k][i] = basis_i(s_k); we need W = M^-1 */
  for(k=0;k < n;k++){
    for(i=0;i < n;i++)
      a[k][i] = ve_basis(spec,i,spec->sk[k]);
    for(i=0;i < n;i++)
      a[k][n+i] = (k == i)?(1.0):(0.0);
  }
  for(k=0;k < n;k++){
    ip = k;
    for(i=k+1;i < n;i++)
      if(fabs(a[i][k]) > fabs(a[ip][k]))
	ip = i;
    if(fabs(a[ip][k]) < 1e-14){
      fprintf(stderr,"ve_solve_weights: singular basis matrix, check sample points\n");
      exit(-1);
    }
    if(ip != k)
      for(j=0;j < 2*n;j++){
	tmp = a[k][j];a[k][j] = a[ip][j];a[ip][j] = tmp;
      }
    fac = a[k][k];
    for(j=0;j < 2*n;j++)
      a[k][j] /= fac;
    for(i=0;i < n;i++){
      if(i == k)continue;
      fac = a[i][k];
      for(j=0;j < 2*n;j++)
	a[i][j] -= fac * a[k][j];
    }
  }
  /* W[i][k] = (M^-1)[i][k] */
  for(i=0;i < n;i++)
    for(k=0;k < n;k++)
      spec->W[i][k] = a[i][n+k];
}

/* 

   amplitude sets for the stress kernel of one source-receiver patch
   pair: evaluate the elastic Green's function at all sampled
   effective moduli, form amplitudes C[iterm][3][3] via W, and return
   the maximum relative held-out residual over the tensor components
   (relative to the largest sampled magnitude of that component; the
   caller decides what residual level to accept)

*/
COMP_PRECISION ve_prony_amplitudes_stress(struct prony_spec *spec,
					  struct med *medium,
					  struct flt *fault,
					  int irec, int isrc,
					  COMP_PRECISION *slip,
					  COMP_PRECISION C[VE_MAX_NP+1][3][3])
{
  COMP_PRECISION smp[VE_MAX_NS][3][3],u[3],model,scl,res,resmax;
  int k,it,i,j,iret;
  for(k=0;k < spec->ns;k++){
    eval_green_at_receiver(fault,irec,isrc,slip,u,smp[k],&iret,
			   GC_STRESS_ONLY,TRUE,medium->full_space,
			   spec->ep[k]);
    if(iret)
      fprintf(stderr,"ve_prony_amplitudes_stress: WARNING: singular %i <- %i sample %i\n",
	      irec,isrc,k);
  }
  for(it=0;it < spec->nterm;it++)
    for(i=0;i < 3;i++)
      for(j=0;j < 3;j++){
	C[it][i][j] = 0.0;
	for(k=0;k < spec->nterm;k++)
	  C[it][i][j] += spec->W[it][k] * smp[k][i][j];
      }
  /* held-out verification */
  resmax = 0.0;
  for(i=0;i < 3;i++)
    for(j=0;j < 3;j++){
      scl = 0.0;
      for(k=0;k < spec->ns;k++)
	if(fabs(smp[k][i][j]) > scl)
	  scl = fabs(smp[k][i][j]);
      if(scl < EPS_COMP_PREC)
	continue;
      for(k=spec->nterm;k < spec->ns;k++){
	model = 0.0;
	for(it=0;it < spec->nterm;it++){
	  model += C[it][i][j] * ve_basis(spec,it,spec->sk[k]);
	}
	res = fabs(model - smp[k][i][j])/scl;
	if(res > resmax)
	  resmax = res;
      }
    }
  return resmax;
}
/* 

   amplitude sets for the displacement kernel at an arbitrary
   observation point x (not a patch): D[iterm][3] displacement
   amplitudes, return value the held-out residual as above; for the
   homogeneous medium the t_M amplitude is expected to vanish (G
   cancels in displacements from dislocations) and the constant term
   carries the permanent (relaxed) deformation

*/
COMP_PRECISION ve_prony_amplitudes_disp(struct prony_spec *spec,
					struct med *medium,
					struct flt *fault,
					COMP_PRECISION *x, int isrc,
					COMP_PRECISION *slip,
					COMP_PRECISION D[VE_MAX_NP+1][3])
{
  COMP_PRECISION up[VE_MAX_NS][3],sm[3][3],model,scl,res,resmax;
  int k,it,i,iret,nsing;
  nsing = 0;
  for(k=0;k < spec->ns;k++){
    eval_green(x,(fault+isrc),slip,up[k],sm,&iret,
	       GC_DISP_ONLY,FALSE,medium->full_space,spec->ep[k]);
    if(iret)
      nsing++;
  }
  if(nsing)
    fprintf(stderr,"ve_prony_amplitudes_disp: WARNING: singular obs point (%g, %g, %g), source %i, %i of %i samples\n",
	    x[INT_X],x[INT_Y],x[INT_Z],isrc,nsing,spec->ns);
  for(it=0;it < spec->nterm;it++)
    for(i=0;i < 3;i++){
      D[it][i] = 0.0;
      for(k=0;k < spec->nterm;k++)
	D[it][i] += spec->W[it][k] * up[k][i];
    }
  resmax = 0.0;
  for(i=0;i < 3;i++){
    scl = 0.0;
    for(k=0;k < spec->ns;k++)
      if(fabs(up[k][i]) > scl)
	scl = fabs(up[k][i]);
    if(scl < EPS_COMP_PREC)
      continue;
    for(k=spec->nterm;k < spec->ns;k++){
      model = 0.0;
      for(it=0;it < spec->nterm;it++){
	model += D[it][i] * ve_basis(spec,it,spec->sk[k]);
      }
      res = fabs(model - up[k][i])/scl;
      if(res > resmax)
	resmax = res;
    }
  }
  return resmax;
}
/* 

   time-domain basis of term iterm for a unit slip step at t = 0

*/
COMP_PRECISION ve_basis_time_step(struct prony_spec *spec, int iterm,
				  COMP_PRECISION t)
{
  if(iterm < spec->np)
    return exp(-t/spec->tau[iterm]);
  else
    return 1.0;
}
/* 

   time-domain basis of term iterm for a unit slip ramp over
   [0, t_ramp] (constant velocity 1/t_ramp, zero after); reduces to
   the step basis in the t_ramp -> 0 limit

*/
COMP_PRECISION ve_basis_time_ramp(struct prony_spec *spec, int iterm,
				  COMP_PRECISION t, COMP_PRECISION t_ramp)
{
  COMP_PRECISION tp;
  if(t_ramp <= 0.0)
    return ve_basis_time_step(spec,iterm,t);
  if(iterm < spec->np){
    tp = spec->tau[iterm];
    if(t <= t_ramp)
      return (tp/t_ramp)*(1.0 - exp(-t/tp));
    else
      return (tp/t_ramp)*(exp(-(t - t_ramp)/tp) - exp(-t/tp));
  }else{			/* constant term follows the slip */
    return (t <= t_ramp)?(t/t_ramp):(1.0);
  }
}

/*

  step 7a of rsf_ve_implementation_plan.md: antiplane elastic plate
  (thickness H, free surface at y = 0) over a Maxwell half space of
  equal rigidity, via the Nur-Mavko / Savage-Prescott image series

  in the Laplace domain the interface reflection coefficient is

    beta(s) = (G - G_bar(s))/(G + G_bar(s)) = a/(s + a),  a = 1/(2 t_M)

  so the s-domain kernel is the elastic antiplane kernel of the
  source plus its images: free-surface images carry +1, each
  interaction with the substrate multiplies by beta(s), giving image
  families at depth shifts 2 n H with weight beta^n; the series is
  geometric in beta and truncated at n_img with the tail bounded by
  the last retained term

  the time-domain step response of the n-th family is the Erlang
  ramp 1 - exp(-a t) sum_{k<n} (a t)^k / k!, i.e. polynomial times
  exponential; the pure-exponential Prony state machinery therefore
  FITS this kernel with rates on a ladder around a, and the held-out
  residual measures the truncation, which is the quantity step 7a
  exists to establish

  limits used as anchors: s -> inf gives beta -> 0, the homogeneous
  half-space kernel (instantaneously the equal-rigidity interface is
  invisible); s -> 0 gives beta -> 1, the plate with a free bottom
  (fluid substrate), evaluable by direct image summation with unit
  weights

  geometry convention as in the 2-D (w = 0) interact elements: y is
  depth (y <= 0 in the plate down to y = -H), antiplane displacement
  and stress from STRIKE-mode slip; images of a segment at depth y_s
  sit at -y_s (free surface) and at y_s -/+ 2 n H and -y_s -/+ 2 n H

  part of interact

*/
/* 

   s-domain layered antiplane kernel sample for one source-receiver
   pair: evaluate the elastic kernel (at the REFERENCE moduli, since
   the plate stays elastic) for the source and its images with
   beta(s)^n weights; mode GC_STRESS_ONLY or GC_DISP_ONLY; x is the
   receiver point (for stress on a patch, pass the patch center and
   resolve outside); n_img image families are summed

*/
void ve_layered_antiplane_sample(struct med *medium, struct flt *fault,
				 int isrc, COMP_PRECISION *slip,
				 COMP_PRECISION *x,
				 COMP_PRECISION beta, int n_img,
				 COMP_PRECISION plate_h,
				 MODE_TYPE mode,
				 COMP_PRECISION *u, COMP_PRECISION sm[3][3])
{
  struct flt img[1];
  COMP_PRECISION us[3],sms[3][3],wfac,ysave;
  int n,iret,is,i,j;
  for(i=0;i<3;i++){
    u[i] = 0.0;
    for(j=0;j<3;j++)
      sm[i][j] = 0.0;
  }
  /* copy the source patch so its depth can be shifted; the elastic
     evaluation uses the medium's reference parameters */
  *img = fault[isrc];
  ysave = fault[isrc].x[INT_Y];
  for(n=0;n <= n_img;n++){
    wfac = 1.0;
    for(is=0;is < n;is++)wfac *= beta;
    if(fabs(wfac) < EPS_COMP_PREC)
      break;
    /* family n: depths ysave - 2 n H and -ysave - 2 n H (the
       free-surface image of the shifted source is generated by the
       half-plane elastic kernel itself, so only the downward shifts
       are summed explicitly) */
    img->x[INT_Y] = ysave - 2.0*(COMP_PRECISION)n*plate_h;
    eval_green(x,img,slip,us,sms,&iret,mode,FALSE,FALSE,
	       medium->elastic);
    for(i=0;i<3;i++){
      u[i] += wfac * us[i];
      for(j=0;j<3;j++)sm[i][j] += wfac * sms[i][j];
    }
  }
}
/* 

   fill a prony_spec for the layered antiplane case: rates on a
   geometric ladder around a = 1/(2 t_M) (the image-series pole), np
   terms plus the constant (relaxed) term, sample points bracketing
   the rates, weights by the same small solve as the homogeneous
   case; the caller checks the held-out residual, which here is a
   FIT quality, not an exactness statement

*/
void ve_spec_layered_antiplane(struct prony_spec *spec,
			       COMP_PRECISION g0, COMP_PRECISION nu0,
			       COMP_PRECISION t_M, COMP_PRECISION rate_b,
			       int np)
{
  COMP_PRECISION a;
  int i;
  if((np < 1)||(np > VE_MAX_NP)){
    fprintf(stderr,"ve_spec_layered_antiplane: np %i out of range\n",np);
    exit(-1);
  }
  spec->g0 = g0;spec->nu0 = nu0;spec->t_M = t_M;
  spec->bulk0 = bulk_mod_from_G_nu(g0,nu0);
  a = rate_b;			/* layered relaxation rate, 1/(2 t_M)
				   for equal rigidities (see
				   ve_layered_gamma_pars) */
  spec->np = np;
  /* 
     geometric ladder of rates spanning [a/15, 1.5a]: the n-th image
     family relaxes on the timescale n/a (Erlang ramp of order n), so
     with of order ten significant images the relaxation spectrum
     extends well below a and a ladder centered on a misses the slow
     tail (which the time-domain check against Nur-Mavko exposes even
     when the held-out residual looks small). NOTE: the ladder floor
     and the image truncation n_img are coupled: images beyond n
     relax on timescales n/a, so raising n_img past roughly 4 times
     the ladder floor ratio (n_img above about 60 for the b/15
     floor here) adds content the fit cannot represent and the
     time-domain accuracy degrades even as the held-out residual
     holds; pair n_img of order 40 to 80 with this ladder
  */
  for(i=0;i < np;i++)
    spec->tau[i] = 1.0/(a * (1.0/15.0) * pow(22.5,
				(COMP_PRECISION)i/
				((np > 1)?((COMP_PRECISION)(np-1)):(1.0))));
  spec->has_const = TRUE;
  spec->nterm = np + 1;
  spec->ns = spec->nterm + VE_NHELD;
  /* fit samples bracketing the ladder, held-out in the interior */
  for(i=0;i < spec->nterm;i++)
    spec->sk[i] = a * (1.0/60.0) * pow(240.0,
				 (COMP_PRECISION)i/
				 ((COMP_PRECISION)spec->nterm - 1.0));
  spec->sk[spec->nterm]   = a * (1.0/25.0);
  spec->sk[spec->nterm+1] = a * 0.80;
  ve_solve_weights(spec);
}
/*

  general rigidity contrast for the layered antiplane problem: with
  plate rigidity g1 (elastic) and substrate rigidity g2 (Maxwell,
  relaxation time t_M), the correspondence-principle interface
  reflection coefficient is

    Gamma(s) = (g1 - g2 s/(s + 1/t_M)) / (g1 + g2 s/(s + 1/t_M))
             = Gamma_0 + (1 - Gamma_0) b/(s + b)

  with the elastic contrast Gamma_0 = (g1 - g2)/(g1 + g2) and the
  layered relaxation rate b = g1/((g1 + g2) t_M); equal rigidities
  give Gamma_0 = 0 and b = 1/(2 t_M), the case treated above

*/
void ve_layered_gamma_pars(COMP_PRECISION g1, COMP_PRECISION g2,
			   COMP_PRECISION t_M,
			   COMP_PRECISION *gamma0, COMP_PRECISION *b)
{
  *gamma0 = (g1 - g2)/(g1 + g2);
  *b = g1/((g1 + g2) * t_M);
}
COMP_PRECISION ve_layered_gamma_of_s(COMP_PRECISION gamma0,
				     COMP_PRECISION b,
				     COMP_PRECISION s)
{
  return gamma0 + (1.0 - gamma0) * b/(s + b);
}
/*

  Rybicki (1971) style elastic image series for the SURFACE
  displacement of a vertical screw fault segment in the plate:
  segment at horizontal position xs spanning depths [c1, c2]
  (0 <= c1 < c2 <= H, depths positive down), unit displacement
  discontinuity in the sense of the antiplane element's
  disp[STRIKE], plate thickness plate_h over a substrate with
  elastic reflection coefficient gam,

    u_z(x, 0) = (1/pi) sum_{n=0}^{n_img} gam^n
                [ atan((x - xs)/(c1 + 2 n H)) -
                  atan((x - xs)/(c2 + 2 n H)) ]

  (the n = 0, c1 = 0 limit atan(x/0+) = sgn(x) pi/2 reproduces the
  homogeneous surface-breaking profile validated against the
  antiplane element). gam = Gamma_0 gives the elastic layered
  solution, gam = 1 the fully relaxed (free bottom) limit, gam = 0
  the homogeneous half space

*/
COMP_PRECISION ve_rybicki_antiplane_uz(COMP_PRECISION x, COMP_PRECISION xs,
				       COMP_PRECISION c1, COMP_PRECISION c2,
				       COMP_PRECISION plate_h,
				       COMP_PRECISION gam, int n_img)
{
  COMP_PRECISION uz,wfac,dx,d1,d2;
  int n;
  uz = 0.0;wfac = 1.0;
  dx = x - xs;
  for(n=0;n <= n_img;n++){
    if(fabs(wfac) < 1e-16)break;
    d1 = c1 + 2.0*(COMP_PRECISION)n*plate_h;
    d2 = c2 + 2.0*(COMP_PRECISION)n*plate_h;
    uz += wfac * (atan2(dx,d1) - atan2(dx,d2));
    wfac *= gam;
  }
  return uz/PI;
}
/*

  Nur and Mavko (1974) style time-domain weight of the n-th image
  family after a step at t = 0: the inverse transform of
  Gamma(s)^n / s, by binomial expansion

    W_n(t) = sum_{m=0}^{n} C(n,m) Gamma_0^{n-m} (1 - Gamma_0)^m
             P(m, b t)

  with P(0,x) = 1 and P(m,x) = 1 - exp(-x) sum_{k<m} x^k/k! the
  regularized lower incomplete gamma (Erlang ramp). W_n(0) =
  Gamma_0^n (elastic) and W_n(inf) = 1 (fully relaxed). for
  Gamma_0 = 0 this collapses to the single term P(n, b t) with
  b = 1/(2 t_M). NOTE on conditioning: for negative Gamma_0 (stiff
  substrate) the binomial sum alternates and loses roughly
  n log10(1+2|Gamma_0|/(1-|Gamma_0|)) digits, so keep n_img
  moderate there; the physically relevant weights fall off as
  |Gamma|^n anyway

*/
COMP_PRECISION ve_nur_mavko_weight(int n, COMP_PRECISION gamma0,
				   COMP_PRECISION b, COMP_PRECISION t)
{
  COMP_PRECISION w,binom,pfac,x,eterm,psum,pm;
  int m,k;
  x = b * t;
  w = 0.0;
  binom = 1.0;			/* C(n,0) */
  for(m=0;m <= n;m++){
    /* pfac = Gamma_0^{n-m} (1-Gamma_0)^m */
    pfac = 1.0;
    for(k=0;k < n-m;k++)
      pfac *= gamma0;
    for(k=0;k < m;k++)
      pfac *= (1.0 - gamma0);
    if(m == 0){
      pm = 1.0;
    }else{
      /* P(m, x) = 1 - exp(-x) sum_{k<m} x^k/k! */
      eterm = exp(-x);psum = 0.0;
      for(k=0;k < m;k++){
	psum += eterm;
	eterm *= x/((COMP_PRECISION)(k+1));
      }
      pm = 1.0 - psum;
      if(pm < 0.0)
	pm = 0.0;	/* guard rounding at x -> 0 */
    }
    w += binom * pfac * pm;
    binom *= ((COMP_PRECISION)(n - m))/((COMP_PRECISION)(m + 1));
  }
  return w;
}
/*

  full time-domain Nur-Mavko surface displacement for the segment of
  ve_rybicki_antiplane_uz: source term plus W_n(t)-weighted images

*/
COMP_PRECISION ve_nur_mavko_uz(COMP_PRECISION x, COMP_PRECISION xs,
			       COMP_PRECISION c1, COMP_PRECISION c2,
			       COMP_PRECISION plate_h,
			       COMP_PRECISION gamma0, COMP_PRECISION b,
			       COMP_PRECISION t, int n_img)
{
  COMP_PRECISION uz,dx,d1,d2;
  int n;
  dx = x - xs;
  uz = 0.0;
  for(n=0;n <= n_img;n++){
    d1 = c1 + 2.0*(COMP_PRECISION)n*plate_h;
    d2 = c2 + 2.0*(COMP_PRECISION)n*plate_h;
    uz += ((n == 0)?(1.0):(ve_nur_mavko_weight(n,gamma0,b,t))) *
      (atan2(dx,d1) - atan2(dx,d2));
  }
  return uz/PI;
}
/*

  time derivative of the Nur-Mavko image weight: with W_n as in
  ve_nur_mavko_weight, 

    Wdot_n(t) = sum_{m=1}^{n} C(n,m) Gamma_0^{n-m} (1 - Gamma_0)^m
                b exp(-b t) (b t)^{m-1} / (m-1)!

  (the Erlang density of order m); the same conditioning caveat for
  negative Gamma_0 applies

*/
COMP_PRECISION ve_nur_mavko_weight_dot(int n, COMP_PRECISION gamma0,
				       COMP_PRECISION b, COMP_PRECISION t)
{
  COMP_PRECISION w,binom,pfac,x,dens;
  int m,k;
  x = b * t;
  w = 0.0;
  binom = (COMP_PRECISION)n;	/* C(n,1) */
  dens = b * exp(-x);		/* order-1 Erlang density */
  for(m=1;m <= n;m++){
    pfac = 1.0;
    for(k=0;k < n-m;k++)pfac *= gamma0;
    for(k=0;k < m;k++)pfac *= (1.0 - gamma0);
    w += binom * pfac * dens;
    binom *= ((COMP_PRECISION)(n - m))/((COMP_PRECISION)(m + 1));
    dens *= x/((COMP_PRECISION)m);
  }
  return w;
}
/*

  time-domain Nur-Mavko surface velocity for the segment of
  ve_rybicki_antiplane_uz (the coseismic term is a step, so only the
  images contribute)

*/
COMP_PRECISION ve_nur_mavko_vz(COMP_PRECISION x, COMP_PRECISION xs,
			       COMP_PRECISION c1, COMP_PRECISION c2,
			       COMP_PRECISION plate_h,
			       COMP_PRECISION gamma0, COMP_PRECISION b,
			       COMP_PRECISION t, int n_img)
{
  COMP_PRECISION vz,dx,d1,d2;
  int n;
  dx = x - xs;
  vz = 0.0;
  for(n=1;n <= n_img;n++){
    d1 = c1 + 2.0*(COMP_PRECISION)n*plate_h;
    d2 = c2 + 2.0*(COMP_PRECISION)n*plate_h;
    vz += ve_nur_mavko_weight_dot(n,gamma0,b,t) *
      (atan2(dx,d1) - atan2(dx,d2));
  }
  return vz/PI;
}
