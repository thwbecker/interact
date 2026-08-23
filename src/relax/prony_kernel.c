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
