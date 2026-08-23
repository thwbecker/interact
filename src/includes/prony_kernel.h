/*

  visco-elastic Prony kernel machinery for interact/rsf_solve

  representation of the step-response interaction kernel of a Maxwell
  (in shear, elastic in bulk) medium as

     R(t) = C_const + sum_{p=1..np} C_p exp(-t/tau_p)

  where the amplitudes C are obtained from samples of the existing
  elastic Green's functions evaluated at the effective moduli of the
  correspondence principle,

     G_bar(s) = G s/(s + 1/t_m),
     nu_bar(s) = (3 K_b - 2 G_bar)/(6 K_b + 2 G_bar),  K_b fixed

  for the homogeneous half or full space the pole set is material
  only (see rsf_ve_design.md): T_m, T_m 3(1-nu)/(1+nu), and
  T_m 3/(2(1+nu)), all simple, and C_const vanishes for stress
  kernels; the constant term is carried anyway so that layered
  kernels (relaxed part nonzero) and displacement kernels (no T_m
  pole, nonzero relaxed part) use the same machinery

  part of interact

*/
#ifndef INT_PRONY_KERNEL_H
#define INT_PRONY_KERNEL_H

#define VE_MAX_NP 6		/* max number of exponential terms */
#define VE_MAX_NS (VE_MAX_NP+1+4)	/* max sample points incl fit + held out */
#define VE_NHELD 2		/* held-out sample points for verification */

struct prony_spec{
  int np;			/* number of exponential terms */
  my_boolean has_const;		/* carry a constant (relaxed) term */
  int nterm;			/* np + (has_const?1:0), amplitude count */
  COMP_PRECISION tau[VE_MAX_NP]; /* relaxation times, same units as T_m */
  int ns;			/* total sample points, nterm + VE_NHELD */
  COMP_PRECISION sk[VE_MAX_NS];	/* Laplace sample points, > 0 */
  /* weights: amplitude[iterm] = sum_k W[iterm][k] * sample[k], over
     the first nterm samples */
  COMP_PRECISION W[VE_MAX_NP+1][VE_MAX_NS];
  /* reference material */
  COMP_PRECISION g0,nu0,bulk0,t_M;
  /* effective elastic parameters at each sample point */
  struct el_par ep[VE_MAX_NS];
};

/* prony_kernel.c; hand-maintained prototypes (src/ve is not scanned
   by the cproto rule) */
void ve_effective_elpar(COMP_PRECISION, COMP_PRECISION, COMP_PRECISION,
			COMP_PRECISION, struct el_par *);
COMP_PRECISION ve_basis(struct prony_spec *, int, COMP_PRECISION);
void ve_spec_homogeneous(struct prony_spec *, COMP_PRECISION,
			 COMP_PRECISION, COMP_PRECISION);
COMP_PRECISION ve_prony_amplitudes_stress(struct prony_spec *,
					  struct med *, struct flt *,
					  int, int, COMP_PRECISION *,
					  COMP_PRECISION [VE_MAX_NP+1][3][3]);
COMP_PRECISION ve_prony_amplitudes_disp(struct prony_spec *,
					struct med *, struct flt *,
					COMP_PRECISION *, int,
					COMP_PRECISION *,
					COMP_PRECISION [VE_MAX_NP+1][3]);
COMP_PRECISION ve_basis_time_step(struct prony_spec *, int,
				  COMP_PRECISION);
COMP_PRECISION ve_basis_time_ramp(struct prony_spec *, int,
				  COMP_PRECISION, COMP_PRECISION);
void ve_solve_weights(struct prony_spec *);
/* prony_layered.c */
void ve_layered_antiplane_sample(struct med *, struct flt *, int,
				 COMP_PRECISION *, COMP_PRECISION *,
				 COMP_PRECISION, int, COMP_PRECISION,
				 MODE_TYPE, COMP_PRECISION *,
				 COMP_PRECISION [3][3]);
void ve_spec_layered_antiplane(struct prony_spec *, COMP_PRECISION,
			       COMP_PRECISION, COMP_PRECISION, int);

#endif
