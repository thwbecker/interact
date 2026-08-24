/* pom.h
 *
 * C port of Kaj M. Johnson's Plate_over_Maxwell MATLAB codes:
 * semi-analytical propagator matrix solutions for postseismic
 * deformation due to slip on a rectangular dislocation in an
 * elastic plate over Maxwell viscoelastic layer(s).
 *
 * All units follow the original codes: km for lengths, years for
 * times, m for slip and displacement, m/yr for velocity.
 */
#ifndef POM_H
#define POM_H

#include <complex.h>

/* mode flags for the layer-over-halfspace solution */
#define POM_DISP 0		/* cumulative postseismic displacement */
#define POM_VEL  1		/* postseismic velocity */

/*
 * elastic plate over viscoelastic layer over viscoelastic halfspace
 *
 * m[10]: fault patch, Okada-style, as in the MATLAB codes:
 *   m[0] length along strike (km)      m[1] down-dip width (km)
 *   m[2] depth to down-dip edge (km)   m[3] dip (deg)   m[4] strike (deg)
 *   m[5] east, m[6] north position of center of down-dip edge (km)
 *   m[7] strike-slip (m, positive left-lateral)
 *   m[8] dip-slip (m, positive reverse)
 *   m[9] tensile (m, positive opening)
 * xs, ys: n station coordinates (km, east and north)
 * H1: elastic plate thickness (km)
 * H2: depth to bottom of viscoelastic layer (km)
 * nu: Poisson's ratio
 * t:  time since earthquake (years), single time
 * tR1, tR2: Maxwell relaxation times of layer and halfspace (years)
 * mode: POM_DISP or POM_VEL
 * nl, nw: number of Gauss-Legendre point sources along strike and dip;
 *   pass 0 to use the defaults of the original codes
 *   (ceil(0.05*L) x ceil(0.05*W) for displacements,
 *    ceil(0.025*L) x ceil(0.05*W) for velocities)
 * ue, un, uz: output arrays of length n (east, north, up)
 *
 * returns 0 on success, nonzero on error
 */
int pom_layer(const double *m, const double *xs, const double *ys, int n,
	      double H1, double H2, double nu, double t,
	      double tR1, double tR2, int mode, int nl, int nw,
	      double *ue, double *un, double *uz);

/*
 * elastic plate over Maxwell viscoelastic halfspace, patch represented
 * by a single point source (Plate_over_Maxwell.m)
 *
 * m[10] as above; mu_lam = mu/lambda (1 for nu = 0.25)
 * t: array of nt times; outputs are ue[i*nt+it] etc.
 *
 * returns 0 on success, nonzero on error
 */
int pom_maxwell(const double *m, const double *xs, const double *ys, int n,
		double H, double mu_lam, const double *t, int nt, double tR,
		double *ue, double *un, double *uz);
int pom_maxwell_z(const double *m, const double *xs, const double *ys, int n,
		  double H, double mu_lam, const double *t, int nt, double tR,
		  double zr, double *ue, double *un, double *uz,
		  double *szz, double *sxz, double *syz);

/* machine-translated coefficient functions (terms_2visc.c) */
void num1_2visc(double pf1, double pf2, double B1, double B2,
		double K, double k, double d1,
		double P0[4][4], double u, double *N);
void num2_2visc(double pf1, double pf2, double B1, double B2,
		double K, double k, double d1,
		double P0[4][4], double u, double *N);
void denom_2visc(double B1, double B2, double K, double k,
		 double d1, double P0[4][4], double u, double *D);
void antiplane_2visc_terms(double pf1, double pf2, double B1,
			   double B2, double k, double d1,
			   double P0[2][2], double u,
			   double *N, double *D);

/* utilities (pom_util.c) */
void pom_source_vectors(double M[3][3], double kj, double g, double lam,
		double mu, double complex F[5][4], double complex f2[5][2]);
void pom_gauleg(double a, double b, int n, double *x, double *w);
int pom_poly_roots(int deg, const double *c, double complex *r);
void pom_momtensor(double strike, double dip, double lam, double mu,
		   double M1[3][3], double M2[3][3], double M3[3][3]);
void pom_prop4(double kj, double dz, double mu, double lam, double g,
	       double P[4][4]);
void pom_prop2(double kj, double dz, double mu, double P[2][2]);

#endif
