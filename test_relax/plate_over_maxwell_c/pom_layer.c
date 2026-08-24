#define _GNU_SOURCE
/* pom_layer.c
 *
 * C port of Plate_over_Maxwell_Layer_over_Halfspace_Displacements.m and
 * Plate_over_Maxwell_Layer_over_Halfspace_Velocities.m
 * (Kaj M. Johnson, Indiana University)
 *
 * Postseismic surface deformation due to uniform slip on a rectangular
 * dislocation in an elastic plate over a Maxwell viscoelastic layer over
 * a Maxwell viscoelastic halfspace; semi-analytical propagator matrix
 * solution, inverted from the Laplace domain by residues over the roots
 * of the (degree 7 after removing s = 0) denominator polynomial, and
 * from the Hankel domain by trapezoidal integration over wavenumber.
 */
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include "pom.h"

#define NK 100			/* wavenumber samples (N in the m-files) */
#define KMIN 0.000001
#define KMAX 0.5		/* 1/km, as in the m-files */
#define NORD 5			/* Fourier orders 0,-1,+1,-2,+2 */
#define RG 3.0e-3		/* gravity correction, as in the m-files */

static const int order_m[NORD] = { 0, -1, 1, -2, 2 };

/* complex polynomial evaluation, ascending coefficients */
static double complex cpolyval(const double complex *c, int n, double complex s)
{
  int i;
  double complex p;

  p = c[n - 1];
  for (i = n - 2; i >= 0; i--)
    p = p * s + c[i];
  return p;
}

/*
 * per-wavenumber ingredients that do not depend on the source term:
 * numerator basis vectors, denominator roots, propagators
 */
struct kterms {
  double N1a[8], N1b[8], N2a[8], N2b[8];	/* spheroidal numerators */
  double complex R[7];		/* nonzero roots of the denominator */
  double a;			/* leading denominator coefficient */
  double NA1[3], NA2[3];	/* antiplane numerator basis (pf=(1,0),(0,1)) */
  double complex R2[2];		/* nonzero antiplane roots */
  double a2;
  double Pzs4[2][4];		/* rows 3,4 of source 4x4 propagator */
  double Pzs2[2][2];		/* source 2x2 propagator */
};

static int build_kterms(double kj, double H1, double H2, double zs,
			double mu, double lam, double g, double bulk,
			double B1, double B2, struct kterms *kt)
{
  int i, j;
  double P4[4][4], Ph4[4][4], Ph2[2][2], Pz4[4][4], Pz2[2][2];
  double D[8], DA[3], NAa[3];
  double d1, cs[8], scale;
  double complex disc;

  d1 = H1 - H2;			/* negative: layer below the plate */
  /* propagator surface -> plate bottom, with gravity correction G */
  pom_prop4(kj, -H1, mu, lam, g, P4);
  for (j = 0; j < 4; j++) {
    Ph4[0][j] = P4[0][j];
    Ph4[1][j] = P4[1][j];
    Ph4[2][j] = P4[2][j];
    Ph4[3][j] = P4[3][j] - RG * P4[1][j];
  }
  pom_prop2(kj, -H1, mu, Ph2);
  /* propagator surface -> source depth */
  pom_prop4(kj, -zs, mu, lam, g, Pz4);
  for (j = 0; j < 4; j++) {
    kt->Pzs4[0][j] = Pz4[2][j];
    kt->Pzs4[1][j] = Pz4[3][j] - RG * Pz4[1][j];
  }
  pom_prop2(kj, -zs, mu, Pz2);
  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      kt->Pzs2[i][j] = Pz2[i][j];

  /* spheroidal numerator bases and denominator */
  num1_2visc(1.0, 0.0, B1, B2, bulk, kj, d1, Ph4, mu, kt->N1a);
  num1_2visc(0.0, 1.0, B1, B2, bulk, kj, d1, Ph4, mu, kt->N1b);
  num2_2visc(1.0, 0.0, B1, B2, bulk, kj, d1, Ph4, mu, kt->N2a);
  num2_2visc(0.0, 1.0, B1, B2, bulk, kj, d1, Ph4, mu, kt->N2b);
  denom_2visc(B1, B2, bulk, kj, d1, Ph4, mu, D);
  kt->a = D[7];
  /* roots of D8 s^7 + ... + D1 = 0, scaled for conditioning */
  scale = 0.0;
  for (i = 0; i < 8; i++)
    if (fabs(D[i]) > scale)
      scale = fabs(D[i]);
  if (scale == 0.0 || D[7] == 0.0)
    return 1;
  for (i = 0; i < 8; i++)
    cs[i] = D[i] / scale;
  if (pom_poly_roots(7, cs, kt->R))
    return 1;

  /* antiplane basis: N depends only on pf2 in the generated terms,
     but keep both basis vectors for fidelity with the m-file */
  antiplane_2visc_terms(1.0, 0.0, B1, B2, kj, d1, Ph2, mu, NAa, DA);
  antiplane_2visc_terms(0.0, 1.0, B1, B2, kj, d1, Ph2, mu, kt->NA2, DA);
  for (i = 0; i < 3; i++)
    kt->NA1[i] = NAa[i];
  kt->a2 = DA[2];
  if (DA[2] == 0.0)
    return 1;
  /* roots of D3 s^2 + D2 s + D1 = 0 */
  disc = (double complex)(DA[1] * DA[1] - 4.0 * DA[2] * DA[0]);
  disc = csqrt(disc);
  kt->R2[0] = (-DA[1] + disc) / (2.0 * DA[2]);
  kt->R2[1] = (-DA[1] - disc) / (2.0 * DA[2]);
  return 0;
}

/* time functions: displacement (e^{st}-1), velocity s e^{st} */
static double complex timefun(double complex s, double t, int mode)
{
  if (mode == POM_VEL)
    return s * cexp(s * t);
  return cexp(s * t) - 1.0;
}

/* spheroidal residue sum for one source vector */
static void spheroidal(const struct kterms *kt, const double complex F[4],
		       double t, int mode,
		       double complex *US, double complex *UR)
{
  int i, l, mm;
  double complex Pzsf[2], Nu1[8], Nu2[8], den, u1, u2, s, tf;

  for (i = 0; i < 2; i++)
    Pzsf[i] = kt->Pzs4[i][0] * F[0] + kt->Pzs4[i][1] * F[1]
      + kt->Pzs4[i][2] * F[2] + kt->Pzs4[i][3] * F[3];
  for (i = 0; i < 8; i++) {
    Nu1[i] = kt->N1a[i] * Pzsf[0] + kt->N1b[i] * Pzsf[1];
    Nu2[i] = kt->N2a[i] * Pzsf[0] + kt->N2b[i] * Pzsf[1];
  }
  u1 = 0.0;
  u2 = 0.0;
  for (l = 0; l < 7; l++) {
    s = kt->R[l];
    den = kt->a * s;
    for (mm = 0; mm < 7; mm++)
      if (mm != l)
	den *= (s - kt->R[mm]);
    tf = timefun(s, t, mode);
    u1 += -cpolyval(Nu1, 8, s) * tf / den;
    u2 += -cpolyval(Nu2, 8, s) * tf / den;
  }
  *US = u1;
  *UR = -u2;
}

/* antiplane residue sum for one source vector */
static double complex toroidal(const struct kterms *kt,
			       const double complex f[2], double t, int mode)
{
  int l;
  double complex Pzsf[2], NN[3], s, den, num, sum;

  Pzsf[0] = kt->Pzs2[0][0] * f[0] + kt->Pzs2[0][1] * f[1];
  Pzsf[1] = kt->Pzs2[1][0] * f[0] + kt->Pzs2[1][1] * f[1];
  for (l = 0; l < 3; l++)
    NN[l] = kt->NA1[l] * Pzsf[0] + kt->NA2[l] * Pzsf[1];
  sum = 0.0;
  for (l = 0; l < 2; l++) {
    s = kt->R2[l];
    den = kt->a2 * s * (s - kt->R2[1 - l]);
    num = -(NN[0] + NN[1] * s + NN[2] * s * s);
    sum += num * timefun(s, t, mode) / den;
  }
  return -sum;
}


/* Bessel J1(x)/x and J2(x)/x with small-argument care */
static double besj2(double x, double b0, double b1)
{
  if (fabs(x) < 1e-6)
    return x * x / 8.0;
  return -b0 + 2.0 / x * b1;
}

static double besj3(double x, double b1, double b2)
{
  if (fabs(x) < 1e-6)
    return x * x * x / 48.0;
  return -b1 + 4.0 / x * b2;
}

int pom_layer(const double *m, const double *xs, const double *ys, int n,
	      double H1, double H2, double nu, double t,
	      double tR1, double tR2, int mode, int nl, int nw,
	      double *ue, double *un, double *uz)
{
  int i, j, iw, il, io, im, err;
  int NL, NW, mech_on[3];
  double mu, lam, g, bulk, B1, B2;
  double L, W, dip, strike, D, dipdir, offset, eoff, noff;
  double sd, cd, ss_, cs_, kj, dk, wgt;
  double slip[3], zs, x0, y0, xr, yr;
  double M1[3][3], M2[3][3], M3[3][3];
  double *xp, *yp, *wL, *wW, *Xc, *Yc, *k;
  double *xpos, *ypos;
  double bb0[NK], bb1[NK], bb2[NK], bb3[NK];
  double (*Mt)[3][3];
  double *Gx[3], *Gy[3], *Gz[3];
  struct kterms kt;
  double complex F[NORD][4], f2[NORD][2];
  double complex(*cUR)[NORD][NK], (*cUS)[NORD][NK], (*cUT)[NORD][NK];

  if (n < 1 || t < 0.0 || tR1 <= 0.0 || tR2 <= 0.0)
    return 1;
  mu = 1.0;
  lam = 2.0 * mu * nu / (1.0 - 2.0 * nu);
  g = lam + 2.0 * mu;
  bulk = lam + 2.0 * mu / 3.0;
  B1 = 1.0 / tR1;
  B2 = 1.0 / tR2;
  L = m[0];
  W = m[1];
  dip = m[3];
  strike = m[4];
  slip[0] = m[7];
  slip[1] = m[8];
  slip[2] = m[9];
  for (im = 0; im < 3; im++)
    mech_on[im] = (slip[im] != 0.0);

  /* defaults of the original codes */
  if (nl <= 0)
    NL = (int)ceil(((mode == POM_VEL) ? 0.025 : 0.05) * L);
  else
    NL = nl;
  if (nw <= 0)
    NW = (int)ceil(0.05 * W);
  else
    NW = nw;
  if (NL < 1)
    NL = 1;
  if (NW < 1)
    NW = 1;

  /* geometry, following the m-files (x and y swapped internally) */
  if (dip <= 90.0 && dip >= -90.0)
    dipdir = (strike + 90.0) * M_PI / 180.0;
  else
    dipdir = (strike - 90.0) * M_PI / 180.0;
  offset = fabs(W * cos(dip * M_PI / 180.0));
  eoff = offset * sin(dipdir);
  noff = offset * cos(dipdir);
  D = m[2] - W * sin(dip * M_PI / 180.0);

  Xc = malloc(sizeof(double) * n);
  Yc = malloc(sizeof(double) * n);
  xp = malloc(sizeof(double) * NL);
  wL = malloc(sizeof(double) * NL);
  yp = malloc(sizeof(double) * NW);
  wW = malloc(sizeof(double) * NW);
  k = malloc(sizeof(double) * NK);
  xpos = malloc(sizeof(double) * NL * NW);
  ypos = malloc(sizeof(double) * NL * NW);
  cUR = malloc(sizeof(*cUR) * 3);
  cUS = malloc(sizeof(*cUS) * 3);
  cUT = malloc(sizeof(*cUT) * 3);
  for (im = 0; im < 3; im++) {
    Gx[im] = calloc(n, sizeof(double));
    Gy[im] = calloc(n, sizeof(double));
    Gz[im] = calloc(n, sizeof(double));
  }
  if (!Xc || !Yc || !xp || !wL || !yp || !wW || !k || !xpos || !ypos
      || !cUR || !cUS || !cUT)
    return 1;

  for (i = 0; i < n; i++) {
    /* station coords relative to fault, with the x/y swap of the m-files */
    Xc[i] = ys[i] - m[6] + noff;
    Yc[i] = xs[i] - m[5] + eoff;
  }
  pom_gauleg(-L / 2.0, L / 2.0, NL, xp, wL);
  pom_gauleg(0.0, W, NW, yp, wW);
  sd = sin(dip * M_PI / 180.0);
  cd = cos(dip * M_PI / 180.0);
  ss_ = sin(strike * M_PI / 180.0);
  cs_ = cos(strike * M_PI / 180.0);
  for (iw = 0; iw < NW; iw++)
    for (il = 0; il < NL; il++) {
      x0 = xp[il];
      y0 = cd * yp[iw];		/* dip rotation about x-axis */
      xr = cs_ * x0 - ss_ * y0;	/* strike rotation about z-axis */
      yr = ss_ * x0 + cs_ * y0;
      xpos[iw * NL + il] = xr;
      ypos[iw * NL + il] = yr;
    }
  dk = (KMAX - KMIN) / (NK - 1);
  for (j = 0; j < NK; j++)
    k[j] = KMIN + dk * j;
  pom_momtensor(strike, dip, lam, mu, M1, M2, M3);
  Mt = malloc(sizeof(double) * 27);
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      Mt[0][i][j] = M1[i][j];
      Mt[1][i][j] = M2[i][j];
      Mt[2][i][j] = M3[i][j];
    }

  err = 0;
  for (iw = 0; iw < NW && !err; iw++) {
    zs = D + sd * yp[iw];
    if (zs < 0.0)
      fprintf(stderr,
	      "pom_layer: warning: fault extends above ground surface\n");
    /* Laplace-domain Bessel coefficients for this source depth */
    for (j = 0; j < NK && !err; j++) {
      kj = k[j];
      if (build_kterms(kj, H1, H2, zs, mu, lam, g, bulk, B1, B2, &kt)) {
	err = 1;
	break;
      }
      for (im = 0; im < 3; im++) {
	if (!mech_on[im])
	  continue;
	pom_source_vectors(Mt[im], kj, g, lam, mu, F, f2);
	for (io = 0; io < NORD; io++) {
	  spheroidal(&kt, F[io], t, mode,
		     &cUS[im][io][j], &cUR[im][io][j]);
	  if (io == 0)
	    cUT[im][io][j] = 0.0;
	  else
	    cUT[im][io][j] = toroidal(&kt, f2[io], t, mode);
	}
      }
    }
    if (err)
      break;
    /* inverse Hankel transform and assembly over along-strike sources */
    for (il = 0; il < NL; il++) {
      for (i = 0; i < n; i++) {
	double X, Y, r, th;
	double complex Uz_t, Ur_t, Uth_t;

	X = Xc[i] - xpos[iw * NL + il];
	Y = Yc[i] - ypos[iw * NL + il];
	r = sqrt(X * X + Y * Y);
	th = atan2(Y, X);
	wgt = wW[iw] * wL[il];
	/* Bessel functions depend on k and r only; evaluate once per
	   station and source, not per mechanism and order */
	for (j = 0; j < NK; j++) {
	  double kr;

	  kr = k[j] * r;
	  bb0[j] = j0(kr);
	  bb1[j] = j1(kr);
	  bb2[j] = besj2(kr, bb0[j], bb1[j]);
	  bb3[j] = besj3(kr, bb1[j], bb2[j]);
	}
	for (im = 0; im < 3; im++) {
	  if (!mech_on[im])
	    continue;
	  Uz_t = 0.0;
	  Ur_t = 0.0;
	  Uth_t = 0.0;
	  for (io = 0; io < NORD; io++) {
	    int M = order_m[io];
	    double complex sz, sr, sth, eim, iz, ir, ith;

	    sz = 0.0;
	    sr = 0.0;
	    sth = 0.0;
	    for (j = 0; j < NK; j++) {
	      double kr, b0, b1, b2, b3, DrJm, bm, bor;
	      double tw;

	      kj = k[j];
	      kr = kj * r;
	      b0 = bb0[j];
	      b1 = bb1[j];
	      b2 = bb2[j];
	      b3 = bb3[j];
	      switch (M) {
	      case 0:
		DrJm = -kj * b1;
		bm = b0;
		bor = 0.0;
		break;
	      case -1:
		DrJm = 0.5 * kj * (b2 - b0);
		bm = b1;
		/* b1/r, finite as r -> 0 */
		bor = (kr < 1e-6) ? 0.5 * kj : b1 / r;
		break;
	      case 1:
		DrJm = 0.5 * kj * (b0 - b2);
		bm = b1;
		bor = (kr < 1e-6) ? 0.5 * kj : b1 / r;
		break;
	      default:		/* |M| = 2 */
		DrJm = 0.5 * kj * (b1 - b3);
		bm = b2;
		bor = (kr < 1e-6) ? 0.0 : b2 / r;
		break;
	      }
	      if (M == 0)
		iz = kj * cUR[im][io][j] * bm;
	      else if (M == -1)
		iz = -kj * cUR[im][io][j] * bm;
	      else
		iz = kj * cUR[im][io][j] * bm;
	      /* the m-files write order -1 with a minus sign on the
	         i*M*UT term, all other orders with a plus sign */
	      if (M == -1)
		ir = cUS[im][io][j] * DrJm
		  - I * (double)M * cUT[im][io][j] * bor;
	      else
		ir = cUS[im][io][j] * DrJm
		  + I * (double)M * cUT[im][io][j] * bor;
	      if (M == -2 || M == 2)
		ith = I * (double)M * cUS[im][io][j] * bor
		  - cUT[im][io][j] * DrJm;
	      else
		ith = 0.0;
	      /* trapezoid over k; uniform spacing */
	      tw = (j == 0 || j == NK - 1) ? 0.5 : 1.0;
	      if (M != 0) {
		eim = cexp(I * (double)M * th);
		iz *= eim;
		ir *= eim;
		ith *= eim;
	      }
	      sz += tw * iz;
	      sr += tw * ir;
	      sth += tw * ith;
	    }
	    sz *= dk / (2.0 * M_PI);
	    sr *= dk / (2.0 * M_PI);
	    sth *= dk / (2.0 * M_PI);
	    Uz_t += sz;
	    Ur_t += creal(sr);
	    Uth_t += creal(sth);
	  }
	  Gx[im][i] += wgt * (creal(Ur_t) * cos(th)
			      + creal(Uth_t) * cos(th + M_PI / 2.0));
	  Gy[im][i] += wgt * (creal(Ur_t) * sin(th)
			      + creal(Uth_t) * sin(th + M_PI / 2.0));
	  Gz[im][i] += wgt * creal(Uz_t);
	}
      }
    }
  }

  if (!err) {
    for (i = 0; i < n; i++) {
      double e, no, up, sgn;

      e = 0.0;
      no = 0.0;
      up = 0.0;
      for (im = 0; im < 3; im++) {
	if (!mech_on[im])
	  continue;
	/* east = Gy, north = Gx, up = -Gz; dip-slip sign flip */
	sgn = (im == 1) ? -1.0 : 1.0;
	e += slip[im] * sgn * Gy[im][i];
	no += slip[im] * sgn * Gx[im][i];
	up += slip[im] * sgn * (-Gz[im][i]);
      }
      ue[i] = e;
      un[i] = no;
      uz[i] = up;
    }
  }

  free(Xc);
  free(Yc);
  free(xp);
  free(wL);
  free(yp);
  free(wW);
  free(k);
  free(xpos);
  free(ypos);
  free(cUR);
  free(cUS);
  free(cUT);
  free(Mt);
  for (im = 0; im < 3; im++) {
    free(Gx[im]);
    free(Gy[im]);
    free(Gz[im]);
  }
  return err;
}
