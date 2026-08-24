#define _GNU_SOURCE
/* pom_maxwell.c
 *
 * C port of Plate_over_Maxwell.m (Kaj M. Johnson, Indiana University):
 * postseismic displacements due to uniform slip on a small rectangular
 * dislocation, represented as a single point source, in an elastic
 * plate over a Maxwell viscoelastic halfspace. Semi-analytical
 * propagator matrix solution; the Laplace-domain denominator is a
 * cubic solved in closed form.
 */
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include "pom.h"

#define NK1 25			/* wavenumber samples (N in the m-file) */
#define KMIN1 0.000001
#define KMAX1 0.25		/* 1/km, as in the m-file */
#define NORD 5

static const int order_m[NORD] = { 0, -1, 1, -2, 2 };

/* cubic roots, replicating the MATLAB Cardano branch choices */
static void cubic_roots(double b2, double b1, double b0,
			double complex *z1, double complex *z2,
			double complex *z3)
{
  double q, r;
  double complex rt, s1, s2;

  q = b1 / 3.0 - b2 * b2 / 9.0;
  r = (b1 * b2 - 3.0 * b0) / 6.0 - b2 * b2 * b2 / 27.0;
  rt = q * q * q + r * r;
  rt = csqrt(rt);
  s1 = cpow(r + rt, 1.0 / 3.0);
  s2 = cpow(r - rt, 1.0 / 3.0);
  *z1 = s1 + s2 - b2 / 3.0;
  *z2 = -(s1 + s2) / 2.0 - b2 / 3.0 + I * sqrt(3.0) * (s1 - s2) / 2.0;
  *z3 = -(s1 + s2) / 2.0 - b2 / 3.0 - I * sqrt(3.0) * (s1 - s2) / 2.0;
}

/* denominator cubic coefficients (roots_right in the m-file) */
static void roots_right(double k, double u, double K, double L, double B,
			double Ph[4][4], double complex *s1,
			double complex *s2, double complex *s3, double *a)
{
  double P31, P32, P33, P34, P41, P42, P43, P44;
  double a3, a2, a1, a0, b2, b1, b0;

  P31 = Ph[2][0];
  P32 = Ph[2][1];
  P33 = Ph[2][2];
  P34 = Ph[2][3];
  P41 = Ph[3][0];
  P42 = Ph[3][1];
  P43 = Ph[3][2];
  P44 = Ph[3][3];
  a3 = (-L * P32 * P41 + L * P31 * P42 - 3 * P32 * P41 * u +
	2 * k * L * P34 * P41 * u + 3 * P31 * P42 * u -
	2 * k * L * P33 * P42 * u + 2 * k * L * P32 * P43 * u -
	2 * k * L * P31 * P44 * u + 2 * k * P33 * P41 * u * u +
	4 * k * P34 * P41 * u * u - 4 * k * P33 * P42 * u * u -
	2 * k * P34 * P42 * u * u - 2 * k * P31 * P43 * u * u +
	4 * k * P32 * P43 * u * u - 4 * k * k * L * P34 * P43 * u * u -
	4 * k * P31 * P44 * u * u + 2 * k * P32 * P44 * u * u +
	4 * k * k * L * P33 * P44 * u * u -
	4 * k * k * P34 * P43 * u * u * u +
	4 * k * k * P33 * P44 * u * u * u);
  a2 = (-2 * B * K * P32 * P41 - 4 * B * L * P32 * P41 +
	2 * B * K * P31 * P42 + 4 * B * L * P31 * P42 -
	12 * B * P32 * P41 * u + 4 * B * k * K * P34 * P41 * u +
	4 * B * k * L * P34 * P41 * u + 12 * B * P31 * P42 * u -
	4 * B * k * K * P33 * P42 * u - 4 * B * k * L * P33 * P42 * u +
	4 * B * k * K * P32 * P43 * u + 4 * B * k * L * P32 * P43 * u -
	4 * B * k * K * P31 * P44 * u - 4 * B * k * L * P31 * P44 * u +
	4 * B * k * P33 * P41 * u * u + 8 * B * k * P34 * P41 * u * u -
	8 * B * k * P33 * P42 * u * u - 4 * B * k * P34 * P42 * u * u -
	4 * B * k * P31 * P43 * u * u + 8 * B * k * P32 * P43 * u * u -
	8 * B * k * k * K * P34 * P43 * u * u -
	8 * B * k * P31 * P44 * u * u + 4 * B * k * P32 * P44 * u * u +
	8 * B * k * k * K * P33 * P44 * u * u);
  a1 = (-8 * B * B * K * P32 * P41 - 4 * B * B * L * P32 * P41 +
	8 * B * B * K * P31 * P42 + 4 * B * B * L * P31 * P42 -
	12 * B * B * P32 * P41 * u + 8 * B * B * k * K * P34 * P41 * u +
	12 * B * B * P31 * P42 * u - 8 * B * B * k * K * P33 * P42 * u +
	8 * B * B * k * K * P32 * P43 * u -
	8 * B * B * k * K * P31 * P44 * u);
  a0 = -8 * B * B * B * K * P32 * P41 + 8 * B * B * B * K * P31 * P42;
  b2 = a2 / a3;
  b1 = a1 / a3;
  b0 = a0 / a3;
  *a = a3;
  cubic_roots(b2, b1, b0, s1, s2, s3);
}

/* V matrix (makeV_right in the m-file); s complex */
static void makeV_right(double k, double K, double L, double complex s,
			double B, double H, double u,
			double complex V[4][2])
{
  V[0][0] = -(s + 2 * B);
  V[0][1] = -s * (2 * B + s) * u
    * (-1 + 4 * B * H * k * K + 2 * H * k * L * s + 2 * H * k * s * u);
  V[1][0] = -(s + 2 * B);
  V[1][1] = (-2 * B - s)
    * (2 * B * K + L * s + 2 * s * u + 4 * B * H * k * K * s * u
       + 2 * H * k * L * s * s * u + 2 * H * k * s * s * u * u);
  V[2][0] = 2 * u * s * k;
  V[2][1] = 4 * H * k * k * s * s * u * u * (2 * B * K + L * s + s * u);
  V[3][0] = 2 * u * s * k;
  V[3][1] = 2 * k * s * u * (2 * B * K + L * s + s * u)
    * (1 + 2 * H * k * s * u);
}

/* one root's contribution to the 2-vector Fhat */
static void root_term(double kj, double K, double L, double B, double H,
		      double u, double Ph4[4][4], double Pzs4[4][4],
		      const double complex F[4], double complex s,
		      double complex sa, double complex sb, double a,
		      double t, double complex out[2])
{
  int i, j, l;
  double complex V[4][2], A[2][2], Adj[2][2], PzsF[2], AdP[2], VA[4][2];
  double complex num[2], den, tf;

  makeV_right(kj, K, L, s, B, H, u, V);
  /* A = Ph4(3:4,:)*V */
  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++) {
      A[i][j] = 0.0;
      for (l = 0; l < 4; l++)
	A[i][j] += Ph4[2 + i][l] * V[l][j];
    }
  Adj[0][0] = A[1][1];
  Adj[1][1] = A[0][0];
  Adj[0][1] = -A[0][1];
  Adj[1][0] = -A[1][0];
  /* PzsF = Pzs4(3:4,:)*F */
  for (i = 0; i < 2; i++) {
    PzsF[i] = 0.0;
    for (l = 0; l < 4; l++)
      PzsF[i] += Pzs4[2 + i][l] * F[l];
  }
  /* AdP = Adj*PzsF */
  for (i = 0; i < 2; i++)
    AdP[i] = Adj[i][0] * PzsF[0] + Adj[i][1] * PzsF[1];
  /* VA = V*AdP as 4-vector, then Num = -Ph4(1:2,:)*V*Adj*PzsF */
  for (i = 0; i < 4; i++)
    VA[i][0] = V[i][0] * AdP[0] + V[i][1] * AdP[1];
  for (i = 0; i < 2; i++) {
    num[i] = 0.0;
    for (l = 0; l < 4; l++)
      num[i] += Ph4[i][l] * VA[l][0];
    num[i] = -num[i];
  }
  den = a * s * (s - sa) * (s - sb);
  tf = cexp(s * t) - 1.0;
  out[0] = num[0] * tf / den;
  out[1] = num[1] * tf / den;
}
/*
 * interior variant: the transient field is continuous across the
 * source (the slip jump lives in the elastic part, which the
 * exp(s t) - 1 residue excludes), and its surface tractions vanish,
 * so the transient state at plate depth zr is the plate propagator
 * applied to the surface transient state, per root:
 *   N(zr) = P4(-zr) (num0, num1, 0, 0)
 * out4 returns all four state components (two displacement-type,
 * two traction-type rows) times tf/den
 */
static void root_term_z(double kj, double K, double L, double B, double H,
			double u, double Ph4[4][4], double Pzs4[4][4],
			double P4z[4][4], int mflip,
			const double complex F[4], double complex s,
			double complex sa, double complex sb, double a,
			double t, double complex out4[4])
{
  int i, j, l;
  double complex V[4][2], A[2][2], Adj[2][2], PzsF[2], AdP[2], VA[4][2];
  double complex num[2], den, tf, N0[4], Nz[4];

  makeV_right(kj, K, L, s, B, H, u, V);
  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++) {
      A[i][j] = 0.0;
      for (l = 0; l < 4; l++)
	A[i][j] += Ph4[2 + i][l] * V[l][j];
    }
  Adj[0][0] = A[1][1];
  Adj[1][1] = A[0][0];
  Adj[0][1] = -A[0][1];
  Adj[1][0] = -A[1][0];
  for (i = 0; i < 2; i++) {
    PzsF[i] = 0.0;
    for (l = 0; l < 4; l++)
      PzsF[i] += Pzs4[2 + i][l] * F[l];
  }
  for (i = 0; i < 2; i++)
    AdP[i] = Adj[i][0] * PzsF[0] + Adj[i][1] * PzsF[1];
  for (i = 0; i < 4; i++)
    VA[i][0] = V[i][0] * AdP[0] + V[i][1] * AdP[1];
  for (i = 0; i < 2; i++) {
    num[i] = 0.0;
    for (l = 0; l < 4; l++)
      num[i] += Ph4[i][l] * VA[l][0];
    num[i] = -num[i];
  }
  (void)mflip;
  N0[0] = num[0];N0[1] = num[1];N0[2] = 0.0;N0[3] = 0.0;
  for (i = 0; i < 4; i++) {
    Nz[i] = 0.0;
    for (l = 0; l < 4; l++)
      Nz[i] += P4z[i][l] * N0[l];
  }
  den = a * s * (s - sa) * (s - sb);
  tf = cexp(s * t) - 1.0;
  for (i = 0; i < 4; i++)
    out4[i] = Nz[i] * tf / den;
}

int pom_maxwell(const double *m, const double *xs, const double *ys, int n,
		double H, double mu_lam, const double *t, int nt, double tR,
		double *ue, double *un, double *uz)
{
  return pom_maxwell_z(m, xs, ys, n, H, mu_lam, t, nt, tR, 0.0,
		       ue, un, uz, NULL, NULL, NULL);
}
/*
 * interior evaluation at receiver depth zr (0 <= zr < H, plate
 * only): transient displacements and, when the pointers are given,
 * the transient traction components on horizontal planes
 * (szz, sxz, syz), all per unit slip as for the surface routine
 */
int pom_maxwell_z(const double *m, const double *xs, const double *ys, int n,
		  double H, double mu_lam, const double *t, int nt, double tR,
		  double zr,
		  double *ue, double *un, double *uz,
		  double *szz, double *sxz, double *syz)
{
  int i, j, io, im, it, mech_on[3];
  double mu, lam, g, bulk, beta;
  double L, W, dip, strike, D, dipdir, offset, eoff, noff;
  double sd, cd, kj, dk, zs, xr, yr, x0, y0;
  double slip[3];
  double M1[3][3], M2[3][3], M3[3][3];
  double (*Mt)[3][3];
  double Ph4[4][4], Ph2[2][2], Pzs4[4][4], Pzs2[2][2];
  double P4z[4][4], P2z[2][2];
  double complex term4[4], ut0z[2];
  double complex *cTS = NULL, *cTR = NULL, *cTT = NULL;
  double *Xc, *Yc, k[NK1];
  double complex *cUR, *cUS, *cUT;	/* [mech][order][j][it] */
  double complex F[NORD][4], f2[NORD][2];
  double *Gx[3], *Gy[3], *Gz[3];
  double *GTx[3], *GTy[3], *GTz[3];
  double *GEh[3], *GGx[3], *GGy[3];

#define CIDX(im,io,j,it) ((((im)*NORD+(io))*NK1+(j))*nt+(it))
  if (n < 1 || nt < 1 || tR <= 0.0)
    return 1;
  mu = 1.0;
  lam = 1.0 / mu_lam;
  bulk = lam + 2.0 * mu / 3.0;
  g = lam + 2.0 * mu;
  beta = 1.0 / tR;
  L = m[0];
  W = m[1];
  dip = m[3];
  strike = m[4];
  slip[0] = -m[7];		/* sign convention of the m-file */
  slip[1] = m[8];
  slip[2] = m[9];
  for (im = 0; im < 3; im++)
    mech_on[im] = (slip[im] != 0.0);

  if (dip <= 90.0 && dip >= -90.0)
    dipdir = (strike + 90.0) * M_PI / 180.0;
  else
    dipdir = (strike - 90.0) * M_PI / 180.0;
  offset = fabs(W * cos(dip * M_PI / 180.0));
  eoff = offset * sin(dipdir);
  noff = offset * cos(dipdir);
  D = m[2] - W * sin(dip * M_PI / 180.0);
  sd = sin(dip * M_PI / 180.0);
  cd = cos(dip * M_PI / 180.0);
  /* single point source at the patch center */
  x0 = 0.0;
  y0 = W / 2.0;
  xr = cos(strike * M_PI / 180.0) * x0 - sin(strike * M_PI / 180.0) * cd * y0;
  yr = sin(strike * M_PI / 180.0) * x0 + cos(strike * M_PI / 180.0) * cd * y0;
  zs = D + sd * y0;
  if (zs < 0.0)
    fprintf(stderr, "pom_maxwell: warning: fault extends above surface\n");

  Xc = malloc(sizeof(double) * n);
  Yc = malloc(sizeof(double) * n);
  cUR = malloc(sizeof(double complex) * 3 * NORD * NK1 * nt);
  cUS = malloc(sizeof(double complex) * 3 * NORD * NK1 * nt);
  cUT = malloc(sizeof(double complex) * 3 * NORD * NK1 * nt);
  cTS = malloc(sizeof(double complex) * 3 * NORD * NK1 * nt);
  cTR = malloc(sizeof(double complex) * 3 * NORD * NK1 * nt);
  cTT = malloc(sizeof(double complex) * 3 * NORD * NK1 * nt);
  Mt = malloc(sizeof(double) * 27);
  for (im = 0; im < 3; im++) {
    Gx[im] = calloc(n * nt, sizeof(double));
    Gy[im] = calloc(n * nt, sizeof(double));
    Gz[im] = calloc(n * nt, sizeof(double));
    GTx[im] = calloc(n * nt, sizeof(double));
    GTy[im] = calloc(n * nt, sizeof(double));
    GTz[im] = calloc(n * nt, sizeof(double));
    GEh[im] = calloc(n * nt, sizeof(double));
    GGx[im] = calloc(n * nt, sizeof(double));
    GGy[im] = calloc(n * nt, sizeof(double));
  }
  if (!Xc || !Yc || !cUR || !cUS || !cUT || !Mt)
    return 1;
  for (i = 0; i < n; i++) {
    Xc[i] = ys[i] - m[6] + noff;
    Yc[i] = xs[i] - m[5] + eoff;
  }
  dk = (KMAX1 - KMIN1) / (NK1 - 1);
  for (j = 0; j < NK1; j++)
    k[j] = KMIN1 + dk * j;
  pom_momtensor(strike, dip, lam, mu, M1, M2, M3);
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      Mt[0][i][j] = M1[i][j];
      Mt[1][i][j] = M2[i][j];
      Mt[2][i][j] = M3[i][j];
    }

  for (j = 0; j < NK1; j++) {
    double complex s1, s2, s3, term[2], acc[2];
    double a, aa, tor_pref;

    kj = k[j];
    pom_prop4(kj, -H, mu, lam, g, Ph4);
    pom_prop2(kj, -H, mu, Ph2);
    pom_prop4(kj, -zs, mu, lam, g, Pzs4);
    pom_prop2(kj, -zs, mu, Pzs2);
    pom_prop4(kj, zr, mu, lam, g, P4z);	/* surface DOWN to receiver
					   depth: the generator acts in
					   depth, and Ph4/Pzs4 use
					   negative arguments because
					   they map upward to the
					   surface */
    pom_prop2(kj, zr, mu, P2z);
    roots_right(kj, mu, bulk, lam, beta, Ph4, &s1, &s2, &s3, &a);
    aa = Ph2[1][0] - mu * fabs(kj) * Ph2[1][1];
    tor_pref = Ph2[0][0] / (Ph2[1][0] * aa) * (aa - Ph2[1][0])
      + mu * fabs(kj) / aa * Ph2[0][1];
    for (im = 0; im < 3; im++) {
      if (!mech_on[im])
	continue;
      pom_source_vectors(Mt[im], kj, g, lam, mu, F, f2);
      for (io = 0; io < NORD; io++) {
	for (it = 0; it < nt; it++) {
	  double complex a4[4];
	  int q;
	  for (q = 0; q < 4; q++)
	    a4[q] = 0.0;
	  root_term_z(kj, bulk, lam, beta, H, mu, Ph4, Pzs4, P4z,
		      (order_m[io] == -1), F[io],
		      s1, s2, s3, a, t[it], term4);
	  for (q = 0; q < 4; q++)
	    a4[q] += term4[q];
	  root_term_z(kj, bulk, lam, beta, H, mu, Ph4, Pzs4, P4z,
		      (order_m[io] == -1), F[io],
		      s2, s1, s3, a, t[it], term4);
	  for (q = 0; q < 4; q++)
	    a4[q] += term4[q];
	  root_term_z(kj, bulk, lam, beta, H, mu, Ph4, Pzs4, P4z,
		      (order_m[io] == -1), F[io],
		      s3, s1, s2, a, t[it], term4);
	  for (q = 0; q < 4; q++)
	    a4[q] += term4[q];
	  cUS[CIDX(im, io, j, it)] = a4[0];
	  cUR[CIDX(im, io, j, it)] = -a4[1];
	  cTS[CIDX(im, io, j, it)] = a4[2];
	  cTR[CIDX(im, io, j, it)] = -a4[3];
	  if (io == 0) {
	    cUT[CIDX(im, io, j, it)] = 0.0;
	    cTT[CIDX(im, io, j, it)] = 0.0;
	  } else {
	    /* toroidal surface transient, then propagate the
	       (u_t, t_t) pair to depth zr */
	    ut0z[0] = tor_pref
	      * (Pzs2[1][0] * f2[io][0] + Pzs2[1][1] * f2[io][1])
	      * (exp(-2.0 * beta * Ph2[1][0] * t[it] / aa) - 1.0);
	    ut0z[1] = 0.0;
	    cUT[CIDX(im, io, j, it)] = P2z[0][0]*ut0z[0] + P2z[0][1]*ut0z[1];
	    cTT[CIDX(im, io, j, it)] = P2z[1][0]*ut0z[0] + P2z[1][1]*ut0z[1];
	  }
	}
      }
    }
  }

  /* inverse Hankel transform and assembly (single point source) */
  for (i = 0; i < n; i++) {
    double X, Y, r, th;

    X = Xc[i] - xr;
    Y = Yc[i] - yr;
    r = sqrt(X * X + Y * Y);
    th = atan2(Y, X);
    for (im = 0; im < 3; im++) {
      if (!mech_on[im])
	continue;
      for (it = 0; it < nt; it++) {
	double complex Uz_t, Ur_t, Uth_t;
	double complex Tz_t, Tr_t, Tth_t;
	double complex Eh_t, Gr_t, Gth_t;

	Uz_t = 0.0;
	Ur_t = 0.0;
	Uth_t = 0.0;
	Tz_t = 0.0;
	Tr_t = 0.0;
	Tth_t = 0.0;
	Eh_t = 0.0;
	Gr_t = 0.0;
	Gth_t = 0.0;
	for (io = 0; io < NORD; io++) {
	  int M = order_m[io];
	  double complex sz, sr, sth, eim, iz, ir, ith;
	  double complex sz2, sr2, sth2, tz, tr, tth;
	  double complex se, sgr, sgth, ee, gr, gth, uzc;

	  sz = 0.0;
	  sr = 0.0;
	  sth = 0.0;
	  sz2 = 0.0;
	  sr2 = 0.0;
	  sth2 = 0.0;
	  se = 0.0;
	  sgr = 0.0;
	  sgth = 0.0;
	  for (j = 0; j < NK1; j++) {
	    double kr, b0, b1, b2, b3, DrJm, bm, bor, tw;

	    kj = k[j];
	    kr = kj * r;
	    b0 = j0(kr);
	    b1 = j1(kr);
	    b2 = (fabs(kr) < 1e-6) ? kr * kr / 8.0 : -b0 + 2.0 / kr * b1;
	    b3 = (fabs(kr) < 1e-6) ? kr * kr * kr / 48.0
	      : -b1 + 4.0 / kr * b2;
	    switch (M) {
	    case 0:
	      DrJm = -kj * b1;
	      bm = b0;
	      bor = 0.0;
	      break;
	    case -1:
	    case 1:
	      DrJm = (M == -1) ? 0.5 * kj * (b2 - b0) : 0.5 * kj * (b0 - b2);
	      bm = b1;
	      bor = (kr < 1e-6) ? 0.5 * kj : b1 / r;
	      break;
	    default:
	      DrJm = 0.5 * kj * (b1 - b3);
	      bm = b2;
	      bor = (kr < 1e-6) ? 0.0 : b2 / r;
	      break;
	    }
	    iz = ((M == -1) ? -kj : kj) * cUR[CIDX(im, io, j, it)] * bm;
	    if (M == -1)
	      ir = cUS[CIDX(im, io, j, it)] * DrJm
		- I * (double)M * cUT[CIDX(im, io, j, it)] * bor;
	    else
	      ir = cUS[CIDX(im, io, j, it)] * DrJm
		+ I * (double)M * cUT[CIDX(im, io, j, it)] * bor;
	    if (M == -2 || M == 2)
	      ith = I * (double)M * cUS[CIDX(im, io, j, it)] * bor
		- cUT[CIDX(im, io, j, it)] * DrJm;
	    else
	      ith = 0.0;
	    /* traction rows: same Hankel-order structure as the
	       paired displacement rows (calibration of the vertical
	       k factor by the finite-difference check) */
	    tz = ((M == -1) ? -kj : kj) * cTR[CIDX(im, io, j, it)] * bm;
	    if (M == -1)
	      tr = cTS[CIDX(im, io, j, it)] * DrJm
		- I * (double)M * cTT[CIDX(im, io, j, it)] * bor;
	    else
	      tr = cTS[CIDX(im, io, j, it)] * DrJm
		+ I * (double)M * cTT[CIDX(im, io, j, it)] * bor;
	    if (M == -2 || M == 2)
	      tth = I * (double)M * cTS[CIDX(im, io, j, it)] * bor
		- cTT[CIDX(im, io, j, it)] * DrJm;
	    else
	      tth = 0.0;
	    /* frame corrections: horizontal divergence of the
	       displacement (spheroidal only) and horizontal gradient
	       of the vertical displacement field */
	    ee = -kj * kj * cUS[CIDX(im, io, j, it)] * bm;
	    uzc = ((M == -1) ? -kj : kj) * cUR[CIDX(im, io, j, it)];
	    gr = uzc * DrJm;
	    gth = I * (double)M * uzc * bor;
	    if (M != 0) {
	      eim = cexp(I * (double)M * th);
	      iz *= eim;
	      ir *= eim;
	      ith *= eim;
	      tz *= eim;
	      tr *= eim;
	      tth *= eim;
	      ee *= eim;
	      gr *= eim;
	      gth *= eim;
	    }
	    tw = (j == 0 || j == NK1 - 1) ? 0.5 : 1.0;
	    sz += tw * iz;
	    sr += tw * ir;
	    sth += tw * ith;
	    sz2 += tw * tz;
	    sr2 += tw * tr;
	    sth2 += tw * tth;
	    se += tw * ee;
	    sgr += tw * gr;
	    sgth += tw * gth;
	  }
	  Uz_t += sz * dk / (2.0 * M_PI);
	  Ur_t += creal(sr * dk / (2.0 * M_PI));
	  Uth_t += creal(sth * dk / (2.0 * M_PI));
	  Tz_t += sz2 * dk / (2.0 * M_PI);
	  Tr_t += creal(sr2 * dk / (2.0 * M_PI));
	  Tth_t += creal(sth2 * dk / (2.0 * M_PI));
	  Eh_t += se * dk / (2.0 * M_PI);
	  Gr_t += creal(sgr * dk / (2.0 * M_PI));
	  Gth_t += creal(sgth * dk / (2.0 * M_PI));
	}
	Gx[im][i * nt + it] = creal(Ur_t) * cos(th)
	  + creal(Uth_t) * cos(th + M_PI / 2.0);
	Gy[im][i * nt + it] = creal(Ur_t) * sin(th)
	  + creal(Uth_t) * sin(th + M_PI / 2.0);
	Gz[im][i * nt + it] = creal(Uz_t);
	GTx[im][i * nt + it] = creal(Tr_t) * cos(th)
	  + creal(Tth_t) * cos(th + M_PI / 2.0);
	GTy[im][i * nt + it] = creal(Tr_t) * sin(th)
	  + creal(Tth_t) * sin(th + M_PI / 2.0);
	GTz[im][i * nt + it] = creal(Tz_t);
	GEh[im][i * nt + it] = creal(Eh_t);
	GGx[im][i * nt + it] = creal(Gr_t) * cos(th)
	  + creal(Gth_t) * cos(th + M_PI / 2.0);
	GGy[im][i * nt + it] = creal(Gr_t) * sin(th)
	  + creal(Gth_t) * sin(th + M_PI / 2.0);
      }
    }
  }

  for (i = 0; i < n; i++)
    for (it = 0; it < nt; it++) {
      double e, no, up, sgn, w;

      e = 0.0;
      no = 0.0;
      up = 0.0;
      w = W * L;		/* patch area */
      for (im = 0; im < 3; im++) {
	if (!mech_on[im])
	  continue;
	sgn = (im == 1) ? -1.0 : 1.0;	/* dip-slip Okada sign */
	no += slip[im] * sgn * w * Gx[im][i * nt + it];
	e += slip[im] * sgn * w * Gy[im][i * nt + it];
	up += slip[im] * sgn * w * (-Gz[im][i * nt + it]);
      }
      ue[i * nt + it] = e;
      un[i * nt + it] = no;
      uz[i * nt + it] = up;
      if (szz && sxz && syz) {
	double te2, tn2, tzz2, eh2, gze, gzn;
	double lam2, mu2;

	/* internal normalized moduli (g = 1 basis, mu_lam ratio) */
	mu2 = 1.0;
	lam2 = 1.0;
	te2 = 0.0;
	tn2 = 0.0;
	tzz2 = 0.0;
	eh2 = 0.0;
	gze = 0.0;
	gzn = 0.0;
	for (im = 0; im < 3; im++) {
	  if (!mech_on[im])
	    continue;
	  sgn = (im == 1) ? -1.0 : 1.0;
	  tn2 += slip[im] * sgn * w * GTx[im][i * nt + it];
	  te2 += slip[im] * sgn * w * GTy[im][i * nt + it];
	  tzz2 += slip[im] * sgn * w * (-GTz[im][i * nt + it]);
	  eh2 += slip[im] * sgn * w * GEh[im][i * nt + it];
	  /* gradient of the OUTPUT vertical displacement (up
	     positive), hence the sign flip mirroring up = -Gz */
	  gzn += slip[im] * sgn * w * (-GGx[im][i * nt + it]);
	  gze += slip[im] * sgn * w * (-GGy[im][i * nt + it]);
	}
	/* the state rows are the physical traction coefficients in
	   the same basis as the displacements (see the generator in
	   pom_prop4: row 3 is sigma_xz-type, row 4 sigma_zz-type),
	   so no constitutive corrections apply; the eh2/gze/gzn
	   fields above are kept because they are the ingredients of
	   the in-plane components sigma_xx, sigma_yy in a later
	   extension. lam2/mu2 reserved for that use */
	(void)lam2;(void)mu2;(void)eh2;(void)gze;(void)gzn;
	/* one global sign maps the storage convention to the
	   physical z-up tensor (calibrated by finite differences of
	   the interior displacements across mechanisms, depths, and
	   times) */
	szz[i * nt + it] = -tzz2;
	sxz[i * nt + it] = -te2;
	syz[i * nt + it] = -tn2;
      }
    }

  free(Xc);
  free(Yc);
  free(cUR);
  free(cUS);
  free(cUT);
  free(cTS);
  free(cTR);
  free(cTT);
  for (im = 0; im < 3; im++) {
    free(GTx[im]);
    free(GTy[im]);
    free(GTz[im]);
    free(GEh[im]);
    free(GGx[im]);
    free(GGy[im]);
  }
  free(Mt);
  for (im = 0; im < 3; im++) {
    free(Gx[im]);
    free(Gy[im]);
    free(Gz[im]);
  }
  return 0;
}
