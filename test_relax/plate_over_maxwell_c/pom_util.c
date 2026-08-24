#define _GNU_SOURCE
/* pom_util.c
 *
 * numerical utilities for the Plate_over_Maxwell C port
 */
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include "pom.h"

/*
 * Gauss-Legendre nodes and weights on [a,b], n >= 1
 * (Newton iteration on Legendre polynomials; replaces the
 * Golub-Welsch eigenvalue construction of the MATLAB codes)
 */
void pom_gauleg(double a, double b, int n, double *x, double *w)
{
  int i, j, m;
  double z, z1, p1, p2, p3, pp, xm, xl;

  m = (n + 1) / 2;
  xm = 0.5 * (b + a);
  xl = 0.5 * (b - a);
  for (i = 0; i < m; i++) {
    z = cos(M_PI * (i + 0.75) / (n + 0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (j = 0; j < n; j++) {
	p3 = p2;
	p2 = p1;
	p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    } while (fabs(z - z1) > 1e-15);
    x[i] = xm - xl * z;
    x[n - 1 - i] = xm + xl * z;
    w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n - 1 - i] = w[i];
  }
}

/*
 * roots of a real polynomial c[0] + c[1] s + ... + c[deg] s^deg
 * by the Aberth-Ehrlich method with complex arithmetic;
 * returns 0 on success. Roots are unordered.
 */
int pom_poly_roots(int deg, const double *c, double complex *r)
{
  int i, j, it, done;
  double rad, cm;
  double complex p, dp, s, corr, sum;
  const int maxit = 200;
  const double tol = 1e-14;

  if (deg < 1 || c[deg] == 0.0)
    return 1;
  /* Cauchy-style bound on root magnitudes */
  cm = 0.0;
  for (i = 0; i < deg; i++) {
    double q = fabs(c[i] / c[deg]);
    if (q > cm)
      cm = q;
  }
  rad = 1.0 + cm;
  /* initial guesses on a circle, non-symmetric angle offset */
  for (i = 0; i < deg; i++)
    r[i] = 0.5 * rad * cexp(I * (2.0 * M_PI * i / deg + 0.4));
  for (it = 0; it < maxit; it++) {
    done = 1;
    for (i = 0; i < deg; i++) {
      s = r[i];
      /* Horner for p and p' */
      p = c[deg];
      dp = 0.0;
      for (j = deg - 1; j >= 0; j--) {
	dp = dp * s + p;
	p = p * s + c[j];
      }
      if (cabs(p) == 0.0)
	continue;
      sum = 0.0;
      for (j = 0; j < deg; j++)
	if (j != i)
	  sum += 1.0 / (s - r[j]);
      corr = (p / dp) / (1.0 - (p / dp) * sum);
      r[i] = s - corr;
      if (cabs(corr) > tol * (1.0 + cabs(s)))
	done = 0;
    }
    if (done)
      break;
  }
  return 0;
}

/*
 * moment tensor basis for unit strike-slip, dip-slip, tensile
 * dislocations (momtensor_inverse in the MATLAB codes)
 */
static void mt_fill(const double *va, const double *vn,
		    double lam, double mu, double M[3][3])
{
  double d;

  d = va[0] * vn[0] + va[1] * vn[1] + va[2] * vn[2];
  M[0][0] = va[0] * vn[0] * 2.0 * mu + lam * d;
  M[1][1] = va[1] * vn[1] * 2.0 * mu + lam * d;
  M[2][2] = va[2] * vn[2] * 2.0 * mu + lam * d;
  M[0][1] = (va[0] * vn[1] + va[1] * vn[0]) * mu;
  M[0][2] = (va[0] * vn[2] + va[2] * vn[0]) * mu;
  M[1][2] = (va[1] * vn[2] + va[2] * vn[1]) * mu;
  M[1][0] = M[0][1];
  M[2][0] = M[0][2];
  M[2][1] = M[1][2];
}

static void mt_clip(double M[3][3])
{
  int i, j;
  double mx;

  mx = 0.0;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (fabs(M[i][j]) > mx)
	mx = fabs(M[i][j]);
  if (mx == 0.0)
    return;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (fabs(M[i][j] / mx) < 1e-13)
	M[i][j] = 0.0;
}

void pom_momtensor(double strike, double dip, double lam, double mu,
		   double M1[3][3], double M2[3][3], double M3[3][3])
{
  double st, di, vs[3], vd[3], vn[3];

  st = strike * M_PI / 180.0;
  di = dip * M_PI / 180.0;
  vs[0] = cos(st);
  vs[1] = sin(st);
  vs[2] = 0.0;
  vd[0] = -cos(di) * sin(st);
  vd[1] = cos(di) * cos(st);
  vd[2] = sin(di);
  /* vn = vd x vs */
  vn[0] = vd[1] * vs[2] - vd[2] * vs[1];
  vn[1] = vd[2] * vs[0] - vd[0] * vs[2];
  vn[2] = vd[0] * vs[1] - vd[1] * vs[0];
  mt_fill(vs, vn, lam, mu, M1);
  mt_fill(vd, vn, lam, mu, M2);
  mt_fill(vn, vn, lam, mu, M3);
  mt_clip(M1);
  mt_clip(M2);
  mt_clip(M3);
}

/* 4x4 matrix product c = a*b */
static void mat4mul(double a[4][4], double b[4][4], double c[4][4])
{
  int i, j, l;
  double s;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) {
      s = 0.0;
      for (l = 0; l < 4; l++)
	s += a[i][l] * b[l][j];
      c[i][j] = s;
    }
}

/*
 * plane-strain 4x4 propagator P(dz) = C3 A^3 + C2 A^2 + C1 A + C0 I
 * with dz = z - z0 (getprop in the MATLAB codes)
 */
void pom_prop4(double kj, double dz, double mu, double lam, double g,
	       double P[4][4])
{
  int i, j;
  double A[4][4], A2[4][4], A3[4][4];
  double C0, C1, C2, C3, kz, sh, ch;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      A[i][j] = 0.0;
  A[0][1] = kj;
  A[0][2] = 1.0 / mu;
  A[1][0] = -kj * lam / g;
  A[1][3] = 1.0 / g;
  A[2][0] = 4.0 * kj * kj * mu * (lam + mu) / g;
  A[2][3] = kj * lam / g;
  A[3][2] = -kj;
  mat4mul(A, A, A2);
  mat4mul(A2, A, A3);
  kz = kj * dz;
  sh = sinh(kz);
  ch = cosh(kz);
  C3 = -(sh - kz * ch) / (2.0 * kj * kj * kj);
  C2 = kz * sh / (2.0 * kj * kj);
  C1 = (3.0 * sh - kz * ch) / (2.0 * kj);
  C0 = (2.0 * ch - kz * sh) / 2.0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      P[i][j] = C3 * A3[i][j] + C2 * A2[i][j] + C1 * A[i][j]
	+ ((i == j) ? C0 : 0.0);
}

/* antiplane 2x2 propagator */
void pom_prop2(double kj, double dz, double mu, double P[2][2])
{
  double ak, sh, ch;

  ak = fabs(kj);
  sh = sinh(dz * ak);
  ch = cosh(dz * ak);
  P[0][0] = ch;
  P[0][1] = sh / (mu * ak);
  P[1][0] = mu * ak * sh;
  P[1][1] = ch;
}

/* source vectors per Fourier order for moment tensor M (getbesselco) */
void pom_source_vectors(double M[3][3], double kj, double g, double lam,
			   double mu, double complex F[5][4],
			   double complex f2[5][2])
{
  double Mxx, Myy, Mzz, Mxy, Mxz, Myz;
  int i, j;

  Mxx = M[0][0];
  Myy = M[1][1];
  Mzz = M[2][2];
  Mxy = M[0][1];
  Mxz = M[0][2];
  Myz = M[1][2];
  for (i = 0; i < 5; i++) {
    for (j = 0; j < 4; j++)
      F[i][j] = 0.0;
    f2[i][0] = 0.0;
    f2[i][1] = 0.0;
  }
  /* order 0 */
  F[0][1] = Mzz / g;
  F[0][2] = kj * ((-Mxx - Myy) / 2.0 + lam * Mzz / g);
  /* order -1 */
  F[1][0] = -0.5 * (Mxz + I * Myz) / mu;
  f2[1][0] = 0.5 * (Myz - I * Mxz) / mu;
  /* order +1 */
  F[2][0] = -0.5 * (-Mxz + I * Myz) / mu;
  f2[2][0] = -0.5 * (Myz + I * Mxz) / mu;
  /* order -2 */
  F[3][2] = kj * (-(Myy - Mxx) / 4.0 + I * Mxy / 2.0);
  f2[3][1] = kj * (I * (Mxx - Myy) / 4.0 - 0.5 * Mxy);
  /* order +2 */
  F[4][2] = kj * (-(Myy - Mxx) / 4.0 - I * Mxy / 2.0);
  f2[4][1] = kj * (-I * (Mxx - Myy) / 4.0 - 0.5 * Mxy);
}
