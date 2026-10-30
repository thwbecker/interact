/*
  2-D in-plane kernels for an OPENING (tensile) segment with constant opening rate (fully
  dynamic extension, fd_plan.md M7 step 1): physical stress (sxx, sxy, syy) in the source frame
  at (xr, yr), time t, from unit opening rate on |x'| <= a, y' = 0 since t = 0 (opening =
  u_y(+) - u_y(-)). Port of fully_dynamic/inplane_open.py: quadrature of the generated point
  field (fd_kernels_open_gen.c, from gen_inplane_open_c.py) between the wave-front crossings,
  with the complement identity for receivers on the segment (the whole-line integral is the
  plane P wave (sxx, sxy, syy) = (-lambda/(2 alpha), 0, -rho alpha/2)).
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd_kernels.h"

void fd_sig_open_ps(double x, double y, double t, double alpha, double beta, double rho, double *s);
void fd_sig_open_p(double x, double y, double t, double alpha, double beta, double rho, double *s);

#define FD_OPEN_NG 12
#define FD_OPEN_MAXBP 64

static int g_gl_init = 0;
static double g_gl_v[FD_OPEN_NG], g_gl_w[FD_OPEN_NG];

/* Gauss-Legendre nodes and weights on [0, 1] */
static void gl_init(void)
{
  int i, j, it;
  double z, z1, p1, p2, p3, pp;
  if(g_gl_init)
    return;
  for(i = 0; i < FD_OPEN_NG; i++){
    z = cos(M_PI*(i + 0.75)/(FD_OPEN_NG + 0.5));
    for(it = 0; it < 100; it++){
      p1 = 1.0;
      p2 = 0.0;
      for(j = 0; j < FD_OPEN_NG; j++){
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*j + 1.0)*z*p2 - j*p3)/(j + 1.0);
      }
      pp = FD_OPEN_NG*(z*p1 - p2)/(z*z - 1.0);
      z1 = z;
      z = z1 - p1/pp;
      if(fabs(z - z1) < 1e-15)
        break;
    }
    g_gl_v[i] = 0.5*(1.0 - z);
    g_gl_w[i] = 1.0/((1.0 - z*z)*pp*pp);      /* 2/((1-z^2) pp^2) on [-1,1], halved on [0,1] */
  }
  g_gl_init = 1;
}

/* the generated point field (the prototype's sign convention) */
static void point_open(double x, double y, double t, double alpha, double beta, double rho, double *s)
{
  double R;
  int k;
  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 0.0;
  R = hypot(x, y);
  if(beta*t > R)
    fd_sig_open_ps(x, y, t, alpha, beta, rho, s);
  else if(alpha*t > R)
    fd_sig_open_p(x, y, t, alpha, beta, rho, s);
  for(k = 0; k < 3; k++)
    if(!isfinite(s[k]))
      s[k] = 0.0;
}

static int cmp_double(const void *a, const void *b)
{
  double x, y;
  x = *(const double *)a;
  y = *(const double *)b;
  return (x > y) - (x < y);
}

/* integral over x' in [x_lo, x_hi] of the point field at (xr - x', yr), breakpoints at the P
   and S front crossings and at bp[0..nbp-1]; cubic map against the endpoint singularities */
static void integrate_S(double x_lo, double x_hi, double xr, double yr, double t, double alpha, double beta, double rho,
    const double *bp, int nbp, double *out)
{
  double B[FD_OPEN_MAXBP], c, h, sgn, lo, hi, L, v, ph, dph, xp, wt, s[3];
  int nb, i, k, q, m;
  gl_init();
  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 0.0;
  if(x_hi <= x_lo)
    return;
  nb = 0;
  B[nb++] = x_lo;
  B[nb++] = x_hi;
  for(m = 0; m < 2; m++){
    c = m ? beta : alpha;
    if(c*t > fabs(yr)){
      h = sqrt((c*t)*(c*t) - yr*yr);
      for(q = 0; q < 2; q++){
        sgn = q ? 1.0 : -1.0;
        B[nb++] = fmin(fmax(xr + sgn*h, x_lo), x_hi);
      }
    }
  }
  for(i = 0; i < nbp && nb < FD_OPEN_MAXBP; i++)
    B[nb++] = fmin(fmax(bp[i], x_lo), x_hi);
  qsort(B, nb, sizeof(double), cmp_double);
  for(i = 0; i + 1 < nb; i++){
    lo = B[i];
    hi = B[i + 1];
    L = hi - lo;
    if(L <= 0.0)
      continue;
    for(q = 0; q < FD_OPEN_NG; q++){
      v = g_gl_v[q];
      ph = 3.0*v*v - 2.0*v*v*v;
      dph = 6.0*v*(1.0 - v);
      xp = lo + L*ph;
      wt = L*g_gl_w[q]*dph;
      point_open(xr - xp, yr, t, alpha, beta, rho, s);
      for(k = 0; k < 3; k++)
        out[k] += s[k]*wt;
    }
  }
}

/* one time, no shift */
static void open_at(double xr, double yr, double t, double a, double alpha, double beta, double rho, double *st)
{
  double bp[FD_OPEN_MAXBP], d, xc, o[3], left[3], right[3], lam, Lc, dl, dr, tot[3];
  int nbp, k, sg;
  lam = rho*(alpha*alpha - 2.0*beta*beta);
  st[0] = 0.0;
  st[1] = 0.0;
  st[2] = 0.0;
  if(t <= 0.0)
    return;
  if(!(fabs(yr) < 0.5*a && fabs(xr) < a)){
    nbp = 0;
    if(fabs(yr) < 2.0*a){                        /* near the source line: refine around the receiver */
      d = fmax(fabs(yr), 1e-3*a);
      xc = fmin(fmax(xr, -a), a);
      bp[nbp++] = xc;
      for(k = 0; k < 8; k++)
        for(sg = -1; sg <= 1; sg += 2)
          bp[nbp++] = xc + sg*d*ldexp(1.0, k);
    }
    integrate_S(-a, a, xr, yr, t, alpha, beta, rho, bp, nbp, o);
    for(k = 0; k < 3; k++)
      st[k] = -o[k];
    return;
  }
  /* on the segment: complement identity */
  Lc = alpha*t*1.000001;
  dl = xr + a;
  dr = a - xr;
  nbp = 0;
  for(k = 0; k < 24; k++)
    bp[nbp++] = -a - dl*ldexp(1.0, k);
  integrate_S(xr - Lc, -a, xr, yr, t, alpha, beta, rho, bp, nbp, left);
  nbp = 0;
  for(k = 0; k < 24; k++)
    bp[nbp++] = a + dr*ldexp(1.0, k);
  integrate_S(a, xr + Lc, xr, yr, t, alpha, beta, rho, bp, nbp, right);
  tot[0] = (alpha*t > fabs(yr)) ? lam/(2.0*alpha) : 0.0;
  tot[1] = 0.0;
  tot[2] = (alpha*t > fabs(yr)) ? rho*alpha/2.0 : 0.0;
  for(k = 0; k < 3; k++)
    st[k] = -(tot[k] - left[k] - right[k]);
}

void fd_seg_ramp_inplane_open(double xr, double yr, double t, double a, double alpha, double beta, double rho,
    double tshift, double st[3])
{
  double lo[3], hi[3];
  int k;
  if(tshift > 0.0){
    open_at(xr, yr, t - tshift, a, alpha, beta, rho, lo);
    open_at(xr, yr, t + tshift, a, alpha, beta, rho, hi);
    for(k = 0; k < 3; k++)
      st[k] = 0.5*(lo[k] + hi[k]);
    return;
  }
  open_at(xr, yr, t, a, alpha, beta, rho, st);
}

void fd_static_open_cs(double x, double y, double a, double G, double nu, double *s);

/* full-plane static stress (sxx, sxy, syy) in the source frame from unit opening on |x'| <= a
   (Crouch and Starfield, physical sign; generated in fd_kernels_open_gen.c); the closed form has
   a removable 0/0 on the normals through the element ends, averaged over a tiny offset there */
void fd_static_inplane_open(double x, double y, double a, double mu, double nu, double st[3])
{
  double lo[3], hi[3], e;
  int k;
  fd_static_open_cs(x, y, a, mu, nu, st);
  if(isfinite(st[0]) && isfinite(st[1]) && isfinite(st[2]))
    return;
  e = 1e-7*a;
  fd_static_open_cs(x - e, y, a, mu, nu, lo);
  fd_static_open_cs(x + e, y, a, mu, nu, hi);
  for(k = 0; k < 3; k++)
    st[k] = 0.5*(lo[k] + hi[k]);
}
