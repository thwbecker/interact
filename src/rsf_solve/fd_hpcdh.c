/* see fd_hpcdh.h; the Python prototype is fully_dynamic/proto/halfplane_cdh.py */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "fd_hpcdh.h"
#include "fd_hpdwn.h"

/* Crouch and Starfield, interact's src/green/tdd_coeff.F */
extern void tdd_coeff_(double *x, double *y, double *cx, double *cy, double *a, double *cosb, double *sinb, int *ihalf,
                       double *pr, double *pr1, double *pr2, double *con, double *cons, double *chi,
                       double *sxxs, double *sxxn, double *syys, double *syyn, double *sxys, double *sxyn,
                       double *uxs, double *uxn, double *uys, double *uyn, int *iret);
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define CDH_NGRID 64      /* arrival times along the segment */
#define CDH_NQ 16         /* Gauss nodes per piece */
#define CDH_NCONT 12      /* continuation steps from the saddle */
#define CDH_TAIL 40.0     /* direct evaluation up to this many first arrivals */

struct cdh_med {
  double alpha, beta, rho, lam, mu;
};

struct cdh_geo {               /* one entry: receiver point, source segment */
  double xr, zr, xs, zs, tx, tz, a;
};

#define CDH_NPIECE 8
struct cdh_ws {                /* one entry's arrival grid and the previous lag's path points */
  double tg[CDH_NGRID + 1];
  double t0g[4][CDH_NGRID + 1];
  int npc[4];                                   /* pieces evaluated at the previous lag, -1 none */
  double complex qp[4][CDH_NPIECE][CDH_NQ];
  double complex dqp[4][CDH_NPIECE][CDH_NQ];    /* their dq/dt, for a linear prediction */
  double tp[4];                                 /* the time of the previous lag */
  unsigned char ok[4][CDH_NPIECE][CDH_NQ];
};

static const double cdh_gx[CDH_NQ] = {
  -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
  -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.0950125098376374,
   0.0950125098376374,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
   0.7554044083550030,  0.8656312023878318,  0.9445750230732326,  0.9894009349916499};
static const double cdh_gw[CDH_NQ] = {
   0.0271524594117541, 0.0622535239386479, 0.0951585116824928, 0.1246289712555339,
   0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.1894506104550685,
   0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
   0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541};

/* complex division without the overflow guards of the C library (__divdc3, half the time of the
   table builder): the moduli here stay far from overflow and underflow */
static inline double complex cdh_div(double complex a, double complex b)
{
  double br, bi, d;
  br = creal(b);
  bi = cimag(b);
  d = 1.0/(br*br + bi*bi);
  return CMPLX((creal(a)*br + cimag(a)*bi)*d, (cimag(a)*br - creal(a)*bi)*d);
}

/* sqrt(1/c^2 - q^2) with a non-negative real part */
static double complex cdh_eta(double complex q, double c)
{
  double complex e;
  e = csqrt(1.0/(c*c) - q*q);
  if(creal(e) < 0.0)
    e = -e;
  return e;
}

/* the up-going incident traction (txz, tzz) at the surface per unit point moment (tangent tx, tz)
   for incident wave c1 (0 P, 1 S), at slowness q, s = 1 (k = i q, omega = i) */
static void cdh_incident(const struct cdh_med *m, double complex q, int c1, double tx, double tz,
                         double complex *txz, double complex *tzz)
{
  double nx, nz, mm[2][2];
  double complex k, w, kb2, nu, Dr[2], Ds[2], g, Gip, u[2], exx, ezz;
  int i, p, r;
  k = I*q;
  w = I;
  kb2 = (w/m->beta)*(w/m->beta);
  nu = cdh_eta(q, c1 ? m->beta : m->alpha);
  nx = tz;
  nz = -tx;
  mm[0][0] = 2.0*tx*nx*m->mu;
  mm[0][1] = (tx*nz + tz*nx)*m->mu;
  mm[1][0] = mm[0][1];
  mm[1][1] = 2.0*tz*nz*m->mu;
  Dr[0] = I*k;
  Dr[1] = nu;
  Ds[0] = -I*k;
  Ds[1] = -nu;
  g = cdh_div(1.0, 2.0*nu);
  for(i = 0; i < 2; i++){
    u[i] = 0.0;
    for(p = 0; p < 2; p++){
      if(c1)
        Gip = -(kb2*(i == p ? 1.0 : 0.0) + Dr[i]*Dr[p])*g/m->rho;      /* w^2 = -1 */
      else
        Gip = (Dr[i]*Dr[p])*g/m->rho;
      for(r = 0; r < 2; r++)
        u[i] += mm[p][r]*Ds[r]*Gip;
    }
  }
  exx = Dr[0]*u[0];
  ezz = Dr[1]*u[1];
  *tzz = m->lam*exx + (m->lam + 2.0*m->mu)*ezz;
  *txz = m->mu*(Dr[1]*u[0] + Dr[0]*u[1]);
}

/* the stresses (sxx, sxz, szz) at the receiver of ray (c1, c2) at slowness q, without the
   exponential: incident c1 reflected as c2 (0 P, 1 S) */
static void cdh_ray(const struct cdh_med *m, double complex q, int c1, int c2, double tx, double tz, double complex sg[3])
{
  double complex k, na, nb, PxzA, PzzA, PxzB, PzzB, det, txz, tzz, A, ux, uz, nu;
  double lam, mu;
  lam = m->lam;
  mu = m->mu;
  k = I*q;
  na = cdh_eta(q, m->alpha);
  nb = cdh_eta(q, m->beta);
  PxzA = mu*(I*k*(-na) + (-na)*I*k);
  PzzA = lam*(I*k*I*k + na*na) + 2.0*mu*na*na;
  PxzB = mu*(I*k*I*k - nb*nb);
  PzzB = lam*(I*k*nb - nb*I*k) + 2.0*mu*(-nb*I*k);
  det = PxzA*PzzB - PxzB*PzzA;
  cdh_incident(m, q, c1, tx, tz, &txz, &tzz);
  if(c2 == 0){
    A = cdh_div(-txz*PzzB + tzz*PxzB, det);
    ux = I*k*A;
    uz = -na*A;
    nu = na;
  }else{
    A = cdh_div(-tzz*PxzA + txz*PzzA, det);
    ux = nb*A;
    uz = I*k*A;
    nu = nb;
  }
  sg[0] = (lam + 2.0*mu)*(I*k*ux) + lam*(-nu*uz);
  sg[1] = mu*(-nu*ux + I*k*uz);
  sg[2] = lam*(I*k*ux) + (lam + 2.0*mu)*(-nu*uz);
}

/* the travel-time function tau(q) = q X + eta1 h1 + eta2 h2 of a ray (source depth h1 at speed
   c1, receiver depth h2 at speed c2, horizontal offset X); the Cagniard path is tau(q) = t */
static double complex cdh_tau(double complex q, double X, double h1, double c1, double h2, double c2)
{
  return q*X + cdh_eta(q, c1)*h1 + cdh_eta(q, c2)*h2;
}

/* d tau / d q (zero at the Snell-law saddle) */
static double complex cdh_dtau(double complex q, double X, double h1, double c1, double h2, double c2)
{
  return X - cdh_div(q*h1, cdh_eta(q, c1)) - cdh_div(q*h2, cdh_eta(q, c2));
}

/* the real stationary point (Snell's law) and the arrival time of a ray: Newton on tau'(q) = 0
   (tau' decreases in q), safeguarded by the bracket */
static void cdh_saddle(double X, double h1, double c1, double h2, double c2, double *qs, double *ts, double *d2)
{
  double b, lo, hi, qn, f, fp, e1, e2;
  int it;
  b = 1.0/fmax(c1, c2)*(1.0 - 1e-12);
  lo = -b;
  hi = b;
  qn = X/sqrt(X*X + (h1 + h2)*(h1 + h2))/fmax(c1, c2);       /* the straight-ray guess */
  for(it = 0; it < 100; it++){
    e1 = sqrt(1.0/(c1*c1) - qn*qn);
    e2 = sqrt(1.0/(c2*c2) - qn*qn);
    f = X - qn*h1/e1 - qn*h2/e2;
    if(fabs(f) < 1e-14*(fabs(X) + h1 + h2))
      break;                                     /* converged at qn */
    if(f > 0.0)
      lo = qn;
    else
      hi = qn;
    if(hi - lo < 1e-15*b){
      qn = 0.5*(lo + hi);
      break;
    }
    fp = -h1/(c1*c1*e1*e1*e1) - h2/(c2*c2*e2*e2*e2);
    qn = qn - f/fp;
    if((qn <= lo) || (qn >= hi))
      qn = 0.5*(lo + hi);
  }
  *qs = qn;
  *ts = creal(cdh_tau(qn, X, h1, c1, h2, c2));
  e1 = sqrt(1.0/(c1*c1) - qn*qn);
  e2 = sqrt(1.0/(c2*c2) - qn*qn);
  *d2 = -h1/(c1*c1*e1*e1*e1) - h2/(c2*c2*e2*e2*e2);
}

/* Newton on tau(q) = t from g (upper half plane); returns 1 if converged */
static int cdh_newton(double t, double X, double h1, double c1, double h2, double c2, double complex *g)
{
  double complex f, dg, e1, e2;
  int it;
  for(it = 0; it < 60; it++){
    e1 = cdh_eta(*g, c1);                       /* tau and tau' share the square roots */
    e2 = cdh_eta(*g, c2);
    f = *g*X + e1*h1 + e2*h2 - t;
    dg = cdh_div(f, X - cdh_div(*g*h1, e1) - cdh_div(*g*h2, e2));
    *g -= dg;
    if(cimag(*g) < 0.0)
      *g = creal(*g) + 1e-30*I;
    if(cabs(dg) < 1e-14*cabs(*g))
      return 1;
  }
  return 0;
}

/* q(t) and dq/dt on the Cagniard path for t above the arrival (returns 0 before it); from the
   guess *qw if warm (the previous lag's point), else by continuation from the saddle */
static int cdh_path(double t, double X, double h1, double c1, double h2, double c2, double complex *q, double complex *dq,
                    const double complex *qw)
{
  double qs, t0, d2, tt;
  double complex g;
  int k;
  g = 0.0;
  if(qw){                                       /* warm: the point has arrived already */
    g = *qw;
    if(cdh_newton(t, X, h1, c1, h2, c2, &g) && (cimag(g) > 0.0)){
      *q = g;
      *dq = cdh_div(1.0, cdh_dtau(g, X, h1, c1, h2, c2));
      return 1;
    }
  }
  cdh_saddle(X, h1, c1, h2, c2, &qs, &t0, &d2);
  if(t <= t0)
    return 0;
  for(k = 1; k <= CDH_NCONT; k++){
    tt = t0 + (t - t0)*((double)k/CDH_NCONT)*((double)k/CDH_NCONT);
    if(k == 1)
      g = qs + I*sqrt(2.0*(tt - t0)/(-d2));
    cdh_newton(tt, X, h1, c1, h2, c2, &g);
  }
  *q = g;
  *dq = cdh_div(1.0, cdh_dtau(g, X, h1, c1, h2, c2));
  return 1;
}

/* the arrival time of ray (c1, c2) from the point tau of the segment */
static double cdh_arrival(const struct cdh_med *m, const struct cdh_geo *e, double tau, int c1, int c2)
{
  double qs, ts, d2;
  cdh_saddle(e->xr - (e->xs + tau*e->tx), e->zs + tau*e->tz, c1 ? m->beta : m->alpha, e->zr, c2 ? m->beta : m->alpha, &qs, &ts, &d2);
  return ts;
}

/* a root of arrival(tau) = t between lo and hi (f(lo) and f(hi) of opposite sign or one of them 0) */
static double cdh_root(const struct cdh_med *m, const struct cdh_geo *e, double t, int c1, int c2, double lo, double hi)
{
  double flo, fhi, x, fx, tol;
  int it, side;
  flo = cdh_arrival(m, e, lo, c1, c2) - t;
  fhi = cdh_arrival(m, e, hi, c1, c2) - t;
  if(flo == 0.0)
    return lo;
  if(fhi == 0.0)
    return hi;
  tol = 1e-13*(e->a + fabs(lo) + fabs(hi));
  side = 0;
  x = 0.5*(lo + hi);
  for(it = 0; it < 100; it++){                  /* Illinois (regula falsi) */
    x = (lo*fhi - hi*flo)/(fhi - flo);
    fx = cdh_arrival(m, e, x, c1, c2) - t;
    if((fx > 0.0) == (fhi > 0.0)){
      hi = x;
      fhi = fx;
      if(side == -1)
        flo *= 0.5;
      side = -1;
    }else{
      lo = x;
      flo = fx;
      if(side == 1)
        fhi *= 0.5;
      side = 1;
    }
    if((hi - lo < tol) || (fx == 0.0))
      break;
  }
  return x;
}

/* the minimum of the arrival time between lo and hi (golden section) */
static double cdh_argmin(const struct cdh_med *m, const struct cdh_geo *e, int c1, int c2, double lo, double hi)
{
  double gr, x1, x2, f1, f2;
  int it;
  gr = 0.5*(sqrt(5.0) - 1.0);
  x1 = hi - gr*(hi - lo);
  x2 = lo + gr*(hi - lo);
  f1 = cdh_arrival(m, e, x1, c1, c2);
  f2 = cdh_arrival(m, e, x2, c1, c2);
  for(it = 0; it < 80; it++){
    if(f1 < f2){
      hi = x2;
      x2 = x1;
      f2 = f1;
      x1 = hi - gr*(hi - lo);
      f1 = cdh_arrival(m, e, x1, c1, c2);
    }else{
      lo = x1;
      x1 = x2;
      f1 = f2;
      x2 = lo + gr*(hi - lo);
      f2 = cdh_arrival(m, e, x2, c1, c2);
    }
  }
  return 0.5*(lo + hi);
}

/* the arrival grid of an entry (once per entry), the warm starts cleared */
static void cdh_ws_init(const struct cdh_med *m, const struct cdh_geo *e, struct cdh_ws *w)
{
  int i, r;
  for(i = 0; i <= CDH_NGRID; i++)
    w->tg[i] = -e->a + 2.0*e->a*i/CDH_NGRID;
  for(r = 0; r < 4; r++){
    for(i = 0; i <= CDH_NGRID; i++)
      w->t0g[r][i] = cdh_arrival(m, e, w->tg[i], r/2, r%2);
    w->npc[r] = -1;
  }
}

/* the breakpoints of the segment for time t: the ends and the points where the arrival time
   equals t, also around an interior minimum of the arrival time that the grid does not resolve */
static int cdh_breakpoints(const struct cdh_med *m, const struct cdh_geo *e, const struct cdh_ws *w, double t, int c1, int c2,
                           double *br)
{
  const double *tg;
  double fg0, fg1, fgm, tm, fmn;
  int i, nb, r;
  r = 2*c1 + c2;
  tg = w->tg;
  nb = 0;
  br[nb++] = -e->a;
  for(i = 0; i < CDH_NGRID; i++){
    fg0 = w->t0g[r][i] - t;
    fg1 = w->t0g[r][i + 1] - t;
    if((fg0 <= 0.0) != (fg1 <= 0.0))
      br[nb++] = cdh_root(m, e, t, c1, c2, tg[i], tg[i + 1]);
    else if(i > 0){
      fgm = w->t0g[r][i - 1] - t;
      if((fg0 < fgm) && (fg0 <= fg1) && (fg0 > 0.0)){
        /* a local minimum above t on the grid: it may dip below t between the grid points */
        tm = cdh_argmin(m, e, c1, c2, tg[i - 1], tg[i + 1]);
        fmn = cdh_arrival(m, e, tm, c1, c2) - t;
        if(fmn < 0.0){
          if(tm < tg[i]){
            br[nb++] = cdh_root(m, e, t, c1, c2, tg[i - 1], tm);
            br[nb++] = cdh_root(m, e, t, c1, c2, tm, tg[i]);
          }else{
            br[nb++] = cdh_root(m, e, t, c1, c2, tg[i], tm);
            br[nb++] = cdh_root(m, e, t, c1, c2, tm, tg[i + 1]);
          }
        }
      }
    }
  }
  br[nb++] = e->a;
  return nb;
}

/* one ray's segment response (sxx, sxz, szz) at time t; w holds the warm starts per piece and node */
static void cdh_seg_ray(const struct cdh_med *m, const struct cdh_geo *e, struct cdh_ws *w, double t, int c1, int c2, double out[3])
{
  double br[2*CDH_NGRID + 4];
  double ta, tb, th, tau, jac, X, h1;
  double complex q, dq, sg[3], qg;
  int i, nb, k, n, r, np_, warm, use;
  r = 2*c1 + c2;
  out[0] = out[1] = out[2] = 0.0;
  nb = cdh_breakpoints(m, e, w, t, c1, c2, br);
  np_ = 0;
  warm = (w->npc[r] == nb - 1);
  for(k = 0; k + 1 < nb; k++){
    ta = br[k];
    tb = br[k + 1];
    if((tb <= ta) || (cdh_arrival(m, e, 0.5*(ta + tb), c1, c2) >= t)){
      if(k < CDH_NPIECE)
        for(n = 0; n < CDH_NQ; n++)
          w->ok[r][k][n] = 0;
      continue;
    }
    for(n = 0; n < CDH_NQ; n++){
      th = 0.5*M_PI*(cdh_gx[n] + 1.0);
      tau = ta + (tb - ta)*0.5*(1.0 - cos(th));
      jac = (tb - ta)*0.5*sin(th)*0.5*M_PI;
      X = e->xr - (e->xs + tau*e->tx);
      h1 = e->zs + tau*e->tz;
      use = warm && (k < CDH_NPIECE) && w->ok[r][k][n];
      if(use)
        qg = w->qp[r][k][n] + w->dqp[r][k][n]*(t - w->tp[r]);    /* predicted from the previous lag */
      if(!cdh_path(t, X, h1, c1 ? m->beta : m->alpha, e->zr, c2 ? m->beta : m->alpha, &q, &dq,
                   use ? &qg : NULL)){
        if(k < CDH_NPIECE)
          w->ok[r][k][n] = 0;
        continue;
      }
      if(k < CDH_NPIECE){
        w->qp[r][k][n] = q;
        w->dqp[r][k][n] = dq;
        w->ok[r][k][n] = 1;
      }
      cdh_ray(m, q, c1, c2, e->tx, e->tz, sg);
      for(i = 0; i < 3; i++)
        out[i] += cdh_gw[n]*jac*cimag(sg[i]*dq)/M_PI;
    }
    np_++;
  }
  w->npc[r] = nb - 1;
  w->tp[r] = t;
}

/* the reflected stresses (sxx, sxz, szz) of the segment at time t: the sum of the four rays */
static void cdh_seg_total(const struct cdh_med *m, const struct cdh_geo *e, struct cdh_ws *w, double t, double out[3])
{
  double r[3];
  int c1, c2, i;
  out[0] = out[1] = out[2] = 0.0;
  for(c1 = 0; c1 < 2; c1++)
    for(c2 = 0; c2 < 2; c2++){
      cdh_seg_ray(m, e, w, t, c1, c2, r);
      for(i = 0; i < 3; i++)
        out[i] += r[i];
    }
}

/* fd_project's convention: receiver tangent (tx, tz) in x, z; y up */
static void cdh_project(const double sgm[3], double tx, double tz, double *ts, double *tn)
{
  double srx, sry, nrx, nry, sxx, sxy, syy, trx, try_;
  srx = tx;
  sry = -tz;
  nrx = -sry;
  nry = srx;
  sxx = sgm[0];
  sxy = -sgm[1];
  syy = sgm[2];
  trx = sxx*nrx + sxy*nry;
  try_ = sxy*nrx + syy*nry;
  *ts = trx*srx + try_*sry;
  *tn = trx*nrx + try_*nry;
}

/* the static reflected tractions of an entry (the tail's rate): Crouch and Starfield's half plane
   minus their full plane, in the x, z frame (their slip is u- - u+, ours the opposite) */
static void cdh_static(const struct cdh_med *m, const struct cdh_geo *e, double rtx, double rtz, double *kst, double *knt)
{
  double x, y, cx, cy, a, cb, sb, nu, pr1, pr2, con, cons, chi, sg[3], o[10], f[10];
  int ih, iret;
  nu = (m->alpha*m->alpha - 2.0*m->beta*m->beta)/(2.0*(m->alpha*m->alpha - m->beta*m->beta));
  pr1 = 1.0 - 2.0*nu;
  pr2 = 2.0*(1.0 - nu);
  con = 1.0/(4.0*M_PI*(1.0 - nu));
  cons = 2.0*m->mu;
  chi = 3.0 - 4.0*nu;
  x = e->xr;
  y = -e->zr;
  cx = e->xs;
  cy = -e->zs;
  a = e->a;
  cb = e->tx;
  sb = -e->tz;
  ih = 1;
  tdd_coeff_(&x, &y, &cx, &cy, &a, &cb, &sb, &ih, &nu, &pr1, &pr2, &con, &cons, &chi,
             o, o + 1, o + 2, o + 3, o + 4, o + 5, o + 6, o + 7, o + 8, o + 9, &iret);
  ih = 0;
  tdd_coeff_(&x, &y, &cx, &cy, &a, &cb, &sb, &ih, &nu, &pr1, &pr2, &con, &cons, &chi,
             f, f + 1, f + 2, f + 3, f + 4, f + 5, f + 6, f + 7, f + 8, f + 9, &iret);
  sg[0] = -(o[0] - f[0]);                      /* sxx */
  sg[1] = (o[4] - f[4]);                       /* sxz = -sxy */
  sg[2] = -(o[2] - f[2]);                      /* szz = syy */
  cdh_project(sg, rtx, rtz, kst, knt);
}

/* cubic Lagrange interpolation at lag l from the anchors ia[k0..k0+np-1] (np 2 to 4) */
static double cdh_lagrange(const int *ia, const double *y, int k0, int np_, int l)
{
  double v, wgt;
  int j, k;
  v = 0.0;
  for(j = 0; j < np_; j++){
    wgt = 1.0;
    for(k = 0; k < np_; k++)
      if(k != j)
        wgt *= (double)(l - ia[k0 + k])/(double)(ia[k0 + j] - ia[k0 + k]);
    v += wgt*y[ia[k0 + j]];
  }
  return v;
}

/* one entry: the lag series (shear, normal) with the late tail from the static rates. Lags up to
   the last arrival over the segment (all rays) plus two are evaluated directly; after it the
   response is smooth and only anchors are evaluated, their spacing growing with the time since
   the last arrival (CDH_SMOOTH of it, at most CDH_STRIDE lags; FD_HPCDH_EVERY=1 evaluates every
   lag), the lags between them by cubic interpolation */
#define CDH_SMOOTH 0.1
#define CDH_STRIDE 32
static void cdh_entry(const struct cdh_med *m, const struct cdh_geo *e, double rtx, double rtz, double dt, int nt,
                      double kst, double knt, double *os, double *on, int every)
{
  struct cdh_ws *w;
  double t0, tc, t1, tl, t, sg[3], s1, n1, s2, n2, bs, bn, cs_, cn;
  int i, r, l, lm, na, k, k0, np_, st;
  int *ia;
  w = (struct cdh_ws *)malloc(sizeof(struct cdh_ws));
  ia = (int *)malloc(sizeof(int)*(size_t)(nt + 1));
  cdh_ws_init(m, e, w);
  t0 = 1e300;
  tl = 0.0;
  for(i = 0; i <= CDH_NGRID; i++){
    if(w->t0g[0][i] < t0)
      t0 = w->t0g[0][i];
    for(r = 0; r < 4; r++)
      if(w->t0g[r][i] > tl)
        tl = w->t0g[r][i];
  }
  tc = CDH_TAIL*t0;
  t1 = 0.75*tc;
  bs = bn = cs_ = cn = 0.0;
  if((nt - 1)*dt > tc){
    cdh_seg_total(m, e, w, t1, sg);
    cdh_project(sg, rtx, rtz, &s1, &n1);
    cdh_seg_total(m, e, w, tc, sg);
    cdh_project(sg, rtx, rtz, &s2, &n2);
    s1 -= kst*t1;
    n1 -= knt*t1;
    s2 -= kst*tc;
    n2 -= knt*tc;
    bs = (s2 - s1)/(1.0/t1 - 1.0/tc);
    bn = (n2 - n1)/(1.0/t1 - 1.0/tc);
    cs_ = s2 + bs/tc;
    cn = n2 + bn/tc;
    for(i = 0; i < 4; i++)
      w->npc[i] = -1;                          /* no warm start from these times */
  }
  /* the last lag evaluated by CdH (before the tail) */
  lm = -1;
  for(l = 0; l < nt; l++)
    if((l*dt > t0) && (l*dt <= tc))
      lm = l;
  /* anchors: every lag up to the last arrival plus two, then growing steps */
  na = 0;
  st = 1;
  for(l = 0; l <= lm; l++){
    t = l*dt;
    if(t <= t0)
      continue;
    if(every || (t <= tl + 2.0*dt) || (na < 4) || (l == lm) || (l >= ia[na - 1] + st)){
      ia[na++] = l;
      cdh_seg_total(m, e, w, t, sg);
      cdh_project(sg, rtx, rtz, os + l, on + l);
    }
    st = (int)floor(CDH_SMOOTH*(t - tl)/dt);
    if(st < 1)
      st = 1;
    if(st > CDH_STRIDE)
      st = CDH_STRIDE;
  }
  /* the lags between anchors */
  for(k = 0; k + 1 < na; k++){
    if(ia[k + 1] - ia[k] < 2)
      continue;
    np_ = (na < 4) ? na : 4;                     /* four anchors, shifted inward at the ends */
    k0 = k - 1;
    if(k0 + np_ > na)
      k0 = na - np_;
    if(k0 < 0)
      k0 = 0;
    for(l = ia[k] + 1; l < ia[k + 1]; l++){
      os[l] = cdh_lagrange(ia, os, k0, np_, l);
      on[l] = cdh_lagrange(ia, on, k0, np_, l);
    }
  }
  for(l = 0; l < nt; l++){
    t = l*dt;
    if(t <= t0){
      os[l] = 0.0;
      on[l] = 0.0;
    }else if(t > tc){
      os[l] = kst*t + cs_ - bs/t;
      on[l] = knt*t + cn - bn/t;
    }
  }
  free(ia);
  free(w);
}

/* the batch (fd_hpcdh.h): per receiver-source entry the exact static (Crouch and Starfield) and the
   lag series (cdh_entry); entries are independent, so OpenMP over them (dynamic schedule, since
   their cost varies with the window), and the result does not depend on threads or MPI ranks */
int fd_hpcdh_batch(int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                   int nr, const double *rx, const double *rz, const double *rtx, const double *rtz,
                   int npr, const double *prx, const double *prz,
                   double alpha, double beta, double rho, double dt, int nt, double *os, double *on, double *ks, double *kn)
{
  struct cdh_med m;
  struct cdh_geo e;
  long p, np;
  int i, j, every;
  every = getenv("FD_HPCDH_EVERY") ? atoi(getenv("FD_HPCDH_EVERY")) : 0;
  m.alpha = alpha;
  m.beta = beta;
  m.rho = rho;
  m.mu = rho*beta*beta;
  m.lam = rho*(alpha*alpha - 2.0*beta*beta);
  np = (long)nr*ns;
#pragma omp parallel for schedule(dynamic, 1) private(p, i, j, e)
  for(p = 0; p < np; p++){
    i = (int)(p/ns);
    j = (int)(p % ns);
    e.xr = rx[i];
    e.zr = rz[i];
    e.xs = sx[j];
    e.zs = sz[j];
    e.tx = stx[j];
    e.tz = stz[j];
    e.a = sa[j];
    cdh_static(&m, &e, rtx[i], rtz[i], ks + (size_t)i*ns + j, kn + (size_t)i*ns + j);
    cdh_entry(&m, &e, rtx[i], rtz[i], dt, nt, ks[(size_t)i*ns + j], kn[(size_t)i*ns + j],
              os + (size_t)p*nt, on + (size_t)p*nt, every);
  }
  return 0;
}
