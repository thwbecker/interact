/*
  windowed dynamic-correction tables, see fd_tables.h. Port of
  shells.PairTable / inplane.W_table / antiplane.W_table_ap of the
  prototype; the test program test_fd_tables compares the result with
  tables dumped by the prototype.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd_kernels.h"
#include "fd_tables.h"

static double g_ctol;
static int g_cexp;
static void ct_compress_pair(struct fd_pair *p, const double *xr, const double *yr, double mu, double nu);
static void ct_free(struct fd_ctab *c);

int fd_segments_parallel(const struct fd_segment *a, const struct fd_segment *b)
{
  double tol = 1e-12;
  return fabs(a->dx - b->dx) < tol * a->dx && (a->s[0] * b->s[0] + a->s[1] * b->s[1]) > 1.0 - tol;
}

/* receiver centre at Toeplitz offset m, in the source frame */
static void offset_point(const struct fd_segment *src, const struct fd_segment *rcv, int m, double *xr, double *yr, double Xg[2])
{
  double d0[2] = {rcv->X0[0] - src->X0[0], rcv->X0[1] - src->X0[1]};
  double base = d0[0] * src->s[0] + d0[1] * src->s[1];
  double perp[2] = {d0[0] - base * src->s[0], d0[1] - base * src->s[1]};
  double ns[2] = {-src->s[1], src->s[0]};
  Xg[0] = src->X0[0] + (base + m * src->dx) * src->s[0] + perp[0];
  Xg[1] = src->X0[1] + (base + m * src->dx) * src->s[1] + perp[1];
  *xr = (Xg[0] - src->X0[0]) * src->s[0] + (Xg[1] - src->X0[1]) * src->s[1];
  *yr = (Xg[0] - src->X0[0]) * ns[0] + (Xg[1] - src->X0[1]) * ns[1];
}

struct fd_pair *fd_pair_build(const struct fd_segment *src, const struct fd_segment *rcv, int mode,
    double dt, double et, double alpha, double beta, double rho,
    double tol, int lmax, int cap_ok, int M, char *errmsg)
{
  struct fd_pair *p;
  int ns = src->n, nr = rcv->n, nm = ns + nr - 1, i, l, lsw, lc, L, idx;
  double mu = rho * beta * beta, nu = (alpha * alpha - 2.0 * beta * beta) / (2.0 * (alpha * alpha - beta * beta));
  double xr, yr, Xg[2], rr, env, dmax = 0.0, rfar, gamma_switch = 10.0;
  double *Ws, *Wn, s0[3], s1[3], d[3], ts, tn, st[3], a2[2], b2[2], ns_src[2], nr_g[2];
  double _Complex *col;
  const double *ss = src->s, *sr = rcv->s;
  if(errmsg)
    errmsg[0] = 0;
  if(!fd_segments_parallel(src, rcv) || src->kind != FD_SEG_FAULT || rcv->kind != FD_SEG_FAULT || getenv("FD_FORCE_DENSE"))
    return fd_pair_build_dense(src, rcv, mode, dt, et, alpha, beta, rho, tol, lmax, cap_ok, M, errmsg);
  if(lmax <= 0){
    /* generous bound: S arrival at the farthest pair times (1 + 2/sqrt(tol)) */
    for(i = -(ns - 1); i <= nr - 1; i++){
      offset_point(src, rcv, i, &xr, &yr, Xg);
      if(hypot(xr, yr) > dmax)
        dmax = hypot(xr, yr);
    }
    lmax = (int)ceil((dmax + src->dx) / beta / dt * (1.0 + 2.0 / sqrt(tol))) + 2;
  }
  p = (struct fd_pair *)calloc(1, sizeof(*p));
  p->src = src;
  p->rcv = rcv;
  p->mode = mode;
  p->dt = dt;
  p->et = et;
  p->nm = nm;
  p->has_normal = (mode == FD_MODE_INPLANE);
  p->m = (int *)malloc(sizeof(int) * nm);
  p->ks_row = (double *)malloc(sizeof(double) * nm);
  p->kn_row = (double *)malloc(sizeof(double) * nm);
  Ws = (double *)malloc(sizeof(double) * (size_t)nm * lmax);
  Wn = (double *)calloc((size_t)nm * lmax, sizeof(double));
  /* receiver normal in the source frame (antiplane traction projection) */
  nr_g[0] = -sr[1];
  nr_g[1] = sr[0];
  ns_src[0] = nr_g[0] * ss[0] + nr_g[1] * ss[1];
  ns_src[1] = nr_g[0] * (-ss[1]) + nr_g[1] * ss[0];
  for(i = 0; i < nm; i++){
    p->m[i] = i - (ns - 1);
    offset_point(src, rcv, p->m[i], &xr, &yr, Xg);
    if(mode == FD_MODE_INPLANE){
      fd_static_inplane(xr, yr, src->a, mu, nu, st);
      fd_project(st[0], st[1], st[2], ss, sr, &p->ks_row[i], &p->kn_row[i]);
      rfar = fmax(hypot(xr + src->a, yr), hypot(xr - src->a, yr));
      lsw = (int)ceil(gamma_switch * rfar / beta / dt - et) + 2;      /* first column using the series */
      if(lsw > lmax)
        lsw = lmax;
      lc = lsw;
      /* closed form for columns 0..lc-1: levels 0, et dt, (et+1) dt, ... */
      fd_seg_ramp_inplane(xr, yr, 0.0, src->a, alpha, beta, rho, s0);
      for(l = 0; l < lc; l++){
        fd_seg_ramp_inplane(xr, yr, (et + l) * dt, src->a, alpha, beta, rho, s1);
        d[0] = s1[0] - s0[0];
        d[1] = s1[1] - s0[1];
        d[2] = s1[2] - s0[2];
        fd_project(d[0], d[1], d[2], ss, sr, &ts, &tn);
        Ws[(size_t)i * lmax + l] = ts - p->ks_row[i] * dt;
        Wn[(size_t)i * lmax + l] = tn - p->kn_row[i] * dt;
        memcpy(s0, s1, sizeof(s0));
      }
      /* series beyond: column l-1 (1-based lag l >= lsw) covers [(et+l-2) dt, (et+l-1) dt] */
      for(l = lsw; l <= lmax; l++){
        fd_series_inplane(xr, yr, (l - 2 + et) * dt, dt, src->a, alpha, beta, rho, d);
        fd_project(d[0], d[1], d[2], ss, sr, &ts, &tn);
        Ws[(size_t)i * lmax + (l - 1)] = ts;
        Wn[(size_t)i * lmax + (l - 1)] = tn;
      }
    }else{
      fd_static_antiplane(xr, yr, src->a, mu, a2);
      p->ks_row[i] = ns_src[0] * a2[0] + ns_src[1] * a2[1];
      p->kn_row[i] = 0.0;
      fd_seg_ramp_antiplane(xr, yr, 0.0, src->a, beta, rho, a2);
      for(l = 0; l < lmax; l++){
        fd_seg_ramp_antiplane(xr, yr, (et + l) * dt, src->a, beta, rho, b2);
        Ws[(size_t)i * lmax + l] = ns_src[0] * (b2[0] - a2[0]) + ns_src[1] * (b2[1] - a2[1]) - p->ks_row[i] * dt;
        a2[0] = b2[0];
        a2[1] = b2[1];
      }
    }
  }
  /* window: last column where |W| exceeds tol times the static interaction envelope */
  L = 0;
  for(i = 0; i < nm; i++){
    offset_point(src, rcv, p->m[i], &xr, &yr, Xg);
    rr = fmax(hypot(xr, yr), src->a);
    env = mu * src->dx * src->dx / (2.0 * M_PI * (mode == FD_MODE_INPLANE ? (1.0 - nu) : 1.0) * rr * rr) * dt;
    for(l = lmax - 1; l >= 0; l--)
      if(fabs(Ws[(size_t)i * lmax + l]) > tol * env || fabs(Wn[(size_t)i * lmax + l]) > tol * env)
        break;
    if(l + 1 > L)
      L = l + 1;
  }
  if(L < 1)
    L = 1;
  if(L >= lmax){
    if(!cap_ok){
      if(errmsg)
        sprintf(errmsg, "window reached lmax = %d; increase lmax", lmax);
      free(Ws);
      free(Wn);
      fd_pair_free(p);
      return NULL;
    }
    L = lmax - 1;
  }
  p->L = L;
  p->M = M > 0 ? M : fd_fft_size_for(nm);
  p->pf = fd_fft_plan_create(p->M, FD_FFT_FORWARD);
  p->pb = fd_fft_plan_create(p->M, FD_FFT_BACKWARD);
  p->Ws_hat = (double _Complex *)malloc(sizeof(double _Complex) * (size_t)L * p->M);
  p->Wn_hat = p->has_normal ? (double _Complex *)malloc(sizeof(double _Complex) * (size_t)L * p->M) : NULL;
  col = (double _Complex *)malloc(sizeof(double _Complex) * p->M);
  for(l = 0; l < L; l++){
    memset(col, 0, sizeof(double _Complex) * p->M);
    for(i = 0; i < nm; i++){
      idx = p->m[i] >= 0 ? p->m[i] : p->M + p->m[i];
      col[idx] = Ws[(size_t)i * lmax + l];
    }
    fd_fft_exec(p->pf, col, p->Ws_hat + (size_t)l * p->M);
    if(p->has_normal){
      memset(col, 0, sizeof(double _Complex) * p->M);
      for(i = 0; i < nm; i++){
        idx = p->m[i] >= 0 ? p->m[i] : p->M + p->m[i];
        col[idx] = Wn[(size_t)i * lmax + l];
      }
      fd_fft_exec(p->pf, col, p->Wn_hat + (size_t)l * p->M);
    }
  }
  free(col);
  free(Ws);
  free(Wn);
  return p;
}


/* ---------------- surface elements (fd_plan.md M7 step 2) ----------------
   Opening sources use the quadrature kernel of fd_kernels_open.c, which is expensive, so the
   stress histories are cached by the receiver's position in the source frame (on a straight
   surface the relative positions repeat) and evaluated only until the window criterion has held
   for FD_OPEN_QUIET lags after the S arrival from the far end of the source. */
#define FD_OPEN_QUIET 32
/* the symmetric time average of the opening kernel (as the prototype's tshift), 1.7e-3 of the
   step rather than the prototype's 1e-4: with surface elements of a whole number of S steps the
   fronts cross element ends, or touch the source line, exactly at sample times, and the
   quadrature is ill-conditioned within beta tshift of such a coincidence. The values are stable
   (to 3e-4 relative) for shifts of 7e-4 to 5e-3 in the cases found (an on-line receiver as an S
   front passes a segment end, a front tangent to the source line); 1e-4 was 2 percent off in
   the second case and 1e-3 met another coincidence (a spike of 1e12 Pa) */
#define FD_OPEN_TSHIFT 1.7e-3
struct open_cache_entry {
  long kx;
  long ky;
  long ka;
  int n;                        /* lags evaluated (values beyond are zero increments) */
  double *s;                    /* 3 x (n + 1) stresses at (et + l - 1) dt, l = 0..n, s[0] at t = 0 */
};
static struct open_cache_entry *g_oc = NULL;
static size_t g_oc_n = 0, g_oc_cap = 0;

void fd_tables_open_cache_clear(void)
{
  size_t i;
  for(i = 0; i < g_oc_cap; i++)
    free(g_oc[i].s);
  free(g_oc);
  g_oc = NULL;
  g_oc_n = 0;
  g_oc_cap = 0;
}
static size_t oc_hash(long kx, long ky, long ka, size_t cap)
{
  unsigned long h;
  h = (unsigned long)kx*1000003UL ^ (unsigned long)ky*998244353UL ^ (unsigned long)ka*19260817UL;
  return (size_t)(h % cap);
}
static struct open_cache_entry *oc_find(long kx, long ky, long ka, int insert)
{
  size_t i, j, cap;
  struct open_cache_entry *old;
  if(g_oc_cap == 0 || (insert && 2*(g_oc_n + 1) > g_oc_cap)){
    old = g_oc;
    cap = g_oc_cap;
    g_oc_cap = cap ? 2*cap : 1024;
    g_oc = (struct open_cache_entry *)calloc(g_oc_cap, sizeof(*g_oc));
    for(i = 0; i < cap; i++)
      if(old[i].s){
        j = oc_hash(old[i].kx, old[i].ky, old[i].ka, g_oc_cap);
        while(g_oc[j].s)
          j = (j + 1) % g_oc_cap;
        g_oc[j] = old[i];
      }
    free(old);
  }
  i = oc_hash(kx, ky, ka, g_oc_cap);
  while(g_oc[i].s){
    if(g_oc[i].kx == kx && g_oc[i].ky == ky && g_oc[i].ka == ka)
      return g_oc + i;
    i = (i + 1) % g_oc_cap;
  }
  if(!insert)
    return NULL;
  g_oc[i].kx = kx;
  g_oc[i].ky = ky;
  g_oc[i].ka = ka;
  g_oc_n++;
  return g_oc + i;
}

/* one dense entry with a surface element involved (in-plane): the source's stress history
   (slip: closed-form ramp and series as below; opening: quadrature, cached), projected on the
   receiver's primary row (shear, or normal for an opening receiver) into W and the normal row
   into Wn; K: static values of both rows. Increments of lags l = 1..lmax as W[l - 1]. */
static void dense_entry_open(const struct fd_segment *src, const struct fd_segment *rcv, double x, double y,
    double dt, double et, double alpha, double beta, double rho, double mu, double nu, double tol, int lmax,
    double *W, double *Wn, double *Kp, double *Kn)
{
  double st[3], s0[3], s1[3], d[3], ts, tn, kps, kpn, rfar, env, rr;
  int l, k, quiet, lS, nl;
  long kx;
  long ky;
  long ka;
  struct open_cache_entry *ce;
  double *hist;
  /* static */
  if(src->kind == FD_SEG_SURF_OPEN)
    fd_static_inplane_open(x, y, src->a, mu, nu, st);
  else
    fd_static_inplane(x, y, src->a, mu, nu, st);
  fd_project(st[0], st[1], st[2], src->s, rcv->s, &kps, &kpn);
  *Kp = (rcv->kind == FD_SEG_SURF_OPEN) ? kpn : kps;
  *Kn = kpn;
  rr = fmax(hypot(x, y), src->a);
  env = mu*src->dx*src->dx/(2.0*M_PI*(1.0 - nu)*rr*rr)*dt;
  rfar = fmax(hypot(x + src->a, y), hypot(x - src->a, y));
  lS = (int)ceil(rfar/beta/dt) + 2;
  /* stress history in the source frame, from the cache for opening sources */
  hist = NULL;
  nl = 0;
  ce = NULL;
  kx = 0;
  ky = 0;
  ka = 0;
  if(src->kind == FD_SEG_SURF_OPEN){
    kx = lround(x*1e3);
    ky = lround(y*1e3);
    ka = lround(src->a*1e3);
    ce = oc_find(kx, ky, ka, 0);
    if(ce){
      hist = ce->s;
      nl = ce->n;
    }
  }
  if(!hist){
    hist = (double *)malloc(sizeof(double)*3*(size_t)(lmax + 1));
    hist[0] = hist[1] = hist[2] = 0.0;
    quiet = 0;
    for(l = 1; l <= lmax; l++){
      if(src->kind == FD_SEG_SURF_OPEN)
        fd_seg_ramp_inplane_open(x, y, (et + l - 1)*dt, src->a, alpha, beta, rho, FD_OPEN_TSHIFT*dt, s1);
      else
        fd_seg_ramp_inplane(x, y, (et + l - 1)*dt, src->a, alpha, beta, rho, s1);
      for(k = 0; k < 3; k++)
        hist[3*l + k] = s1[k];
      nl = l;
      /* stop once the increments minus the static slope stay below the window criterion */
      for(k = 0; k < 3; k++){
        d[k] = s1[k] - hist[3*(l - 1) + k];
        s0[k] = st[k]*dt;
      }
      fd_project(d[0] - s0[0], d[1] - s0[1], d[2] - s0[2], src->s, rcv->s, &ts, &tn);
      if(l > lS && fabs(ts) <= tol*env && fabs(tn) <= tol*env)
        quiet++;
      else
        quiet = 0;
      if(quiet >= FD_OPEN_QUIET)
        break;
    }
    if(src->kind == FD_SEG_SURF_OPEN){
      ce = oc_find(kx, ky, ka, 1);
      ce->s = hist;
      ce->n = nl;
    }
  }
  for(l = 1; l <= lmax; l++){
    if(l <= nl){
      for(k = 0; k < 3; k++)
        d[k] = hist[3*l + k] - hist[3*(l - 1) + k];
    }else{
      for(k = 0; k < 3; k++)
        d[k] = st[k]*dt;                       /* the static rate: zero increment after the window */
    }
    fd_project(d[0], d[1], d[2], src->s, rcv->s, &ts, &tn);
    W[l - 1] = ((rcv->kind == FD_SEG_SURF_OPEN) ? tn : ts) - (*Kp)*dt;
    Wn[l - 1] = tn - kpn*dt;
  }
  if(src->kind != FD_SEG_SURF_OPEN)
    free(hist);
}

/* dense pair (fd_plan.md N1): the same columns as fd_pair_build for every (receiver, source)
   element pair, receiver centre in the source element's frame */
struct fd_pair *fd_pair_build_dense(const struct fd_segment *src, const struct fd_segment *rcv, int mode,
    double dt, double et, double alpha, double beta, double rho,
    double tol, int lmax, int cap_ok, int M, char *errmsg)
{
  struct fd_pair *p;
  int ns = src->n, nr = rcv->n, i, j, l, lsw, L;
  size_t np = (size_t)ns*nr, e;
  double mu = rho*beta*beta, nu = (alpha*alpha - 2.0*beta*beta)/(2.0*(alpha*alpha - beta*beta));
  double rfar, gamma_switch = 10.0, s0[3], s1[3], d[3], st[3], ts, tn, a2[2], b2[2], ns_src[2], nr_g[2], dmax = 0.0;
  double *Ws, *Wn, *xr, *yr;
  const double *ss = src->s, *sr = rcv->s;
  double Xr[2];
  double Xs[2];
  double dd[2];
  double x;
  double y;
  double rr;
  double env;
  if(errmsg)
    errmsg[0] = 0;
  xr = (double *)malloc(sizeof(double)*np);
  yr = (double *)malloc(sizeof(double)*np);
  for(i = 0; i < nr; i++)
    for(j = 0; j < ns; j++){
      Xr[0] = rcv->X0[0] + i*rcv->dx*sr[0];
      Xr[1] = rcv->X0[1] + i*rcv->dx*sr[1];
      Xs[0] = src->X0[0] + j*src->dx*ss[0];
      Xs[1] = src->X0[1] + j*src->dx*ss[1];
      dd[0] = Xr[0] - Xs[0];
      dd[1] = Xr[1] - Xs[1];
      e = (size_t)i*ns + j;
      xr[e] = dd[0]*ss[0] + dd[1]*ss[1];
      yr[e] = -dd[0]*ss[1] + dd[1]*ss[0];
      if(hypot(xr[e], yr[e]) > dmax)
        dmax = hypot(xr[e], yr[e]);
    }
  if(lmax <= 0)
    lmax = (int)ceil((dmax + src->dx)/beta/dt*(1.0 + 2.0/sqrt(tol))) + 2;
  p = (struct fd_pair *)calloc(1, sizeof(*p));
  p->src = src;
  p->rcv = rcv;
  p->mode = mode;
  p->dt = dt;
  p->et = et;
  p->dense = 1;
  p->has_normal = (mode == FD_MODE_INPLANE);     /* (surface receivers do not need it; kept for uniformity) */
  p->Kd_s = (double *)malloc(sizeof(double)*np);
  p->Kd_n = (double *)calloc(np, sizeof(double));
  Ws = (double *)calloc(np*lmax, sizeof(double));
  Wn = (double *)calloc(np*lmax, sizeof(double));
  nr_g[0] = -sr[1];
  nr_g[1] = sr[0];
  ns_src[0] = nr_g[0]*ss[0] + nr_g[1]*ss[1];
  ns_src[1] = nr_g[0]*(-ss[1]) + nr_g[1]*ss[0];
  for(e = 0; e < np; e++){
    x = xr[e];
    y = yr[e];
    if(mode == FD_MODE_INPLANE && (src->kind == FD_SEG_SURF_OPEN || rcv->kind == FD_SEG_SURF_OPEN)){
      /* surface elements (M7 step 2): opening sources and/or normal-traction receivers */
      dense_entry_open(src, rcv, x, y, dt, et, alpha, beta, rho, mu, nu, tol, lmax, Ws + e*lmax, Wn + e*lmax, &p->Kd_s[e], &p->Kd_n[e]);
      continue;
    }
    if(mode == FD_MODE_INPLANE){
      fd_static_inplane(x, y, src->a, mu, nu, st);
      fd_project(st[0], st[1], st[2], ss, sr, &p->Kd_s[e], &p->Kd_n[e]);
      rfar = fmax(hypot(x + src->a, y), hypot(x - src->a, y));
      lsw = (int)ceil(gamma_switch*rfar/beta/dt - et) + 2;
      if(lsw > lmax)
        lsw = lmax;
      fd_seg_ramp_inplane(x, y, 0.0, src->a, alpha, beta, rho, s0);
      for(l = 0; l < lsw; l++){
        fd_seg_ramp_inplane(x, y, (et + l)*dt, src->a, alpha, beta, rho, s1);
        d[0] = s1[0] - s0[0];
        d[1] = s1[1] - s0[1];
        d[2] = s1[2] - s0[2];
        fd_project(d[0], d[1], d[2], ss, sr, &ts, &tn);
        Ws[e*lmax + l] = ts - p->Kd_s[e]*dt;
        Wn[e*lmax + l] = tn - p->Kd_n[e]*dt;
        memcpy(s0, s1, sizeof(s0));
      }
      for(l = lsw; l <= lmax; l++){
        fd_series_inplane(x, y, (l - 2 + et)*dt, dt, src->a, alpha, beta, rho, d);
        fd_project(d[0], d[1], d[2], ss, sr, &ts, &tn);
        Ws[e*lmax + (l - 1)] = ts;
        Wn[e*lmax + (l - 1)] = tn;
      }
    }else{
      fd_static_antiplane(x, y, src->a, mu, a2);
      p->Kd_s[e] = ns_src[0]*a2[0] + ns_src[1]*a2[1];
      fd_seg_ramp_antiplane(x, y, 0.0, src->a, beta, rho, a2);
      for(l = 0; l < lmax; l++){
        fd_seg_ramp_antiplane(x, y, (et + l)*dt, src->a, beta, rho, b2);
        Ws[e*lmax + l] = ns_src[0]*(b2[0] - a2[0]) + ns_src[1]*(b2[1] - a2[1]) - p->Kd_s[e]*dt;
        a2[0] = b2[0];
        a2[1] = b2[1];
      }
    }
  }
  L = 0;
  for(e = 0; e < np; e++){
    rr = fmax(hypot(xr[e], yr[e]), src->a);
    env = mu*src->dx*src->dx/(2.0*M_PI*(mode == FD_MODE_INPLANE ? (1.0 - nu) : 1.0)*rr*rr)*dt;
    for(l = lmax - 1; l >= 0; l--)
      if(fabs(Ws[e*lmax + l]) > tol*env || fabs(Wn[e*lmax + l]) > tol*env)
        break;
    if(l + 1 > L)
      L = l + 1;
  }
  if(L < 1)
    L = 1;
  if(L >= lmax){
    if(!cap_ok){
      if(errmsg)
        sprintf(errmsg, "window reached lmax = %d; increase lmax", lmax);
      free(Ws);
      free(Wn);
      free(xr);
      free(yr);
      fd_pair_free(p);
      return NULL;
    }
    L = lmax - 1;
  }
  p->L = L;
  p->M = M;
  p->Wd_s = (double *)malloc(sizeof(double)*(size_t)L*np);
  p->Wd_n = p->has_normal ? (double *)malloc(sizeof(double)*(size_t)L*np) : NULL;
  for(l = 0; l < L; l++)
    for(e = 0; e < np; e++){
      p->Wd_s[(size_t)l*np + e] = Ws[e*lmax + l];
      if(p->Wd_n)
        p->Wd_n[(size_t)l*np + e] = Wn[e*lmax + l];
    }
  free(Ws);
  free(Wn);
  if(g_ctol > 0.0 && !getenv("FD_CT_OFF"))
    ct_compress_pair(p, xr, yr, mu, nu);
  free(xr);
  free(yr);
  return p;
}

void fd_pair_free(struct fd_pair *p)
{
  if(!p)
    return;
  free(p->m);
  free(p->ks_row);
  free(p->kn_row);
  free(p->Ws_hat);
  free(p->Wn_hat);
  free(p->Wd_s);
  free(p->Wd_n);
  free(p->Kd_s);
  free(p->Kd_n);
  if(p->compressed){
    ct_free(&p->cs);
    ct_free(&p->cn);
  }
  if(!p->dense){
    fd_fft_plan_destroy(p->pf);
    fd_fft_plan_destroy(p->pb);
  }
  free(p);
}

static int row_index(const struct fd_pair *p, int i, int j)
{
  return (i - j) + (p->src->n - 1);
}

void fd_pair_static_dense(const struct fd_pair *p, double *Ks, double *Kn)
{
  int i, j, ns = p->src->n, nr = p->rcv->n;
  if(p->dense){
    memcpy(Ks, p->Kd_s, sizeof(double)*(size_t)ns*nr);
    if(Kn)
      memcpy(Kn, p->Kd_n, sizeof(double)*(size_t)ns*nr);
    return;
  }
  for(i = 0; i < nr; i++)
    for(j = 0; j < ns; j++){
      Ks[i * ns + j] = p->ks_row[row_index(p, i, j)];
      if(Kn)
        Kn[i * ns + j] = p->kn_row[row_index(p, i, j)];
    }
}

void fd_pair_W_dense(const struct fd_pair *p, int l1, double *Ws, double *Wn)
{
  int i, j, k, ns = p->src->n, nr = p->rcv->n, m, idx;
  size_t np;
  size_t e;
  if(p->dense && p->compressed){
    np = (size_t)ns*nr;
    for(e = 0; e < np; e++){
      Ws[e] = fd_ctab_value(&p->cs, e, l1);
      if(Wn)
        Wn[e] = p->has_normal ? fd_ctab_value(&p->cn, e, l1) : 0.0;
    }
    return;
  }
  if(p->dense){
    np = (size_t)ns*nr;
    memcpy(Ws, p->Wd_s + (size_t)(l1 - 1)*np, sizeof(double)*np);
    if(Wn){
      if(p->Wd_n)
        memcpy(Wn, p->Wd_n + (size_t)(l1 - 1)*np, sizeof(double)*np);
      else
        memset(Wn, 0, sizeof(double)*np);
    }
    return;
  }
  double _Complex *c = (double _Complex *)malloc(sizeof(double _Complex) * p->M);
  double _Complex *cn = p->has_normal ? (double _Complex *)malloc(sizeof(double _Complex) * p->M) : NULL;
  fd_fft_exec(p->pb, p->Ws_hat + (size_t)(l1 - 1) * p->M, c);
  if(cn)
    fd_fft_exec(p->pb, p->Wn_hat + (size_t)(l1 - 1) * p->M, cn);
  for(k = 0; k < p->M; k++){
    c[k] /= p->M;
    if(cn)
      cn[k] /= p->M;
  }
  for(i = 0; i < nr; i++)
    for(j = 0; j < ns; j++){
      m = i - j;
      idx = m >= 0 ? m : p->M + m;
      Ws[i * ns + j] = creal(c[idx]);
      if(Wn)
        Wn[i * ns + j] = cn ? creal(cn[idx]) : 0.0;
    }
  free(c);
  free(cn);
}

/* ---------------- compressed dense tables (fd_plan.md D2) ----------------
   The dense lag tables of non-parallel pairs hold, per (receiver, source) element pair, the
   increments W of all lags of the window. Between and after the wavefront arrivals they are
   smooth in the lag (a constant before the P arrival, smooth between P and S, the slowly
   decaying 2-D tail after S), so runs of lags are replaced by least-squares quadratics, each
   verified at every lag against the full increments to ctol times the entry's envelope (the
   normalisation of the window criterion), and the history of a run uses prefix moments of the
   source's history ring, so its cost does not depend on the run's length. */
static double g_ctol = 1e-7;                    /* default: 1e-5 of the default window tolerance of rsf_fd */

static int g_cexp = 1;                          /* lags 1..g_cexp stay explicit (FD_CT_EXPLICIT) */

void fd_tables_set_compress(double ctol)
{
  g_ctol = ctol;
  if(getenv("FD_CT_EXPLICIT"))
    g_cexp = atoi(getenv("FD_CT_EXPLICIT"));
}

void fd_tables_set_compress_explicit(int nexp)
{
  g_cexp = getenv("FD_CT_EXPLICIT") ? atoi(getenv("FD_CT_EXPLICIT")) : nexp;
}

double fd_tables_get_compress(void)
{
  return g_ctol;
}

/* least-squares quadratic through w[a..b] (index l is lag l + 1) in a local basis, converted to
   1, y, y^2 with y = (l + 1)/FD_CT_LNORM; 0 if a point misses by more than tol */
static int ct_fit(const double *w, int a, int b, double tol, double *cf)
{
  double lm, h, t, A[3][4], f, piv, al, be, y, r;
  double e[3];
  int l, i, j, k, ip;
  lm = 0.5*(a + b);
  h = 0.5*(b - a);
  for(i = 0; i < 3; i++)
    for(j = 0; j < 4; j++)
      A[i][j] = 0.0;
  for(l = a; l <= b; l++){
    t = (l - lm)/h;
    A[0][0] += 1.0;
    A[0][1] += t;
    A[0][2] += t*t;
    A[1][2] += t*t*t;
    A[2][2] += t*t*t*t;
    A[0][3] += w[l];
    A[1][3] += t*w[l];
    A[2][3] += t*t*w[l];
  }
  A[1][0] = A[0][1];
  A[1][1] = A[0][2];
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];
  for(k = 0; k < 3; k++){
    ip = k;
    for(i = k + 1; i < 3; i++)
      if(fabs(A[i][k]) > fabs(A[ip][k]))
        ip = i;
    for(j = 0; j < 4; j++){
      f = A[k][j];
      A[k][j] = A[ip][j];
      A[ip][j] = f;
    }
    piv = A[k][k];
    if(piv == 0.0)
      return 0;
    for(i = k + 1; i < 3; i++){
      f = A[i][k]/piv;
      for(j = k; j < 4; j++)
        A[i][j] -= f*A[k][j];
    }
  }
  for(k = 2; k >= 0; k--){
    f = A[k][3];
    for(j = k + 1; j < 3; j++)
      f -= A[k][j]*e[j];
    e[k] = f/A[k][k];
  }
  al = FD_CT_LNORM/h;
  be = -(lm + 1.0)/h;
  cf[0] = e[0] + e[1]*be + e[2]*be*be;
  cf[1] = e[1]*al + 2.0*e[2]*al*be;
  cf[2] = e[2]*al*al;
  for(l = a; l <= b; l++){
    y = (l + 1)/FD_CT_LNORM;
    r = cf[0] + cf[1]*y + cf[2]*y*y - w[l];
    if(fabs(r) > tol)
      return 0;
  }
  return 1;
}

#define FD_CT_MINRUN 6
static void ct_init(struct fd_ctab *c, size_t ne)
{
  memset(c, 0, sizeof(*c));
  c->ne = ne;
  c->hi = (int *)calloc(ne + 1, sizeof(int));
  c->nm = (unsigned short *)calloc(ne + 1, sizeof(unsigned short));
  c->woff = (size_t *)calloc(ne + 1, sizeof(size_t));
  c->moff = (size_t *)calloc(ne + 1, sizeof(size_t));
  c->nwcap = 16*ne + 16;
  c->W = (double *)malloc(sizeof(double)*c->nwcap);
  c->nruncap = 2*ne + 16;
  c->run = (struct fd_crun *)malloc(sizeof(struct fd_crun)*c->nruncap);
}
static void ct_free(struct fd_ctab *c)
{
  free(c->hi);
  free(c->nm);
  free(c->woff);
  free(c->moff);
  free(c->W);
  free(c->run);
  memset(c, 0, sizeof(*c));
}
static void ct_push_w(struct fd_ctab *c, double v)
{
  if(c->nw >= c->nwcap){
    c->nwcap *= 2;
    c->W = (double *)realloc(c->W, sizeof(double)*c->nwcap);
  }
  c->W[c->nw] = v;
  c->nw++;
}
static void ct_push_run(struct fd_ctab *c, int a, int b, const double *cf)
{
  if(c->nrun >= c->nruncap){
    c->nruncap *= 2;
    c->run = (struct fd_crun *)realloc(c->run, sizeof(struct fd_crun)*c->nruncap);
  }
  c->run[c->nrun].a = a;
  c->run[c->nrun].b = b;
  c->run[c->nrun].c[0] = cf[0];
  c->run[c->nrun].c[1] = cf[1];
  c->run[c->nrun].c[2] = cf[2];
  c->nrun++;
}
/* the longest run from index a (at most to hi) on one quadratic: doubling, then bisection;
   returns its last index, or a - 1 if not even FD_CT_MINRUN lags fit */
static int ct_longest(const double *w, int a, int hi, double tol, double *cf)
{
  int ok, bad, mid, len;
  double ct[3];
  if(a + FD_CT_MINRUN - 1 > hi || !ct_fit(w, a, a + FD_CT_MINRUN - 1, tol, cf))
    return a - 1;
  ok = a + FD_CT_MINRUN - 1;
  len = FD_CT_MINRUN;
  bad = -1;
  while(bad < 0){
    len *= 2;
    mid = a + len - 1;
    if(mid > hi)
      mid = hi;
    if(mid == ok)
      break;
    if(ct_fit(w, a, mid, tol, ct)){
      ok = mid;
      memcpy(cf, ct, sizeof(ct));
      if(mid == hi)
        break;
    }else
      bad = mid;
  }
  while(bad > ok + 1){
    mid = (ok + bad)/2;
    if(ct_fit(w, a, mid, tol, ct)){
      ok = mid;
      memcpy(cf, ct, sizeof(ct));
    }else
      bad = mid;
  }
  return ok;
}
/* entry e from the full window w[0..L-1] (index l is lag l + 1); tol: absolute tolerance */
static void ct_entry(struct fd_ctab *c, size_t e, const double *w, int L, double tol)
{
  int hi, pos, b, n;
  double cf[3];
  hi = L - 1;
  while(hi >= 0 && fabs(w[hi]) <= tol)
    hi--;
  c->hi[e] = hi + 1;
  c->woff[e] = c->nw;
  c->moff[e] = c->nrun;
  c->nm[e] = 0;
  if(hi < 0)
    return;
  c->nlag_win += (size_t)(hi + 1);
  for(pos = 0; pos < g_cexp && pos <= hi; pos++)  /* the first lags stay explicit */
    ct_push_w(c, w[pos]);
  if(pos < 1)
    pos = 1;
  n = 0;
  while(pos <= hi){
    b = (n < 65535) ? ct_longest(w, pos, hi, tol, cf) : pos - 1;
    if(b >= pos){
      ct_push_run(c, pos + 1, b + 1, cf);
      n++;
      pos = b + 1;
    }else{
      ct_push_w(c, w[pos]);
      pos++;
    }
  }
  c->nm[e] = (unsigned short)n;
}

double fd_ctab_value(const struct fd_ctab *c, size_t e, int l)
{
  int k, pos, u;
  const struct fd_crun *r;
  double y;
  if(l < 1 || l > c->hi[e])
    return 0.0;
  r = c->run + c->moff[e];
  u = 0;
  pos = 1;
  for(k = 0; k < c->nm[e]; k++){
    if(l < r[k].a)
      return c->W[c->woff[e] + u + (l - pos)];
    u += r[k].a - pos;
    if(l <= r[k].b){
      y = l/FD_CT_LNORM;
      return r[k].c[0] + r[k].c[1]*y + r[k].c[2]*y*y;
    }
    pos = r[k].b + 1;
  }
  return c->W[c->woff[e] + u + (l - pos)];
}

/* compress a dense pair built with full windows (Wd_s, Wd_n: L x np, lag-major) */
static void ct_compress_pair(struct fd_pair *p, const double *xr, const double *yr, double mu, double nu)
{
  size_t np, e;
  int l, c, L;
  double *w, rr, env;
  const double *Wsrc;
  struct fd_ctab *ct;
  np = (size_t)p->src->n*p->rcv->n;
  L = p->L;
  w = (double *)malloc(sizeof(double)*(L + 1));
  for(c = 0; c < (p->Wd_n ? 2 : 1); c++){
    ct = c ? &p->cn : &p->cs;
    Wsrc = c ? p->Wd_n : p->Wd_s;
    ct_init(ct, np);
    for(e = 0; e < np; e++){
      rr = fmax(hypot(xr[e], yr[e]), p->src->a);
      env = mu*p->src->dx*p->src->dx/(2.0*M_PI*(p->mode == FD_MODE_INPLANE ? (1.0 - nu) : 1.0)*rr*rr)*p->dt;
      for(l = 0; l < L; l++)
        w[l] = Wsrc[(size_t)l*np + e];
      ct_entry(ct, e, w, L, g_ctol*env);
    }
  }
  free(w);
  free(p->Wd_s);
  free(p->Wd_n);
  p->Wd_s = NULL;
  p->Wd_n = NULL;
  p->compressed = 1;
}

double fd_pair_dense_bytes(const struct fd_pair *p, double *win)
{
  double b;
  size_t np;
  int c;
  const struct fd_ctab *ct;
  np = (size_t)p->src->n*p->rcv->n;
  if(win)
    *win = 8.0*(double)p->L*np*(p->has_normal ? 2 : 1);
  if(!p->compressed)
    return 8.0*(double)p->L*np*(p->has_normal ? 2 : 1);
  b = 0.0;
  for(c = 0; c < (p->has_normal ? 2 : 1); c++){
    ct = c ? &p->cn : &p->cs;
    b += 8.0*ct->nw + (double)sizeof(struct fd_crun)*ct->nrun + (double)ct->ne*(sizeof(int) + sizeof(unsigned short) + 2.0*sizeof(size_t));
  }
  return b;
}

void fd_ctab_prefix(const double *ring, int Lbuf, long pos0, int lag0, int lav, int L, int n, double *P, const double **vl)
{
  int l, k, j;
  double y, yk;
  double *Pk, *Pm;
  const double *v;
  if(lav > L)
    lav = L;
  for(l = 0; l <= L; l++)
    vl[l] = (l >= lag0 && l <= lav) ? ring + (size_t)((pos0 - l) % Lbuf)*n : NULL;
  for(k = 0; k < 3; k++){
    Pk = P + (size_t)k*(L + 1)*n;
    for(j = 0; j < n; j++)
      Pk[j] = 0.0;
    for(l = 1; l <= L; l++){
      Pm = Pk + (size_t)l*n;
      v = vl[l];
      if(v){
        y = l/FD_CT_LNORM;
        yk = (k == 0) ? 1.0 : (k == 1 ? y : y*y);
        for(j = 0; j < n; j++)
          Pm[j] = Pm[j - n] + yk*v[j];
      }else
        memcpy(Pm, Pm - n, sizeof(double)*n);
    }
  }
}

void fd_ctab_hist(const struct fd_ctab *c, int nr, int ns, const double *P, const double *const *vl, int lag0, int lav, int L, double *out)
{
  int i, j, k, l, pos, hi, u, ma, mb;
  size_t e, L1;
  double a, b;
  const double *w, *P0, *P1, *P2;
  const struct fd_crun *r;
  L1 = (size_t)L + 1;
  P0 = P;
  P1 = P + L1*ns;
  P2 = P + 2*L1*ns;
  if(lav > L)
    lav = L;
  for(i = 0; i < nr; i++){
    a = 0.0;
    for(j = 0; j < ns; j++){
      e = (size_t)i*ns + j;
      hi = c->hi[e];
      if(hi < 1)
        continue;
      w = c->W + c->woff[e];
      r = c->run + c->moff[e];
      u = 0;
      pos = 1;
      for(k = 0; k < c->nm[e]; k++){
        ma = r[k].a;
        mb = r[k].b;
        for(l = pos; l < ma; l++){
          if(l >= lag0 && l <= lav)
            a += w[u]*vl[l][j];
          u++;
        }
        if(ma <= lav && mb >= lag0){
          b = r[k].c[0]*(P0[(size_t)mb*ns + j] - P0[(size_t)(ma - 1)*ns + j]);
          b += r[k].c[1]*(P1[(size_t)mb*ns + j] - P1[(size_t)(ma - 1)*ns + j]);
          b += r[k].c[2]*(P2[(size_t)mb*ns + j] - P2[(size_t)(ma - 1)*ns + j]);
          a += b;
        }
        pos = mb + 1;
      }
      for(l = pos; l <= hi; l++){
        if(l >= lag0 && l <= lav)
          a += w[u]*vl[l][j];
        u++;
      }
    }
    out[i] += a;
  }
}
