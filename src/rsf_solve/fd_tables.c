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
  if(errmsg) errmsg[0] = 0;
  if(!fd_segments_parallel(src, rcv)){ if(errmsg) sprintf(errmsg, "segments not parallel with equal elements"); return NULL; }
  if(lmax <= 0){
    /* generous bound: S arrival at the farthest pair times (1 + 2/sqrt(tol)) */
    for(i = -(ns - 1); i <= nr - 1; i++){ offset_point(src, rcv, i, &xr, &yr, Xg); if(hypot(xr, yr) > dmax) dmax = hypot(xr, yr); }
    lmax = (int)ceil((dmax + src->dx) / beta / dt * (1.0 + 2.0 / sqrt(tol))) + 2;
  }
  p = (struct fd_pair *)calloc(1, sizeof(*p));
  p->src = src; p->rcv = rcv; p->mode = mode; p->dt = dt; p->et = et; p->nm = nm;
  p->has_normal = (mode == FD_MODE_INPLANE);
  p->m = (int *)malloc(sizeof(int) * nm);
  p->ks_row = (double *)malloc(sizeof(double) * nm); p->kn_row = (double *)malloc(sizeof(double) * nm);
  Ws = (double *)malloc(sizeof(double) * (size_t)nm * lmax); Wn = (double *)calloc((size_t)nm * lmax, sizeof(double));
  /* receiver normal in the source frame (antiplane traction projection) */
  nr_g[0] = -sr[1]; nr_g[1] = sr[0];
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
      if(lsw > lmax) lsw = lmax;
      lc = lsw;
      /* closed form for columns 0..lc-1: levels 0, et dt, (et+1) dt, ... */
      fd_seg_ramp_inplane(xr, yr, 0.0, src->a, alpha, beta, rho, s0);
      for(l = 0; l < lc; l++){
        fd_seg_ramp_inplane(xr, yr, (et + l) * dt, src->a, alpha, beta, rho, s1);
        d[0] = s1[0] - s0[0]; d[1] = s1[1] - s0[1]; d[2] = s1[2] - s0[2];
        fd_project(d[0], d[1], d[2], ss, sr, &ts, &tn);
        Ws[(size_t)i * lmax + l] = ts - p->ks_row[i] * dt;
        Wn[(size_t)i * lmax + l] = tn - p->kn_row[i] * dt;
        memcpy(s0, s1, sizeof(s0));
      }
      /* series beyond: column l-1 (1-based lag l >= lsw) covers [(et+l-2) dt, (et+l-1) dt] */
      for(l = lsw; l <= lmax; l++){
        fd_series_inplane(xr, yr, (l - 2 + et) * dt, dt, src->a, alpha, beta, rho, d);
        fd_project(d[0], d[1], d[2], ss, sr, &ts, &tn);
        Ws[(size_t)i * lmax + (l - 1)] = ts; Wn[(size_t)i * lmax + (l - 1)] = tn;
      }
    }else{
      fd_static_antiplane(xr, yr, src->a, mu, a2);
      p->ks_row[i] = ns_src[0] * a2[0] + ns_src[1] * a2[1]; p->kn_row[i] = 0.0;
      fd_seg_ramp_antiplane(xr, yr, 0.0, src->a, beta, rho, a2);
      for(l = 0; l < lmax; l++){
        fd_seg_ramp_antiplane(xr, yr, (et + l) * dt, src->a, beta, rho, b2);
        Ws[(size_t)i * lmax + l] = ns_src[0] * (b2[0] - a2[0]) + ns_src[1] * (b2[1] - a2[1]) - p->ks_row[i] * dt;
        a2[0] = b2[0]; a2[1] = b2[1];
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
      if(fabs(Ws[(size_t)i * lmax + l]) > tol * env || fabs(Wn[(size_t)i * lmax + l]) > tol * env) break;
    if(l + 1 > L) L = l + 1;
  }
  if(L < 1) L = 1;
  if(L >= lmax){
    if(!cap_ok){ if(errmsg) sprintf(errmsg, "window reached lmax = %d; increase lmax", lmax); free(Ws); free(Wn); fd_pair_free(p); return NULL; }
    L = lmax - 1;
  }
  p->L = L;
  p->M = M > 0 ? M : fd_fft_size_for(nm);
  p->pf = fd_fft_plan_create(p->M, FD_FFT_FORWARD); p->pb = fd_fft_plan_create(p->M, FD_FFT_BACKWARD);
  p->Ws_hat = (double _Complex *)malloc(sizeof(double _Complex) * (size_t)L * p->M);
  p->Wn_hat = p->has_normal ? (double _Complex *)malloc(sizeof(double _Complex) * (size_t)L * p->M) : NULL;
  col = (double _Complex *)malloc(sizeof(double _Complex) * p->M);
  for(l = 0; l < L; l++){
    memset(col, 0, sizeof(double _Complex) * p->M);
    for(i = 0; i < nm; i++){ idx = p->m[i] >= 0 ? p->m[i] : p->M + p->m[i]; col[idx] = Ws[(size_t)i * lmax + l]; }
    fd_fft_exec(p->pf, col, p->Ws_hat + (size_t)l * p->M);
    if(p->has_normal){
      memset(col, 0, sizeof(double _Complex) * p->M);
      for(i = 0; i < nm; i++){ idx = p->m[i] >= 0 ? p->m[i] : p->M + p->m[i]; col[idx] = Wn[(size_t)i * lmax + l]; }
      fd_fft_exec(p->pf, col, p->Wn_hat + (size_t)l * p->M);
    }
  }
  free(col); free(Ws); free(Wn);
  return p;
}

void fd_pair_free(struct fd_pair *p)
{
  if(!p) return;
  free(p->m); free(p->ks_row); free(p->kn_row); free(p->Ws_hat); free(p->Wn_hat);
  fd_fft_plan_destroy(p->pf); fd_fft_plan_destroy(p->pb); free(p);
}

static int row_index(const struct fd_pair *p, int i, int j){ return (i - j) + (p->src->n - 1); }

void fd_pair_static_dense(const struct fd_pair *p, double *Ks, double *Kn)
{
  int i, j, ns = p->src->n, nr = p->rcv->n;
  for(i = 0; i < nr; i++) for(j = 0; j < ns; j++){
    Ks[i * ns + j] = p->ks_row[row_index(p, i, j)];
    if(Kn) Kn[i * ns + j] = p->kn_row[row_index(p, i, j)];
  }
}

void fd_pair_W_dense(const struct fd_pair *p, int l1, double *Ws, double *Wn)
{
  int i, j, k, ns = p->src->n, nr = p->rcv->n, m, idx;
  double _Complex *c = (double _Complex *)malloc(sizeof(double _Complex) * p->M);
  double _Complex *cn = p->has_normal ? (double _Complex *)malloc(sizeof(double _Complex) * p->M) : NULL;
  fd_fft_exec(p->pb, p->Ws_hat + (size_t)(l1 - 1) * p->M, c);
  if(cn) fd_fft_exec(p->pb, p->Wn_hat + (size_t)(l1 - 1) * p->M, cn);
  for(k = 0; k < p->M; k++){ c[k] /= p->M; if(cn) cn[k] /= p->M; }
  for(i = 0; i < nr; i++) for(j = 0; j < ns; j++){
    m = i - j; idx = m >= 0 ? m : p->M + m;
    Ws[i * ns + j] = creal(c[idx]);
    if(Wn) Wn[i * ns + j] = cn ? creal(cn[idx]) : 0.0;
  }
  free(c); free(cn);
}
