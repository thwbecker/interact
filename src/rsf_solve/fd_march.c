/* fully dynamic marcher core, see fd_march.h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd_march.h"

/* ---------------- system ---------------- */
struct fd_sys *fd_sys_create(struct fd_seglist *segs, struct fd_pair **pairs, int nsrc, const int *src_of, const int *src_rev,
                             int mode, double dt, double alpha, double beta, double rho, int nsign)
{
  struct fd_sys *s = (struct fd_sys *)calloc(1, sizeof(*s));
  int k, l, i, j, ns, nr, jj;
  double *Ks, *Kn;
  s->segs = segs; s->pairs = pairs; s->nseg = segs->nseg; s->mode = mode; s->dt = dt;
  s->nsrc = nsrc > 0 ? nsrc : segs->nseg;
  s->src_of = (int *)malloc(sizeof(int) * (size_t)s->nsrc); s->src_rev = (int *)malloc(sizeof(int) * (size_t)s->nsrc);
  for(k = 0; k < s->nsrc; k++){ s->src_of[k] = src_of ? src_of[k] : k; s->src_rev[k] = src_rev ? src_rev[k] : 0; }
  s->alpha = alpha; s->beta = beta; s->rho = rho; s->nsign = nsign;
  s->eta = rho * beta * beta / (2.0 * beta);
  s->has_normal = (mode == FD_MODE_INPLANE);
  s->off = (int *)malloc(sizeof(int) * (size_t)(s->nseg + 1));
  s->off[0] = 0;
  for(k = 0; k < s->nseg; k++) s->off[k + 1] = s->off[k] + segs->seg[k].n;
  s->N = s->off[s->nseg];
  s->M = pairs[0]->M; s->L = 0;
  for(k = 0; k < s->nsrc * s->nseg; k++) if(pairs[k]->L > s->L) s->L = pairs[k]->L;
  s->Lbuf = s->L + 1;
  s->hist = (double _Complex **)malloc(sizeof(double _Complex *) * (size_t)s->nsrc);
  s->hk = (int *)calloc((size_t)s->nsrc, sizeof(int));
  for(k = 0; k < s->nsrc; k++) s->hist[k] = (double _Complex *)calloc((size_t)s->Lbuf * s->M, sizeof(double _Complex));
  s->pf = fd_fft_plan_create(s->M, FD_FFT_FORWARD); s->pb = fd_fft_plan_create(s->M, FD_FFT_BACKWARD);
  s->w1 = (double _Complex *)malloc(sizeof(double _Complex) * s->M);
  s->w2 = (double _Complex *)malloc(sizeof(double _Complex) * s->M);
  s->w3 = (double _Complex *)malloc(sizeof(double _Complex) * s->M);
  /* dense static operators in segment order */
  s->Ks = (double *)calloc((size_t)s->N * s->N, sizeof(double));
  s->Kn = (double *)calloc((size_t)s->N * s->N, sizeof(double));
  for(k = 0; k < s->nsrc; k++)
    for(l = 0; l < s->nseg; l++){
      struct fd_pair *p = pairs[k * s->nseg + l];     /* source k, receiver l */
      int r = s->src_of[k], rev = s->src_rev[k];
      ns = p->src->n; nr = p->rcv->n;
      Ks = (double *)malloc(sizeof(double) * ns * nr); Kn = (double *)malloc(sizeof(double) * ns * nr);
      fd_pair_static_dense(p, Ks, Kn);
      for(i = 0; i < nr; i++)
        for(j = 0; j < ns; j++){
          jj = rev ? ns - 1 - j : j;                 /* column of the real element feeding source element j */
          s->Ks[(size_t)(s->off[l] + i) * s->N + s->off[r] + jj] += Ks[i * ns + j];
          s->Kn[(size_t)(s->off[l] + i) * s->N + s->off[r] + jj] += Kn[i * ns + j];
        }
      free(Ks); free(Kn);
    }
  s->nops = 0;
  return s;
}

void fd_sys_free(struct fd_sys *s)
{
  int k;
  if(!s) return;
  for(k = 0; k < s->nsrc; k++) free(s->hist[k]);
  free(s->hist); free(s->hk); free(s->off); free(s->src_of); free(s->src_rev); free(s->w1); free(s->w2); free(s->w3); free(s->Ks); free(s->Kn);
  for(k = 0; k < s->nops; k++){ free(s->ops[k].Ks); free(s->ops[k].Kn); free(s->ops[k].ab); }
  fd_fft_plan_destroy(s->pf); fd_fft_plan_destroy(s->pb);
  free(s);
}

void fd_sys_reset_history(struct fd_sys *s)
{
  int k;
  for(k = 0; k < s->nsrc; k++){ memset(s->hist[k], 0, sizeof(double _Complex) * (size_t)s->Lbuf * s->M); s->hk[k] = 0; }
}

void fd_sys_push(struct fd_sys *s, const double *V)
{
  int k, i, r, n;
  for(k = 0; k < s->nsrc; k++){
    r = s->src_of[k]; n = s->segs->seg[r].n;
    memset(s->w1, 0, sizeof(double _Complex) * s->M);
    for(i = 0; i < n; i++) s->w1[i] = V[s->off[r] + (s->src_rev[k] ? n - 1 - i : i)];
    fd_fft_exec(s->pf, s->w1, s->hist[k] + (size_t)(s->hk[k] % s->Lbuf) * s->M);
    s->hk[k]++;
  }
}

/* sum_{l = lag0}^{L} W[l] V^{(next - l)} for all receivers (shear, normal);
   the ring buffer is walked in at most two contiguous reversed segments */
void fd_sys_history(struct fd_sys *s, int lag0, double *ts, double *tn)
{
  int a, b, l, lmax, pos, nseg, q, i, m, hk;
  memset(ts, 0, sizeof(double) * s->N); if(tn) memset(tn, 0, sizeof(double) * s->N);
  for(a = 0; a < s->nsrc; a++){          /* source */
    hk = s->hk[a];
    for(b = 0; b < s->nseg; b++){        /* receiver */
      struct fd_pair *p = s->pairs[a * s->nseg + b];
      lmax = p->L < hk ? p->L : hk;
      if(lmax < lag0) continue;
      memset(s->w1, 0, sizeof(double _Complex) * s->M);
      if(p->has_normal) memset(s->w2, 0, sizeof(double _Complex) * s->M);
      l = lag0;
      while(l <= lmax){
        pos = (hk - l) % s->Lbuf;
        nseg = lmax - l + 1 < pos + 1 ? lmax - l + 1 : pos + 1;
        for(q = 0; q < nseg; q++){
          const double _Complex *Vh = s->hist[a] + (size_t)(pos - q) * s->M;
          const double _Complex *Ws = p->Ws_hat + (size_t)(l - 1 + q) * s->M;
          for(m = 0; m < s->M; m++) s->w1[m] += Ws[m] * Vh[m];
          if(p->has_normal){
            const double _Complex *Wn = p->Wn_hat + (size_t)(l - 1 + q) * s->M;
            for(m = 0; m < s->M; m++) s->w2[m] += Wn[m] * Vh[m];
          }
        }
        l += nseg;
      }
      fd_fft_exec(s->pb, s->w1, s->w3);
      for(i = 0; i < p->rcv->n; i++) ts[s->off[b] + i] += creal(s->w3[i]) / s->M;
      if(p->has_normal && tn){
        fd_fft_exec(s->pb, s->w2, s->w3);
        for(i = 0; i < p->rcv->n; i++) tn[s->off[b] + i] += creal(s->w3[i]) / s->M;
      }
    }
  }
}

void fd_sys_static(struct fd_sys *s, const double *delta, double *ts, double *tn)
{
  int i, j;
  for(i = 0; i < s->N; i++){
    double a = 0.0, b = 0.0;
    const double *rs = s->Ks + (size_t)i * s->N, *rn = s->Kn + (size_t)i * s->N;
    for(j = 0; j < s->N; j++){ a += rs[j] * delta[j]; b += rn[j] * delta[j]; }
    ts[i] = a; if(tn) tn[i] = b;
  }
}

/* K_k = sum_{l<=k} W[l] + k Ks dt, dense in segment order, plus its band */
struct fd_ops *fd_sys_ops(struct fd_sys *s, int k)
{
  int q, a, b, l, i, j, ns, nr, bw;
  double *Ws, *Wn, thr, amax = 0.0;
  struct fd_ops *op;
  for(q = 0; q < s->nops; q++) if(s->ops[q].k == k) return s->ops + q;
  if(s->nops >= 8) return NULL;
  op = s->ops + s->nops;
  op->k = k;
  op->Ks = (double *)malloc(sizeof(double) * (size_t)s->N * s->N);
  op->Kn = (double *)malloc(sizeof(double) * (size_t)s->N * s->N);
  for(i = 0; i < (size_t)s->N * s->N; i++){ op->Ks[i] = s->Ks[i] * s->dt * k; op->Kn[i] = s->Kn[i] * s->dt * k; }
  for(a = 0; a < s->nsrc; a++)
    for(b = 0; b < s->nseg; b++){
      struct fd_pair *p = s->pairs[a * s->nseg + b];
      int r = s->src_of[a], rev = s->src_rev[a], jj;
      ns = p->src->n; nr = p->rcv->n;
      Ws = (double *)malloc(sizeof(double) * ns * nr); Wn = (double *)malloc(sizeof(double) * ns * nr);
      for(l = 1; l <= (k < p->L ? k : p->L); l++){
        fd_pair_W_dense(p, l, Ws, Wn);
        for(i = 0; i < nr; i++)
          for(j = 0; j < ns; j++){
            jj = rev ? ns - 1 - j : j;
            op->Ks[(size_t)(s->off[b] + i) * s->N + s->off[r] + jj] += Ws[i * ns + j];
            op->Kn[(size_t)(s->off[b] + i) * s->N + s->off[r] + jj] += Wn[i * ns + j];
          }
      }
      free(Ws); free(Wn);
    }
  /* prune FFT round-trip noise so that the band is the physical one */
  for(i = 0; i < (size_t)s->N * s->N; i++) if(fabs(op->Ks[i]) > amax) amax = fabs(op->Ks[i]);
  thr = 1e-13 * amax;
  for(i = 0; i < (size_t)s->N * s->N; i++){ if(fabs(op->Ks[i]) <= thr) op->Ks[i] = 0.0; if(fabs(op->Kn[i]) <= thr) op->Kn[i] = 0.0; }
  bw = fd_band_width(s->N, op->Ks, 0.0);
  op->bw = bw;
  op->ab = (double *)malloc(sizeof(double) * (size_t)(3 * bw + 1) * s->N);
  fd_band_from_dense(s->N, bw, op->Ks, op->ab);
  s->nops++;
  return op;
}

/* ---------------- slip weakening ---------------- */
static void dense_mult(int N, const double *A, const double *x, double *y)
{
  int i, j;
  for(i = 0; i < N; i++){ double a = 0.0; const double *r = A + (size_t)i * N; for(j = 0; j < N; j++) a += r[j] * x[j]; y[i] = a; }
}

/* exact coupled solve with an active-set loop (stepper.SlipWeakening._solve_coupled) */
int fd_slipweak_solve(struct fd_slipweak *f, struct fd_sys *s, struct fd_ops *op, const double *H, const double *Hn,
                      const double *sigma0, int normal_coupling, double t, double dt, const double *delta,
                      const double *Vguess, double *V, double *sigma)
{
  int N = s->N, bw = op->bw, ldab = 3 * bw + 1, i, j, it, changed, max_iter = 12;
  double *sig = (double *)malloc(sizeof(double) * N), *mu_slip = (double *)malloc(sizeof(double) * N);
  double *slope0 = (double *)malloc(sizeof(double) * N), *str0 = (double *)malloc(sizeof(double) * N);
  double *ab = (double *)malloc(sizeof(double) * (size_t)ldab * N), *rhs = (double *)calloc(N, sizeof(double));
  double *Vn = (double *)malloc(sizeof(double) * N), *tau = (double *)malloc(sizeof(double) * N), *tmp = (double *)malloc(sizeof(double) * N);
  int *active = (int *)calloc(N, sizeof(int)), *sat = (int *)calloc(N, sizeof(int));
  double muf, fr, strength0, slope;
  /* normal stress with the guessed slip rate */
  if(normal_coupling){
    dense_mult(N, op->Kn, Vguess, tmp);
    for(i = 0; i < N; i++) sig[i] = sigma0[i] + s->nsign * (Hn[i] + tmp[i]);
  }else for(i = 0; i < N; i++) sig[i] = sigma0[i];
  for(i = 0; i < N; i++){
    double r = delta[i] / f->d0; if(r > 1.0) r = 1.0;
    if(f->dugdale){ mu_slip[i] = delta[i] < f->d0 ? f->mu_s[i] : f->mu_d[i]; slope0[i] = 0.0; }
    else{ mu_slip[i] = f->mu_s[i] - (f->mu_s[i] - f->mu_d[i]) * r; slope0[i] = delta[i] < f->d0 ? sig[i] * (f->mu_s[i] - f->mu_d[i]) * dt / f->d0 : 0.0; }
    if(f->force_time){
      fr = (t - f->force_time[i]) / (f->force_t0 > 1e-12 ? f->force_t0 : 1e-12);
      if(fr < 0.0) fr = 0.0;
      if(fr > 1.0) fr = 1.0;
      muf = f->mu_s[i] - (f->mu_s[i] - f->mu_d[i]) * fr;
      if(muf < mu_slip[i]){ mu_slip[i] = muf; slope0[i] = 0.0; }
    }
    str0[i] = sig[i] * mu_slip[i];
  }
  /* initial active set */
  dense_mult(N, op->Ks, Vguess, tmp);
  for(i = 0; i < N; i++){
    active[i] = (H[i] + tmp[i] > str0[i]) || (Vguess[i] > 0.0);
    if(f->locked && f->locked[i]) active[i] = 0;
  }
  for(i = 0; i < N; i++) V[i] = 0.0;
  for(it = 0; it < max_iter; it++){
    memcpy(ab, op->ab, sizeof(double) * (size_t)ldab * N);
    for(i = 0; i < N; i++){
      strength0 = sat[i] ? sig[i] * f->mu_d[i] : str0[i];
      slope = sat[i] ? 0.0 : slope0[i];
      if(active[i]){
        ab[(2 * bw) + (size_t)i * ldab] += slope;           /* diagonal (i,i) */
        rhs[i] = strength0 - H[i];
      }else{
        for(j = (i - bw > 0 ? i - bw : 0); j <= (i + bw < N - 1 ? i + bw : N - 1); j++)
          ab[(2 * bw + i - j) + (size_t)j * ldab] = (j == i) ? 1.0 : 0.0;
        rhs[i] = 0.0;
      }
    }
    if(fd_band_solve(N, bw, ab, rhs, Vn) != 0){ fprintf(stderr, "fd_slipweak_solve: band solve failed\n"); break; }
    for(i = 0; i < N; i++) if(!active[i]) Vn[i] = 0.0;
    dense_mult(N, op->Ks, Vn, tmp);
    changed = 0;
    for(i = 0; i < N; i++){
      int neg, over, cross, uncross;
      strength0 = sat[i] ? sig[i] * f->mu_d[i] : str0[i];
      tau[i] = H[i] + tmp[i];
      neg = active[i] && Vn[i] < 0.0;
      over = (!active[i]) && tau[i] > strength0 && !(f->locked && f->locked[i]);
      cross = active[i] && !sat[i] && delta[i] < f->d0 && delta[i] + Vn[i] * dt > f->d0;
      uncross = active[i] && sat[i] && delta[i] < f->d0 && delta[i] + Vn[i] * dt < f->d0;
      if(neg || over || cross || uncross) changed = 1;
      V[i] = Vn[i];
      if(neg) active[i] = 0;
      if(over) active[i] = 1;
      if(cross) sat[i] = 1;
      if(uncross) sat[i] = 0;
    }
    if(!changed) break;
  }
  for(i = 0; i < N; i++) if(V[i] < 0.0) V[i] = 0.0;
  if(normal_coupling){
    dense_mult(N, op->Kn, V, tmp);
    for(i = 0; i < N; i++) sigma[i] = sigma0[i] + s->nsign * (Hn[i] + tmp[i]);
  }else for(i = 0; i < N; i++) sigma[i] = sigma0[i];
  free(sig); free(mu_slip); free(slope0); free(str0); free(ab); free(rhs); free(Vn); free(tau); free(tmp); free(active); free(sat);
  return it;
}

/* ---------------- state and step ---------------- */
struct fd_state *fd_state_create(int N)
{
  struct fd_state *st = (struct fd_state *)calloc(1, sizeof(*st));
  st->N = N; st->k = 1;
  st->delta = calloc(N, sizeof(double)); st->V = calloc(N, sizeof(double)); st->tau = calloc(N, sizeof(double));
  st->sigma = calloc(N, sizeof(double)); st->tau0 = calloc(N, sizeof(double)); st->taudot = calloc(N, sizeof(double));
  st->sigma0 = calloc(N, sizeof(double)); st->Hn_dot = calloc(N, sizeof(double)); st->psi = calloc(N, sizeof(double));
  st->H = calloc(N, sizeof(double)); st->Hn = calloc(N, sizeof(double)); st->hs = calloc(N, sizeof(double));
  st->hn = calloc(N, sizeof(double)); st->work = calloc(N, sizeof(double));
  return st;
}
void fd_state_free(struct fd_state *st)
{
  if(!st) return;
  free(st->delta); free(st->V); free(st->tau); free(st->sigma); free(st->tau0); free(st->taudot); free(st->sigma0);
  free(st->Hn_dot); free(st->psi); free(st->H); free(st->Hn); free(st->hs); free(st->hn); free(st->work); free(st);
}

int fd_step_slipweak(struct fd_sys *s, struct fd_state *st, struct fd_slipweak *f, int k)
{
  int i, q, N = s->N;
  double dt = s->dt * k, t1 = st->t + dt;
  struct fd_ops *op = fd_sys_ops(s, k);
  if(!op) return -1;
  /* known parts: static + loading + history beyond the current step */
  fd_sys_static(s, st->delta, st->H, st->Hn);
  fd_sys_history(s, k + 1, st->hs, st->hn);
  for(i = 0; i < N; i++){
    st->H[i] += st->tau0[i] + st->taudot[i] * t1 + st->hs[i];
    st->Hn[i] += st->Hn_dot[i] * t1 + st->hn[i];
  }
  fd_slipweak_solve(f, s, op, st->H, st->Hn, st->sigma0, st->normal_coupling, t1, dt, st->delta, st->V, st->work, st->sigma);
  memcpy(st->V, st->work, sizeof(double) * N);
  dense_mult(N, op->Ks, st->V, st->work);
  for(i = 0; i < N; i++) st->tau[i] = st->H[i] + st->work[i];
  for(q = 0; q < k; q++) fd_sys_push(s, st->V);
  for(i = 0; i < N; i++) st->delta[i] += st->V[i] * dt;
  st->t = t1; st->n++; st->k = k;
  return 0;
}

/* ---------------- rate and state (aging law) ---------------- */
double fd_rsf_tau(const struct fd_rsf *f, int i, double V, double psi, double sigma)
{
  return f->a[i] * sigma * asinh(V / (2.0 * f->v0) * exp(psi / f->a[i]));
}
double fd_rsf_dtau(const struct fd_rsf *f, int i, double V, double psi, double sigma)
{
  double e = exp(psi / f->a[i]), z = V / (2.0 * f->v0) * e;
  return f->a[i] * sigma * e / (2.0 * f->v0) / sqrt(1.0 + z * z);
}
/* exact aging-law solution over dt at the mean slip rate (Lapusta et al. 2000):
   theta = dc/v0 exp((psi - f0)/b); theta' = dc/V + (theta - dc/V) exp(-V dt/dc) */
double fd_rsf_psi_advance(const struct fd_rsf *f, int i, double psi, double Va, double Vb, double dt)
{
  double V = 0.5 * (fabs(Va) + fabs(Vb)), th, thn;
  if(V < 1e-30) V = 1e-30;
  th = f->dc[i] / f->v0 * exp((psi - f->f0) / f->b[i]);
  thn = f->dc[i] / V + (th - f->dc[i] / V) * exp(-V * dt / f->dc[i]);
  return f->f0 + f->b[i] * log(f->v0 * thn / f->dc[i]);
}

/* elementwise bracketed Newton in u = ln V on Hloc + k11 V = tau_f(V): g is
   strictly decreasing in u since k11 < 0 */
static double rsf_local(const struct fd_rsf *f, int i, double Hloc, double k11, double sigma, double psi, double Vg)
{
  double u, lo = -700.0, hi = 50.0, g, dg, un, V;
  int it;
  if(Hloc <= 0.0) return 0.0;
  u = log(Vg > 1e-300 ? Vg : 1e-300);
  if(u < -700.0) u = -700.0;
  if(u > 50.0) u = 50.0;
  for(it = 0; it < 200; it++){
    V = exp(u);
    g = Hloc + k11 * V - fd_rsf_tau(f, i, V, psi, sigma);
    if(g > 0.0){ if(u > lo) lo = u; } else { if(u < hi) hi = u; }
    dg = (k11 - fd_rsf_dtau(f, i, V, psi, sigma)) * V;
    un = (dg != 0.0 && isfinite(dg)) ? u - g / dg : 0.5 * (lo + hi);
    if(!isfinite(un) || un <= lo || un >= hi) un = 0.5 * (lo + hi);
    if(fabs(un - u) < f->tol){ u = un; break; }
    u = un;
  }
  return exp(u);
}

int fd_rsf_solve(const struct fd_rsf *f, struct fd_sys *s, struct fd_ops *op, const double *H, const double *sigma,
                 const double *psi, const double *Vg, double *V)
{
  int N = s->N, bw = op->bw, ldab = 3 * bw + 1, i, j, it, converged = 0;
  double *u = (double *)malloc(sizeof(double) * N), *g = (double *)malloc(sizeof(double) * N), *du = (double *)malloc(sizeof(double) * N);
  double *ab = (double *)malloc(sizeof(double) * (size_t)ldab * N), *tmp = (double *)malloc(sizeof(double) * N), *k11 = (double *)malloc(sizeof(double) * N);
  double dmax, Vi;
  for(i = 0; i < N; i++){ u[i] = log(Vg[i] > 1e-30 ? Vg[i] : 1e-30); if(u[i] < -70.0) u[i] = -70.0; if(u[i] > 10.0) u[i] = 10.0; }
  for(it = 0; it < f->max_iter; it++){
    for(i = 0; i < N; i++) V[i] = exp(u[i]);
    dense_mult(N, op->Ks, V, tmp);
    for(i = 0; i < N; i++) g[i] = H[i] + tmp[i] - fd_rsf_tau(f, i, V[i], psi[i], sigma[i]);
    /* Jacobian in u: (K - diag(dtau/dV)) diag(V), band storage: column j scaled by V_j */
    for(j = 0; j < N; j++)
      for(i = (j - bw > 0 ? j - bw : 0); i <= (j + bw < N - 1 ? j + bw : N - 1); i++)
        ab[(2 * bw + i - j) + (size_t)j * ldab] = op->ab[(2 * bw + i - j) + (size_t)j * ldab] * V[j];
    for(i = 0; i < N; i++) ab[(2 * bw) + (size_t)i * ldab] -= fd_rsf_dtau(f, i, V[i], psi[i], sigma[i]) * V[i];
    for(i = 0; i < N; i++){ if(!isfinite(g[i]) || !isfinite(ab[(2 * bw) + (size_t)i * ldab])) break; }
    if(i < N) break;
    if(fd_band_solve(N, bw, ab, g, du) != 0) break;
    dmax = 0.0;
    for(i = 0; i < N; i++){
      double d = -du[i];
      if(!isfinite(d)){ dmax = INFINITY; break; }
      if(d > 2.0) d = 2.0;
      if(d < -2.0) d = -2.0;
      u[i] += d; if(u[i] < -70.0) u[i] = -70.0; if(u[i] > 10.0) u[i] = 10.0;
      if(fabs(d) > dmax) dmax = fabs(d);
    }
    if(!isfinite(dmax)) break;
    if(dmax < f->tol){ converged = 1; break; }
  }
  if(converged){
    for(i = 0; i < N; i++) V[i] = exp(u[i]);
    free(u); free(g); free(du); free(ab); free(tmp); free(k11);
    return it;
  }
  /* fallback: fixed point of the elementwise solve with lagged neighbours */
  for(i = 0; i < N; i++){ k11[i] = op->Ks[(size_t)i * N + i]; V[i] = exp(u[i]); }
  for(it = 0; it < 60; it++){
    dense_mult(N, op->Ks, V, tmp);
    dmax = 0.0;
    for(i = 0; i < N; i++){
      Vi = rsf_local(f, i, H[i] + tmp[i] - k11[i] * V[i], k11[i], sigma[i], psi[i], V[i]);
      if(fabs(Vi - V[i]) / (V[i] > 1e-30 ? V[i] : 1e-30) > dmax) dmax = fabs(Vi - V[i]) / (V[i] > 1e-30 ? V[i] : 1e-30);
      g[i] = Vi;
    }
    memcpy(V, g, sizeof(double) * N);
    if(dmax < f->tol) break;
  }
  free(u); free(g); free(du); free(ab); free(tmp); free(k11);
  return -(it + 1);
}

/* one step of k dt_min: predictor-corrector in psi (exact aging update at the
   mean slip rate), coupled solve for V twice (stepper.RateState.solve) */
int fd_step_rsf(struct fd_sys *s, struct fd_state *st, const struct fd_rsf *f, int k)
{
  int i, q, N = s->N;
  double dt = s->dt * k, t1 = st->t + dt, *psi_s = st->work, *Vs, *sig;
  struct fd_ops *op = fd_sys_ops(s, k);
  if(!op) return -1;
  Vs = (double *)malloc(sizeof(double) * N); sig = (double *)malloc(sizeof(double) * N);
  fd_sys_static(s, st->delta, st->H, st->Hn);
  fd_sys_history(s, k + 1, st->hs, st->hn);
  for(i = 0; i < N; i++){
    st->H[i] += st->tau0[i] + st->taudot[i] * t1 + st->hs[i];
    st->Hn[i] += st->Hn_dot[i] * t1 + st->hn[i];
  }
  /* predictor */
  if(st->normal_coupling){ dense_mult(N, op->Kn, st->V, sig); for(i = 0; i < N; i++) sig[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]); }
  else for(i = 0; i < N; i++) sig[i] = st->sigma0[i];
  for(i = 0; i < N; i++) psi_s[i] = fd_rsf_psi_advance(f, i, st->psi[i], st->V[i], st->V[i], dt);
  fd_rsf_solve(f, s, op, st->H, sig, psi_s, st->V, Vs);
  /* corrector */
  if(st->normal_coupling){ dense_mult(N, op->Kn, Vs, sig); for(i = 0; i < N; i++) sig[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]); }
  for(i = 0; i < N; i++) psi_s[i] = fd_rsf_psi_advance(f, i, st->psi[i], st->V[i], Vs[i], dt);
  fd_rsf_solve(f, s, op, st->H, sig, psi_s, Vs, st->V);
  memcpy(st->psi, psi_s, sizeof(double) * N);
  dense_mult(N, op->Ks, st->V, Vs);
  for(i = 0; i < N; i++) st->tau[i] = st->H[i] + Vs[i];
  if(st->normal_coupling){ dense_mult(N, op->Kn, st->V, sig); for(i = 0; i < N; i++) st->sigma[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]); }
  else for(i = 0; i < N; i++) st->sigma[i] = st->sigma0[i];
  for(q = 0; q < k; q++) fd_sys_push(s, st->V);
  for(i = 0; i < N; i++) st->delta[i] += st->V[i] * dt;
  st->t = t1; st->n++; st->k = k;
  free(Vs); free(sig);
  return 0;
}
