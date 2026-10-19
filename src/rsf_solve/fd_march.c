/* fully dynamic marcher core, see fd_march.h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd_march.h"
#include <time.h>
double fd_prof[6];   /* history, static, solve, push, ops, other: seconds (FD_PROFILE) */
static double fd_now(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }


/* in-place radix-2 FFT along the first index of a [F][M] array (rows of M
   complex values), batched over the M columns so that every butterfly is a
   contiguous, vectorisable loop; sign -1 forward, +1 backward (unnormalised) */
static void fft_rows(double _Complex *a, int F, int M, int sign)
{
  int i, j, bit, len, half, h, m;
  double _Complex *tmp = (double _Complex *)malloc(sizeof(double _Complex) * (size_t)M), w, wl;
  for(i = 1, j = 0; i < F; i++){                 /* bit reversal of rows */
    for(bit = F >> 1; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if(i < j){
      memcpy(tmp, a + (size_t)i * M, sizeof(double _Complex) * M);
      memcpy(a + (size_t)i * M, a + (size_t)j * M, sizeof(double _Complex) * M);
      memcpy(a + (size_t)j * M, tmp, sizeof(double _Complex) * M);
    }
  }
  for(len = 2; len <= F; len <<= 1){
    half = len >> 1;
    wl = cexp(sign * 2.0 * M_PI * I / len);
    for(i = 0; i < F; i += len){
      w = 1.0;
      for(h = 0; h < half; h++){
        double _Complex *u = a + (size_t)(i + h) * M, *v = a + (size_t)(i + h + half) * M;
        for(m = 0; m < M; m++){
          double _Complex t = w * v[m];
          v[m] = u[m] - t; u[m] = u[m] + t;
        }
        w *= wl;
      }
    }
  }
  free(tmp);
}


/* spectrum of input window j of source k (inputs (j-1)B .. (j+1)B-1), [f][m] */
static void form_window(struct fd_sys *s, int k, long j)
{
  long idx;
  int rr, sl = (int)(j % s->nslot[k]);
  double _Complex *Xs = s->X[k] + (size_t)sl * s->M * s->F;
  for(rr = 0; rr < s->F; rr++){
    idx = (j - 1) * s->B + rr;
    if(idx >= 0) memcpy(Xs + (size_t)rr * s->M, s->hist[k] + (size_t)(idx % s->Lbuf) * s->M, sizeof(double _Complex) * s->M);
    else memset(Xs + (size_t)rr * s->M, 0, sizeof(double _Complex) * s->M);
  }
  fft_rows(Xs, s->F, s->M, -1);
  s->Xj[k][sl] = j;
}

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
  if(s->conv){
    int np = s->nsrc * s->nseg;
    for(k = 0; k < np; k++){ free(s->Hs[k]); if(s->Hn) free(s->Hn[k]); free(s->Ys[k]); if(s->Yn) free(s->Yn[k]); }
    for(k = 0; k < s->nsrc; k++){ free(s->X[k]); free(s->Xj[k]); }
    free(s->Xj);
    free(s->Hs); free(s->Hn); free(s->Ys); free(s->Yn); free(s->X); free(s->P); free(s->nslot); free(s->cblk);
    free(s->tw1); free(s->tw2); free(s->tz); free(s->tzn);
    fd_fft_plan_destroy(s->tf); fd_fft_plan_destroy(s->tb);
  }
  if(s->Kc){
    for(k = 0; k < s->nsrc * s->nseg; k++){ free(s->Acs[k]); free(s->Acn[k]); }
    free(s->Acs); free(s->Acn); free(s->Jc);
  }
  free(s->hist); free(s->hk); free(s->off); free(s->src_of); free(s->src_rev); free(s->w1); free(s->w2); free(s->w3); free(s->Ks); free(s->Kn);
  for(k = 0; k < s->nops; k++){ free(s->ops[k].Ks); free(s->ops[k].Kn); free(s->ops[k].ab); }
  fd_fft_plan_destroy(s->pf); fd_fft_plan_destroy(s->pb);
  free(s);
}

void fd_sys_reset_history(struct fd_sys *s)
{
  int k;
  s->coarse_from = -1;
  for(k = 0; k < s->nsrc; k++){ memset(s->hist[k], 0, sizeof(double _Complex) * (size_t)s->Lbuf * s->M); s->hk[k] = 0; }
  if(s->conv){
    for(k = 0; k < s->nsrc; k++){
      int q;
      memset(s->X[k], 0, sizeof(double _Complex) * (size_t)s->nslot[k] * s->M * s->F);
      for(q = 0; q < s->nslot[k]; q++) s->Xj[k][q] = -1;
    }
    for(k = 0; k < s->nsrc * s->nseg; k++) s->cblk[k] = -1;
  }
}

void fd_sys_push(struct fd_sys *s, const double *V)
{
  fd_sys_push_k(s, V, 1);
}

/* push the same slice k times (a step of k fine substeps): one spatial FFT
   per source, copied into k ring slots; windows formed at block boundaries */
void fd_sys_push_k(struct fd_sys *s, const double *V, int kk)
{
  int k, i, r, n, q;
  if(s->Kc){
    if((kk == s->Kc) && (s->hk[0] % s->Kc == 0)){ if(s->coarse_from < 0) s->coarse_from = s->hk[0]; }
    else s->coarse_from = -1;
  }
  for(k = 0; k < s->nsrc; k++){
    double _Complex *first = NULL;
    r = s->src_of[k]; n = s->segs->seg[r].n;
    for(q = 0; q < kk; q++){
      double _Complex *slot = s->hist[k] + (size_t)(s->hk[k] % s->Lbuf) * s->M;
      if(!first){
        memset(s->w1, 0, sizeof(double _Complex) * s->M);
        for(i = 0; i < n; i++) s->w1[i] = V[s->off[r] + (s->src_rev[k] ? n - 1 - i : i)];
        fd_fft_exec(s->pf, s->w1, slot);
        first = slot;
      }else
        memcpy(slot, first, sizeof(double _Complex) * s->M);
      s->hk[k]++;
      if((s->conv == 1) && (s->hk[k] % s->B == 0))
        form_window(s, k, s->hk[k] / s->B - 1);      /* complete now; the ring will not keep its inputs */
    }
  }
}

int fd_sys_set_coarse(struct fd_sys *s, int K)
{
  int np = s->nsrc * s->nseg, k, j, l, m, comp;
  if(K < 2 || s->conv == 1) return 1;
  s->Kc = K; s->coarse_from = -1;
  s->Jc = calloc((size_t)np, sizeof(int));
  s->Acs = calloc((size_t)np, sizeof(double _Complex *)); s->Acn = calloc((size_t)np, sizeof(double _Complex *));
  for(k = 0; k < np; k++){
    struct fd_pair *p = s->pairs[k];
    s->Jc[k] = (p->L - 1) / K;                     /* blocks j = 1..Jc with jK + 1 <= L */
    for(comp = 0; comp < (p->has_normal ? 2 : 1); comp++){
      const double _Complex *W = comp ? p->Wn_hat : p->Ws_hat;
      double _Complex *A = calloc((size_t)(s->Jc[k] + 1) * p->M, sizeof(double _Complex));
      for(j = 1; j <= s->Jc[k]; j++)
        for(l = j * K + 1; l <= (j + 1) * K && l <= p->L; l++)
          for(m = 0; m < p->M; m++) A[(size_t)j * p->M + m] += W[(size_t)(l - 1) * p->M + m];
      if(comp) s->Acn[k] = A; else s->Acs[k] = A;
    }
  }
  return 0;
}

double fd_sys_bytes(const struct fd_sys *s)
{
  double b = 16.0 * s->nsrc * (double)s->Lbuf * s->M;
  int a, k;
  for(k = 0; k < s->nsrc * s->nseg; k++){
    const struct fd_pair *p = s->pairs[k];
    int lt = (s->conv == 1) ? (p->L < s->D ? p->L : s->D) : p->L;
    b += 16.0 * lt * p->M * (p->has_normal ? 2 : 1);
    if(s->conv) b += 16.0 * (double)s->P[k] * s->M * s->F * (p->has_normal ? 2 : 1) + 16.0 * s->M * s->B * (p->has_normal ? 2 : 1);
  }
  if(s->conv) for(a = 0; a < s->nsrc; a++) b += 16.0 * (double)s->nslot[a] * s->M * s->F;
  return b;
}

int fd_sys_set_partitioned(struct fd_sys *s, int B, int kmax)
{
  return fd_sys_set_conv(s, 1, B, kmax);
}

/* mode 1: partitioned only (tables beyond lag D and the ring beyond D + 2B
   + kmax are released: about twice the memory of the direct sum); mode 2:
   hybrid, direct representation kept as well (about three times the memory)
   and used for steps with 2k > B, where a far block would serve fewer than
   two outputs */
int fd_sys_set_conv(struct fd_sys *s, int mode, int B, int kmax)
{
  int np = s->nsrc * s->nseg, k, p, q, a, b, Lmax, best = 0;
  double cost, cbest = 1e300;
  if(s->conv) return 0;
  if(kmax < 1) kmax = 1;
  if(B <= 0){
    /* minimise direct lags D(B) + far work 2 L / B over powers of two */
    for(b = 16; b <= 4096; b *= 2){
      int qq = 1 + (kmax + b - 1) / b;
      cost = (double)(qq * b - 1) + 2.0 * s->L / b + 4.0 * log2(2.0 * b);
      if(cost < cbest){ cbest = cost; best = b; }
    }
    B = best;
  }
  if(B < 2 || (B & (B - 1))) return 1;
  s->B = B; s->F = 2 * B; s->kmax = kmax;
  s->q0 = 1 + (kmax + B - 1) / B;           /* (q0 - 1) B >= kmax */
  s->D = s->q0 * B - 1;
  if(s->D >= s->L) return 0;               /* window shorter than the direct part: nothing to gain */
  s->conv = mode;
  s->tf = fd_fft_plan_create(s->F, FD_FFT_FORWARD); s->tb = fd_fft_plan_create(s->F, FD_FFT_BACKWARD);
  s->tw1 = malloc(sizeof(double _Complex) * s->F); s->tw2 = malloc(sizeof(double _Complex) * s->F);
  s->tz = malloc(sizeof(double _Complex) * (size_t)s->M * s->F); s->tzn = malloc(sizeof(double _Complex) * (size_t)s->M * s->F);
  s->P = calloc((size_t)np, sizeof(int)); s->cblk = malloc(sizeof(long) * (size_t)np);
  s->Hs = calloc((size_t)np, sizeof(double _Complex *)); s->Hn = s->has_normal ? calloc((size_t)np, sizeof(double _Complex *)) : NULL;
  s->Ys = calloc((size_t)np, sizeof(double _Complex *)); s->Yn = s->has_normal ? calloc((size_t)np, sizeof(double _Complex *)) : NULL;
  for(k = 0; k < np; k++){
    struct fd_pair *pr = s->pairs[k];
    s->cblk[k] = -1;
    s->P[k] = pr->L > s->D ? (pr->L - s->D + B - 1) / B : 0;
    s->Hs[k] = calloc((size_t)s->P[k] * s->M * s->F + 1, sizeof(double _Complex));
    if(s->has_normal) s->Hn[k] = calloc((size_t)s->P[k] * s->M * s->F + 1, sizeof(double _Complex));
    s->Ys[k] = calloc((size_t)s->M * B, sizeof(double _Complex));
    if(s->has_normal) s->Yn[k] = calloc((size_t)s->M * B, sizeof(double _Complex));
    for(p = 0; p < s->P[k]; p++){
      int lp = s->D + 1 + p * B;
      {
        int comp;
        for(comp = 0; comp < (s->has_normal ? 2 : 1); comp++){
          const double _Complex *W = comp ? pr->Wn_hat : pr->Ws_hat;
          double _Complex *H = (comp ? s->Hn[k] : s->Hs[k]) + (size_t)p * s->F * s->M;   /* [f][m] */
          for(q = 0; q < s->F; q++){
            if(q < B && lp + q <= pr->L) memcpy(H + (size_t)q * s->M, W + (size_t)(lp + q - 1) * pr->M, sizeof(double _Complex) * s->M);
            else memset(H + (size_t)q * s->M, 0, sizeof(double _Complex) * s->M);
          }
          fft_rows(H, s->F, s->M, -1);
        }
      }
    }
  }
  /* windows: per source, enough slots for the deepest partition plus slack */
  s->X = malloc(sizeof(double _Complex *) * (size_t)s->nsrc); s->nslot = malloc(sizeof(int) * (size_t)s->nsrc);
  s->Xj = malloc(sizeof(long *) * (size_t)s->nsrc);
  for(a = 0; a < s->nsrc; a++){
    int pm = 0;
    for(b = 0; b < s->nseg; b++) if(s->P[a * s->nseg + b] > pm) pm = s->P[a * s->nseg + b];
    s->nslot[a] = pm + s->q0 + 2;
    s->X[a] = calloc((size_t)s->nslot[a] * s->M * s->F, sizeof(double _Complex));
    s->Xj[a] = malloc(sizeof(long) * (size_t)s->nslot[a]);
  }
  /* time-domain tables are only needed up to lag D now; the ring only up to
     D + 2B + kmax back (partitioned-only mode) */
  if(mode == 2){ fd_sys_reset_history(s); return 0; }
  for(k = 0; k < np; k++){
    struct fd_pair *pr = s->pairs[k];
    if(pr->L > s->D){
      pr->Ws_hat = realloc(pr->Ws_hat, sizeof(double _Complex) * (size_t)s->D * pr->M);
      if(pr->has_normal) pr->Wn_hat = realloc(pr->Wn_hat, sizeof(double _Complex) * (size_t)s->D * pr->M);
    }
  }
  Lmax = s->D + 2 * B + kmax + 1;
  if(Lmax < s->Lbuf){
    for(a = 0; a < s->nsrc; a++){ free(s->hist[a]); s->hist[a] = calloc((size_t)Lmax * s->M, sizeof(double _Complex)); }
    s->Lbuf = Lmax;
  }
  fd_sys_reset_history(s);
  return 0;
}

/* far part (lags > D) for output block c of pair k (source a) */
static void far_block(struct fd_sys *s, int k, int a, long c)
{
  int p, m, f, comp, nc = s->has_normal ? 2 : 1;
  for(comp = 0; comp < nc; comp++){
    double _Complex *Z = comp ? s->tzn : s->tz;
    memset(Z, 0, sizeof(double _Complex) * (size_t)s->M * s->F);
    for(p = 0; p < s->P[k]; p++){
      long j = c - s->q0 - p;
      const double _Complex *H, *Xw;
      if(j < 0) break;                       /* windows before the start are zero */
      H = (comp ? s->Hn[k] : s->Hs[k]) + (size_t)p * s->M * s->F;
      if(s->Xj[a][j % s->nslot[a]] != j) form_window(s, a, j);   /* lazy (hybrid mode) */
      Xw = s->X[a] + (size_t)(j % s->nslot[a]) * s->M * s->F;
      for(m = 0; m < s->M * s->F; m++) Z[m] += H[m] * Xw[m];
    }
    fft_rows(Z, s->F, s->M, +1);
    {
      double _Complex *Y = comp ? s->Yn[k] : s->Ys[k];                    /* [i][m] */
      for(f = 0; f < s->B * s->M; f++) Y[f] = Z[(size_t)s->B * s->M + f] / s->F;
    }
  }
  s->cblk[k] = c;
}

/* history term at the end of a step of k = lag0 - 1 fine substeps that
   have not been pushed yet: sum_{l = lag0}^{L} W[l] V^{(next + lag0 - 1 - l)},
   i.e. lag l counted from the end of the step (the current step's own
   lags 1..k enter through the cumulative kernel K_k). Before 2026-10-19
   this used V^{(next - l)}, a misalignment by k slices (one step for
   k = 1); see fd_c_design.md section 13. The ring buffer is walked in at
   most two contiguous reversed segments. */
void fd_sys_history(struct fd_sys *s, int lag0, double *ts, double *tn)
{
  int a, b, l, lmax, pos, nseg, q, i, m, hk, sh = lag0 - 1, use_far;
  memset(ts, 0, sizeof(double) * s->N); if(tn) memset(tn, 0, sizeof(double) * s->N);
  for(a = 0; a < s->nsrc; a++){          /* source */
    hk = s->hk[a];
    for(b = 0; b < s->nseg; b++){        /* receiver */
      struct fd_pair *p = s->pairs[a * s->nseg + b];
      if(s->Kc && (sh == s->Kc) && (s->coarse_from >= 0) && (hk % s->Kc == 0) &&
         ((long)hk - s->coarse_from >= (long)p->L - s->Kc) && !getenv("FD_COARSE_OFF")){
        /* coarse path: block j of Kc slices ends Kc*j slices before the end of the step */
        int kk = a * s->nseg + b, j;
        memset(s->w1, 0, sizeof(double _Complex) * s->M);
        if(p->has_normal) memset(s->w2, 0, sizeof(double _Complex) * s->M);
        for(j = 1; j <= s->Jc[kk]; j++){
          const double _Complex *Vh = s->hist[a] + (size_t)((hk - j * s->Kc) % s->Lbuf) * s->M;
          const double _Complex *As = s->Acs[kk] + (size_t)j * s->M;
          for(m = 0; m < s->M; m++) s->w1[m] += As[m] * Vh[m];
          if(p->has_normal){
            const double _Complex *An = s->Acn[kk] + (size_t)j * s->M;
            for(m = 0; m < s->M; m++) s->w2[m] += An[m] * Vh[m];
          }
        }
        fd_fft_exec(s->pb, s->w1, s->w3);
        for(i = 0; i < p->rcv->n; i++) ts[s->off[b] + i] += creal(s->w3[i]) / s->M;
        if(p->has_normal && tn){
          fd_fft_exec(s->pb, s->w2, s->w3);
          for(i = 0; i < p->rcv->n; i++) tn[s->off[b] + i] += creal(s->w3[i]) / s->M;
        }
        continue;
      }
      lmax = p->L < hk + sh ? p->L : hk + sh;
      use_far = s->conv && s->P[a * s->nseg + b] > 0 && (s->conv == 1 || 2 * sh <= s->B);
      if(use_far && lmax > s->D) lmax = s->D;
      if(lmax < lag0) continue;
      memset(s->w1, 0, sizeof(double _Complex) * s->M);
      if(p->has_normal) memset(s->w2, 0, sizeof(double _Complex) * s->M);
      l = lag0;
      while(l <= lmax){
        pos = (hk + sh - l) % s->Lbuf;
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
      if(use_far){
        long T = (long)hk + sh, cb = T / s->B;
        int kk = a * s->nseg + b, it = (int)(T - cb * s->B);
        if(s->cblk[kk] != cb) far_block(s, kk, a, cb);
        for(m = 0; m < s->M; m++) s->w1[m] += s->Ys[kk][(size_t)it * s->M + m];
        if(p->has_normal) for(m = 0; m < s->M; m++) s->w2[m] += s->Yn[kk][(size_t)it * s->M + m];
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
  if(s->nops >= 16) return NULL;
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
  fd_band_mult(N, op->bw, op->ab, Vguess, tmp);
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
    fd_band_mult(N, op->bw, op->ab, Vn, tmp);
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
  int i, N = s->N;
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
  fd_band_mult(N, op->bw, op->ab, st->V, st->work);
  for(i = 0; i < N; i++) st->tau[i] = st->H[i] + st->work[i];
  fd_sys_push_k(s, st->V, k);
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
    fd_band_mult(N, op->bw, op->ab, V, tmp);
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
    fd_band_mult(N, op->bw, op->ab, V, tmp);
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
  int i, N = s->N;
  double dt = s->dt * k, t1 = st->t + dt, *psi_s = st->work, *Vs, *sig;
  struct fd_ops *op = fd_sys_ops(s, k);
  if(!op) return -1;
  double tp0 = fd_now(), tp1;
  Vs = (double *)calloc((size_t)N, sizeof(double)); sig = (double *)calloc((size_t)N, sizeof(double));
  for(i = 0; i < N; i++) sig[i] = st->sigma0[i];
  fd_sys_static(s, st->delta, st->H, st->Hn);
  tp1 = fd_now(); fd_prof[1] += tp1 - tp0; tp0 = tp1;
  fd_sys_history(s, k + 1, st->hs, st->hn);
  tp1 = fd_now(); fd_prof[0] += tp1 - tp0; tp0 = tp1;
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
  tp1 = fd_now(); fd_prof[2] += tp1 - tp0; tp0 = tp1;
  memcpy(st->psi, psi_s, sizeof(double) * N);
  fd_band_mult(N, op->bw, op->ab, st->V, Vs);
  for(i = 0; i < N; i++) st->tau[i] = st->H[i] + Vs[i];
  if(st->normal_coupling){ dense_mult(N, op->Kn, st->V, sig); for(i = 0; i < N; i++) st->sigma[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]); }
  else for(i = 0; i < N; i++) st->sigma[i] = st->sigma0[i];
  tp1 = fd_now(); fd_prof[5] += tp1 - tp0; tp0 = tp1;
  fd_sys_push_k(s, st->V, k);
  tp1 = fd_now(); fd_prof[3] += tp1 - tp0; tp0 = tp1;
  for(i = 0; i < N; i++) st->delta[i] += st->V[i] * dt;
  st->t = t1; st->n++; st->k = k;
  free(Vs); free(sig);
  return 0;
}
