/* fully dynamic marcher core, see fd_march.h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd_march.h"
#include <time.h>
double fd_prof[6];   /* history, static, solve, push, ops, other: seconds (FD_PROFILE) */
static double fd_now(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9 * ts.tv_nsec;
}


/* in-place radix-2 FFT along the first index of a [F][M] array (rows of M
   complex values), batched over the M columns so that every butterfly is a
   contiguous, vectorisable loop; sign -1 forward, +1 backward (unnormalised) */
static void fft_rows(double _Complex *a, int F, int M, int sign)
{
  int i, j, bit, len, half, h, m;
  double _Complex *tmp = (double _Complex *)malloc(sizeof(double _Complex) * (size_t)M), w, wl;
  _Complex double *u;
  _Complex double *v;
  _Complex double t;
  for(i = 1, j = 0; i < F; i++){                 /* bit reversal of rows */
    for(bit = F >> 1; j & bit; bit >>= 1)
      j ^= bit;
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
        u = a + (size_t)(i + h) * M;
        v = a + (size_t)(i + h + half) * M;
        for(m = 0; m < M; m++){
          t = w * v[m];
          v[m] = u[m] - t;
          u[m] = u[m] + t;
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
    if(idx >= 0)
      memcpy(Xs + (size_t)rr * s->M, s->hist[k] + (size_t)(idx % s->Lbuf) * s->M, sizeof(double _Complex) * s->M);
    else
      memset(Xs + (size_t)rr * s->M, 0, sizeof(double _Complex) * s->M);
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
  int nmax;
  struct fd_pair *p;
  int r;
  int rev;
  s->segs = segs;
  s->pairs = pairs;
  s->nseg = segs->nseg;
  s->mode = mode;
  s->dt = dt;
  s->nsrc = nsrc > 0 ? nsrc : segs->nseg;
  s->src_of = (int *)malloc(sizeof(int) * (size_t)s->nsrc);
  s->src_rev = (int *)malloc(sizeof(int) * (size_t)s->nsrc);
  for(k = 0; k < s->nsrc; k++){
    s->src_of[k] = src_of ? src_of[k] : k;
    s->src_rev[k] = src_rev ? src_rev[k] : 0;
  }
  s->alpha = alpha;
  s->beta = beta;
  s->rho = rho;
  s->nsign = nsign;
  s->eta = rho * beta * beta / (2.0 * beta);
  s->has_normal = (mode == FD_MODE_INPLANE);
  s->off = (int *)malloc(sizeof(int) * (size_t)(s->nseg + 1));
  s->off[0] = 0;
  for(k = 0; k < s->nseg; k++)
    s->off[k + 1] = s->off[k] + segs->seg[k].n;
  s->N = s->off[s->nseg];
  /* surface segments (M7 step 2b) must follow all fault segments */
  s->Nf = s->N;
  for(k = 0; k < s->nseg; k++)
    if(segs->seg[k].kind != FD_SEG_FAULT){
      s->Nf = s->off[k];
      break;
    }
  for(; k < s->nseg; k++)
    if(segs->seg[k].kind == FD_SEG_FAULT){
      free(s->off);
      free(s);
      return NULL;
    }
  s->Nu = s->N - s->Nf;
  s->etai = (double *)malloc(sizeof(double)*(s->N > 0 ? s->N : 1));
  for(k = 0; k < s->nseg; k++)
    for(i = s->off[k]; i < s->off[k + 1]; i++)
      s->etai[i] = (segs->seg[k].kind == FD_SEG_SURF_OPEN) ? 0.5*rho*alpha : rho*beta*beta/(2.0*beta);
  s->M = 0;
  s->L = 0;
  s->ndense = 0;
  for(k = 0; k < s->nsrc * s->nseg; k++){
    if(pairs[k]->L > s->L)
      s->L = pairs[k]->L;
    if(pairs[k]->dense)
      s->ndense++;
    else if(pairs[k]->M > s->M)
      s->M = pairs[k]->M;      /* Toeplitz pairs share one M */
  }
  if(s->M == 0){                                         /* only dense pairs: the spectral ring is still kept */
    nmax = 0;
    for(k = 0; k < s->nseg; k++)
      if(segs->seg[k].n > nmax)
        nmax = segs->seg[k].n;
    s->M = fd_fft_size_for(nmax);
  }
  s->Lbuf = s->L + 1;
  /* real-space ring of the source slices for dense (non-parallel) pairs (fd_plan.md N1) */
  s->histr = (double **)calloc((size_t)s->nsrc, sizeof(double *));
  if(s->ndense)
    for(k = 0; k < s->nsrc; k++)
      s->histr[k] = (double *)calloc((size_t)s->Lbuf * segs->seg[s->src_of[k]].n, sizeof(double));
  s->hist = (double _Complex **)malloc(sizeof(double _Complex *) * (size_t)s->nsrc);
  s->hk = (int *)calloc((size_t)s->nsrc, sizeof(int));
  s->spec = (int *)calloc((size_t)s->nsrc, sizeof(int));
  for(k = 0; k < s->nsrc; k++){
    for(l = 0; l < s->nseg; l++)
      if(!pairs[k*s->nseg + l]->dense)
        s->spec[k] = 1;
    s->hist[k] = s->spec[k] ? (double _Complex *)calloc((size_t)s->Lbuf * s->M, sizeof(double _Complex)) : NULL;
  }
  s->pf = fd_fft_plan_create(s->M, FD_FFT_FORWARD);
  s->pb = fd_fft_plan_create(s->M, FD_FFT_BACKWARD);
  s->w1 = (double _Complex *)malloc(sizeof(double _Complex) * s->M);
  s->w2 = (double _Complex *)malloc(sizeof(double _Complex) * s->M);
  s->w3 = (double _Complex *)malloc(sizeof(double _Complex) * s->M);
  /* dense static operators in segment order */
  s->Ks = (double *)calloc((size_t)s->N * s->N, sizeof(double));
  s->Kn = (double *)calloc((size_t)s->N * s->N, sizeof(double));
  for(k = 0; k < s->nsrc; k++)
    for(l = 0; l < s->nseg; l++){
      p = pairs[k * s->nseg + l];     /* source k, receiver l */
      r = s->src_of[k];
      rev = s->src_rev[k];
      ns = p->src->n;
      nr = p->rcv->n;
      Ks = (double *)malloc(sizeof(double) * ns * nr);
      Kn = (double *)malloc(sizeof(double) * ns * nr);
      fd_pair_static_dense(p, Ks, Kn);
      for(i = 0; i < nr; i++)
        for(j = 0; j < ns; j++){
          jj = rev ? ns - 1 - j : j;                 /* column of the real element feeding source element j */
          s->Ks[(size_t)(s->off[l] + i) * s->N + s->off[r] + jj] += Ks[i * ns + j];
          s->Kn[(size_t)(s->off[l] + i) * s->N + s->off[r] + jj] += Kn[i * ns + j];
        }
      free(Ks);
      free(Kn);
    }
  s->nops = 0;
  return s;
}

void fd_sys_free(struct fd_sys *s)
{
  int k;
  int np;
  if(!s)
    return;
  for(k = 0; k < s->nsrc; k++){
    free(s->hist[k]);
    free(s->histr[k]);
  }
  free(s->histr);
  free(s->ctP);
  free(s->ctV);
  free(s->etai);
  free(s->spec);
  if(s->conv){
    np = s->nsrc * s->nseg;
    for(k = 0; k < np; k++){
      free(s->Hs[k]);
      if(s->Hn)
        free(s->Hn[k]);
      free(s->Ys[k]);
      if(s->Yn)
        free(s->Yn[k]);
    }
    for(k = 0; k < s->nsrc; k++){
      free(s->X[k]);
      free(s->Xj[k]);
    }
    free(s->Xj);
    free(s->Hs);
    free(s->Hn);
    free(s->Ys);
    free(s->Yn);
    free(s->X);
    free(s->P);
    free(s->nslot);
    free(s->cblk);
    free(s->tw1);
    free(s->tw2);
    free(s->tz);
    free(s->tzn);
    fd_fft_plan_destroy(s->tf);
    fd_fft_plan_destroy(s->tb);
  }
  if(s->Kc){
    for(k = 0; k < s->nsrc * s->nseg; k++){
      free(s->Acs[k]);
      free(s->Acn[k]);
      free(s->Bcs[k]);
      free(s->Bcn[k]);
    }
    free(s->Acs);
    free(s->Acn);
    free(s->Bcs);
    free(s->Bcn);
    free(s->Jc);
  }
  free(s->hist);
  free(s->hk);
  free(s->off);
  free(s->src_of);
  free(s->src_rev);
  free(s->w1);
  free(s->w2);
  free(s->w3);
  free(s->Ks);
  free(s->Kn);
  for(k = 0; k < s->nops; k++){
    free(s->ops[k].Ks);
    free(s->ops[k].Kn);
    free(s->ops[k].ab);
    free(s->ops[k].Cs);
    free(s->ops[k].Cn);
    if(s->ops[k].el){
      free(s->ops[k].el->Ks);
      free(s->ops[k].el->Kn);
      free(s->ops[k].el->ab);
      free(s->ops[k].el);
    }
    free(s->ops[k].lu);
    free(s->ops[k].piv);
    free(s->ops[k].P);
    free(s->ops[k].Bfs);
    free(s->ops[k].Knfs);
  }
  fd_fft_plan_destroy(s->pf);
  fd_fft_plan_destroy(s->pb);
  free(s);
}

void fd_sys_reset_history(struct fd_sys *s)
{
  int k;
  int q;
  s->coarse_from = -1;
  for(k = 0; k < s->nsrc; k++){
    if(s->hist[k])
      memset(s->hist[k], 0, sizeof(double _Complex) * (size_t)s->Lbuf * s->M);
    s->hk[k] = 0;
    if(s->histr[k])
      memset(s->histr[k], 0, sizeof(double) * (size_t)s->Lbuf * s->segs->seg[s->src_of[k]].n);
  }
  if(s->conv){
    for(k = 0; k < s->nsrc; k++){
      memset(s->X[k], 0, sizeof(double _Complex) * (size_t)s->nslot[k] * s->M * s->F);
      for(q = 0; q < s->nslot[k]; q++)
        s->Xj[k][q] = -1;
    }
    for(k = 0; k < s->nsrc * s->nseg; k++)
      s->cblk[k] = -1;
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
  _Complex double *first;
  _Complex double *slot;
  double *rs;
  if(s->Kc){
    if((kk == s->Kc) && (s->hk[0] % s->Kc == 0)){
      if(s->coarse_from < 0)
        s->coarse_from = s->hk[0];
    } else
      s->coarse_from = -1;
  }
  for(k = 0; k < s->nsrc; k++){
    first = NULL;
    r = s->src_of[k];
    n = s->segs->seg[r].n;
    for(q = 0; q < kk; q++){
      if(s->histr[k]){
        rs = s->histr[k] + (size_t)(s->hk[k] % s->Lbuf) * n;
        for(i = 0; i < n; i++)
          rs[i] = V[s->off[r] + (s->src_rev[k] ? n - 1 - i : i)];
      }
      if(!s->spec[k]){                           /* dense pairs only: no spectral ring */
        s->hk[k]++;
        continue;
      }
      slot = s->hist[k] + (size_t)(s->hk[k] % s->Lbuf) * s->M;
      if(!first){
        memset(s->w1, 0, sizeof(double _Complex) * s->M);
        for(i = 0; i < n; i++)
          s->w1[i] = V[s->off[r] + (s->src_rev[k] ? n - 1 - i : i)];
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

void fd_sys_push_lin(struct fd_sys *s, const double *Vold, const double *V, int kk)
{
  int k, i, r, n, q;
  double w;
  _Complex double *slot;
  int m;
  double *rs;
  int ii;
  s->lin = 1;
  if(s->Kc){                                     /* aligned steps of Kc substeps: a block is one step */
    if((kk == s->Kc) && (s->hk[0] % s->Kc == 0)){
      if(s->coarse_from < 0)
        s->coarse_from = s->hk[0];
    } else
      s->coarse_from = -1;
  }
  for(k = 0; k < s->nsrc; k++){
    r = s->src_of[k];
    n = s->segs->seg[r].n;
    if(s->spec[k]){
      memset(s->w1, 0, sizeof(double _Complex) * s->M);
      memset(s->w3, 0, sizeof(double _Complex) * s->M);
      for(i = 0; i < n; i++){
        s->w1[i] = Vold[s->off[r] + (s->src_rev[k] ? n - 1 - i : i)];
        s->w3[i] = V[s->off[r] + (s->src_rev[k] ? n - 1 - i : i)];
      }
      fd_fft_exec(s->pf, s->w1, s->w2);          /* w2 = F(Vold) */
      fd_fft_exec(s->pf, s->w3, s->w1);          /* w1 = F(V) */
    }
    for(q = 1; q <= kk; q++){
      w = (q - 0.5) / kk;
      if(s->spec[k]){
        slot = s->hist[k] + (size_t)(s->hk[k] % s->Lbuf) * s->M;
        for(m = 0; m < s->M; m++)
          slot[m] = (1.0 - w) * s->w2[m] + w * s->w1[m];
      }
      if(s->histr[k]){
        rs = s->histr[k] + (size_t)(s->hk[k] % s->Lbuf) * n;
        for(i = 0; i < n; i++){
          ii = s->off[r] + (s->src_rev[k] ? n - 1 - i : i);
          rs[i] = (1.0 - w) * Vold[ii] + w * V[ii];
        }
      }
      s->hk[k]++;
      if((s->conv == 1) && (s->hk[k] % s->B == 0))
        form_window(s, k, s->hk[k] / s->B - 1);
    }
  }
}

int fd_sys_set_coarse(struct fd_sys *s, int K)
{
  int np = s->nsrc * s->nseg, k, j, l, m, comp;
  struct fd_pair *p;
  const _Complex double *W;
  _Complex double *A;
  _Complex double *Bm;
  double w;
  if(K < 2 || s->conv == 1)
    return 1;
  s->Kc = K;
  s->coarse_from = -1;
  s->Jc = calloc((size_t)np, sizeof(int));
  s->Acs = calloc((size_t)np, sizeof(double _Complex *));
  s->Acn = calloc((size_t)np, sizeof(double _Complex *));
  s->Bcs = calloc((size_t)np, sizeof(double _Complex *));
  s->Bcn = calloc((size_t)np, sizeof(double _Complex *));
  for(k = 0; k < np; k++){
    p = s->pairs[k];
    if(p->dense)
      continue;                         /* dense pairs: direct sum (fd_sys_history) */
    s->Jc[k] = (p->L - 1) / K;                     /* blocks j = 1..Jc with jK + 1 <= L */
    for(comp = 0; comp < (p->has_normal ? 2 : 1); comp++){
      W = comp ? p->Wn_hat : p->Ws_hat;
      A = calloc((size_t)(s->Jc[k] + 1) * p->M, sizeof(double _Complex));
      Bm = calloc((size_t)(s->Jc[k] + 1) * p->M, sizeof(double _Complex));
      for(j = 1; j <= s->Jc[k]; j++)
        for(l = j * K + 1; l <= (j + 1) * K && l <= p->L; l++){
          w = ((j + 1) * K - l + 0.5) / K;
          for(m = 0; m < p->M; m++){
            A[(size_t)j * p->M + m] += W[(size_t)(l - 1) * p->M + m];
            Bm[(size_t)j * p->M + m] += w * W[(size_t)(l - 1) * p->M + m];
          }
        }
      if(comp){
        s->Acn[k] = A;
        s->Bcn[k] = Bm;
      } else {
        s->Acs[k] = A;
        s->Bcs[k] = Bm;
      }
    }
  }
  return 0;
}

double fd_sys_bytes(const struct fd_sys *s)
{
  double b = 16.0 * s->nsrc * (double)s->Lbuf * s->M;
  int a, k;
  const struct fd_pair *p;
  int lt;
  for(k = 0; k < s->nsrc * s->nseg; k++){
    p = s->pairs[k];
    lt = (s->conv == 1) ? (p->L < s->D ? p->L : s->D) : p->L;
    if(p->dense){
      b += fd_pair_dense_bytes(p, NULL);       /* compressed (D2) or full windows */
      continue;
    }
    b += 16.0 * lt * p->M * (p->has_normal ? 2 : 1);
    if(s->conv)
      b += 16.0 * (double)s->P[k] * s->M * s->F * (p->has_normal ? 2 : 1) + 16.0 * s->M * s->B * (p->has_normal ? 2 : 1);
  }
  if(s->conv)
    for(a = 0; a < s->nsrc; a++)
      b += 16.0 * (double)s->nslot[a] * s->M * s->F;
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
  int qq;
  struct fd_pair *pr;
  int lp;
  int comp;
  const _Complex double *W;
  _Complex double *H;
  int pm;
  struct fd_pair *pr_2;
  if(s->ndense)
    return 0;                          /* partitioned convolution needs Toeplitz pairs: off */
  int np = s->nsrc * s->nseg, k, p, q, a, b, Lmax, best = 0;
  double cost, cbest = 1e300;
  if(s->conv)
    return 0;
  if(kmax < 1)
    kmax = 1;
  if(B <= 0){
    /* minimise direct lags D(B) + far work 2 L / B over powers of two */
    for(b = 16; b <= 4096; b *= 2){
      qq = 1 + (kmax + b - 1) / b;
      cost = (double)(qq * b - 1) + 2.0 * s->L / b + 4.0 * log2(2.0 * b);
      if(cost < cbest){
        cbest = cost;
        best = b;
      }
    }
    B = best;
  }
  if(B < 2 || (B & (B - 1)))
    return 1;
  s->B = B;
  s->F = 2 * B;
  s->kmax = kmax;
  s->q0 = 1 + (kmax + B - 1) / B;           /* (q0 - 1) B >= kmax */
  s->D = s->q0 * B - 1;
  if(s->D >= s->L)
    return 0;               /* window shorter than the direct part: nothing to gain */
  s->conv = mode;
  s->tf = fd_fft_plan_create(s->F, FD_FFT_FORWARD);
  s->tb = fd_fft_plan_create(s->F, FD_FFT_BACKWARD);
  s->tw1 = malloc(sizeof(double _Complex) * s->F);
  s->tw2 = malloc(sizeof(double _Complex) * s->F);
  s->tz = malloc(sizeof(double _Complex) * (size_t)s->M * s->F);
  s->tzn = malloc(sizeof(double _Complex) * (size_t)s->M * s->F);
  s->P = calloc((size_t)np, sizeof(int));
  s->cblk = malloc(sizeof(long) * (size_t)np);
  s->Hs = calloc((size_t)np, sizeof(double _Complex *));
  s->Hn = s->has_normal ? calloc((size_t)np, sizeof(double _Complex *)) : NULL;
  s->Ys = calloc((size_t)np, sizeof(double _Complex *));
  s->Yn = s->has_normal ? calloc((size_t)np, sizeof(double _Complex *)) : NULL;
  for(k = 0; k < np; k++){
    pr = s->pairs[k];
    s->cblk[k] = -1;
    s->P[k] = pr->L > s->D ? (pr->L - s->D + B - 1) / B : 0;
    s->Hs[k] = calloc((size_t)s->P[k] * s->M * s->F + 1, sizeof(double _Complex));
    if(s->has_normal)
      s->Hn[k] = calloc((size_t)s->P[k] * s->M * s->F + 1, sizeof(double _Complex));
    s->Ys[k] = calloc((size_t)s->M * B, sizeof(double _Complex));
    if(s->has_normal)
      s->Yn[k] = calloc((size_t)s->M * B, sizeof(double _Complex));
    for(p = 0; p < s->P[k]; p++){
      lp = s->D + 1 + p * B;
      {
        for(comp = 0; comp < (s->has_normal ? 2 : 1); comp++){
          W = comp ? pr->Wn_hat : pr->Ws_hat;
          H = (comp ? s->Hn[k] : s->Hs[k]) + (size_t)p * s->F * s->M;   /* [f][m] */
          for(q = 0; q < s->F; q++){
            if(q < B && lp + q <= pr->L)
              memcpy(H + (size_t)q * s->M, W + (size_t)(lp + q - 1) * pr->M, sizeof(double _Complex) * s->M);
            else
              memset(H + (size_t)q * s->M, 0, sizeof(double _Complex) * s->M);
          }
          fft_rows(H, s->F, s->M, -1);
        }
      }
    }
  }
  /* windows: per source, enough slots for the deepest partition plus slack */
  s->X = malloc(sizeof(double _Complex *) * (size_t)s->nsrc);
  s->nslot = malloc(sizeof(int) * (size_t)s->nsrc);
  s->Xj = malloc(sizeof(long *) * (size_t)s->nsrc);
  for(a = 0; a < s->nsrc; a++){
    pm = 0;
    for(b = 0; b < s->nseg; b++)
      if(s->P[a * s->nseg + b] > pm)
        pm = s->P[a * s->nseg + b];
    s->nslot[a] = pm + s->q0 + 2;
    s->X[a] = calloc((size_t)s->nslot[a] * s->M * s->F, sizeof(double _Complex));
    s->Xj[a] = malloc(sizeof(long) * (size_t)s->nslot[a]);
  }
  /* time-domain tables are only needed up to lag D now; the ring only up to
     D + 2B + kmax back (partitioned-only mode) */
  if(mode == 2){
    fd_sys_reset_history(s);
    return 0;
  }
  for(k = 0; k < np; k++){
    pr_2 = s->pairs[k];
    if(pr_2->L > s->D){
      pr_2->Ws_hat = realloc(pr_2->Ws_hat, sizeof(double _Complex) * (size_t)s->D * pr_2->M);
      if(pr_2->has_normal)
        pr_2->Wn_hat = realloc(pr_2->Wn_hat, sizeof(double _Complex) * (size_t)s->D * pr_2->M);
    }
  }
  Lmax = s->D + 2 * B + kmax + 1;
  if(Lmax < s->Lbuf){
    for(a = 0; a < s->nsrc; a++){
      free(s->hist[a]);
      s->hist[a] = calloc((size_t)Lmax * s->M, sizeof(double _Complex));
    }
    s->Lbuf = Lmax;
  }
  fd_sys_reset_history(s);
  return 0;
}

/* far part (lags > D) for output block c of pair k (source a) */
static void far_block(struct fd_sys *s, int k, int a, long c)
{
  int p, m, f, comp, nc = s->has_normal ? 2 : 1;
  _Complex double *Z;
  long j;
  const _Complex double *H;
  const _Complex double *Xw;
  _Complex double *Y;
  for(comp = 0; comp < nc; comp++){
    Z = comp ? s->tzn : s->tz;
    memset(Z, 0, sizeof(double _Complex) * (size_t)s->M * s->F);
    for(p = 0; p < s->P[k]; p++){
      j = c - s->q0 - p;

      if(j < 0)
        break;                       /* windows before the start are zero */
      H = (comp ? s->Hn[k] : s->Hs[k]) + (size_t)p * s->M * s->F;
      if(s->Xj[a][j % s->nslot[a]] != j)
        form_window(s, a, j);   /* lazy (hybrid mode) */
      Xw = s->X[a] + (size_t)(j % s->nslot[a]) * s->M * s->F;
      for(m = 0; m < s->M * s->F; m++)
        Z[m] += H[m] * Xw[m];
    }
    fft_rows(Z, s->F, s->M, +1);
    {
      Y = comp ? s->Yn[k] : s->Ys[k];                    /* [i][m] */
      for(f = 0; f < s->B * s->M; f++)
        Y[f] = Z[(size_t)s->B * s->M + f] / s->F;
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
  struct fd_pair *p;
  int nr;
  int nsr;
  int j;
  size_t np2;
  const double *v;
  const double *Ws;
  const double *Wn;
  double acs;
  double acn;
  int kk;
  int j_2;
  const _Complex double *Vh;
  const _Complex double *As;
  const _Complex double *An;
  const _Complex double *V1;
  const _Complex double *VK;
  const _Complex double *Bs;
  double w1;
  double fac;
  const _Complex double *An_2;
  const _Complex double *Bn;
  const _Complex double *Vh_2;
  const _Complex double *Ws_2;
  const _Complex double *Wn_2;
  long T;
  long cb;
  int kk_2;
  int it;
  int have_pref;
  int Lp;
  size_t need;
  memset(ts, 0, sizeof(double) * s->N);
  if(tn)
    memset(tn, 0, sizeof(double) * s->N);
  for(a = 0; a < s->nsrc; a++){          /* source */
    hk = s->hk[a];
    have_pref = 0;
    for(b = 0; b < s->nseg; b++){        /* receiver */
      p = s->pairs[a * s->nseg + b];
      if(p->dense && p->compressed){
        /* compressed dense pair (fd_plan.md D2): prefix moments of the source's history, built
           once per source and call */
        nsr = p->src->n;
        lmax = p->L < hk + sh ? p->L : hk + sh;
        if(!have_pref){
          Lp = 0;
          for(q = 0; q < s->nseg; q++)
            if(s->pairs[a*s->nseg + q]->L > Lp)
              Lp = s->pairs[a*s->nseg + q]->L;
          need = (size_t)3*(Lp + 1)*nsr;
          if(need > s->ctP_cap){
            free(s->ctP);
            s->ctP = (double *)malloc(sizeof(double)*need);
            s->ctP_cap = need;
          }
          if(Lp + 1 > s->ctV_cap){
            free(s->ctV);
            s->ctV = (const double **)malloc(sizeof(double *)*(Lp + 1));
            s->ctV_cap = Lp + 1;
          }
          fd_ctab_prefix(s->histr[a], s->Lbuf, (long)hk + sh, lag0, hk + sh, Lp, nsr, s->ctP, s->ctV);
          have_pref = 1;
        }
        fd_ctab_hist(&p->cs, p->rcv->n, nsr, s->ctP, s->ctV, lag0, lmax, Lp, ts + s->off[b]);
        if(p->has_normal && tn)
          fd_ctab_hist(&p->cn, p->rcv->n, nsr, s->ctP, s->ctV, lag0, lmax, Lp, tn + s->off[b]);
        continue;
      }
      if(p->dense){
        /* non-parallel pair: direct sum over the real-space ring, W[l] dense (receiver x source) */
        nr = p->rcv->n;
        nsr = p->src->n;
        np2 = (size_t)nr * nsr;
        lmax = p->L < hk + sh ? p->L : hk + sh;
        for(l = lag0; l <= lmax; l++){
          v = s->histr[a] + (size_t)((hk + sh - l) % s->Lbuf) * nsr;
          Ws = p->Wd_s + (size_t)(l - 1) * np2;
          Wn = p->Wd_n ? p->Wd_n + (size_t)(l - 1) * np2 : NULL;
          for(i = 0; i < nr; i++){
            acs = 0.0;
            acn = 0.0;
            for(j = 0; j < nsr; j++){
              acs += Ws[(size_t)i * nsr + j] * v[j];
              if(Wn)
                acn += Wn[(size_t)i * nsr + j] * v[j];
            }
            ts[s->off[b] + i] += acs;
            if(Wn && tn)
              tn[s->off[b] + i] += acn;
          }
        }
        continue;
      }
      if(s->Kc && (sh == s->Kc) && (s->coarse_from >= 0) && (hk % s->Kc == 0) &&
          ((long)hk - s->coarse_from >= (long)p->L - s->Kc) && !getenv("FD_COARSE_OFF")){
        /* coarse path: block j of Kc slices ends Kc*j slices before the end of the step */
        kk = a * s->nseg + b;
        memset(s->w1, 0, sizeof(double _Complex) * s->M);
        if(p->has_normal)
          memset(s->w2, 0, sizeof(double _Complex) * s->M);
        for(j_2 = 1; j_2 <= s->Jc[kk]; j_2++){
          Vh = s->hist[a] + (size_t)((hk - j_2 * s->Kc) % s->Lbuf) * s->M;
          As = s->Acs[kk] + (size_t)j_2 * s->M;
          for(m = 0; m < s->M; m++)
            s->w1[m] += As[m] * Vh[m];
          if(p->has_normal){
            An = s->Acn[kk] + (size_t)j_2 * s->M;
            for(m = 0; m < s->M; m++)
              s->w2[m] += An[m] * Vh[m];
          }
          if(s->lin){
            /* rate linear within the block: slices s_q = Vo + w_q (Vn - Vo), w_q = (q - 1/2)/Kc;
               sum = A s_1 + (B - w_1 A)(s_Kc - s_1) Kc/(Kc - 1) */
            V1 = Vh;
            VK = s->hist[a] + (size_t)((hk - (j_2 - 1) * s->Kc - 1) % s->Lbuf) * s->M;
            Bs = s->Bcs[kk] + (size_t)j_2 * s->M;
            w1 = 0.5 / s->Kc;
            fac = (double)s->Kc / (s->Kc - 1);
            for(m = 0; m < s->M; m++)
              s->w1[m] += (Bs[m] - w1 * As[m]) * (VK[m] - V1[m]) * fac;
            if(p->has_normal){
              An_2 = s->Acn[kk] + (size_t)j_2 * s->M;
              Bn = s->Bcn[kk] + (size_t)j_2 * s->M;
              for(m = 0; m < s->M; m++)
                s->w2[m] += (Bn[m] - w1 * An_2[m]) * (VK[m] - V1[m]) * fac;
            }
          }
        }
        fd_fft_exec(s->pb, s->w1, s->w3);
        for(i = 0; i < p->rcv->n; i++)
          ts[s->off[b] + i] += creal(s->w3[i]) / s->M;
        if(p->has_normal && tn){
          fd_fft_exec(s->pb, s->w2, s->w3);
          for(i = 0; i < p->rcv->n; i++)
            tn[s->off[b] + i] += creal(s->w3[i]) / s->M;
        }
        continue;
      }
      lmax = p->L < hk + sh ? p->L : hk + sh;
      use_far = s->conv && s->P[a * s->nseg + b] > 0 && (s->conv == 1 || 2 * sh <= s->B);
      if(use_far && lmax > s->D)
        lmax = s->D;
      if(lmax < lag0)
        continue;
      memset(s->w1, 0, sizeof(double _Complex) * s->M);
      if(p->has_normal)
        memset(s->w2, 0, sizeof(double _Complex) * s->M);
      l = lag0;
      while(l <= lmax){
        pos = (hk + sh - l) % s->Lbuf;
        nseg = lmax - l + 1 < pos + 1 ? lmax - l + 1 : pos + 1;
        for(q = 0; q < nseg; q++){
          Vh_2 = s->hist[a] + (size_t)(pos - q) * s->M;
          Ws_2 = p->Ws_hat + (size_t)(l - 1 + q) * s->M;
          for(m = 0; m < s->M; m++)
            s->w1[m] += Ws_2[m] * Vh_2[m];
          if(p->has_normal){
            Wn_2 = p->Wn_hat + (size_t)(l - 1 + q) * s->M;
            for(m = 0; m < s->M; m++)
              s->w2[m] += Wn_2[m] * Vh_2[m];
          }
        }
        l += nseg;
      }
      if(use_far){
        T = (long)hk + sh;
        cb = T / s->B;
        kk_2 = a * s->nseg + b;
        it = (int)(T - cb * s->B);
        if(s->cblk[kk_2] != cb)
          far_block(s, kk_2, a, cb);
        for(m = 0; m < s->M; m++)
          s->w1[m] += s->Ys[kk_2][(size_t)it * s->M + m];
        if(p->has_normal)
          for(m = 0; m < s->M; m++)
            s->w2[m] += s->Yn[kk_2][(size_t)it * s->M + m];
      }
      fd_fft_exec(s->pb, s->w1, s->w3);
      for(i = 0; i < p->rcv->n; i++)
        ts[s->off[b] + i] += creal(s->w3[i]) / s->M;
      if(p->has_normal && tn){
        fd_fft_exec(s->pb, s->w2, s->w3);
        for(i = 0; i < p->rcv->n; i++)
          tn[s->off[b] + i] += creal(s->w3[i]) / s->M;
      }
    }
  }
}

void fd_sys_static(struct fd_sys *s, const double *delta, double *ts, double *tn)
{
  int i, j;
  double a;
  double b;
  const double *rs;
  const double *rn;
  for(i = 0; i < s->N; i++){
    a = 0.0;
    b = 0.0;
    rs = s->Ks + (size_t)i * s->N;
    rn = s->Kn + (size_t)i * s->N;
    for(j = 0; j < s->N; j++){
      a += rs[j] * delta[j];
      b += rn[j] * delta[j];
    }
    ts[i] = a;
    if(tn)
      tn[i] = b;
  }
}

/* ---------------- virtual free surface: elimination of the surface unknowns (M7 step 2b) ----------------
   The surface rows say H_s + B_ss V_s + B_sf V_f = 0 (zero traction), so V_s = -(q + P V_f) with
   q = B_ss^-1 H_s and P = B_ss^-1 B_sf; the fault rows become tau_f = (H_f - B_fs q) +
   (B_ff - B_fs P) V_f and the normal rows (Hn_f - Kn_fs q) + (Kn_ff - Kn_fs P) V_f. The friction
   solvers then work on the fault alone with op->el. B_ss is dense and factorised once per k. */
extern void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
/* the trailing size_t is the hidden length of the character argument (gfortran ABI) */
extern void dgetrs_(char *tr, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info, size_t ltr);

static int surf_eliminate(struct fd_sys *s, struct fd_ops *op)
{
  int N, Nf, Nu, i, j, m, info, bw;
  char tr;
  double *Ae, *Ke, amax, thr, acc;
  struct fd_ops *el;
  N = s->N;
  Nf = s->Nf;
  Nu = s->Nu;
  op->lu = (double *)malloc(sizeof(double)*(size_t)Nu*Nu);
  op->piv = (int *)malloc(sizeof(int)*(size_t)Nu);
  op->P = (double *)malloc(sizeof(double)*(size_t)Nu*(Nf > 0 ? Nf : 1));
  op->Bfs = (double *)malloc(sizeof(double)*(size_t)(Nf > 0 ? Nf : 1)*Nu);
  op->Knfs = (double *)malloc(sizeof(double)*(size_t)(Nf > 0 ? Nf : 1)*Nu);
  for(j = 0; j < Nu; j++)
    for(i = 0; i < Nu; i++)
      op->lu[(size_t)j*Nu + i] = op->Ks[(size_t)(Nf + i)*N + Nf + j];
  for(j = 0; j < Nf; j++)
    for(i = 0; i < Nu; i++)
      op->P[(size_t)j*Nu + i] = op->Ks[(size_t)(Nf + i)*N + j];
  for(i = 0; i < Nf; i++)
    for(j = 0; j < Nu; j++){
      op->Bfs[(size_t)i*Nu + j] = op->Ks[(size_t)i*N + Nf + j];
      op->Knfs[(size_t)i*Nu + j] = op->Kn[(size_t)i*N + Nf + j];
    }
  dgetrf_(&Nu, &Nu, op->lu, &Nu, op->piv, &info);
  if(info != 0){
    fprintf(stderr, "fd_march: surface operator singular (k %d, dgetrf %d)\n", op->k, info);
    return 1;
  }
  tr = 'N';
  if(Nf > 0)
    dgetrs_(&tr, &Nu, &Nf, op->lu, &Nu, op->piv, op->P, &Nu, &info, 1);
  el = (struct fd_ops *)calloc(1, sizeof(struct fd_ops));
  el->k = op->k;
  el->nl = op->nl;
  Ae = (double *)malloc(sizeof(double)*(size_t)(Nf > 0 ? Nf : 1)*Nf);
  Ke = (double *)malloc(sizeof(double)*(size_t)(Nf > 0 ? Nf : 1)*Nf);
  amax = 0.0;
  for(i = 0; i < Nf; i++)
    for(j = 0; j < Nf; j++){
      acc = op->Ks[(size_t)i*N + j];
      for(m = 0; m < Nu; m++)
        acc -= op->Bfs[(size_t)i*Nu + m]*op->P[(size_t)j*Nu + m];
      Ae[(size_t)i*Nf + j] = acc;
      if(fabs(acc) > amax)
        amax = fabs(acc);
      acc = op->Kn[(size_t)i*N + j];
      for(m = 0; m < Nu; m++)
        acc -= op->Knfs[(size_t)i*Nu + m]*op->P[(size_t)j*Nu + m];
      Ke[(size_t)i*Nf + j] = acc;
    }
  thr = 1e-13*amax;                             /* as for the full operator: prune round-off for the band */
  for(i = 0; i < Nf*Nf; i++)
    if(fabs(Ae[i]) <= thr)
      Ae[i] = 0.0;
  el->Ks = Ae;
  el->Kn = Ke;
  bw = fd_band_width(Nf, Ae, 0.0);
  el->bw = bw;
  el->ab = (double *)malloc(sizeof(double)*(size_t)(3*bw + 1)*(Nf > 0 ? Nf : 1));
  fd_band_from_dense(Nf, bw, Ae, el->ab);
  op->el = el;
  return 0;
}

/* the fault's known parts with the surface eliminated: He = H_f - B_fs q, Hne = Hn_f - Kn_fs q,
   q = B_ss^-1 H_s (returned) */
static void surf_reduce(const struct fd_sys *s, const struct fd_ops *op, const double *H, const double *Hn, double *He, double *Hne, double *q)
{
  int i, m, Nf, Nu, one, info;
  char tr;
  double acc;
  Nf = s->Nf;
  Nu = s->Nu;
  for(m = 0; m < Nu; m++)
    q[m] = H[Nf + m];
  tr = 'N';
  one = 1;
  dgetrs_(&tr, &Nu, &one, (double *)op->lu, &Nu, (int *)op->piv, q, &Nu, &info, 1);
  for(i = 0; i < Nf; i++){
    acc = H[i];
    for(m = 0; m < Nu; m++)
      acc -= op->Bfs[(size_t)i*Nu + m]*q[m];
    He[i] = acc;
    if(Hn && Hne){
      acc = Hn[i];
      for(m = 0; m < Nu; m++)
        acc -= op->Knfs[(size_t)i*Nu + m]*q[m];
      Hne[i] = acc;
    }
  }
}
/* the surface rates from the fault rates: V_s = -(q + P V_f) */
static void surf_recover(const struct fd_sys *s, const struct fd_ops *op, const double *q, double *V)
{
  int j, m, Nf, Nu;
  double acc;
  Nf = s->Nf;
  Nu = s->Nu;
  for(m = 0; m < Nu; m++){
    acc = q[m];
    for(j = 0; j < Nf; j++)
      acc += op->P[(size_t)j*Nu + m]*V[j];
    V[Nf + m] = -acc;
  }
}

/* K_k = sum_{l<=k} W[l] + k Ks dt, dense in segment order, plus its band */
struct fd_ops *fd_sys_ops(struct fd_sys *s, int k)
{
  int q, a, b, l, i, j, ns, nr, bw;
  double *Ws, *Wn, thr, amax = 0.0;
  struct fd_ops *op;
  struct fd_pair *p;
  int r;
  int rev;
  int jj;
  for(q = 0; q < s->nops; q++)
    if(s->ops[q].k == k && !s->ops[q].nl)
      return s->ops + q;
  if(s->nops >= 16)
    return NULL;
  op = s->ops + s->nops;
  memset(op, 0, sizeof(*op));
  op->k = k;
  op->Ks = (double *)malloc(sizeof(double) * (size_t)s->N * s->N);
  op->Kn = (double *)malloc(sizeof(double) * (size_t)s->N * s->N);
  for(i = 0; i < s->N * s->N; i++){
    op->Ks[i] = s->Ks[i] * s->dt * k;
    op->Kn[i] = s->Kn[i] * s->dt * k;
  }
  for(a = 0; a < s->nsrc; a++)
    for(b = 0; b < s->nseg; b++){
      p = s->pairs[a * s->nseg + b];
      r = s->src_of[a];
      rev = s->src_rev[a];
      ns = p->src->n;
      nr = p->rcv->n;
      Ws = (double *)malloc(sizeof(double) * ns * nr);
      Wn = (double *)malloc(sizeof(double) * ns * nr);
      for(l = 1; l <= (k < p->L ? k : p->L); l++){
        fd_pair_W_dense(p, l, Ws, Wn);
        for(i = 0; i < nr; i++)
          for(j = 0; j < ns; j++){
            jj = rev ? ns - 1 - j : j;
            op->Ks[(size_t)(s->off[b] + i) * s->N + s->off[r] + jj] += Ws[i * ns + j];
            op->Kn[(size_t)(s->off[b] + i) * s->N + s->off[r] + jj] += Wn[i * ns + j];
          }
      }
      free(Ws);
      free(Wn);
    }
  /* prune FFT round-trip noise so that the band is the physical one */
  for(i = 0; i < s->N * s->N; i++)
    if(fabs(op->Ks[i]) > amax)
      amax = fabs(op->Ks[i]);
  thr = 1e-13 * amax;
  for(i = 0; i < s->N * s->N; i++){
    if(fabs(op->Ks[i]) <= thr)
      op->Ks[i] = 0.0;
    if(fabs(op->Kn[i]) <= thr)
      op->Kn[i] = 0.0;
  }
  bw = fd_band_width(s->N, op->Ks, 0.0);
  op->bw = bw;
  op->ab = (double *)malloc(sizeof(double) * (size_t)(3 * bw + 1) * s->N);
  fd_band_from_dense(s->N, bw, op->Ks, op->ab);
  if(s->Nu > 0 && surf_eliminate(s, op))
    return NULL;
  s->nops++;
  return op;
}

/* B_k and C_k for second-order marching (see fd_ops): with lag l = k - j + 1 and
   w(l) = (k - l + 1/2)/k, B_k = sum_l w(l) K[l] + (w(1) - 1) eta I and
   C_k = sum_l (1 - w(l)) K[l] + (1 - w(1)) eta I, K[l] = W[l] + Ks dt; eta on the shear diagonal only */
struct fd_ops *fd_sys_ops_nl(struct fd_sys *s, int k)
{
  int q, a, b, l, i, j, ns, nr, bw;
  double *Ws, *Wn, thr, amax = 0.0, w1 = (k - 0.5) / k;
  size_t NN = (size_t)s->N * s->N, ii;
  struct fd_ops *op;
  struct fd_pair *p;
  int r;
  int rev;
  int jj;
  double w;
  size_t idx;
  for(q = 0; q < s->nops; q++)
    if(s->ops[q].k == k && s->ops[q].nl)
      return s->ops + q;
  if(s->nops >= 16)
    return NULL;
  op = s->ops + s->nops;
  memset(op, 0, sizeof(*op));
  op->k = k;
  op->nl = 1;
  op->Ks = (double *)malloc(sizeof(double) * NN);
  op->Kn = (double *)malloc(sizeof(double) * NN);
  op->Cs = (double *)malloc(sizeof(double) * NN);
  op->Cn = (double *)malloc(sizeof(double) * NN);
  for(ii = 0; ii < NN; ii++){                     /* static part: sum_l w(l) = sum_l (1 - w(l)) = k/2 */
    op->Ks[ii] = op->Cs[ii] = s->Ks[ii] * s->dt * 0.5 * k;
    op->Kn[ii] = op->Cn[ii] = s->Kn[ii] * s->dt * 0.5 * k;
  }
  for(a = 0; a < s->nsrc; a++)
    for(b = 0; b < s->nseg; b++){
      p = s->pairs[a * s->nseg + b];
      r = s->src_of[a];
      rev = s->src_rev[a];
      ns = p->src->n;
      nr = p->rcv->n;
      Ws = (double *)malloc(sizeof(double) * ns * nr);
      Wn = (double *)malloc(sizeof(double) * ns * nr);
      for(l = 1; l <= (k < p->L ? k : p->L); l++){
        w = (k - l + 0.5) / k;
        fd_pair_W_dense(p, l, Ws, Wn);
        for(i = 0; i < nr; i++)
          for(j = 0; j < ns; j++){
            jj = rev ? ns - 1 - j : j;
            idx = (size_t)(s->off[b] + i) * s->N + s->off[r] + jj;
            op->Ks[idx] += w * Ws[i * ns + j];
            op->Cs[idx] += (1.0 - w) * Ws[i * ns + j];
            op->Kn[idx] += w * Wn[i * ns + j];
            op->Cn[idx] += (1.0 - w) * Wn[i * ns + j];
          }
      }
      free(Ws);
      free(Wn);
    }
  for(i = 0; i < s->N; i++){
    op->Ks[(size_t)i * s->N + i] += (w1 - 1.0) * s->etai[i];
    op->Cs[(size_t)i * s->N + i] += (1.0 - w1) * s->etai[i];
  }
  for(ii = 0; ii < NN; ii++)
    if(fabs(op->Ks[ii]) > amax)
      amax = fabs(op->Ks[ii]);
  thr = 1e-13 * amax;
  for(ii = 0; ii < NN; ii++){
    if(fabs(op->Ks[ii]) <= thr)
      op->Ks[ii] = 0.0;
    if(fabs(op->Kn[ii]) <= thr)
      op->Kn[ii] = 0.0;
    if(fabs(op->Cs[ii]) <= thr)
      op->Cs[ii] = 0.0;
    if(fabs(op->Cn[ii]) <= thr)
      op->Cn[ii] = 0.0;
  }
  bw = fd_band_width(s->N, op->Ks, 0.0);
  op->bw = bw;
  op->ab = (double *)malloc(sizeof(double) * (size_t)(3 * bw + 1) * s->N);
  fd_band_from_dense(s->N, bw, op->Ks, op->ab);
  if(s->Nu > 0 && surf_eliminate(s, op))
    return NULL;
  s->nops++;
  return op;
}

/* ---------------- slip weakening ---------------- */
static void dense_mult(int N, const double *A, const double *x, double *y)
{
  int i, j;
  double a;
  const double *r;
  for(i = 0; i < N; i++){
    a = 0.0;
    r = A + (size_t)i * N;
    for(j = 0; j < N; j++)
      a += r[j] * x[j];
    y[i] = a;
  }
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
  double r;
  int neg;
  int over;
  int cross;
  int uncross;
  /* normal stress with the guessed slip rate */
  if(normal_coupling){
    dense_mult(N, op->Kn, Vguess, tmp);
    for(i = 0; i < N; i++)
      sig[i] = sigma0[i] + s->nsign * (Hn[i] + tmp[i]);
  }else
    for(i = 0; i < N; i++)
      sig[i] = sigma0[i];
  for(i = 0; i < N; i++){
    r = delta[i] / f->d0;
    if(r > 1.0)
      r = 1.0;
    if(f->dugdale){
      mu_slip[i] = delta[i] < f->d0 ? f->mu_s[i] : f->mu_d[i];
      slope0[i] = 0.0;
    } else{
      mu_slip[i] = f->mu_s[i] - (f->mu_s[i] - f->mu_d[i]) * r;
      slope0[i] = delta[i] < f->d0 ? sig[i] * (f->mu_s[i] - f->mu_d[i]) * dt / f->d0 : 0.0;
    }
    if(f->force_time){
      fr = (t - f->force_time[i]) / (f->force_t0 > 1e-12 ? f->force_t0 : 1e-12);
      if(fr < 0.0)
        fr = 0.0;
      if(fr > 1.0)
        fr = 1.0;
      muf = f->mu_s[i] - (f->mu_s[i] - f->mu_d[i]) * fr;
      if(muf < mu_slip[i]){
        mu_slip[i] = muf;
        slope0[i] = 0.0;
      }
    }
    str0[i] = sig[i] * mu_slip[i];
  }
  /* initial active set */
  fd_band_mult(N, op->bw, op->ab, Vguess, tmp);
  for(i = 0; i < N; i++){
    active[i] = (H[i] + tmp[i] > str0[i]) || (Vguess[i] > 0.0);
    if(f->locked && f->locked[i])
      active[i] = 0;
  }
  for(i = 0; i < N; i++)
    V[i] = 0.0;
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
    if(fd_band_solve(N, bw, ab, rhs, Vn) != 0){
      fprintf(stderr, "fd_slipweak_solve: band solve failed\n");
      break;
    }
    for(i = 0; i < N; i++)
      if(!active[i])
        Vn[i] = 0.0;
    fd_band_mult(N, op->bw, op->ab, Vn, tmp);
    changed = 0;
    for(i = 0; i < N; i++){
      strength0 = sat[i] ? sig[i] * f->mu_d[i] : str0[i];
      tau[i] = H[i] + tmp[i];
      neg = active[i] && Vn[i] < 0.0;
      over = (!active[i]) && tau[i] > strength0 && !(f->locked && f->locked[i]);
      cross = active[i] && !sat[i] && delta[i] < f->d0 && delta[i] + Vn[i] * dt > f->d0;
      uncross = active[i] && sat[i] && delta[i] < f->d0 && delta[i] + Vn[i] * dt < f->d0;
      if(neg || over || cross || uncross)
        changed = 1;
      V[i] = Vn[i];
      if(neg)
        active[i] = 0;
      if(over)
        active[i] = 1;
      if(cross)
        sat[i] = 1;
      if(uncross)
        sat[i] = 0;
    }
    if(!changed)
      break;
  }
  for(i = 0; i < N; i++)
    if(V[i] < 0.0)
      V[i] = 0.0;
  if(normal_coupling){
    dense_mult(N, op->Kn, V, tmp);
    for(i = 0; i < N; i++)
      sigma[i] = sigma0[i] + s->nsign * (Hn[i] + tmp[i]);
  }else
    for(i = 0; i < N; i++)
      sigma[i] = sigma0[i];
  free(sig);
  free(mu_slip);
  free(slope0);
  free(str0);
  free(ab);
  free(rhs);
  free(Vn);
  free(tau);
  free(tmp);
  free(active);
  free(sat);
  return it;
}

/* ---------------- state and step ---------------- */
struct fd_state *fd_state_create(int N)
{
  struct fd_state *st = (struct fd_state *)calloc(1, sizeof(*st));
  st->N = N;
  st->k = 1;
  st->delta = calloc(N, sizeof(double));
  st->V = calloc(N, sizeof(double));
  st->tau = calloc(N, sizeof(double));
  st->sigma = calloc(N, sizeof(double));
  st->tau0 = calloc(N, sizeof(double));
  st->taudot = calloc(N, sizeof(double));
  st->sigma0 = calloc(N, sizeof(double));
  st->Hn_dot = calloc(N, sizeof(double));
  st->psi = calloc(N, sizeof(double));
  st->H = calloc(N, sizeof(double));
  st->Hn = calloc(N, sizeof(double));
  st->hs = calloc(N, sizeof(double));
  st->hn = calloc(N, sizeof(double));
  st->work = calloc(N, sizeof(double));
  return st;
}
void fd_state_free(struct fd_state *st)
{
  if(!st)
    return;
  free(st->delta);
  free(st->V);
  free(st->tau);
  free(st->sigma);
  free(st->tau0);
  free(st->taudot);
  free(st->sigma0);
  free(st->Hn_dot);
  free(st->psi);
  free(st->H);
  free(st->Hn);
  free(st->hs);
  free(st->hn);
  free(st->work);
  free(st);
}

/* the friction solves on the fault with the surface unknowns eliminated (M7 step 2b); without a
   surface they are the plain solves */
static void surf_solve_sw(struct fd_slipweak *f, struct fd_sys *s, struct fd_ops *op, const double *H, const double *Hn,
    const double *sigma0, int normal_coupling, double t, double dt, const double *delta, const double *Vguess, double *V, double *sigma)
{
  double *He, *Hne, *q;
  int N;
  if(s->Nu == 0){
    fd_slipweak_solve(f, s, op, H, Hn, sigma0, normal_coupling, t, dt, delta, Vguess, V, sigma);
    return;
  }
  He = (double *)malloc(sizeof(double)*(s->Nf + 1));
  Hne = (double *)malloc(sizeof(double)*(s->Nf + 1));
  q = (double *)malloc(sizeof(double)*s->Nu);
  surf_reduce(s, op, H, Hn, He, Hne, q);
  N = s->N;
  s->N = s->Nf;
  fd_slipweak_solve(f, s, op->el, He, Hne, sigma0, normal_coupling, t, dt, delta, Vguess, V, sigma);
  s->N = N;
  surf_recover(s, op, q, V);
  free(He);
  free(Hne);
  free(q);
}
static void surf_solve_rsf(const struct fd_rsf *f, struct fd_sys *s, struct fd_ops *op, const double *H, const double *sigma,
    const double *psi, const double *Vg, double *V)
{
  double *He, *q;
  int N;
  if(s->Nu == 0){
    fd_rsf_solve(f, s, op, H, sigma, psi, Vg, V);
    return;
  }
  He = (double *)malloc(sizeof(double)*(s->Nf + 1));
  q = (double *)malloc(sizeof(double)*s->Nu);
  surf_reduce(s, op, H, NULL, He, NULL, q);
  N = s->N;
  s->N = s->Nf;
  fd_rsf_solve(f, s, op->el, He, sigma, psi, Vg, V);
  s->N = N;
  surf_recover(s, op, q, V);
  free(He);
  free(q);
}

static int fd_step_slipweak_nl(struct fd_sys *s, struct fd_state *st, struct fd_slipweak *f, int k)
{
  int i, N = s->N;
  double dt = s->dt * k, t1 = st->t + dt, *Vold, *d1;
  struct fd_ops *op = fd_sys_ops_nl(s, k);
  if(!op)
    return -1;
  Vold = (double *)malloc(sizeof(double) * N);
  d1 = (double *)calloc((size_t)N, sizeof(double));      /* calloc: silences -Wmaybe-uninitialized */
  memcpy(Vold, st->V, sizeof(double) * N);
  fd_sys_static(s, st->delta, st->H, st->Hn);
  fd_sys_history(s, k + 1, st->hs, st->hn);
  dense_mult(N, op->Cs, Vold, st->work);
  for(i = 0; i < N; i++){
    st->H[i] += st->tau0[i] + st->taudot[i] * t1 + st->hs[i] + st->work[i];
    st->Hn[i] += st->Hn_dot[i] * t1 + st->hn[i];
    d1[i] = st->delta[i] + 0.5 * dt * Vold[i];     /* slip at t1 = d1 + (dt/2) V */
  }
  if(st->normal_coupling){                         /* normal stress: Hn + Cn V_old + Bn V */
    dense_mult(N, op->Cn, Vold, st->work);
    for(i = 0; i < N; i++)
      st->Hn[i] += st->work[i];
  }
  surf_solve_sw(f, s, op, st->H, st->Hn, st->sigma0, st->normal_coupling, t1, 0.5 * dt, d1, st->V, st->work, st->sigma);
  memcpy(st->V, st->work, sizeof(double) * N);
  fd_band_mult(N, op->bw, op->ab, st->V, st->work);
  for(i = 0; i < N; i++)
    st->tau[i] = st->H[i] + st->work[i];
  fd_sys_push_lin(s, Vold, st->V, k);
  for(i = 0; i < N; i++)
    st->delta[i] += 0.5 * (Vold[i] + st->V[i]) * dt;
  st->t = t1;
  st->n++;
  st->k = k;
  free(Vold);
  free(d1);
  return 0;
}

int fd_step_slipweak(struct fd_sys *s, struct fd_state *st, struct fd_slipweak *f, int k)
{
  int i, N = s->N;
  double dt = s->dt * k, t1 = st->t + dt;
  struct fd_ops *op;
  if(st->nl)
    return fd_step_slipweak_nl(s, st, f, k);
  op = fd_sys_ops(s, k);
  if(!op)
    return -1;
  /* known parts: static + loading + history beyond the current step */
  fd_sys_static(s, st->delta, st->H, st->Hn);
  fd_sys_history(s, k + 1, st->hs, st->hn);
  for(i = 0; i < N; i++){
    st->H[i] += st->tau0[i] + st->taudot[i] * t1 + st->hs[i];
    st->Hn[i] += st->Hn_dot[i] * t1 + st->hn[i];
  }
  surf_solve_sw(f, s, op, st->H, st->Hn, st->sigma0, st->normal_coupling, t1, dt, st->delta, st->V, st->work, st->sigma);
  memcpy(st->V, st->work, sizeof(double) * N);
  fd_band_mult(N, op->bw, op->ab, st->V, st->work);
  for(i = 0; i < N; i++)
    st->tau[i] = st->H[i] + st->work[i];
  fd_sys_push_k(s, st->V, k);
  for(i = 0; i < N; i++)
    st->delta[i] += st->V[i] * dt;
  st->t = t1;
  st->n++;
  st->k = k;
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
  if(V < 1e-30)
    V = 1e-30;
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
  if(Hloc <= 0.0)
    return 0.0;
  u = log(Vg > 1e-300 ? Vg : 1e-300);
  if(u < -700.0)
    u = -700.0;
  if(u > 50.0)
    u = 50.0;
  for(it = 0; it < 200; it++){
    V = exp(u);
    g = Hloc + k11 * V - fd_rsf_tau(f, i, V, psi, sigma);
    if(g > 0.0){
      if(u > lo)
        lo = u;
    } else {
      if(u < hi)
        hi = u;
    }
    dg = (k11 - fd_rsf_dtau(f, i, V, psi, sigma)) * V;
    un = (dg != 0.0 && isfinite(dg)) ? u - g / dg : 0.5 * (lo + hi);
    if(!isfinite(un) || un <= lo || un >= hi)
      un = 0.5 * (lo + hi);
    if(fabs(un - u) < f->tol){
      u = un;
      break;
    }
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
  double d;
  for(i = 0; i < N; i++){
    u[i] = log(Vg[i] > 1e-30 ? Vg[i] : 1e-30);
    if(u[i] < -70.0)
      u[i] = -70.0;
    if(u[i] > 10.0)
      u[i] = 10.0;
  }
  for(it = 0; it < f->max_iter; it++){
    for(i = 0; i < N; i++)
      V[i] = exp(u[i]);
    fd_band_mult(N, op->bw, op->ab, V, tmp);
    for(i = 0; i < N; i++)
      g[i] = H[i] + tmp[i] - fd_rsf_tau(f, i, V[i], psi[i], sigma[i]);
    /* Jacobian in u: (K - diag(dtau/dV)) diag(V), band storage: column j scaled by V_j */
    for(j = 0; j < N; j++)
      for(i = (j - bw > 0 ? j - bw : 0); i <= (j + bw < N - 1 ? j + bw : N - 1); i++)
        ab[(2 * bw + i - j) + (size_t)j * ldab] = op->ab[(2 * bw + i - j) + (size_t)j * ldab] * V[j];
    for(i = 0; i < N; i++)
      ab[(2 * bw) + (size_t)i * ldab] -= fd_rsf_dtau(f, i, V[i], psi[i], sigma[i]) * V[i];
    for(i = 0; i < N; i++){
      if(!isfinite(g[i]) || !isfinite(ab[(2 * bw) + (size_t)i * ldab]))
        break;
    }
    if(i < N)
      break;
    if(fd_band_solve(N, bw, ab, g, du) != 0)
      break;
    dmax = 0.0;
    for(i = 0; i < N; i++){
      d = -du[i];
      if(!isfinite(d)){
        dmax = INFINITY;
        break;
      }
      if(d > 2.0)
        d = 2.0;
      if(d < -2.0)
        d = -2.0;
      u[i] += d;
      if(u[i] < -70.0)
        u[i] = -70.0;
      if(u[i] > 10.0)
        u[i] = 10.0;
      if(fabs(d) > dmax)
        dmax = fabs(d);
    }
    if(!isfinite(dmax))
      break;
    if(dmax < f->tol){
      converged = 1;
      break;
    }
  }
  if(converged){
    for(i = 0; i < N; i++)
      V[i] = exp(u[i]);
    free(u);
    free(g);
    free(du);
    free(ab);
    free(tmp);
    free(k11);
    return it;
  }
  /* fallback: fixed point of the elementwise solve with lagged neighbours */
  for(i = 0; i < N; i++){
    k11[i] = op->Ks[(size_t)i * N + i];
    V[i] = exp(u[i]);
  }
  for(it = 0; it < 60; it++){
    fd_band_mult(N, op->bw, op->ab, V, tmp);
    dmax = 0.0;
    for(i = 0; i < N; i++){
      Vi = rsf_local(f, i, H[i] + tmp[i] - k11[i] * V[i], k11[i], sigma[i], psi[i], V[i]);
      if(fabs(Vi - V[i]) / (V[i] > 1e-30 ? V[i] : 1e-30) > dmax)
        dmax = fabs(Vi - V[i]) / (V[i] > 1e-30 ? V[i] : 1e-30);
      g[i] = Vi;
    }
    memcpy(V, g, sizeof(double) * N);
    if(dmax < f->tol)
      break;
  }
  free(u);
  free(g);
  free(du);
  free(ab);
  free(tmp);
  free(k11);
  return -(it + 1);
}

/* one step of k dt_min: predictor-corrector in psi (exact aging update at the
   mean slip rate), coupled solve for V twice (stepper.RateState.solve) */
static int fd_step_rsf_nl(struct fd_sys *s, struct fd_state *st, const struct fd_rsf *f, int k)
{
  int i, N = s->N;
  double dt = s->dt * k, t1 = st->t + dt, *psi_s = st->work, *Vs, *sig, *Vold;
  struct fd_ops *op = fd_sys_ops_nl(s, k);
  double tp0 = fd_now(), tp1;
  if(!op)
    return -1;
  Vs = (double *)calloc((size_t)N, sizeof(double));
  sig = (double *)calloc((size_t)N, sizeof(double));
  Vold = (double *)malloc(sizeof(double) * N);
  memcpy(Vold, st->V, sizeof(double) * N);
  for(i = 0; i < N; i++)
    sig[i] = st->sigma0[i];
  fd_sys_static(s, st->delta, st->H, st->Hn);
  tp1 = fd_now();
  fd_prof[1] += tp1 - tp0;
  tp0 = tp1;
  fd_sys_history(s, k + 1, st->hs, st->hn);
  tp1 = fd_now();
  fd_prof[0] += tp1 - tp0;
  tp0 = tp1;
  dense_mult(N, op->Cs, Vold, Vs);
  for(i = 0; i < N; i++){
    st->H[i] += st->tau0[i] + st->taudot[i] * t1 + st->hs[i] + Vs[i];
    st->Hn[i] += st->Hn_dot[i] * t1 + st->hn[i];
  }
  if(st->normal_coupling){                         /* normal stress: Hn + Cn V_old + Bn V */
    dense_mult(N, op->Cn, Vold, Vs);
    for(i = 0; i < N; i++)
      st->Hn[i] += Vs[i];
    dense_mult(N, op->Kn, Vold, sig);
    for(i = 0; i < N; i++)
      sig[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]);
  }
  /* predictor, corrector: state over the step at the mean rate */
  for(i = 0; i < s->Nf; i++)
    psi_s[i] = fd_rsf_psi_advance(f, i, st->psi[i], Vold[i], Vold[i], dt);
  surf_solve_rsf(f, s, op, st->H, sig, psi_s, Vold, Vs);
  if(st->normal_coupling){
    dense_mult(N, op->Kn, Vs, sig);
    for(i = 0; i < N; i++)
      sig[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]);
  }
  for(i = 0; i < s->Nf; i++)
    psi_s[i] = fd_rsf_psi_advance(f, i, st->psi[i], Vold[i], Vs[i], dt);
  surf_solve_rsf(f, s, op, st->H, sig, psi_s, Vs, st->V);
  tp1 = fd_now();
  fd_prof[2] += tp1 - tp0;
  tp0 = tp1;
  memcpy(st->psi, psi_s, sizeof(double) * N);
  fd_band_mult(N, op->bw, op->ab, st->V, Vs);
  for(i = 0; i < N; i++)
    st->tau[i] = st->H[i] + Vs[i];
  if(st->normal_coupling){
    dense_mult(N, op->Kn, st->V, sig);
    for(i = 0; i < N; i++)
      st->sigma[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]);
  } else
    for(i = 0; i < N; i++)
      st->sigma[i] = st->sigma0[i];
  fd_sys_push_lin(s, Vold, st->V, k);
  tp1 = fd_now();
  fd_prof[3] += tp1 - tp0;
  tp0 = tp1;
  for(i = 0; i < N; i++)
    st->delta[i] += 0.5 * (Vold[i] + st->V[i]) * dt;
  st->t = t1;
  st->n++;
  st->k = k;
  free(Vs);
  free(sig);
  free(Vold);
  return 0;
}

int fd_step_rsf(struct fd_sys *s, struct fd_state *st, const struct fd_rsf *f, int k)
{
  int i, N = s->N;
  double dt = s->dt * k, t1 = st->t + dt, *psi_s = st->work, *Vs, *sig;
  struct fd_ops *op;
  if(st->nl)
    return fd_step_rsf_nl(s, st, f, k);
  op = fd_sys_ops(s, k);
  if(!op)
    return -1;
  double tp0 = fd_now(), tp1;
  Vs = (double *)calloc((size_t)N, sizeof(double));
  sig = (double *)calloc((size_t)N, sizeof(double));
  for(i = 0; i < N; i++)
    sig[i] = st->sigma0[i];
  fd_sys_static(s, st->delta, st->H, st->Hn);
  tp1 = fd_now();
  fd_prof[1] += tp1 - tp0;
  tp0 = tp1;
  fd_sys_history(s, k + 1, st->hs, st->hn);
  tp1 = fd_now();
  fd_prof[0] += tp1 - tp0;
  tp0 = tp1;
  for(i = 0; i < N; i++){
    st->H[i] += st->tau0[i] + st->taudot[i] * t1 + st->hs[i];
    st->Hn[i] += st->Hn_dot[i] * t1 + st->hn[i];
  }
  /* predictor */
  if(st->normal_coupling){
    dense_mult(N, op->Kn, st->V, sig);
    for(i = 0; i < N; i++)
      sig[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]);
  } else
    for(i = 0; i < N; i++)
      sig[i] = st->sigma0[i];
  for(i = 0; i < s->Nf; i++)
    psi_s[i] = fd_rsf_psi_advance(f, i, st->psi[i], st->V[i], st->V[i], dt);
  surf_solve_rsf(f, s, op, st->H, sig, psi_s, st->V, Vs);
  /* corrector */
  if(st->normal_coupling){
    dense_mult(N, op->Kn, Vs, sig);
    for(i = 0; i < N; i++)
      sig[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]);
  }
  for(i = 0; i < s->Nf; i++)
    psi_s[i] = fd_rsf_psi_advance(f, i, st->psi[i], st->V[i], Vs[i], dt);
  surf_solve_rsf(f, s, op, st->H, sig, psi_s, Vs, st->V);
  tp1 = fd_now();
  fd_prof[2] += tp1 - tp0;
  tp0 = tp1;
  memcpy(st->psi, psi_s, sizeof(double) * N);
  fd_band_mult(N, op->bw, op->ab, st->V, Vs);
  for(i = 0; i < N; i++)
    st->tau[i] = st->H[i] + Vs[i];
  if(st->normal_coupling){
    dense_mult(N, op->Kn, st->V, sig);
    for(i = 0; i < N; i++)
      st->sigma[i] = st->sigma0[i] + s->nsign * (st->Hn[i] + sig[i]);
  } else
    for(i = 0; i < N; i++)
      st->sigma[i] = st->sigma0[i];
  tp1 = fd_now();
  fd_prof[5] += tp1 - tp0;
  tp0 = tp1;
  fd_sys_push_k(s, st->V, k);
  tp1 = fd_now();
  fd_prof[3] += tp1 - tp0;
  tp0 = tp1;
  for(i = 0; i < N; i++)
    st->delta[i] += st->V[i] * dt;
  st->t = t1;
  st->n++;
  st->k = k;
  free(Vs);
  free(sig);
  return 0;
}
