/* see fd_hpdwn.h */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "fd_hpdwn.h"
#include "fd_fft.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* complex square root with non-negative real part */
static double complex csqrt_pos(double complex z)
{
  double complex r;
  r = csqrt(z);
  if(creal(r) < 0.0)
    r = -r;
  return r;
}

/* traction (txz, tzz) at z = 0 of the upgoing incident field of unit slip on source j, per unit
   of the common factor exp(i k (x - xs)), wavenumber k, frequency w */
static void incident_traction(double k, double complex w, double alpha, double beta, double rho,
                              double tx, double tz, double a, double zs, double complex *txz, double complex *tzz)
{
  double lam, mu, nx, nz, m[2][2];
  double complex kb2, nu, Dr[2], Ds[2], g, Gip, u[2], acc, exx, ezz, sxz, szz, cc, seg, ph;
  int c, i, p, q;
  lam = rho*(alpha*alpha - 2.0*beta*beta);
  mu = rho*beta*beta;
  nx = tz;                                       /* normal as in halfplane_dwn.py */
  nz = -tx;
  m[0][0] = 2.0*tx*nx*mu;
  m[0][1] = (tx*nz + tz*nx)*mu;
  m[1][0] = m[0][1];
  m[1][1] = 2.0*tz*nz*mu;
  kb2 = (w/beta)*(w/beta);
  *txz = 0.0;
  *tzz = 0.0;
  for(c = 0; c < 2; c++){                        /* 0: P (alpha), 1: S (beta) */
    nu = csqrt_pos(k*k - (c ? kb2 : (w/alpha)*(w/alpha)));
    Dr[0] = I*k;
    Dr[1] = nu;                                   /* upgoing side: Zsign = -1, d/dz -> +nu */
    Ds[0] = -I*k;
    Ds[1] = -nu;
    g = 1.0/(2.0*nu);
    for(i = 0; i < 2; i++){
      acc = 0.0;
      for(p = 0; p < 2; p++){
        if(c)
          Gip = (kb2*(i == p ? 1.0 : 0.0) + Dr[i]*Dr[p])*g/(rho*w*w);
        else
          Gip = -(Dr[i]*Dr[p])*g/(rho*w*w);
        for(q = 0; q < 2; q++)
          acc += m[p][q]*Ds[q]*Gip;
      }
      u[i] = acc;
    }
    exx = Dr[0]*u[0];
    ezz = Dr[1]*u[1];
    szz = lam*exx + (lam + 2.0*mu)*ezz;
    sxz = mu*(Dr[1]*u[0] + Dr[0]*u[1]);
    cc = -I*k*tx - nu*tz;
    if(cabs(cc*a) < 1e-8)
      seg = 2.0*a;
    else
      seg = 2.0*csinh(cc*a)/cc;
    ph = cexp(-nu*zs);
    *txz += sxz*seg*ph;
    *tzz += szz*seg*ph;
  }
}

int fd_hpdwn_reflected(int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                       int nr, const double *rx, const double *rz, int nt, const double *t,
                       double alpha, double beta, double rho, double T, int nw, double wIT, double L, int M,
                       double *oxx, double *oxz, double *ozz)
{
  double lam, mu, wI, dw, k, fac;
  double complex w, na, nb, PxzA, PzzA, PxzB, PzzB, det, txz, tzz, A, B, ea, eb, uxA, uzA, uxB, uzB, ph, Bw, e;
  double complex *acc;
  double complex S[3];
  int jw, m, i, j, l, c;
  size_t np, o;
  lam = rho*(alpha*alpha - 2.0*beta*beta);
  mu = rho*beta*beta;
  wI = wIT/T;
  dw = 2.0*M_PI/T;
  np = (size_t)nr*ns;
  acc = (double complex *)calloc(3*np*nt, sizeof(double complex));
  if(!acc)
    return 1;
  for(jw = 0; jw <= nw/2; jw++){
    w = jw*dw + I*wI;
    Bw = -1.0/(w*w);                               /* the ramp */
    fac = (jw == 0) ? 1.0 : 2.0;
    for(m = -M; m <= M; m++){
      k = 2.0*M_PI*m/L;
      na = csqrt_pos(k*k - (w/alpha)*(w/alpha));
      nb = csqrt_pos(k*k - (w/beta)*(w/beta));
      PxzA = mu*(I*k*(-na) + (-na)*I*k);
      PzzA = lam*(I*k*I*k + na*na) + 2.0*mu*na*na;
      PxzB = mu*(I*k*I*k - nb*nb);
      PzzB = lam*(I*k*nb - nb*I*k) + 2.0*mu*(-nb*I*k);
      det = PxzA*PzzB - PxzB*PzzA;
      for(j = 0; j < ns; j++){
        incident_traction(k, w, alpha, beta, rho, stx[j], stz[j], sa[j], sz[j], &txz, &tzz);
        A = (-txz*PzzB + tzz*PxzB)/det;
        B = (-tzz*PxzA + txz*PzzA)/det;
        uxA = I*k*A;
        uzA = -na*A;
        uxB = nb*B;
        uzB = I*k*B;
        for(i = 0; i < nr; i++){
          ea = cexp(-na*rz[i]);
          eb = cexp(-nb*rz[i]);
          ph = cexp(I*k*(rx[i] - sx[j]))*Bw*fac/L;
          o = ((size_t)i*ns + j)*3;
          /* sxx, sxz, szz at the receiver, then their time samples */
            S[0] = (((lam + 2.0*mu)*(I*k*uxA) + lam*(-na*uzA))*ea + ((lam + 2.0*mu)*(I*k*uxB) + lam*(-nb*uzB))*eb)*ph;
            S[1] = (mu*((-na*uxA + I*k*uzA)*ea + (-nb*uxB + I*k*uzB)*eb))*ph;
            S[2] = ((lam*(I*k*uxA) + (lam + 2.0*mu)*(-na*uzA))*ea + (lam*(I*k*uxB) + (lam + 2.0*mu)*(-nb*uzB))*eb)*ph;
            for(l = 0; l < nt; l++){
              e = cexp(-I*(jw*dw)*t[l]);
              for(c = 0; c < 3; c++)
                acc[(o + c)*nt + l] += S[c]*e;
            }
        }
      }
    }
  }
  for(i = 0; i < nr; i++)
    for(j = 0; j < ns; j++)
      for(l = 0; l < nt; l++){
        o = ((size_t)i*ns + j);
        fac = exp(wI*t[l])*dw/(2.0*M_PI);
        oxx[o*nt + l] = creal(acc[(o*3 + 0)*nt + l])*fac;
        oxz[o*nt + l] = creal(acc[(o*3 + 1)*nt + l])*fac;
        ozz[o*nt + l] = creal(acc[(o*3 + 2)*nt + l])*fac;
      }
  free(acc);
  return 0;
}

extern void zgemm_(char *ta, char *tb, int *m, int *n, int *k, double complex *alpha, double complex *a, int *lda,
                   double complex *b, int *ldb, double complex *beta, double complex *c, int *ldc, size_t lta, size_t ltb);

/* the per-wavenumber source amplitudes A (P) and B (S) of the reflected waves for source j,
   including exp(-i k xs), the ramp transform and 1 / L */
static void source_amplitudes(double k, double complex w, double alpha, double beta, double rho, double L,
                              double sxj, double szj, double tx, double tz, double a, double complex *A, double complex *B)
{
  double lam, mu;
  double complex na, nb, PxzA, PzzA, PxzB, PzzB, det, txz, tzz, f;
  lam = rho*(alpha*alpha - 2.0*beta*beta);
  mu = rho*beta*beta;
  na = csqrt_pos(k*k - (w/alpha)*(w/alpha));
  nb = csqrt_pos(k*k - (w/beta)*(w/beta));
  PxzA = mu*(I*k*(-na) + (-na)*I*k);
  PzzA = lam*(I*k*I*k + na*na) + 2.0*mu*na*na;
  PxzB = mu*(I*k*I*k - nb*nb);
  PzzB = lam*(I*k*nb - nb*I*k) + 2.0*mu*(-nb*I*k);
  det = PxzA*PzzB - PxzB*PzzA;
  incident_traction(k, w, alpha, beta, rho, tx, tz, a, szj, &txz, &tzz);
  f = cexp(-I*k*sxj)*(-1.0/(w*w))/L;
  *A = f*(-txz*PzzB + tzz*PxzB)/det;
  *B = f*(-tzz*PxzA + txz*PzzA)/det;
}

/* the projected tractions (shear, normal) at receiver i per unit A (P) and per unit B (S),
   including exp(i k xr - nu zr) */
static void receiver_factors(double k, double complex w, double alpha, double beta, double rho,
                             double xr, double zr, double tx, double tz,
                             double complex *sA, double complex *nA, double complex *sB, double complex *nB)
{
  double lam, mu, srx, sry, nrx, nry;
  double complex na, nb, ea, eb, uxA, uzA, uxB, uzB, sxx, sxz, szz, sxy, syy, trx, try_;
  int c;
  lam = rho*(alpha*alpha - 2.0*beta*beta);
  mu = rho*beta*beta;
  na = csqrt_pos(k*k - (w/alpha)*(w/alpha));
  nb = csqrt_pos(k*k - (w/beta)*(w/beta));
  ea = cexp(I*k*xr - na*zr);
  eb = cexp(I*k*xr - nb*zr);
  srx = tx;                                       /* y-up frame */
  sry = -tz;
  nrx = -sry;
  nry = srx;
  for(c = 0; c < 2; c++){
    if(c == 0){
      uxA = I*k;
      uzA = -na;
      sxx = ((lam + 2.0*mu)*(I*k*uxA) + lam*(-na*uzA))*ea;
      sxz = mu*(-na*uxA + I*k*uzA)*ea;
      szz = (lam*(I*k*uxA) + (lam + 2.0*mu)*(-na*uzA))*ea;
    }else{
      uxB = nb;
      uzB = I*k;
      sxx = ((lam + 2.0*mu)*(I*k*uxB) + lam*(-nb*uzB))*eb;
      sxz = mu*(-nb*uxB + I*k*uzB)*eb;
      szz = (lam*(I*k*uxB) + (lam + 2.0*mu)*(-nb*uzB))*eb;
    }
    sxy = -sxz;
    syy = szz;
    trx = sxx*nrx + sxy*nry;
    try_ = sxy*nrx + syy*nry;
    if(c == 0){
      *sA = trx*srx + try_*sry;
      *nA = trx*nrx + try_*nry;
    }else{
      *sB = trx*srx + try_*sry;
      *nB = trx*nrx + try_*nry;
    }
  }
}

/* the reflected tractions of all pairs at one complex frequency: Gs, Gn (nr x ns, column-major) */
static void spectra_at(double complex w, int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                       int nr, const double *rx, const double *rz, const double *rtx, const double *rtz,
                       double alpha, double beta, double rho, double L, int M,
                       double complex *R, double complex *S, double complex *Gs, double complex *Gn)
{
  int nk, m, i, j, two_nk;
  double k;
  double complex one, zero;
  double complex *G;
  char tn;
  nk = 2*M + 1;
  two_nk = 2*nk;
  /* R: receivers x (A wavenumbers, B wavenumbers), shear then normal (2 nr rows); S: (2 nk) x ns */
  for(m = 0; m < nk; m++){
    k = 2.0*M_PI*(m - M)/L;
    for(j = 0; j < ns; j++)
      source_amplitudes(k, w, alpha, beta, rho, L, sx[j], sz[j], stx[j], stz[j], sa[j],
                        &S[(size_t)j*two_nk + m], &S[(size_t)j*two_nk + nk + m]);
    for(i = 0; i < nr; i++)
      receiver_factors(k, w, alpha, beta, rho, rx[i], rz[i], rtx[i], rtz[i],
                       &R[(size_t)m*2*nr + i], &R[(size_t)m*2*nr + nr + i],
                       &R[(size_t)(nk + m)*2*nr + i], &R[(size_t)(nk + m)*2*nr + nr + i]);
  }
  one = 1.0;
  zero = 0.0;
  tn = 'N';
  i = 2*nr;
  /* (2 nr x 2 nk) (2 nk x ns) = 2 nr x ns: shear rows then normal rows */
  G = (double complex *)malloc(sizeof(double complex)*(size_t)2*nr*ns);
  zgemm_(&tn, &tn, &i, &ns, &two_nk, &one, R, &i, S, &two_nk, &zero, G, &i, 1, 1);
  for(j = 0; j < ns; j++)
    for(m = 0; m < nr; m++){
      Gs[(size_t)j*nr + m] = G[(size_t)j*2*nr + m];
      Gn[(size_t)j*nr + m] = G[(size_t)j*2*nr + nr + m];
    }
  free(G);
}

/* one block of receivers: the spectra of its pairs at all frequencies, then the lags by FFT */
static int batch_block(int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                       int nr, const double *rx, const double *rz, const double *rtx, const double *rtz,
                       double alpha, double beta, double rho, double dt, int nt, int q, int N, double T, double L, int M,
                       double zsum, double *os, double *on)
{
  int nw, jw, i, j, l, Mw;
  double wI, dw, fac;
  double complex w, *R, *S, *Gs, *Gn, *specs, *specn, *x, *y;
  fd_fft_plan plan;
  size_t np, p;
  wI = (getenv("FD_HPDWN_WIT") ? atof(getenv("FD_HPDWN_WIT")) : 18.0)/T;
  dw = 2.0*M_PI/T;
  nw = N/2 + 1;
  np = (size_t)nr*ns;
  R = (double complex *)malloc(sizeof(double complex)*(size_t)2*nr*2*(2*M + 1));
  S = (double complex *)malloc(sizeof(double complex)*(size_t)2*(2*M + 1)*ns);
  Gs = (double complex *)malloc(sizeof(double complex)*np);
  Gn = (double complex *)malloc(sizeof(double complex)*np);
  specs = (double complex *)malloc(sizeof(double complex)*np*nw);
  specn = (double complex *)malloc(sizeof(double complex)*np*nw);
  x = (double complex *)malloc(sizeof(double complex)*N);
  y = (double complex *)malloc(sizeof(double complex)*N);
  if(!R || !S || !Gs || !Gn || !specs || !specn || !x || !y)
    return 1;
  for(jw = 0; jw < nw; jw++){
    w = jw*dw + I*wI;
    /* wavenumbers to the S-wave number of this frequency plus the evanescent decay over the
       shortest image path (exp(-k zsum) below 1e-17), not beyond the M of the block */
    Mw = (int)ceil((1.2*jw*dw/beta + 40.0/zsum)*L/(2.0*M_PI));
    if(Mw > M)
      Mw = M;
    spectra_at(w, ns, sx, sz, stx, stz, sa, nr, rx, rz, rtx, rtz, alpha, beta, rho, L, Mw, R, S, Gs, Gn);
    fac = (jw == 0) ? 1.0 : 2.0;
    for(p = 0; p < np; p++){
      specs[p*nw + jw] = fac*Gs[p];
      specn[p*nw + jw] = fac*Gn[p];
    }
  }
  /* the lags: Re of sum_j X_j exp(-2 pi i j m / N) at m = l q, times exp(omega_I t) dw / (2 pi) */
  plan = fd_fft_plan_create(N, FD_FFT_FORWARD);
  for(p = 0; p < np; p++){
    i = (int)(p % nr);
    j = (int)(p / nr);
    for(l = 0; l < N; l++)
      x[l] = (l < nw) ? specs[p*nw + l] : 0.0;
    fd_fft_exec(plan, x, y);
    for(l = 0; l < nt; l++)
      os[((size_t)i*ns + j)*nt + l] = creal(y[l*q])*exp(wI*l*dt)*dw/(2.0*M_PI);
    for(l = 0; l < N; l++)
      x[l] = (l < nw) ? specn[p*nw + l] : 0.0;
    fd_fft_exec(plan, x, y);
    for(l = 0; l < nt; l++)
      on[((size_t)i*ns + j)*nt + l] = creal(y[l*q])*exp(wI*l*dt)*dw/(2.0*M_PI);
  }
  fd_fft_plan_destroy(plan);
  free(R);
  free(S);
  free(Gs);
  free(Gn);
  free(specs);
  free(specn);
  free(x);
  free(y);
  return 0;
}

int fd_hpdwn_batch(int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                   int nr, const double *rx, const double *rz, const double *rtx, const double *rtz,
                   double alpha, double beta, double rho, double dt, int nt, double *os, double *on, double *ks, double *kn)
{
  int N, M, i, j, q, nb, i0, n;
  double T, L, xmin, xmax, h, zsum, amin, mem;
  double complex w, *R, *S, *Gs, *Gn;
  /* the time grid is the lag grid refined q times (the reflected fronts need frequencies well
     above the Nyquist frequency of the lags; FD_HPDWN_Q, default 4: errors at the fronts of
     about 1.4 percent of the largest value, 0.7 at 8, 0.4 at 16, the cost growing with q^2) */
  q = getenv("FD_HPDWN_Q") ? atoi(getenv("FD_HPDWN_Q")) : 4;
  if(q < 1)
    q = 1;
  /* the window: twice the lags. The complex frequency damps the wrap-around by exp(-omega_I T)
     but amplifies the truncation errors of the damped signal by exp(omega_I t) at the last lag
     (exp(omega_I T / 2)); FD_HPDWN_WIT sets omega_I T (default 18, which passed the checks;
     a window of 1.25 times the lags with it failed) */
  N = fd_fft_size_for((int)ceil(2.0*nt*q));
  T = N*dt/q;
  xmin = 1e300;
  xmax = -1e300;
  amin = 1e300;
  zsum = 1e300;
  for(j = 0; j < ns; j++){
    xmin = fmin(xmin, sx[j] - sa[j]);
    xmax = fmax(xmax, sx[j] + sa[j]);
    amin = fmin(amin, sa[j]);
  }
  for(i = 0; i < nr; i++){
    xmin = fmin(xmin, rx[i]);
    xmax = fmax(xmax, rx[i]);
  }
  for(j = 0; j < ns; j++)
    for(i = 0; i < nr; i++)
      zsum = fmin(zsum, fabs(sz[j]) + fabs(rz[i]));
  L = alpha*T + 2.0*(xmax - xmin);
  h = 0.25*fmin(amin, zsum);
  M = (int)ceil(L/h);
  /* receivers in blocks, the spectra of a block within FD_HPDWN_MEM bytes (default 2e9) */
  mem = getenv("FD_HPDWN_MEM") ? atof(getenv("FD_HPDWN_MEM")) : 2e9;
  nb = (int)(mem/(2.0*16.0*ns*(N/2 + 1)));
  if(nb < 1)
    nb = 1;
  for(i0 = 0; i0 < nr; i0 += nb){
    n = (nr - i0 < nb) ? nr - i0 : nb;
    if(batch_block(ns, sx, sz, stx, stz, sa, n, rx + i0, rz + i0, rtx + i0, rtz + i0, alpha, beta, rho, dt, nt, q, N, T, L, M,
                   zsum, os + (size_t)i0*ns*nt, on + (size_t)i0*ns*nt))
      return 1;
  }
  /* static: the transfer function at a small frequency (period of 4,000 km), per unit slip */
  if(ks || kn){
    L = fmax(4e6, 100.0*(xmax - xmin));
    w = 2.0*M_PI*1e-4 + I*1e-5;
    M = (int)ceil(L/fmin(h, 100.0));
    if(M > 400000)
      M = 400000;
    R = (double complex *)malloc(sizeof(double complex)*(size_t)2*nr*2*(2*M + 1));
    S = (double complex *)malloc(sizeof(double complex)*(size_t)2*(2*M + 1)*ns);
    Gs = (double complex *)malloc(sizeof(double complex)*(size_t)nr*ns);
    Gn = (double complex *)malloc(sizeof(double complex)*(size_t)nr*ns);
    if(!R || !S || !Gs || !Gn)
      return 1;
    spectra_at(w, ns, sx, sz, stx, stz, sa, nr, rx, rz, rtx, rtz, alpha, beta, rho, L, M, R, S, Gs, Gn);
    for(i = 0; i < nr; i++)
      for(j = 0; j < ns; j++){
        /* spectra_at includes the ramp transform -1 / w^2 */
        if(ks)
          ks[(size_t)i*ns + j] = creal(Gs[(size_t)j*nr + i]*(-(w*w)));
        if(kn)
          kn[(size_t)i*ns + j] = creal(Gn[(size_t)j*nr + i]*(-(w*w)));
      }
    free(R);
    free(S);
    free(Gs);
    free(Gn);
  }
  return 0;
}
