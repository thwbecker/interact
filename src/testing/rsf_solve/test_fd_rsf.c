/* test of fd_rsf_solve (fd_march.c): the coupled rate-and-state solve of the dynamic phase with
   slip in both directions (fd_plan.md, 2026-11-18). A banded system of N elements, K with
   negative diagonal (the instantaneous response) and smaller neighbours, H mixing positive and
   negative stresses. Checks: the residual H + K V - tau_f(V) is small; each V has the sign of the
   element's stress at rest; V(-H) = -V(H) (the friction law is odd); for H > 0 everywhere every
   V > 0. TFR_MAXIT=1 forces the elementwise fallback. Prints PASS or FAIL. usage: test_fd_rsf [N] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd_march.h"

static void band_from_dense(int N, int bw, const double *K, double *ab)
{
  int i, j, ldab;
  ldab = 3 * bw + 1;
  memset(ab, 0, sizeof(double) * (size_t)ldab * N);
  for(j = 0; j < N; j++)
    for(i = (j - bw > 0 ? j - bw : 0); i <= (j + bw < N - 1 ? j + bw : N - 1); i++)
      ab[(2 * bw + i - j) + (size_t)j * ldab] = K[(size_t)i * N + j];
}

static double check(const struct fd_rsf *f, int N, const double *K, const double *H, const double *sig,
                    const double *psi, const double *V, int *badsign)
{
  int i, j;
  double r, rest, worst;
  worst = 0.0;
  *badsign = 0;
  for(i = 0; i < N; i++){
    r = H[i];
    rest = H[i];
    for(j = 0; j < N; j++){
      r += K[(size_t)i * N + j] * V[j];
      if(j != i)
        rest += K[(size_t)i * N + j] * V[j];
    }
    r -= fd_rsf_tau(f, i, V[i], psi[i], sig[i]);
    if(fabs(r) / (fabs(H[i]) + 1.0) > worst)
      worst = fabs(r) / (fabs(H[i]) + 1.0);
    if((rest > 0.0 && V[i] < 0.0) || (rest < 0.0 && V[i] > 0.0))
      (*badsign)++;
  }
  return worst;
}

int main(int argc, char **argv)
{
  int N, bw, i, j, it1, it2, it3, bad1, bad2, bad3, fail;
  double *K, *H, *Hm, *sig, *psi, *Vg, *V, *Vm, *Vp, *a, *b, *dc, r1, r2, r3, sym;
  struct fd_rsf f;
  struct fd_sys s;
  struct fd_ops op;
  N = argc > 1 ? atoi(argv[1]) : 12;
  bw = 2;
  K = calloc((size_t)N * N, sizeof(double));
  H = malloc(sizeof(double) * N);
  Hm = malloc(sizeof(double) * N);
  sig = malloc(sizeof(double) * N);
  psi = malloc(sizeof(double) * N);
  Vg = malloc(sizeof(double) * N);
  V = malloc(sizeof(double) * N);
  Vm = malloc(sizeof(double) * N);
  Vp = malloc(sizeof(double) * N);
  a = malloc(sizeof(double) * N);
  b = malloc(sizeof(double) * N);
  dc = malloc(sizeof(double) * N);
  for(i = 0; i < N; i++){
    a[i] = 0.015;
    b[i] = 0.02;
    dc[i] = 0.01;
    sig[i] = 50e6;
    psi[i] = 0.6 + 0.02 * log(1e-6 / 1e-3);      /* steady state near 1 mm/s */
    Vg[i] = 1e-9;
    for(j = (i - bw > 0 ? i - bw : 0); j <= (i + bw < N - 1 ? i + bw : N - 1); j++)
      K[(size_t)i * N + j] = (i == j) ? -2.0e7 : 3.0e6 / (1 + abs(i - j));
  }
  f.N = N;
  f.a = a;
  f.b = b;
  f.dc = dc;
  f.v0 = 1e-6;
  f.f0 = 0.6;
  f.tol = 1e-12;
  f.max_iter = (getenv("TFR_MAXIT") ? atoi(getenv("TFR_MAXIT")) : 50);
  /* stresses about the strength at 1 mm/s, every third element reversed */
  for(i = 0; i < N; i++)
    H[i] = ((i % 3) == 2) ? -1.05 * fd_rsf_tau(&f, i, 1e-3, psi[i], sig[i]) : 1.02 * fd_rsf_tau(&f, i, 1e-3, psi[i], sig[i]);
  memset(&s, 0, sizeof(s));
  s.N = N;
  memset(&op, 0, sizeof(op));
  op.bw = bw;
  op.Ks = K;
  op.ab = malloc(sizeof(double) * (size_t)(3 * bw + 1) * N);
  band_from_dense(N, bw, K, op.ab);
  /* mixed signs */
  it1 = fd_rsf_solve(&f, &s, &op, H, sig, psi, Vg, V);
  r1 = check(&f, N, K, H, sig, psi, V, &bad1);
  /* odd symmetry */
  for(i = 0; i < N; i++)
    Hm[i] = -H[i];
  it2 = fd_rsf_solve(&f, &s, &op, Hm, sig, psi, Vg, Vm);
  r2 = check(&f, N, K, Hm, sig, psi, Vm, &bad2);
  sym = 0.0;
  for(i = 0; i < N; i++)
    if(fabs(Vm[i] + V[i]) / (fabs(V[i]) + 1e-30) > sym)
      sym = fabs(Vm[i] + V[i]) / (fabs(V[i]) + 1e-30);
  /* all positive */
  for(i = 0; i < N; i++)
    H[i] = fabs(H[i]);
  it3 = fd_rsf_solve(&f, &s, &op, H, sig, psi, Vg, Vp);
  r3 = check(&f, N, K, H, sig, psi, Vp, &bad3);
  for(i = 0; i < N; i++)
    if(Vp[i] <= 0.0)
      bad3++;
  printf("test_fd_rsf: N %d: mixed signs: %d iterations, residual %.1e, %d wrong signs; V:", N, it1, r1, bad1);
  for(i = 0; i < N && i < 6; i++)
    printf(" %.3e", V[i]);
  printf("\n");
  printf("test_fd_rsf: reversed stresses: %d iterations, residual %.1e, %d wrong signs, |V(-H) + V(H)| / |V| %.1e\n", it2, r2, bad2, sym);
  printf("test_fd_rsf: all positive: %d iterations, residual %.1e, %d wrong signs\n", it3, r3, bad3);
  fail = (r1 > 1e-6) || (r2 > 1e-6) || (r3 > 1e-6) || bad1 || bad2 || bad3 || (sym > 1e-8);
  if(!getenv("TFR_MAXIT") && ((it1 < 0) || (it2 < 0) || (it3 < 0)))
    fail = 1;                                    /* Newton expected to converge */
  printf("test_fd_rsf: %s\n", fail ? "FAIL" : "PASS");
  return fail;
}
