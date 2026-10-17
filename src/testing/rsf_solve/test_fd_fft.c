/*
  test of fd_fft against a direct DFT and the round trip; exit code 0
  on success. build: make rsf_testing; run: bin/test_fd_fft
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "fd_fft.h"

static double rel_err(const double _Complex *a, const double _Complex *b, int n)
{
  double num = 0, den = 0; int i;
  for(i = 0; i < n; i++){ num += cabs(a[i] - b[i]) * cabs(a[i] - b[i]); den += cabs(b[i]) * cabs(b[i]); }
  return sqrt(num / (den > 0 ? den : 1));
}

int main(void)
{
  int sizes[] = {1, 2, 8, 64, 1024, 4096}, ns = 6, k, i, j, n, fail = 0;
  double _Complex *x, *f, *g, *d, s;
  fd_fft_plan pf, pb;
  printf("fd_fft backend: %s\n", fd_fft_backend());
  for(k = 0; k < ns; k++){
    n = fd_fft_size_for(sizes[k]);
    x = malloc(sizeof(*x) * n); f = malloc(sizeof(*f) * n); g = malloc(sizeof(*g) * n); d = malloc(sizeof(*d) * n);
    srand(12345 + n);
    for(i = 0; i < n; i++) x[i] = (rand() / (double)RAND_MAX - 0.5) + I * (rand() / (double)RAND_MAX - 0.5);
    pf = fd_fft_plan_create(n, FD_FFT_FORWARD); pb = fd_fft_plan_create(n, FD_FFT_BACKWARD);
    if(!pf || !pb){ printf("n=%d: plan creation failed\n", n); fail = 1; continue; }
    fd_fft_exec(pf, x, f);
    for(i = 0; i < n && n <= 1024; i++){        /* direct DFT, O(n^2), skipped for the largest size */
      s = 0;
      for(j = 0; j < n; j++) s += x[j] * cexp(-2.0 * M_PI * I * (double)i * (double)j / (double)n);
      d[i] = s;
    }
    fd_fft_exec(pb, f, g);
    for(i = 0; i < n; i++) g[i] /= n;
    printf("n=%5d: round trip rel err %.2e", n, rel_err(g, x, n));
    if(rel_err(g, x, n) > 1e-12) fail = 1;
    if(n <= 1024){
      printf(", vs direct DFT %.2e", rel_err(f, d, n));
      if(rel_err(f, d, n) > 1e-10) fail = 1;
    }
    printf("\n");
    fd_fft_plan_destroy(pf); fd_fft_plan_destroy(pb);
    free(x); free(f); free(g); free(d);
  }
  printf(fail ? "test_fd_fft: FAIL\n" : "test_fd_fft: PASS\n");
  return fail;
}
