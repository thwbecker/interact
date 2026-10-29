/*
  complex FFT for the fully dynamic tables, see fd_fft.h

  built-in backend: iterative radix-2 decimation in time with
  precomputed twiddles and a bit-reversal table, O(n log n), n a power
  of two; adequate for the table sizes of the prototype (n up to 2^14),
  replaceable by FFTW with -DUSE_FFTW without touching callers
*/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fd_fft.h"

#ifdef USE_FFTW
#include <fftw3.h>
struct fd_fft_plan_s {
  int n, sign;
  fftw_plan plan;
  fftw_complex *in, *out;
};

int fd_fft_size_for(int n)
{
  return n < 1 ? 1 : n;
}

fd_fft_plan fd_fft_plan_create(int n, int sign)
{
  fd_fft_plan p = (fd_fft_plan)malloc(sizeof(*p));
  p->n = n;
  p->sign = sign;
  p->in  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
  p->out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);
  p->plan = fftw_plan_dft_1d(n, p->in, p->out, sign < 0 ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  return p;
}
void fd_fft_exec(fd_fft_plan p, const double _Complex *in, double _Complex *out)
{
  memcpy(p->in, in, sizeof(fftw_complex) * p->n);
  fftw_execute(p->plan);
  memcpy(out, p->out, sizeof(fftw_complex) * p->n);
}
void fd_fft_plan_destroy(fd_fft_plan p)
{
  if(!p)
    return;
  fftw_destroy_plan(p->plan);
  fftw_free(p->in);
  fftw_free(p->out);
  free(p);
}
const char *fd_fft_backend(void)
{
  return "fftw3";
}

#else  /* built-in radix-2 */

struct fd_fft_plan_s {
  int n, log2n, sign;
  double _Complex *tw;
  int *rev;
};

int fd_fft_size_for(int n)
{
  int m = 1;
  while(m < n)
    m <<= 1;
  return m;
}

fd_fft_plan fd_fft_plan_create(int n, int sign)
{
  int i, j, bit, l2 = 0;
  fd_fft_plan p;
  if(n < 1 || (n & (n - 1)) != 0)
    return NULL;        /* must be a power of two */
  p = (fd_fft_plan)malloc(sizeof(*p));
  p->n = n;
  p->sign = sign < 0 ? -1 : 1;
  while((1 << l2) < n)
    l2++;
  p->log2n = l2;
  p->tw = (double _Complex *)malloc(sizeof(double _Complex) * (n / 2 > 0 ? n / 2 : 1));
  for(i = 0; i < n / 2; i++)
    p->tw[i] = cexp(p->sign * 2.0 * M_PI * I * (double)i / (double)n);
  p->rev = (int *)malloc(sizeof(int) * n);
  for(i = 0; i < n; i++){
    for(j = 0, bit = 0; bit < l2; bit++)
      if(i & (1 << bit))
        j |= 1 << (l2 - 1 - bit);
    p->rev[i] = j;
  }
  return p;
}

void fd_fft_exec(fd_fft_plan p, const double _Complex *in, double _Complex *out)
{
  int n = p->n, len, half, step, i, j;
  double _Complex u, v, w;
  for(i = 0; i < n; i++)
    out[p->rev[i]] = in[i];
  for(len = 2; len <= n; len <<= 1){
    half = len >> 1;
    step = n / len;
    for(i = 0; i < n; i += len){
      for(j = 0; j < half; j++){
        w = p->tw[j * step];
        u = out[i + j];
        v = out[i + j + half] * w;
        out[i + j] = u + v;
        out[i + j + half] = u - v;
      }
    }
  }
}

void fd_fft_plan_destroy(fd_fft_plan p)
{
  if(!p)
    return;
  free(p->tw);
  free(p->rev);
  free(p);
}
const char *fd_fft_backend(void)
{
  return "builtin radix-2";
}
#endif
