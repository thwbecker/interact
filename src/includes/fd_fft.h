/*
  complex FFT for the fully dynamic tables (src/rsf_solve/fd_fft.c)

  plan-based interface modelled on FFTW so that the backend can be
  switched at compile time: with -DUSE_FFTW the plan wraps
  fftw_plan_dft_1d (any n); otherwise a built-in iterative radix-2
  transform is used and n must be a power of two (fd_fft_size_for
  rounds up). Transforms are unnormalized, as in FFTW: forward with
  sign -1, backward with sign +1, backward(forward(x)) = n x.
*/
#ifndef FD_FFT_H
#define FD_FFT_H
#include <complex.h>

#define FD_FFT_FORWARD  (-1)
#define FD_FFT_BACKWARD (+1)

typedef struct fd_fft_plan_s *fd_fft_plan;

int         fd_fft_size_for(int n);                 /* smallest allowed size >= n */
fd_fft_plan fd_fft_plan_create(int n, int sign);    /* n as returned by fd_fft_size_for */
void        fd_fft_exec(fd_fft_plan p, const double _Complex *in, double _Complex *out);
void        fd_fft_plan_destroy(fd_fft_plan p);
const char *fd_fft_backend(void);
#endif
