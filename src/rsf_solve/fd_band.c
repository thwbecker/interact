/*
  banded linear algebra for the first-lag system (fd_march.h). LAPACK
  general band storage, column-major: entry (i, j) of A at
  ab[(2 bw + i - j) + j ldab], ldab = 3 bw + 1 (kl = ku = bw, with the
  extra kl rows that dgbsv needs for pivoting).
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd_march.h"

extern void dgbsv_(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);

/* the half bandwidth of the dense N x N matrix A (row major): the largest |i - j| with
   |A_ij| > thr. O(N^2), once per operator assembly */
int fd_band_width(int N, const double *A, double thr)
{
  int i, j, bw = 0;
  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      if(fabs(A[(size_t)i * N + j]) > thr && abs(i - j) > bw)
        bw = abs(i - j);
  return bw;
}

/* copy the band |i - j| <= bw of the dense A (row major) into LAPACK band storage ab
   (ldab = 3 bw + 1, zeroed first; the rows 0 .. bw - 1 stay free for dgbsv's fill-in) */
void fd_band_from_dense(int N, int bw, const double *A, double *ab)
{
  int i, j, ldab = 3 * bw + 1;
  memset(ab, 0, sizeof(double) * (size_t)ldab * N);
  for(j = 0; j < N; j++)
    for(i = (j - bw > 0 ? j - bw : 0); i <= (j + bw < N - 1 ? j + bw : N - 1); i++)
      ab[(2 * bw + i - j) + (size_t)j * ldab] = A[(size_t)i * N + j];
}

/* solve A x = rhs for the band matrix ab by LAPACK dgbsv (partial pivoting) on a copy, so ab
   is kept for the next right-hand side; returns dgbsv's info (0: success) */
int fd_band_solve(int N, int bw, const double *ab, const double *rhs, double *x)
{
  int ldab = 3 * bw + 1, kl = bw, ku = bw, nrhs = 1, info, n = N, ldb = N;
  double *abc = (double *)malloc(sizeof(double) * (size_t)ldab * N);
  int *ipiv = (int *)malloc(sizeof(int) * N);
  memcpy(abc, ab, sizeof(double) * (size_t)ldab * N);
  memcpy(x, rhs, sizeof(double) * N);
  dgbsv_(&n, &kl, &ku, &nrhs, abc, &ldab, ipiv, x, &ldb, &info);
  free(abc);
  free(ipiv);
  return info;
}

/* y = A x for the band matrix ab (the rows bw .. 3 bw of the storage) */
void fd_band_mult(int N, int bw, const double *ab, const double *x, double *y)
{
  int i, j, ldab = 3 * bw + 1;
  double s;
  for(i = 0; i < N; i++){
    s = 0.0;
    for(j = (i - bw > 0 ? i - bw : 0); j <= (i + bw < N - 1 ? i + bw : N - 1); j++)
      s += ab[(2 * bw + i - j) + (size_t)j * ldab] * x[j];
    y[i] = s;
  }
}
