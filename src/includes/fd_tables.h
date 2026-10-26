/*
  windowed dynamic-correction tables of the fully dynamic extension
  (fd_tables.c), one per (source segment, receiver segment) pair of
  parallel, equal-element segments (Toeplitz in the element offset).

  W[l], l = 0..L-1: column 0 is the current-step kernel S(et dt) - Ks dt,
  column l is S((et+l) dt) - S((et+l-1) dt) - Ks dt (prototype shells.py,
  inplane.W_table). Stored as the FFT along the padded offset axis
  (length M, power of two >= ns + nr - 1): offset m >= 0 at index m,
  m < 0 at index M + m.
*/
#ifndef FD_TABLES_H
#define FD_TABLES_H
#include <complex.h>
#include "fd_fft.h"

#define FD_MODE_INPLANE   0
#define FD_MODE_ANTIPLANE 1

struct fd_segment {
  double X0[2];		/* centre of element 0 */
  double s[2];		/* unit tangent */
  double dx, a;		/* element length, half-length */
  int n;		/* number of elements */
};

struct fd_pair {
  const struct fd_segment *src, *rcv;
  int mode, L, M, has_normal, nm;
  double dt, et;
  int *m;			/* offsets, nm = ns + nr - 1 entries, from -(ns-1) to nr-1 */
  double *ks_row, *kn_row;	/* static traction per offset (stress-from-slip) */
  double _Complex *Ws_hat, *Wn_hat;	/* L x M */
  fd_fft_plan pf, pb;
  /* non-parallel segments (or unequal elements, fd_plan.md N1): no Toeplitz structure; dense
     tables W[l][i][j] (receiver i, source j), l = 0 .. L-1 as the columns above, and the
     static blocks; Ws_hat/Wn_hat, m, ks_row, kn_row are then unused */
  int dense;
  double *Wd_s, *Wd_n, *Kd_s, *Kd_n;
};

/* build the table; lmax caps the window (cap_ok: accept the cap silently),
   lmax <= 0 uses the generous bound of the prototype; M <= 0 chooses the
   smallest power of two >= ns + nr - 1. Returns NULL on failure. */
struct fd_pair *fd_pair_build(const struct fd_segment *src, const struct fd_segment *rcv, int mode,
                              double dt, double et, double alpha, double beta, double rho,
                              double tol, int lmax, int cap_ok, int M, char *errmsg);
void fd_pair_free(struct fd_pair *p);
/* dense W[l] (l = 1..L in the prototype's 1-based count) and static blocks,
   receiver x source row-major, for tests and for the cumulative kernels */
void fd_pair_W_dense(const struct fd_pair *p, int l1, double *Ws, double *Wn);
void fd_pair_static_dense(const struct fd_pair *p, double *Ks, double *Kn);
int  fd_segments_parallel(const struct fd_segment *a, const struct fd_segment *b);
/* the dense variant, built for any pair (also used by fd_pair_build when the segments are not parallel) */
struct fd_pair *fd_pair_build_dense(const struct fd_segment *src, const struct fd_segment *rcv, int mode,
                                    double dt, double et, double alpha, double beta, double rho,
                                    double tol, int lmax, int cap_ok, int M, char *errmsg);
#endif
