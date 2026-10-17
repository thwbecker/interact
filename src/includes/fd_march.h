/*
  fully dynamic marcher (fd_march.c): fixed-step, end-of-step collocation,
  variable substeps k dt_min, history convolution on the fine grid,
  coupled first-lag system solved exactly (banded). Port of the
  prototype's shells.System / stepper.Stepper; PETSc-free so that it can
  be unit tested (src/testing/rsf_solve/fd_dynamic.c) and used by
  rsf_fd.c.

  Element ordering: segment order (elements of segment 0, then 1, ...);
  off[k] is the first index of segment k. rsf_fd.c maps between this and
  interact's global order with the segment list's gidx.

  Normal stress: sigma_eff = sigma0 + nsign (Kn delta + history_n + K1n V),
  nsign = +1 for interact's convention (In compression positive, see
  fd_c_design.md section 9), -1 for the prototype's (Kn tension positive).
*/
#ifndef FD_MARCH_H
#define FD_MARCH_H
#include <complex.h>
#include "fd_tables.h"
#include "fd_segments.h"

struct fd_ops {                 /* operators for one substep multiple k */
  int k, bw;
  double *Ks, *Kn;              /* dense N x N cumulative kernels K_k (shear, normal) */
  double *ab;                   /* LAPACK band storage of Ks, ldab = 3 bw + 1, column-major */
};

struct fd_sys {
  int N, nseg, M, L, Lbuf, mode, has_normal, nsign;
  double dt, eta, alpha, beta, rho;
  struct fd_seglist *segs;
  struct fd_pair **pairs;       /* pairs[src*nseg + rcv] */
  int *off;
  /* history: per source segment a ring of Lbuf FFT'd slices */
  double _Complex **hist;
  int *hk;
  fd_fft_plan pf, pb;
  double _Complex *w1, *w2, *w3;
  /* static operators (dense, segment order) */
  double *Ks, *Kn;
  /* cumulative kernels per k (cached) */
  int nops;
  struct fd_ops ops[8];
};

struct fd_sys *fd_sys_create(struct fd_seglist *segs, struct fd_pair **pairs, int mode, double dt,
                             double alpha, double beta, double rho, int nsign);
void fd_sys_free(struct fd_sys *s);
void fd_sys_reset_history(struct fd_sys *s);
void fd_sys_push(struct fd_sys *s, const double *V);                /* one fine-grid slice */
void fd_sys_history(struct fd_sys *s, int lag0, double *ts, double *tn); /* sum over lags >= lag0 */
struct fd_ops *fd_sys_ops(struct fd_sys *s, int k);                /* K_k, built on first use */
void fd_sys_static(struct fd_sys *s, const double *delta, double *ts, double *tn);

/* banded solve (fd_band.c): dense N x N to LAPACK band storage, and
   solve with dgbsv on a copy of the band */
int  fd_band_width(int N, const double *A, double thr);
void fd_band_from_dense(int N, int bw, const double *A, double *ab);  /* ab: (3 bw + 1) x N */
int  fd_band_solve(int N, int bw, const double *ab, const double *rhs, double *x);
void fd_band_mult(int N, int bw, const double *ab, const double *x, double *y); /* y = A x from the band */

/* friction: slip weakening (linear or Dugdale, optional forced weakening) */
struct fd_slipweak {
  int N, dugdale;
  double *mu_s, *mu_d, d0;
  double *force_time, force_t0;   /* NULL: no forcing */
  int *locked;                    /* NULL: none */
};
int fd_slipweak_solve(struct fd_slipweak *f, struct fd_sys *s, struct fd_ops *op, const double *H, const double *Hn,
                      const double *sigma0, int normal_coupling, double t, double dt, const double *delta,
                      const double *Vguess, double *V, double *sigma);

/* marcher state */
struct fd_state {
  int N, k, n;
  double t;
  double *delta, *V, *tau, *sigma, *tau0, *taudot, *sigma0, *Hn_dot, *psi;
  double *H, *Hn, *hs, *hn, *work;
  int normal_coupling;
};
struct fd_state *fd_state_create(int N);
void fd_state_free(struct fd_state *st);
/* one step of k dt_min with slip weakening */
int fd_step_slipweak(struct fd_sys *s, struct fd_state *st, struct fd_slipweak *f, int k);
#endif
