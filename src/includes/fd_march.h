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
  double *Ks, *Kn;              /* dense N x N cumulative kernels K_k (shear, normal); for nl: B_k */
  double *ab;                   /* LAPACK band storage of Ks, ldab = 3 bw + 1, column-major */
  /* second-order marching (nl = 1, Noda et al. 2020): tau = H + C_k V_old + B_k V with the
     slip rate linear across the k substeps; Ks/Kn/ab hold B_k, Cs/Cn the C_k operators */
  int nl;
  double *Cs, *Cn;
};

struct fd_sys {
  int N, nseg, nsrc, M, L, Lbuf, mode, has_normal, nsign;
  double dt, eta, alpha, beta, rho;
  struct fd_seglist *segs;      /* receivers (real segments) */
  struct fd_pair **pairs;       /* pairs[src*nseg + rcv], src over nsrc sources */
  int *src_of, *src_rev;        /* source s takes the slip of real segment src_of[s], reversed if src_rev[s] (image) */
  int *off;
  /* history: per source segment a ring of Lbuf FFT'd slices */
  double _Complex **hist;
  double **histr;               /* real-space ring per source, only with dense pairs (ndense > 0) */
  int ndense;
  int *hk;
  fd_fft_plan pf, pb;
  double _Complex *w1, *w2, *w3;
  /* static operators (dense, segment order) */
  double *Ks, *Kn;
  /* cumulative kernels per k (cached) */
  int nops;
  struct fd_ops ops[16];
  /* M6: uniformly partitioned convolution over the lag axis (conv = 1).
     Lags lag0..D are summed directly from the ring; lags D+1..L in
     partitions of B lags by overlap-save FFTs of length F = 2B along time,
     per spatial mode, once per block of B output times. */
  int conv, B, F, D, q0, kmax;
  int *P;                       /* partitions per pair */
  double _Complex **Hs, **Hn;   /* kernel spectra per pair, [p][m][f] */
  double _Complex **X;          /* window spectra per source (frequency-domain delay line), [slot][m][f] */
  int *nslot;                   /* slots per source */
  long **Xj;                    /* window index held in each slot, -1 none */
  double _Complex **Ys, **Yn;   /* far part for the current output block per pair, [m][i] */
  long *cblk;                   /* block index held in Ys/Yn per pair, -1 none */
  /* coarse history for long substeps: when every slice in a pair's window
     comes from steps of exactly Kc substeps that started at multiples of Kc,
     the input is constant on aligned Kc blocks and the history equals a sum
     over blocks with the block-aggregated kernel Ac[j] = sum_{l=jKc+1}^{(j+1)Kc} W[l]
     (exact, Kc times fewer terms) */
  int Kc, *Jc;                  /* block length (0: off), blocks per pair */
  double _Complex **Acs, **Acn; /* aggregated kernels per pair, [j][m], j = 1..Jc */
  double _Complex **Bcs, **Bcn; /* first moments for a linear rate within a block (second-order
                                   marching): sum_l W[l] w(l), w = ((j+1)Kc - l + 1/2)/Kc */
  int lin;                      /* history slices come from fd_sys_push_lin */
  long coarse_from;             /* fine index where the current run of aligned Kc steps began, -1 none */
  fd_fft_plan tf, tb;           /* time FFTs of length F */
  double _Complex *tw1, *tw2, *tz, *tzn;
};

/* nsrc sources: the first nseg are the real segments themselves (src_of = identity, not
   reversed); further ones are images (src_of[s] = real segment, src_rev[s] = 1 when the
   image's element order is the reverse of the real one, as for a vertical antiplane
   fault mirrored at the free surface). src_of/src_rev may be NULL when nsrc == nseg. */
struct fd_sys *fd_sys_create(struct fd_seglist *segs, struct fd_pair **pairs, int nsrc, const int *src_of, const int *src_rev,
                             int mode, double dt, double alpha, double beta, double rho, int nsign);
void fd_sys_free(struct fd_sys *s);
/* switch to the partitioned convolution (M6); call once after fd_sys_create,
   before any push. B: partition length (power of two, 0 = automatic),
   kmax: largest substep multiple the caller will use. Returns 0 on success. */
int fd_sys_set_partitioned(struct fd_sys *s, int B, int kmax);
int fd_sys_set_conv(struct fd_sys *s, int mode, int B, int kmax);   /* mode 1 partitioned, 2 hybrid */
/* enable the coarse history for steps of exactly K substeps (call before
   fd_sys_set_conv with mode 1, which releases the full tables) */
int fd_sys_set_coarse(struct fd_sys *s, int K);
double fd_sys_bytes(const struct fd_sys *s);
void fd_sys_reset_history(struct fd_sys *s);
void fd_sys_push(struct fd_sys *s, const double *V);
void fd_sys_push_k(struct fd_sys *s, const double *V, int k);     /* same slice k times */                /* one fine-grid slice */
void fd_sys_history(struct fd_sys *s, int lag0, double *ts, double *tn); /* sum over lags >= lag0 */
struct fd_ops *fd_sys_ops(struct fd_sys *s, int k);                /* K_k, built on first use */
struct fd_ops *fd_sys_ops_nl(struct fd_sys *s, int k);             /* B_k, C_k for second-order marching */
/* push k substep slices with the slip rate linear from Vold to V (second-order marching);
   suspends the coarse history for this step */
void fd_sys_push_lin(struct fd_sys *s, const double *Vold, const double *V, int kk);
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

/* friction: rate and state, aging law, regularized (rsf_solve's form:
   V = 2 v0 exp(-psi/a) sinh(tau/(a sigma))) */
struct fd_rsf {
  int N;
  double *a, *b, *dc;
  double v0, f0, tol;
  int max_iter;
};
double fd_rsf_tau(const struct fd_rsf *f, int i, double V, double psi, double sigma);
double fd_rsf_dtau(const struct fd_rsf *f, int i, double V, double psi, double sigma);
double fd_rsf_psi_advance(const struct fd_rsf *f, int i, double psi, double Va, double Vb, double dt);
/* coupled solve H + K V = tau_f(V, psi, sigma): banded Newton in log V with a bracketed
   elementwise fallback; returns the number of Newton iterations (negative if the fallback ran) */
int fd_rsf_solve(const struct fd_rsf *f, struct fd_sys *s, struct fd_ops *op, const double *H, const double *sigma,
                 const double *psi, const double *Vg, double *V);

/* marcher state */
struct fd_state {
  int N, k, n;
  double t;
  double *delta, *V, *tau, *sigma, *tau0, *taudot, *sigma0, *Hn_dot, *psi;
  double *H, *Hn, *hs, *hn, *work;
  int normal_coupling;
  int nl;                       /* 1: second-order marching (fd_sys_ops_nl, fd_sys_push_lin) */
};
struct fd_state *fd_state_create(int N);
void fd_state_free(struct fd_state *st);
/* one step of k dt_min with slip weakening / rate and state */
int fd_step_slipweak(struct fd_sys *s, struct fd_state *st, struct fd_slipweak *f, int k);
int fd_step_rsf(struct fd_sys *s, struct fd_state *st, const struct fd_rsf *f, int k);
#endif
