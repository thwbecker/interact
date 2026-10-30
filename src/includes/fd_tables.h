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

/* segment kinds (fd_plan.md M7 step 2): a fault segment has slip as unknown and its shear
   traction as equation (friction); a virtual free surface element (in-plane) is represented by
   two segments of the same geometry, one with slip and zero shear traction (FD_SEG_SURF_SLIP),
   one with opening and zero normal traction (FD_SEG_SURF_OPEN); antiplane surfaces need the
   first only */
#define FD_SEG_FAULT     0
#define FD_SEG_SURF_SLIP 1
#define FD_SEG_SURF_OPEN 2
struct fd_segment {
  double X0[2];		/* centre of element 0 */
  double s[2];		/* unit tangent */
  double dx, a;		/* element length, half-length */
  int n;		/* number of elements */
  int kind;		/* FD_SEG_* */
};

/* compressed dense lag tables (fd_plan.md D2): per entry the increments of lags 1..hi (the
   window trimmed where they fall below the compression tolerance), runs of lags on a quadratic
   kept as coefficients of 1, y, y^2 (y = l/FD_CT_LNORM), the other lags explicitly */
#define FD_CT_LNORM 4096.0
struct fd_crun {
  int a, b;             /* lags a..b (1-based) */
  double c[3];
};
struct fd_ctab {
  size_t ne;
  int *hi;              /* per entry: last lag kept (0: none) */
  unsigned short *nm;   /* per entry: number of runs */
  size_t *woff, *moff;  /* per entry: first explicit value in W, first run in run */
  double *W;
  size_t nw, nwcap;
  struct fd_crun *run;
  size_t nrun, nruncap;
  size_t nlag_win;      /* statistics: lags 1..hi summed over the entries */
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
  /* compressed dense tables (fd_plan.md D2): Wd_s, Wd_n are then NULL */
  int compressed;
  struct fd_ctab cs, cn;
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
/* compression of the dense tables (fd_plan.md D2): ctol > 0 compresses with that tolerance
   relative to the window's envelope of each entry; ctol <= 0 keeps the full windows */
/* the cache of opening-source stress histories (surface elements), release after building */
void fd_tables_open_cache_clear(void);
void fd_tables_set_compress(double ctol);
/* lags 1..nexp of every entry stay explicit (default 1) */
void fd_tables_set_compress_explicit(int nexp);
double fd_tables_get_compress(void);
/* increment of a compressed table at lag l (1-based) for entry e */
double fd_ctab_value(const struct fd_ctab *c, size_t e, int l);
/* bytes of the lag data of a dense pair (compressed or not) and of its full windows */
double fd_pair_dense_bytes(const struct fd_pair *p, double *win);
/* prefix moments of one source's history for the compressed tables: P[k][m] = sum of
   y_l^k v(l) over lags lag0 <= l <= min(m, lav), m = 0..L; v(l) = ring + ((pos0 - l) mod Lbuf) n */
void fd_ctab_prefix(const double *ring, int Lbuf, long pos0, int lag0, int lav, int L, int n, double *P, const double **vl);
/* history of a compressed dense pair: out[i] += sum over source j and lags lag0..lav, receiver i
   of nr, source of ns elements */
void fd_ctab_hist(const struct fd_ctab *c, int nr, int ns, const double *P, const double *const *vl, int lag0, int lav, int L, double *out);
#endif
