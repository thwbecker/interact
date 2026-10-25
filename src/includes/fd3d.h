/*
  fully dynamic extension, 3-D: one planar fault of equal square elements on a
  regular grid in a whole space, slip along one fixed in-plane direction
  (rsf_solve's fixed rake). Port of fully_dynamic/tada_quad.py and fd3d.py.

  Kernels: stress response to a constant slip rate on a quadrantal fault
  (Tada 2006, GJI 164, Appendix A); a square element is four quadrants.
  Frame: fault plane x3 = 0, slip along x1, element grid axes x1 (index i,
  nx elements) and x2 (index j, nz elements). The traction computed is
  sigma_31 (shear along the slip direction).

  Tables: W[l](di, dj) = K[l] - Ks dt, K[l] = S(l dt) - S((l-1) dt), stored
  as 2-D FFTs of size Mx x Mz (powers of two >= 2n - 1); W vanishes exactly
  beyond the lag at which the S wave from an element's far corner has passed.
*/
#ifndef FD3D_H
#define FD3D_H
#include <complex.h>
#include "fd_fft.h"

/* storage precision of the tables and the history ring: double by default
   (the tests compare with the prototype at 1e-8); -DFD3D_SINGLE halves the
   memory again at about 1e-5 relative differences in slip-weakening runs */
#ifdef FD3D_SINGLE
typedef float _Complex fd3d_cplx;
#else
typedef double _Complex fd3d_cplx;
#endif

/* quadrant Green's function L_31/1 (dimensionless) and the element ramp stress */
double fd3d_L31(double x1, double x2, double x3, double t, double cL, double cT);
double fd3d_square_ramp(double x1, double x2, double x3, double t, double a, double cL, double cT, double rho);
/* full stress tensor (11, 22, 33, 12, 31, 32), comp 1 slip along x1 or 3 opening, rectangle a1 x a2 */
void fd3d_rect_ramp_tensor(double x1, double x2, double x3, double t, double a1, double a2,
                           double cL, double cT, double rho, int comp, double s[6]);

struct fd3d_tab {
  int nx, nz, Mx, Mz, L;
  double dx, dt, cL, cT, rho;
  int Mh;                       /* stored spectrum: kz = 0..Mz/2 (Hermitian half), Mx (Mz/2 + 1) values */
  fd3d_cplx *W_hat;             /* L x (Mx Mh), lag l at index l - 1 */
  double _Complex *Ks_hat;      /* static stress per unit slip, Mx Mz */
  double k1[3][3];              /* first-lag kernel on the 3 x 3 stencil, [di + 1][dj + 1] */
  double ks_self;
  int slipcomp;                 /* 1: slip along grid axis 1 (default), 2: along axis 2 */
  double *ks_off;               /* static stress per unit slip by offset, (2 nx - 1) x (2 nz - 1), [di + nx - 1][dj + nz - 1] */
  fd_fft_plan px, pz, ix, iz;   /* forward and backward 1-D plans */
  double _Complex *work;        /* Mx Mz */
};
struct fd3d_tab *fd3d_tab_create(int nx, int nz, double dx, double dt, double cL, double cT, double rho, int lmax);
/* comp 2: slip along grid axis 2 (e.g. dip slip with axis 1 along strike), traction s32 */
struct fd3d_tab *fd3d_tab_create_c(int nx, int nz, double dx, double dt, double cL, double cT, double rho, int lmax, int comp);
void fd3d_tab_free(struct fd3d_tab *t);
void fd3d_fft2(const struct fd3d_tab *t, double _Complex *a, int sign);   /* in place, unnormalised */

/* marcher state (scalar slip) */
struct fd3d {
  struct fd3d_tab *t;
  int N, nhist, hk;
  /* free surface by a mirror image (vertical fault, slip along strike, element row 0
     adjacent to the surface): the tables cover 2 nzp rows, image row nzp - 1 - j carries
     the slip of row j. nzp = t->nz and mirror = 0 without it */
  int nzp, mirror;
  /* optional external static operator (e.g. interact's exact half-space Is):
     out = Ks delta in element order; NULL: the tables' own static operator */
  void (*static_fn)(void *ctx, const double *delta, double *out);
  void *static_ctx;
  double time;
  double *delta, *V, *tau, *tau0, *sigma, *H, *work;
  fd3d_cplx *hist;              /* ring of L slices of FFT(V), half spectrum */
  double _Complex *acc;
  /* coarse history for aligned steps of Kc substeps (0: off) */
  int Kc;
  /* friction: 0 slip weakening, 1 rate and state (aging law) */
  int law;
  /* 0: piecewise-constant slip rate, end-of-step collocation (first order);
     1: second-order marching after Noda, Sato and Kurihara (2020): midpoint slip
     rates in the history, tau = H + (K1 + eta I)(V_old + V)/2 - eta V, trapezoid slip */
  int next;                     /* extension rows above the top row (surface-breaking faults), 0: none */
  int nl;
  /* operators for steps of k = 2^ik substeps (fd3d_step_k), built on first use */
  double _Complex *kop_B[8], *kop_C[8];
  double kop_self[8], kop_self_top[8];
  int generic;                  /* test: k = 1 also through the iterative k-step path */
  double kdamp;                 /* relaxation of the long-step iteration (default 0.5) */
  double *Vold;                 /* slip rates at the start of the current step */
  const double *vpres;          /* tests: prescribed fault slip rates for the next step (NULL: friction) */
  double *thnew;                /* internal: state at the end of a step being iterated (layer) */
  struct fd3d_layer *layer;     /* free-surface layer of opening elements (NULL: none) */
  double mu_s, mu_d, d0;                        /* slip weakening (uniform) */
  double *a, *b, *Lc, *theta, f0, V0;           /* rate and state */
  int *locked;                                  /* optional: no slip */
};
struct fd3d *fd3d_create(struct fd3d_tab *t);
/* mirror = 1: the tables were built for 2 nzp rows, see struct fd3d */
struct fd3d *fd3d_create2(struct fd3d_tab *t, int nzp, int mirror);
/* next > 0: the fault plane continues above the top row by next rows tied to row 0 (table nzp + next rows) */
struct fd3d *fd3d_create3(struct fd3d_tab *t, int nzp, int mirror, int next);
void fd3d_free(struct fd3d *f);
/* one fine step with an optional extra load (Pa, N values, may be NULL); returns Jacobi iterations */
int fd3d_step(struct fd3d *f, const double *load);
/* free-surface layer for a vertical strike-slip fault with its mirror image (fd_plan.md,
   stage A): opening-only rectangles on the surface plane z = 0, zero normal traction at
   their centres. Dense real-space tables in this first version (small problems, tests).
   Geometry: fault plane y = 0, strike along x, fault element (i, j) centred at
   (x0f + (i + 1/2) dx, 0, -(j + 1/2) dx); layer element (p, q) centred at
   (xs0 + (p + 1/2) hx, ys0 + (q + 1/2) hy, 0). */
struct fd3d_layer {
  int ns, nxs, nys, L, hk, nhist;
  double hx, hy, xs0, ys0, x0f, eta_p;
  double *Wss, *Kss, *Wfs, *Kfs, *Wsf, *Ksf;   /* [l][rcv][src]: layer<-layer, layer<-fault (with image), fault<-layer */
  double *K1ss, *K1fs, *K1sf;                    /* first lags */
  double *V, *Vold, *delta, *H, *tau;            /* layer opening rates, opening, known part, normal traction */
  double *hs, *hf;                               /* history rings of midpoint rates: layer (L x ns), fault (L x N) */
  double *Bss_lu; int *piv;                      /* factorised layer self operator */
  /* FFT-structured version (fft = 1): the layer is a narrow strip on a lattice of hx = dx/m
     along strike; layer-layer tables 2-D Toeplitz (Mx2 x My2), fault<->layer tables Toeplitz
     along strike (Mx) with dense (fault row) x (layer row) blocks per wavenumber */
  int fft, m, Mx, Mh, Mx2, My2, Mh2, nx, nz;
  double c0;
  fd_fft_plan pf, pb, pf2x, pb2x, pf2y, pb2y;
  double _Complex *Wss2, *Kss2, *K1ss2;          /* [l][Mx2*Mh2], half spectra in y */
  double _Complex *Wfs_h, *Kfs_h, *K1fs_h;       /* [l][k][q][j], k < Mh */
  double _Complex *Wsf_h, *Ksf_h, *K1sf_h;       /* [l][k][j][q] */
  double _Complex *hs2, *hsx, *hfx;              /* rings: layer 2-D, layer rows, fault rows (spectral) */
  double self_ss;                                /* self coefficient of the first-lag layer operator */
  long n_steps, n_outer, n_inner;                /* iteration statistics */
  /* full surface (surf = 1, dipping faults, no mirror): three unknowns per surface element
     (slip along x, along y, opening), zero traction (sxz, syz, szz); fault receivers carry
     the shear traction along the slip and the normal traction (tension positive) */
  int surf, fcomp;                               /* fcomp: 1 strike slip, 2 dip slip */
  double dip, ytop, ztop;                        /* fault: dip, top edge line (y, z) */
  double _Complex *T2[5][3];                     /* surface-surface: (0,0) (0,1) (1,0) (1,1) (2,2); [W lags | Ks | K1]; [k][q'][q] */
  double *yc, *wy, *self_c;                      /* row centres and widths; first-lag self coefficients [comp][row] */
  double _Complex *Tfs[3][3];                    /* fault -> surface traction c: [k][q][j] */
  double _Complex *Tsf[2][3][3];                 /* surface comp c -> fault traction r (0 shear, 1 normal): [k][j][q] */
  double self_sh, self_op;                       /* first-lag self coefficients: surface shear, opening */
  double *V3, *Vold3, *delta3, *H3;              /* surface unknowns, 3 x ns (component major) */
  double _Complex *h2_3, *hsx3;                  /* rings: 3 components */
  double *sig0, *tn;                             /* fault: normal stress at the phase start, normal traction change */
  void (*static_n_fn)(void *ctx, const double *delta, double *out);   /* exact static normal traction (e.g. -In delta) */
  int normal_coupling;                           /* 1: the fault's normal stress follows tn (default 1 in tests) */
  void *static_n_ctx;
};
struct fd3d_layer *fd3d_layer_create(struct fd3d *f, double x0f, int nxs, int nys, double hx, double hy,
                                     double xs0, double ys0, int L);
void fd3d_layer_free(struct fd3d_layer *ly);
/* FFT version: strip of nys rows of width hx = dx/m centred on the trace (y from -nys hx/2),
   nxs columns starting at xs0 (on the lattice relative to the fault's x0f), window L */
struct fd3d_layer *fd3d_layer_create_fft(struct fd3d *f, double x0f, int m, int nxs, int nys, double xs0, int L);
/* full virtual surface for a planar fault of any dip (no mirror): fault element (i, j) centred at
   (x0f + (i + 1/2) dx, ytop - (j + 1/2) dx cos(dip), ztop - (j + 1/2) dx sin(dip)), slip fcomp
   (1 strike, 2 dip, in the fault frame e1 = x, e2 = up dip; fsign -1 reverses the marcher's slip
   direction); surface lattice dx/m starting at (xs0, ys0) with nxs x nys elements */
struct fd3d_layer *fd3d_surf_create(struct fd3d *f, double x0f, double ytop, double ztop, double dip_deg, int fcomp,
                                    int m, int nxs, int nys, double xs0, double ys0, int L, double fsign);
/* the same with rows of given centres yc[q] and widths wy[q] (graded rows, e.g. fine next to the trace) */
struct fd3d_layer *fd3d_surf_create_rows(struct fd3d *f, double x0f, double ytop, double ztop, double dip_deg, int fcomp,
                                         int m, int nxs, int nys, double xs0, const double *yc, const double *wy, int L, double fsign);
int fd3d_step_layer(struct fd3d *f, const double *load);

/* one step of k = 1, 2, 4, ..., 64 fine substeps (k = 1: fd3d_step); returns iterations,
   -1 on bad arguments, -2 if the iteration did not converge (state unchanged: retry with smaller k) */
int fd3d_step_k(struct fd3d *f, const double *load, int k);
void fd3d_push(struct fd3d *f, const double *V);   /* one slice into the history ring */
void fd3d_reset(struct fd3d *f);                   /* empty history, zero slip and time */
#endif
