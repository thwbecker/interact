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
  double *ks_off;               /* static stress per unit slip by offset, (2 nx - 1) x (2 nz - 1), [di + nx - 1][dj + nz - 1] */
  fd_fft_plan px, pz, ix, iz;   /* forward and backward 1-D plans */
  double _Complex *work;        /* Mx Mz */
};
struct fd3d_tab *fd3d_tab_create(int nx, int nz, double dx, double dt, double cL, double cT, double rho, int lmax);
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
  int nl;
  /* operators for steps of k = 2^ik substeps (fd3d_step_k), built on first use */
  double _Complex *kop_B[8], *kop_C[8];
  double kop_self[8], kop_self_top[8];
  int generic;
  double kdamp;                 /* relaxation of the long-step iteration (default 0.5) */                  /* test: k = 1 also through the iterative k-step path */
  double mu_s, mu_d, d0;                        /* slip weakening (uniform) */
  double *a, *b, *Lc, *theta, f0, V0;           /* rate and state */
  int *locked;                                  /* optional: no slip */
};
struct fd3d *fd3d_create(struct fd3d_tab *t);
/* mirror = 1: the tables were built for 2 nzp rows, see struct fd3d */
struct fd3d *fd3d_create2(struct fd3d_tab *t, int nzp, int mirror);
void fd3d_free(struct fd3d *f);
/* one fine step with an optional extra load (Pa, N values, may be NULL); returns Jacobi iterations */
int fd3d_step(struct fd3d *f, const double *load);
/* one step of k = 1, 2, 4, ..., 64 fine substeps (k = 1: fd3d_step); returns iterations,
   -1 on bad arguments, -2 if the iteration did not converge (state unchanged: retry with smaller k) */
int fd3d_step_k(struct fd3d *f, const double *load, int k);
void fd3d_push(struct fd3d *f, const double *V);   /* one slice into the history ring */
void fd3d_reset(struct fd3d *f);                   /* empty history, zero slip and time */
#endif
