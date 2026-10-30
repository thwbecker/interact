/* 2-D elastodynamic kernels of the fully dynamic extension (fd_kernels.c) */
#ifndef FD_KERNELS_H
#define FD_KERNELS_H
void fd_I_inplane(double x1, double x2, double t, double alpha, double beta, double I[3]);
void fd_seg_ramp_inplane(double xr, double yr, double t, double a, double alpha, double beta, double rho, double s[3]);
/* opening (tensile) segment, unit opening rate since t = 0: physical (sxx, sxy, syy) in the
   source frame (fd_kernels_open.c, port of inplane_open.py); tshift > 0 averages t -+ tshift */
void fd_seg_ramp_inplane_open(double xr, double yr, double t, double a, double alpha, double beta, double rho,
                              double tshift, double s[3]);
/* full-plane static stress (sxx, sxy, syy) of unit opening on |x'| <= a, source frame */
void fd_static_inplane_open(double x, double y, double a, double mu, double nu, double st[3]);
void fd_I_antiplane(double x1, double x2, double t, double beta, double I[2]);
void fd_seg_ramp_antiplane(double xr, double yr, double t, double a, double beta, double rho, double s[2]);
void fd_static_inplane(double xr, double yr, double a, double mu, double nu, double s[3]);
void fd_static_antiplane(double xr, double yr, double a, double mu, double s[2]);
void fd_series_inplane(double xr, double yr, double t, double dt, double a, double alpha, double beta, double rho, double s[3]);
void fd_project(double sxx, double sxy, double syy, const double ss[2], const double sr[2], double *ts, double *tn);
#endif
