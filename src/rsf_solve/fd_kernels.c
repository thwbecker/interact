/*
  fully dynamic extension of rsf_solve: 2-D elastodynamic kernels

  Pure C, no PETSc dependency, so that the tables can be unit tested
  against the Python prototype (fully_dynamic/inplane.py, antiplane.py,
  ref_ando2007.py, inplane_series.py).

  Conventions (prototype): source segment along its own x axis with
  half-length a, receiver point (xr, yr) in the source frame, unit slip
  rate on |x| <= a since t = 0. Returned stresses are physical
  stress-from-slip (stress drop convention of interact's Is needs a
  global sign flip at assembly). In-plane: (sxx, sxy, syy); antiplane:
  (s31, s32).

  - fd_I_inplane:   Tada and Madariaga (2001) eqs A4-A6 = Ando et al.
                    (2007) eqs 4-6, I11, I12, I22
  - fd_I_antiplane: Ando et al. (2007) eqs 30-31
  - fd_seg_ramp_*:  -(mu/(2 beta)) [I(xr + a) - I(xr - a)]
  - fd_static_*:    static limits (edge and screw dislocation differences)
  - fd_series_inplane: large-time expansion of the lag kernel (t^-1, t^-3,
                    t^-5), generated coefficients in fd_series_coeffs.c
  - fd_project:     source-frame stress to shear and normal traction on
                    a receiver with tangent sr
*/
#include <math.h>
#include <stdlib.h>
#include "fd_kernels.h"

static double heav(double x){ return x > 0.0 ? 1.0 : 0.0; }
static double sgn(double x){ return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0); }

void fd_I_inplane(double x1, double x2, double t, double alpha, double beta, double I[3])
{
  double r = hypot(x1, x2), p = beta / alpha, HL, HT, gL, gT, gL3, gT3, q, q2, a1, a2, Hpl, arg, arc, ax;
  if(r <= 0.0){ I[0] = I[1] = I[2] = 0.0; return; }
  HL = heav(t - r / alpha); HT = heav(t - r / beta);
  gL = alpha * t / r; gL = gL * gL - 1.0; gL = gL > 0.0 ? sqrt(gL) : 0.0;
  gT = beta * t / r;  gT = gT * gT - 1.0; gT = gT > 0.0 ? sqrt(gT) : 0.0;
  gL3 = gL * gL * gL; gT3 = gT * gT * gT;
  q  = 2.0 * (3.0 * x1 * x1 - x2 * x2) / (3.0 * r * r);
  q2 = 2.0 * (3.0 * x2 * x2 - x1 * x1) / (3.0 * r * r);
  a1 = x1 * x1 / (r * r); a2 = x2 * x2 / (r * r);
  /* A4 */
  I[0] = -(1.0 / M_PI) * HL * (2.0 * x2 / r) * p * (q * p * p * gL3 + (1.0 - 2.0 * a2 * p * p) * gL)
         + (1.0 / M_PI) * HT * (2.0 * x2 / r) * (q * gT3 + (1.0 - 2.0 * a2) * gT);
  /* A5 */
  I[2] =  (1.0 / M_PI) * HL * (2.0 * x2 / r) * p * (q * p * p * gL3 + (2.0 * a1 * p * p - 1.0) * gL)
         - (1.0 / M_PI) * HT * (2.0 * x2 / r) * (q * gT3 + (2.0 * a1 - 1.0) * gT);
  /* A6 */
  Hpl = heav(x1) * heav(t - fabs(x2) / beta);
  arg = beta * t; arg = arg * arg - x2 * x2; arg = fabs(x1) / sqrt(arg > 1e-300 ? arg : 1e-300);
  if(arg > 1.0) arg = 1.0;
  arc = acos(arg);
  ax = fabs(x1);
  I[1] = Hpl
         - (1.0 / M_PI) * sgn(x1) * HL * (2.0 * ax / r) * p * (q2 * p * p * gL3 + 2.0 * a2 * p * p * gL)
         + (1.0 / M_PI) * sgn(x1) * HT * ((2.0 * ax / r) * (q2 * gT3 + 2.0 * a2 * gT) - arc);
}

void fd_seg_ramp_inplane(double xr, double yr, double t, double a, double alpha, double beta, double rho, double s[3])
{
  double mu = rho * beta * beta, Ip[3], Im[3], f = -(mu / (2.0 * beta));
  int k;
  fd_I_inplane(xr + a, yr, t, alpha, beta, Ip);
  fd_I_inplane(xr - a, yr, t, alpha, beta, Im);
  for(k = 0; k < 3; k++) s[k] = f * (Ip[k] - Im[k]);
}

void fd_I_antiplane(double x1, double x2, double t, double beta, double I[2])
{
  double r = hypot(x1, x2), Hs, gb, Hpl, arg, arc;
  if(r <= 0.0){ I[0] = I[1] = 0.0; return; }
  Hs = heav(t - r / beta);
  gb = beta * t / r; gb = gb * gb - 1.0; gb = gb > 0.0 ? sqrt(gb) : 0.0;
  I[0] = -(x2 / (M_PI * r)) * Hs * gb;                               /* eq. 30 */
  Hpl = heav(x1) * heav(t - fabs(x2) / beta);
  arg = beta * t; arg = arg * arg - x2 * x2; arg = fabs(x1) / sqrt(arg > 1e-300 ? arg : 1e-300);
  if(arg > 1.0) arg = 1.0;
  arc = acos(arg);
  I[1] = Hpl + (1.0 / M_PI) * sgn(x1) * Hs * ((fabs(x1) / r) * gb - arc);  /* eq. 31 */
}

void fd_seg_ramp_antiplane(double xr, double yr, double t, double a, double beta, double rho, double s[2])
{
  double mu = rho * beta * beta, Ip[2], Im[2], f = -(mu / (2.0 * beta));
  fd_I_antiplane(xr + a, yr, t, beta, Ip);
  fd_I_antiplane(xr - a, yr, t, beta, Im);
  s[0] = f * (Ip[0] - Im[0]); s[1] = f * (Ip[1] - Im[1]);
}

/* static limits */
static void edge(double x, double y, double mu, double nu, double s[3])
{
  double D = mu / (2.0 * M_PI * (1.0 - nu)), r2 = x * x + y * y;
  s[0] = -D * y * (3.0 * x * x + y * y) / (r2 * r2);
  s[1] =  D * x * (x * x - y * y) / (r2 * r2);
  s[2] =  D * y * (x * x - y * y) / (r2 * r2);
}
void fd_static_inplane(double xr, double yr, double a, double mu, double nu, double s[3])
{
  double p1[3], p2[3]; int k;
  edge(xr - a, yr, mu, nu, p1); edge(xr + a, yr, mu, nu, p2);
  for(k = 0; k < 3; k++) s[k] = p1[k] - p2[k];
}
void fd_static_antiplane(double xr, double yr, double a, double mu, double s[2])
{
  double x, y, r2;
  x = xr - a; y = yr; r2 = x * x + y * y;
  s[0] = -mu * y / (2.0 * M_PI * r2); s[1] = mu * x / (2.0 * M_PI * r2);
  x = xr + a; r2 = x * x + y * y;
  s[0] -= -mu * y / (2.0 * M_PI * r2); s[1] -= mu * x / (2.0 * M_PI * r2);
}

/* large-time series; coefficient functions in fd_series_coeffs.c */
double fd_c_11_m1(double, double, double, double); double fd_c_11_m3(double, double, double, double); double fd_c_11_m5(double, double, double, double);
double fd_c_22_m1(double, double, double, double); double fd_c_22_m3(double, double, double, double); double fd_c_22_m5(double, double, double, double);
double fd_c_12p_m1(double, double, double, double); double fd_c_12p_m3(double, double, double, double); double fd_c_12p_m5(double, double, double, double);

static double cm(int name, int k, double x1, double x2, double alpha, double beta)
{
  double v, X;
  if(name == 1){                          /* I12: odd extension of the non-constant part */
    X = fabs(x1);
    v = k == 1 ? fd_c_12p_m1(X, x2, alpha, beta) : (k == 3 ? fd_c_12p_m3(X, x2, alpha, beta) : fd_c_12p_m5(X, x2, alpha, beta));
    return x1 >= 0.0 ? v : -v;
  }
  if(name == 0)
    return k == 1 ? fd_c_11_m1(x1, x2, alpha, beta) : (k == 3 ? fd_c_11_m3(x1, x2, alpha, beta) : fd_c_11_m5(x1, x2, alpha, beta));
  return k == 1 ? fd_c_22_m1(x1, x2, alpha, beta) : (k == 3 ? fd_c_22_m3(x1, x2, alpha, beta) : fd_c_22_m5(x1, x2, alpha, beta));
}

/* dynamic correction over [t, t + dt] at large beta t / r, physical (sxx, sxy, syy):
   -(mu/(2 beta)) sum_k [c_k(x+a) - c_k(x-a)] [(t+dt)^-k - t^-k], k = 1, 3, 5 */
void fd_series_inplane(double xr, double yr, double t, double dt, double a, double alpha, double beta, double rho, double s[3])
{
  double mu = rho * beta * beta, f = -(mu / (2.0 * beta)), acc, dc;
  int name, k;
  for(name = 0; name < 3; name++){
    acc = 0.0;
    for(k = 1; k <= 5; k += 2){
      dc = cm(name, k, xr + a, yr, alpha, beta) - cm(name, k, xr - a, yr, alpha, beta);
      acc += dc * (pow(t + dt, -(double)k) - pow(t, -(double)k));
    }
    s[name] = f * acc;
  }
}

/* rotate a source-frame tensor to global and return shear and normal
   traction on a receiver with tangent sr (normal rot90(sr) = (-sr_y, sr_x));
   source frame axes in global: ss and rot90(ss) */
void fd_project(double sxx, double sxy, double syy, const double ss[2], const double sr[2], double *ts, double *tn)
{
  double ns[2] = {-ss[1], ss[0]}, nr[2] = {-sr[1], sr[0]};
  double Q[2][2] = {{ss[0], ns[0]}, {ss[1], ns[1]}};        /* columns: source axes */
  double S[2][2] = {{sxx, sxy}, {sxy, syy}}, T[2][2], G[2][2], tr[2];
  int i, j, k;
  for(i = 0; i < 2; i++) for(j = 0; j < 2; j++){ T[i][j] = 0.0; for(k = 0; k < 2; k++) T[i][j] += Q[i][k] * S[k][j]; }
  for(i = 0; i < 2; i++) for(j = 0; j < 2; j++){ G[i][j] = 0.0; for(k = 0; k < 2; k++) G[i][j] += T[i][k] * Q[j][k]; }
  tr[0] = G[0][0] * nr[0] + G[0][1] * nr[1];
  tr[1] = G[1][0] * nr[0] + G[1][1] * nr[1];
  *ts = tr[0] * sr[0] + tr[1] * sr[1];
  *tn = tr[0] * nr[0] + tr[1] * nr[1];
}
