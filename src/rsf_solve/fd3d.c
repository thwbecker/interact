/* 3-D fully dynamic extension: kernels, tables, marcher; see fd3d.h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd3d.h"

/* ---------------- Tada (2006) Appendix A ---------------- */
static double H(double x){ return x > 0.0 ? 1.0 : 0.0; }
static double sgn(double x){ return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0); }
static double sq(double a){ return a > 0.0 ? sqrt(a) : 0.0; }

static double g1(double x1, double x2, double x3, double t, double c)
{
  double chi1 = sqrt(x1*x1 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3), s = sq(t*t - chi1*chi1/(c*c));
  double A = 2.0*H(x2)*H(t - chi1/c) - sgn(x2)*H(t - R/c);
  return s*s*s*A + (x2*(2.0*R*R + chi1*chi1)/(2.0*R*R*R)*t*t*t - 3.0*x2*chi1*chi1/(2.0*c*c*R)*t)*H(t - R/c);
}
static double g2(double x1, double x2, double x3, double t, double c)
{
  double chi1 = sqrt(x1*x1 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double A = 2.0*H(x2)*H(t - chi1/c) - sgn(x2)*H(t - R/c);
  return sq(t*t - chi1*chi1/(c*c))*A + x2/R*t*H(t - R/c);
}
static double g3(double x1, double x2, double x3, double t, double c)
{
  double chi1 = sqrt(x1*x1 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double A = 2.0*H(x2)*H(t - chi1/c) - sgn(x2)*H(t - R/c);
  return atan(sq(t*t - chi1*chi1/(c*c))/(x1/c))*A + atan(x2/x1)*H(t - R/c);
}
static double g4(double R, double t, double c){ return (t*t*t - R*R/(c*c)*t)*H(t - R/c); }

double fd3d_L31(double x1, double x2, double x3, double t, double cL, double cT)
{
  double p = cT/cL, chi1 = sqrt(x1*x1 + x3*x3), chi2 = sqrt(x2*x2 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double c3 = cT*cT*cT, chi1_2 = chi1*chi1, chi1_4 = chi1_2*chi1_2, chi1_6 = chi1_4*chi1_2, R5 = R*R*R*R*R;
  double v = H(x1)*H(x2)*H(t - fabs(x3)/cT)
    - (g3(x1, x2, x3, t, cT) + g3(x2, x1, x3, t, cT))/(2.0*M_PI)
    - 2.0*c3*x1*(3.0*x3*x3 - x1*x1)/(3.0*M_PI*chi1_6)*(g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT));
  if(x3 != 0.0)
    v += - c3*x1*x2*x3*x3/(M_PI*R5*chi1_2)*(g4(R, t, cL) - g4(R, t, cT))
         - 2.0*cT*x1*x3*x3/(M_PI*chi1_4)*(p*p*g2(x1, x2, x3, t, cL) - g2(x1, x2, x3, t, cT));
  v += cT*x2/(2.0*M_PI*chi2*chi2)*g2(x2, x1, x3, t, cT);
  return v;
}

/* full quadrant Green's functions (Tada 2006, A18-A27): L_ij/1 (slip along x1) and
   L_ij/3 (opening), returned as s[6] = (11, 22, 33, 12, 31, 32) */
static double g5(double R, double t, double c){ return (t - R/c)*H(t - R/c); }
static void L1_all(double x1, double x2, double x3, double t, double cL, double cT, double s[6])
{
  double p = cT/cL, chi1 = sqrt(x1*x1 + x3*x3), chi2 = sqrt(x2*x2 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double c3 = cT*cT*cT, chi1_2 = chi1*chi1, R5 = R*R*R*R*R, pi = M_PI;
  double dg1 = g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT), dg4 = g4(R, t, cL) - g4(R, t, cT);
  double g2L = g2(x1, x2, x3, t, cL), g2T = g2(x1, x2, x3, t, cT), g2sT = g2(x2, x1, x3, t, cT);
  s[0] = -2.0*c3*x3*(3.0*x1*x1 - x3*x3)/(3.0*pi*chi1_2*chi1_2*chi1_2)*dg1 - c3*x1*x1*x2*x3/(pi*R5*chi1_2)*dg4
         - cT*x3/(pi*chi1_2)*((1.0 - 2.0*p*p + 2.0*p*p*x1*x1/chi1_2)*g2L + (1.0 - 2.0*x1*x1/chi1_2)*g2T);
  s[3] = c3*x1*x3/(pi*R5)*dg4 - cT*x3/(2.0*pi*chi2*chi2)*g2sT;
  s[1] = c3*x2*x3/(pi*R5)*dg4 - (1.0 - 2.0*p*p)*cT*x3/(pi*chi1_2)*g2L;
  s[4] = fd3d_L31(x1, x2, x3, t, cL, cT);
  s[5] = c3/(pi*R5)*(x3*x3 - R*R/3.0)*dg4 + cT/(pi*R)*(2.0*p*p/3.0*g5(R, t, cL) - g5(R, t, cT)/6.0);
  s[2] = 2.0*c3*x3*(3.0*x1*x1 - x3*x3)/(3.0*pi*chi1_2*chi1_2*chi1_2)*dg1 - c3*x2*x3*x3*x3/(pi*R5*chi1_2)*dg4
         - cT*x3/(pi*chi1_2)*((1.0 - 2.0*p*p + 2.0*p*p*x3*x3/chi1_2)*g2L + (1.0 - 2.0*x3*x3/chi1_2)*g2T);
}
static double L11_3(double x1, double x2, double x3, double t, double cL, double cT)
{
  double p = cT/cL, chi1 = sqrt(x1*x1 + x3*x3), chi2 = sqrt(x2*x2 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double c3 = cT*cT*cT, chi1_2 = chi1*chi1, R5 = R*R*R*R*R, pi = M_PI;
  double dg1 = g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT), dg4 = g4(R, t, cL) - g4(R, t, cT);
  return (1.0 - 2.0*p*p)/p*H(x1)*H(x2)*H(t - fabs(x3)/cL)
         - (1.0 - 2.0*p*p)/(2.0*pi*p)*(g3(x1, x2, x3, t, cL) + g3(x2, x1, x3, t, cL))
         - 2.0*c3*x1*(3.0*x3*x3 - x1*x1)/(3.0*pi*chi1_2*chi1_2*chi1_2)*dg1 - c3*x1*x2*x3*x3/(pi*R5*chi1_2)*dg4
         - 2.0*cT*x1*x3*x3/(pi*chi1_2*chi1_2)*(p*p*g2(x1, x2, x3, t, cL) - g2(x1, x2, x3, t, cT))
         + (1.0 - 2.0*p*p)*cT*x2/(pi*chi2*chi2)*g2(x2, x1, x3, t, cL);
}
static void L3_all(double x1, double x2, double x3, double t, double cL, double cT, double s[6])
{
  double p = cT/cL, chi1 = sqrt(x1*x1 + x3*x3), chi2 = sqrt(x2*x2 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double c3 = cT*cT*cT, R5 = R*R*R*R*R, pi = M_PI, dg4 = g4(R, t, cL) - g4(R, t, cT), a[6], b[6];
  s[0] = L11_3(x1, x2, x3, t, cL, cT);
  s[1] = L11_3(x2, x1, x3, t, cL, cT);
  s[3] = c3/(pi*R5)*(x3*x3 - R*R/3.0)*dg4 + cT/(pi*R)*((1.0 - 4.0*p*p/3.0)*g5(R, t, cL) - 2.0/3.0*g5(R, t, cT));
  L1_all(x1, x2, x3, t, cL, cT, a); s[4] = a[2];          /* L_31/3 = L_33/1 */
  L1_all(x2, x1, x3, t, cL, cT, b); s[5] = b[2];          /* L_32/3 by symmetry */
  s[2] = 1.0/p*H(x1)*H(x2)*H(t - fabs(x3)/cL) - 1.0/(2.0*pi*p)*(g3(x1, x2, x3, t, cL) + g3(x2, x1, x3, t, cL))
         + 2.0*c3*x1*(3.0*x3*x3 - x1*x1)/(3.0*pi*pow(chi1, 6))*(g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT))
         + 2.0*c3*x2*(3.0*x3*x3 - x2*x2)/(3.0*pi*pow(chi2, 6))*(g1(x2, x1, x3, t, cL) - g1(x2, x1, x3, t, cT))
         + c3*x1*x2*x3*x3/(pi*R5)*(1.0/(chi1*chi1) + 1.0/(chi2*chi2))*dg4
         + cT*x1/(pi*chi1*chi1)*((2.0*(1.0 - p*p) + 2.0*p*p*x3*x3/(chi1*chi1))*g2(x1, x2, x3, t, cL) - 2.0*x3*x3/(chi1*chi1)*g2(x1, x2, x3, t, cT))
         + cT*x2/(pi*chi2*chi2)*((2.0*(1.0 - p*p) + 2.0*p*p*x3*x3/(chi2*chi2))*g2(x2, x1, x3, t, cL) - 2.0*x3*x3/(chi2*chi2)*g2(x2, x1, x3, t, cT));
}
static void rect_tensor0(double x1, double x2, double x3, double t, double a1, double a2, double cL, double cT, double rho, int comp, double s[6])
{
  double mu = rho*cT*cT, q[6]; int k, i1, i2;
  for(k = 0; k < 6; k++) s[k] = 0.0;
  for(i1 = 0; i1 < 2; i1++) for(i2 = 0; i2 < 2; i2++){
    double c1 = i1 ? a1 : -a1, c2 = i2 ? a2 : -a2, sg = (i1 == i2) ? 1.0 : -1.0;
    if(comp == 3) L3_all(x1 - c1, x2 - c2, x3, t, cL, cT, q);
    else if(comp == 2){
      /* slip along x2: L_ij/2(x1, x2) = L_i'j'/1(x2, x1) with the indices 1 <-> 2 swapped */
      double r[6];
      L1_all(x2 - c2, x1 - c1, x3, t, cL, cT, r);
      q[0] = r[1]; q[1] = r[0]; q[2] = r[2]; q[3] = r[3]; q[4] = r[5]; q[5] = r[4];
    }else L1_all(x1 - c1, x2 - c2, x3, t, cL, cT, q);
    for(k = 0; k < 6; k++) s[k] += sg*q[k];
  }
  for(k = 0; k < 6; k++) s[k] *= -mu/(2.0*cT);
}
/* full stress tensor (11, 22, 33, 12, 31, 32) at (x1, x2, x3) from unit rate since t = 0 on the
   rectangle |xi1| < a1, |xi2| < a2: comp 1 slip along x1, comp 2 slip along x2, comp 3 opening. Receivers on the
   extension of an element edge hit removable 0/0: averaged over a tiny symmetric offset */
void fd3d_rect_ramp_tensor(double x1, double x2, double x3, double t, double a1, double a2,
                           double cL, double cT, double rho, int comp, double s[6])
{
  double e = 1e-9*(a1 > a2 ? a1 : a2), lo[6], hi[6]; int k;
  /* in the source plane the quadrant fields are singular along the shared edge line and cancel
     only between quadrants: a larger offset (1e-3 a) avoids the cancellation, off the plane the
     0/0 is removable and 1e-7 a suffices (same rule as tada_quad.py) */
  double d = (fabs(x3) < e ? 1e-3 : 1e-7)*(a1 > a2 ? a1 : a2);
  if((fabs(fabs(x1) - a1) < e) || (fabs(fabs(x2) - a2) < e)){
    rect_tensor0(x1 - d, x2 - d, x3, t, a1, a2, cL, cT, rho, comp, lo);
    rect_tensor0(x1 + d, x2 + d, x3, t, a1, a2, cL, cT, rho, comp, hi);
    for(k = 0; k < 6; k++) s[k] = 0.5*(lo[k] + hi[k]);
  }else
    rect_tensor0(x1, x2, x3, t, a1, a2, cL, cT, rho, comp, s);
}

/* sigma_31 at (x1, x2, x3) from unit slip rate along x1 on the square |xi| < a since t = 0 */
double fd3d_square_ramp(double x1, double x2, double x3, double t, double a, double cL, double cT, double rho)
{
  double mu = rho*cT*cT, s = 0.0;
  s += fd3d_L31(x1 + a, x2 + a, x3, t, cL, cT);
  s -= fd3d_L31(x1 - a, x2 + a, x3, t, cL, cT);
  s -= fd3d_L31(x1 + a, x2 - a, x3, t, cL, cT);
  s += fd3d_L31(x1 - a, x2 - a, x3, t, cL, cT);
  return -mu/(2.0*cT)*s;
}

/* ---------------- 2-D FFT from the 1-D plans ---------------- */
void fd3d_fft2(const struct fd3d_tab *t, double _Complex *a, int sign)
{
  int i, j, Mx = t->Mx, Mz = t->Mz;
  double _Complex *col = (double _Complex *)malloc(sizeof(double _Complex)*(Mx > Mz ? Mx : Mz));
  double _Complex *out = (double _Complex *)malloc(sizeof(double _Complex)*(Mx > Mz ? Mx : Mz));
  for(i = 0; i < Mx; i++){                       /* rows: index j contiguous */
    fd_fft_exec(sign < 0 ? t->pz : t->iz, a + (size_t)i*Mz, out);
    memcpy(a + (size_t)i*Mz, out, sizeof(double _Complex)*Mz);
  }
  for(j = 0; j < Mz; j++){
    for(i = 0; i < Mx; i++) col[i] = a[(size_t)i*Mz + j];
    fd_fft_exec(sign < 0 ? t->px : t->ix, col, out);
    for(i = 0; i < Mx; i++) a[(size_t)i*Mz + j] = out[i];
  }
  free(col); free(out);
}

/* store the Hermitian half (kz = 0..Mz/2) of a full spectrum */
static void to_half(const struct fd3d_tab *t, const double _Complex *full, fd3d_cplx *half)
{
  int i, j, Mz = t->Mz, Mh = t->Mh;
  for(i = 0; i < t->Mx; i++) for(j = 0; j < Mh; j++) half[(size_t)i*Mh + j] = (fd3d_cplx)full[(size_t)i*Mz + j];
}
/* rebuild the full spectrum of a real field from its half (double accumulator) */
static void from_half(const struct fd3d_tab *t, double _Complex *full)
{
  int i, j, Mx = t->Mx, Mz = t->Mz;
  for(i = 0; i < Mx; i++) for(j = Mz/2 + 1; j < Mz; j++)
    full[(size_t)i*Mz + j] = conj(full[(size_t)((Mx - i) % Mx)*Mz + (Mz - j)]);
}

static void put_offsets(const struct fd3d_tab *t, const double *A, double _Complex *dst)
{
  /* A: (2 nx - 1) x (2 nz - 1), offset (di, dj) at [di + nx - 1][dj + nz - 1] -> FFT order */
  int i, j, nx = t->nx, nz = t->nz;
  memset(dst, 0, sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  for(i = -(nx - 1); i <= nx - 1; i++)
    for(j = -(nz - 1); j <= nz - 1; j++)
      dst[(size_t)((i + t->Mx) % t->Mx)*t->Mz + (j + t->Mz) % t->Mz] = A[(size_t)(i + nx - 1)*(2*nz - 1) + (j + nz - 1)];
}

/* in-plane traction along the slip for slip along grid axis 1 (comp 1, s31) or axis 2 (comp 2, s32) */
static double tab_kernel(int comp, double x1, double x2, double t, double a, double cL, double cT, double rho)
{
  double s[6];
  if(comp == 1) return fd3d_square_ramp(x1, x2, 0.0, t, a, cL, cT, rho);
  fd3d_rect_ramp_tensor(x1, x2, 0.0, t, a, a, cL, cT, rho, 2, s);
  return s[5];
}
static int tab_comp = 1;
struct fd3d_tab *fd3d_tab_create_c(int nx, int nz, double dx, double dt, double cL, double cT, double rho, int lmax, int comp)
{
  struct fd3d_tab *t;
  tab_comp = comp;
  t = fd3d_tab_create(nx, nz, dx, dt, cL, cT, rho, lmax);
  tab_comp = 1;
  if(t) t->slipcomp = comp;
  return t;
}
struct fd3d_tab *fd3d_tab_create(int nx, int nz, double dx, double dt, double cL, double cT, double rho, int lmax)
{
  struct fd3d_tab *t = (struct fd3d_tab *)calloc(1, sizeof(*t));
  int ox = 2*nx - 1, oz = 2*nz - 1, i, j, l, Lmax = 0;
  double a = dx/2, *S0, *S1, *A, rmax, T1, T2;
  t->nx = nx; t->nz = nz; t->dx = dx; t->dt = dt; t->cL = cL; t->cT = cT; t->rho = rho;
  t->Mx = fd_fft_size_for(ox); t->Mz = fd_fft_size_for(oz);
  t->px = fd_fft_plan_create(t->Mx, FD_FFT_FORWARD); t->ix = fd_fft_plan_create(t->Mx, FD_FFT_BACKWARD);
  t->pz = fd_fft_plan_create(t->Mz, FD_FFT_FORWARD); t->iz = fd_fft_plan_create(t->Mz, FD_FFT_BACKWARD);
  t->work = (double _Complex *)malloc(sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  for(i = -(nx - 1); i <= nx - 1; i++)
    for(j = -(nz - 1); j <= nz - 1; j++){
      rmax = sqrt((fabs(i*dx) + a)*(fabs(i*dx) + a) + (fabs(j*dx) + a)*(fabs(j*dx) + a));
      l = (int)ceil(rmax/cT/dt) + 1;
      if(l > Lmax) Lmax = l;
    }
  t->L = Lmax + 1;
  if(lmax > 0 && t->L > lmax) t->L = lmax;
  S0 = (double *)calloc((size_t)ox*oz, sizeof(double)); S1 = (double *)malloc(sizeof(double)*(size_t)ox*oz);
  A = (double *)malloc(sizeof(double)*(size_t)ox*oz);
  /* static: late-time slope */
  rmax = sqrt(2.0)*(nx > nz ? nx : nz)*dx + dx;
  T1 = 1.5*rmax/cT + 10.0*dt; T2 = 2.0*T1;
  for(i = 0; i < ox; i++)
    for(j = 0; j < oz; j++){
      double x1 = (i - (nx - 1))*dx, x2 = (j - (nz - 1))*dx;
      A[(size_t)i*oz + j] = (tab_kernel(tab_comp, x1, x2, T2, a, cL, cT, rho) - tab_kernel(tab_comp, x1, x2, T1, a, cL, cT, rho))/(T2 - T1);
    }
  t->ks_self = A[(size_t)(nx - 1)*oz + (nz - 1)];
  t->ks_off = (double *)malloc(sizeof(double)*(size_t)ox*oz);
  memcpy(t->ks_off, A, sizeof(double)*(size_t)ox*oz);
  t->Ks_hat = (double _Complex *)malloc(sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  put_offsets(t, A, t->Ks_hat); fd3d_fft2(t, t->Ks_hat, -1);
  t->Mh = t->Mz/2 + 1;
  t->W_hat = (fd3d_cplx *)malloc(sizeof(fd3d_cplx)*(size_t)t->L*t->Mx*t->Mh);
  for(l = 1; l <= t->L; l++){
    double *Ks = (double *)malloc(sizeof(double)*(size_t)ox*oz);
    /* recompute the static rows cheaply from the stored slope array: keep A as Ks */
    memcpy(Ks, A, sizeof(double)*(size_t)ox*oz);
    for(i = 0; i < ox; i++)
      for(j = 0; j < oz; j++){
        double x1 = (i - (nx - 1))*dx, x2 = (j - (nz - 1))*dx;
        S1[(size_t)i*oz + j] = tab_kernel(tab_comp, x1, x2, l*dt, a, cL, cT, rho);
      }
    if(l == 1){
      int di, dj;
      for(di = -1; di <= 1; di++) for(dj = -1; dj <= 1; dj++)
        t->k1[di + 1][dj + 1] = (nx > 1 || di == 0) && (nz > 1 || dj == 0) ?
          S1[(size_t)(di + nx - 1)*oz + (dj + nz - 1)] : 0.0;
    }
    for(i = 0; i < ox*oz; i++) Ks[i] = S1[i] - S0[i] - Ks[i]*dt;      /* W[l] */
    put_offsets(t, Ks, t->work);
    fd3d_fft2(t, t->work, -1);
    to_half(t, t->work, t->W_hat + (size_t)(l - 1)*t->Mx*t->Mh);
    memcpy(S0, S1, sizeof(double)*(size_t)ox*oz);
    free(Ks);
  }
  free(S0); free(S1); free(A);
  return t;
}

void fd3d_tab_free(struct fd3d_tab *t)
{
  if(!t) return;
  free(t->W_hat); free(t->Ks_hat); free(t->work); free(t->ks_off);
  fd_fft_plan_destroy(t->px); fd_fft_plan_destroy(t->ix); fd_fft_plan_destroy(t->pz); fd_fft_plan_destroy(t->iz);
  free(t);
}

/* ---------------- marcher ---------------- */
struct fd3d *fd3d_create(struct fd3d_tab *t)
{
  return fd3d_create2(t, t->nz, 0);
}

struct fd3d *fd3d_create3(struct fd3d_tab *t, int nzp, int mirror, int next);
struct fd3d *fd3d_create2(struct fd3d_tab *t, int nzp, int mirror){ return fd3d_create3(t, nzp, mirror, 0); }
/* next > 0: the fault plane continues above its top row by next rows whose slip rate copies row 0
   (a surface-breaking fault without an edge in the surface; the table must have nzp + next rows) */
struct fd3d *fd3d_create3(struct fd3d_tab *t, int nzp, int mirror, int next)
{
  struct fd3d *f = (struct fd3d *)calloc(1, sizeof(*f));
  int N = t->nx*nzp;
  if(mirror && t->nz != 2*nzp){ free(f); return NULL; }
  if(next && (mirror || t->nz != nzp + next)){ free(f); return NULL; }
  f->nzp = nzp; f->mirror = mirror; f->next = next;
  f->kdamp = getenv("FD3D_KDAMP") ? atof(getenv("FD3D_KDAMP")) : 0.5;
  size_t M = (size_t)t->Mx*t->Mz;
  f->t = t; f->N = N;
  f->delta = calloc(N, sizeof(double)); f->V = calloc(N, sizeof(double)); f->tau = calloc(N, sizeof(double));
  f->tau0 = calloc(N, sizeof(double)); f->sigma = calloc(N, sizeof(double)); f->H = calloc(N, sizeof(double));
  f->work = calloc(N, sizeof(double));
  f->hist = (fd3d_cplx *)calloc((size_t)t->L*t->Mx*t->Mh, sizeof(fd3d_cplx));
  f->acc = (double _Complex *)malloc(sizeof(double _Complex)*M);
  return f;
}

void fd3d_free(struct fd3d *f)
{
  if(!f) return;
  free(f->delta); free(f->V); free(f->tau); free(f->tau0); free(f->sigma); free(f->H); free(f->work);
  { int q; for(q = 0; q < 8; q++){ free(f->kop_B[q]); free(f->kop_C[q]); } }
  free(f->Vold); free(f->dpath); free(f->coh);
  if(f->layer) fd3d_layer_free(f->layer);
  free(f->hist); free(f->acc); free(f);
}

static void field_hat(const struct fd3d *f, const double *x, double _Complex *out)
{
  const struct fd3d_tab *t = f->t;
  int i, j;
  int nzp = f->nzp;
  memset(out, 0, sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  for(i = 0; i < t->nx; i++) for(j = 0; j < nzp; j++){
    if(f->mirror){
      out[(size_t)i*t->Mz + nzp + j] = x[i*nzp + j];
      out[(size_t)i*t->Mz + nzp - 1 - j] = x[i*nzp + j];      /* image */
    }else if(f->next){
      out[(size_t)i*t->Mz + f->next + j] = x[i*nzp + j];
      if(j == 0){ int r; for(r = 0; r < f->next; r++) out[(size_t)i*t->Mz + r] = x[i*nzp]; }   /* extension rows copy row 0 */
    }else
      out[(size_t)i*t->Mz + j] = x[i*nzp + j];
  }
  fd3d_fft2(t, out, -1);
}

/* tau0 + Ks delta + history (lags >= 2) */
static void known_k(struct fd3d *f, int k);
static void known(struct fd3d *f){ known_k(f, 1); }

/* tau0 + Ks delta + history for a step of k substeps: the most recent stored
   slice is k + 1 lags old at the end of the step */
static void known_k(struct fd3d *f, int k)
{
  const struct fd3d_tab *t = f->t;
  size_t M = (size_t)t->Mx*t->Mz, m;
  int l, i, j, lmax = t->L < f->nhist + k ? t->L : f->nhist + k;
  int row0 = f->mirror ? f->nzp : f->next, nzp = f->nzp;
  if(f->static_fn)
    memset(f->acc, 0, sizeof(double _Complex)*M);
  else{
    field_hat(f, f->delta, t->work);
    for(m = 0; m < M; m++) f->acc[m] = t->Ks_hat[m]*t->work[m];
  }
  {
    size_t Mhh = (size_t)t->Mx*t->Mh; int ii, jj;
    double _Complex *ah = (double _Complex *)calloc(Mhh, sizeof(double _Complex));
    for(l = k + 1; l <= lmax; l++){
      const fd3d_cplx *W = t->W_hat + (size_t)(l - 1)*Mhh;
      const fd3d_cplx *Vh = f->hist + (size_t)((f->hk - (l - k) + 2*t->L) % t->L)*Mhh;
      for(m = 0; m < Mhh; m++) ah[m] += (double _Complex)(W[m]*Vh[m]);
    }
    for(ii = 0; ii < t->Mx; ii++) for(jj = 0; jj < t->Mh; jj++) f->acc[(size_t)ii*t->Mz + jj] += ah[(size_t)ii*t->Mh + jj];
    free(ah);
  }
  from_half(t, f->acc);
  fd3d_fft2(t, f->acc, +1);
  for(i = 0; i < t->nx; i++) for(j = 0; j < nzp; j++)
    f->H[i*nzp + j] = f->tau0[i*nzp + j] + creal(f->acc[(size_t)i*t->Mz + row0 + j])/M;
  if(f->static_fn){
    f->static_fn(f->static_ctx, f->delta, f->work);
    for(i = 0; i < f->N; i++) f->H[i] += f->work[i];
  }
}

/* first-lag self coefficient of row j (with the mirror, the top row also feels its own image) */
static double kdiag(const struct fd3d *f, int j)
{
  /* the top row also sees its own image (mirror) or its own copy in the extension row above */
  return f->t->k1[1][1] + (((f->mirror || f->next) && j == 0) ? f->t->k1[1][2] : 0.0);
}

/* first-lag traction from neighbours (excluding self) */
static void apply_k1(const struct fd3d *f, const double *V, double *out)
{
  const struct fd3d_tab *t = f->t;
  int i, j, di, dj, nx = t->nx, nz = f->nzp;
  for(i = 0; i < nx; i++)
    for(j = 0; j < nz; j++){
      double s = 0.0;
      for(di = -1; di <= 1; di++)
        for(dj = -1; dj <= 1; dj++){
          int si = i - di, sj = j - dj;               /* receiver - source = (di, dj) */
          if((di == 0 && dj == 0) || si < 0 || si >= nx || sj >= nz) continue;
          if(sj < 0){
            if(!f->mirror && !f->next) continue;
            if(di == 0) continue;                     /* own image or copy: in kdiag */
            sj = f->mirror ? -1 - sj : 0;             /* image of row -1 - sj, or the extension's copy of row 0 */
          }
          s += t->k1[di + 1][dj + 1]*V[si*nz + sj];
        }
      out[i*nz + j] = s;
    }
}

static double rs_strength(const struct fd3d *f, int i, double v, double theta)
{
  return f->a[i]*f->sigma[i]*asinh(v/(2.0*f->V0)*exp((f->f0 + f->b[i]*log(f->V0*theta/f->Lc[i]))/f->a[i]));
}
static double rs_dstrength(const struct fd3d *f, int i, double v, double theta)
{
  double e = exp((f->f0 + f->b[i]*log(f->V0*theta/f->Lc[i]))/f->a[i])/(2.0*f->V0), z = v*e;
  return f->a[i]*f->sigma[i]*e/sqrt(1.0 + z*z);
}
static double theta_adv(const struct fd3d *f, int i, double th, double va, double vb, double dt)
{
  double v = 0.5*(fabs(va) + fabs(vb));
  if(v < 1e-30) v = 1e-30;
  return f->Lc[i]/v + (th - f->Lc[i]/v)*exp(-v*dt/f->Lc[i]);
}

/* local solve at element i: Hp + k V = strength(V) sign(V), scalar */
static double local_sw(const struct fd3d *f, int i, double Hp, double k, double str0, double slope)
{
  double A = -k - slope;
  if(f->locked && f->locked[i]) return 0.0;
  if(fabs(Hp) <= str0) return 0.0;
  return sgn(Hp)*(fabs(Hp) - str0)/A;
}
static double local_rs(const struct fd3d *f, int i, double Hp, double k, double theta, double vguess)
{
  double h = fabs(Hp), A = -k, lo = -80.0, hi, u, un, v, F, g, dg;
  int it;
  if(f->locked && f->locked[i]) return 0.0;
  if(h <= 0.0) return 0.0;
  hi = log(h/A) + 1.0;
  u = log(fabs(vguess) > 1e-30 ? fabs(vguess) : 1e-30);
  if(u > hi - 0.5) u = hi - 0.5;
  if(u < lo + 1.0) u = lo + 1.0;
  for(it = 0; it < 100; it++){
    v = exp(u);
    F = rs_strength(f, i, v, theta);
    g = h - F - A*v;                                   /* decreasing in v */
    dg = -(rs_dstrength(f, i, v, theta) + A)*v;
    if(g > 0) lo = u; else hi = u;
    un = (dg != 0.0) ? u - g/dg : 0.5*(lo + hi);
    if(!isfinite(un) || un <= lo || un >= hi) un = 0.5*(lo + hi);
    if(fabs(un - u) < 1e-13){ u = un; break; }
    u = un;
  }
  return sgn(Hp)*exp(u);
}

/* push one fine-step slice of slip rate into the history ring */
void fd3d_push(struct fd3d *f, const double *V)
{
  struct fd3d_tab *t = f->t;
  field_hat(f, V, t->work);
  to_half(t, t->work, f->hist + (size_t)(f->hk % t->L)*t->Mx*t->Mh);
  f->hk = (f->hk + 1) % t->L;
  if(f->nhist < t->L) f->nhist++;
}

/* empty history, zero slip and time (start of a dynamic phase) */
void fd3d_reset(struct fd3d *f)
{
  struct fd3d_tab *t = f->t;
  struct fd3d_layer *ly = f->layer;
  if(ly && ly->surf){                             /* the surface starts from rest at each phase */
    memset(ly->V3, 0, sizeof(double)*3*ly->ns); memset(ly->Vold3, 0, sizeof(double)*3*ly->ns); memset(ly->delta3, 0, sizeof(double)*3*ly->ns);
    memset(ly->hsx3, 0, sizeof(double _Complex)*3*(size_t)ly->L*ly->Mh*ly->nys);
    memset(ly->hfx, 0, sizeof(double _Complex)*(size_t)ly->L*ly->Mh*ly->nz);
    memset(ly->tn, 0, sizeof(double)*f->N);
    ly->hk = 0; ly->nhist = 0;
  }else if(ly){                                    /* the layer starts from rest at each phase */
    memset(ly->V, 0, sizeof(double)*ly->ns); memset(ly->Vold, 0, sizeof(double)*ly->ns); memset(ly->delta, 0, sizeof(double)*ly->ns);
    if(ly->fft){
      memset(ly->hs2, 0, sizeof(double _Complex)*(size_t)ly->L*ly->Mx2*ly->Mh2);
      memset(ly->hsx, 0, sizeof(double _Complex)*(size_t)ly->L*ly->Mh*ly->nys);
      memset(ly->hfx, 0, sizeof(double _Complex)*(size_t)ly->L*ly->Mh*ly->nz);
    }else{
      memset(ly->hs, 0, sizeof(double)*(size_t)ly->L*ly->ns); memset(ly->hf, 0, sizeof(double)*(size_t)ly->L*f->N);
    }
    ly->hk = 0; ly->nhist = 0;
  }
  memset(f->hist, 0, sizeof(fd3d_cplx)*(size_t)t->L*t->Mx*t->Mh);
  f->hk = 0; f->nhist = 0; f->time = 0.0;
  memset(f->delta, 0, sizeof(double)*f->N);
}

/* one fine step, split in three parts so that the free-surface layer can iterate the
   fault solve with an extra known traction (step_solve's Hx) before finishing */
static void step_begin(struct fd3d *f, const double *load)
{
  struct fd3d_tab *t = f->t;
  int N = f->N, i;
  double eta = 0.5*t->rho*t->cT;
  known(f);
  if(load) for(i = 0; i < N; i++) f->H[i] += load[i];
  if(!f->Vold) f->Vold = (double *)malloc(sizeof(double)*N);
  memcpy(f->Vold, f->V, sizeof(double)*N);
  /* NL: neighbour and self terms are halved (the step's midpoint rate is (V_old + V)/2),
     the V_old half goes into H, the self term becomes kdiag/2 - eta/2 */
  if(f->nl){
    apply_k1(f, f->Vold, f->work);
    for(i = 0; i < N; i++) f->H[i] += 0.5*(f->work[i] + (kdiag(f, i % f->nzp) + eta)*f->Vold[i]);
  }
}

static int step_solve_g(struct fd3d *f, const double *Hx, const double *guess, const char *mask);
static int step_solve(struct fd3d *f, const double *Hx){ return step_solve_g(f, Hx, NULL, NULL); }
/* guess: starting iterate for the Jacobi sweeps (NULL: V_old); the fixed point does not depend on it.
   mask: if given, only these cells are updated (the others keep the guess, as in a converged
   solve whose known part did not change there) */
static int step_solve_g(struct fd3d *f, const double *Hx, const double *guess, const char *mask)
{
  struct fd3d_tab *t = f->t;
  int N = f->N, i, it, pass, itot = 0;
  double dt = t->dt, *Vn = (double *)malloc(sizeof(double)*N), *Vold = f->Vold;
  double *th_p = NULL, *Vs = NULL, eta = 0.5*t->rho*t->cT, hw = f->nl ? 0.5 : 1.0;
  double *H = (double *)malloc(sizeof(double)*N);
  for(i = 0; i < N; i++) H[i] = f->H[i] + (Hx ? Hx[i] : 0.0);
  if(f->vpres){                                   /* prescribed slip rates (tests) */
    memcpy(f->V, f->vpres, sizeof(double)*N);
    free(Vn); free(H);
    return 0;
  }
  if(f->law == 1){ th_p = (double *)malloc(sizeof(double)*N); Vs = (double *)malloc(sizeof(double)*N); }
  for(pass = 0; pass < (f->law == 1 ? 2 : 1); pass++){
    double *V = (double *)malloc(sizeof(double)*N);
    memcpy(V, pass == 0 ? (guess ? guess : Vold) : Vs, sizeof(double)*N);
    if(f->law == 1)
      for(i = 0; i < N; i++) th_p[i] = theta_adv(f, i, f->theta[i], Vold[i], pass == 0 ? Vold[i] : Vs[i], dt);
    for(it = 0; it < 400; it++){
      double ch = 0.0, vm = 0.0;
      apply_k1(f, V, f->work);
      for(i = 0; i < N; i++){
        double Hp;
        if(mask && !mask[i]){ Vn[i] = V[i]; continue; }
        Hp = H[i] + hw*f->work[i];
        double kd = f->nl ? 0.5*kdiag(f, i % f->nzp) - 0.5*eta : kdiag(f, i % f->nzp);
        if(f->law == 0){
          /* slip at the collocation time: delta + dt V (first order) or delta + dt (V_old + V)/2 */
          /* slip distance: net slip, or (path = 1) the accumulated slip path (TPV10: backward motion
             adds to it); strength c + mu max(sigma, 0) (TPV10: tension acts as zero normal stress) */
          double dm = f->path ? f->dpath[i] + (f->nl ? 0.5*dt*fabs(Vold[i]) : 0.0) : fabs(f->delta[i] + (f->nl ? 0.5*dt*Vold[i] : 0.0));
          double sg = f->sigma[i] > 0.0 ? f->sigma[i] : 0.0, coh = f->coh ? f->coh[i] : 0.0;
          double mu = dm < f->d0 ? f->mu_s - (f->mu_s - f->mu_d)*dm/f->d0 : f->mu_d;
          double slope = dm < f->d0 ? sg*(f->mu_s - f->mu_d)*hw*dt/f->d0 : 0.0;
          Vn[i] = local_sw(f, i, Hp, kd, coh + sg*mu, slope);
          if(slope > 0.0){
            /* the strength falls linearly with |V| only until the slip distance reaches d0 (at
               |V| = vc); beyond, it stays at the dynamic level. If the step weakens faster than
               the radiation term can balance (-kd <= slope; the linearised solve is then ill posed),
               or (sw_saturate = 1) the solution passes vc, take the saturated branch (exact for the
               piecewise-linear strength; the prototype and the earlier tests use the linearised
               solve, so the second case is optional) */
            double vc = (f->d0 - dm)/(hw*dt), sd = coh + sg*f->mu_d;
            if(-kd - slope <= 0.0 || (f->sw_saturate && fabs(Vn[i]) > vc)){
              double vs = (fabs(Hp) > sd) ? (fabs(Hp) - sd)/(-kd) : 0.0;
              if(fabs(Hp) > coh + sg*mu) Vn[i] = sgn(Hp)*(vs > vc ? vs : vc);
            }
          }
        }else
          Vn[i] = local_rs(f, i, Hp, kd, th_p[i], V[i]);
        if(fabs(Vn[i] - V[i]) > ch) ch = fabs(Vn[i] - V[i]);
        if(fabs(Vn[i]) > vm) vm = fabs(Vn[i]);
      }
      memcpy(V, Vn, sizeof(double)*N);
      if(ch <= 1e-13*(vm > 1.0 ? vm : 1.0)) break;     /* Jacobi converged (same rule as the prototype) */
    }
    itot += it + 1;
    if(f->law == 1 && pass == 0) memcpy(Vs, V, sizeof(double)*N);
    else{
      memcpy(f->V, V, sizeof(double)*N);
      if(f->law == 1) memcpy(f->thnew ? f->thnew : f->theta, th_p, sizeof(double)*N);
    }
    free(V);
  }
  free(Vn); free(th_p); free(Vs); free(H);
  return itot;
}

static void step_end(struct fd3d *f)
{
  struct fd3d_tab *t = f->t;
  int N = f->N, i;
  double dt = t->dt, eta = 0.5*t->rho*t->cT, *Vold = f->Vold;
  apply_k1(f, f->V, f->work);
  if(f->nl){
    double *Vm = (double *)calloc((size_t)N, sizeof(double));   /* calloc: silences -Wmaybe-uninitialized */
    for(i = 0; i < N; i++){
      f->tau[i] = f->H[i] + 0.5*f->work[i] + (0.5*kdiag(f, i % f->nzp) - 0.5*eta)*f->V[i];
      Vm[i] = 0.5*(Vold[i] + f->V[i]);                    /* midpoint rate of the step */
    }
    fd3d_push(f, Vm);
    for(i = 0; i < N; i++) f->delta[i] += Vm[i]*dt;
    if(f->path) for(i = 0; i < N; i++) f->dpath[i] += 0.5*(fabs(Vold[i]) + fabs(f->V[i]))*dt;
    free(Vm);
  }else{
    for(i = 0; i < N; i++) f->tau[i] = f->H[i] + f->work[i] + kdiag(f, i % f->nzp)*f->V[i];
    fd3d_push(f, f->V);
    for(i = 0; i < N; i++) f->delta[i] += f->V[i]*dt;
    if(f->path) for(i = 0; i < N; i++) f->dpath[i] += fabs(f->V[i])*dt;
  }
  f->time += dt;
}

static int fd3d_step_layer_fft(struct fd3d *f, const double *load);
static int fd3d_step_surf(struct fd3d *f, const double *load);
int fd3d_step(struct fd3d *f, const double *load)
{
  int it;
  if(f->layer && f->layer->surf) return fd3d_step_surf(f, load);
  if(f->layer && f->layer->fft) return fd3d_step_layer_fft(f, load);
  if(f->layer) return fd3d_step_layer(f, load);
  step_begin(f, load);
  it = step_solve(f, NULL);
  step_end(f);
  return it;
}

/* ---------------- variable substeps (k = 2, 4, ..., 64) ----------------
   A step of k fine substeps: the current step's own contribution is the operator
   K_k = sum_{l<=k} K[l] (first order, constant rate) or, for second-order marching,
   B_k V + C_k V_old with the rate linear across the step (weights w(l) = (k-l+1/2)/k,
   eta on the diagonal as in 2-D). These operators reach well beyond the 3 x 3
   stencil, so the implicit solve iterates: each sweep applies the full operator by
   FFT and solves locally with its self coefficient. Meant for slow phases (rate
   and state), where the rate dependence of friction dominates the self term. */

static int ilog2i(int k){ int i = 0; while((1 << i) < k) i++; return i; }

static void fd3d_kops(struct fd3d *f, int k)
{
  struct fd3d_tab *t = f->t;
  int ik = ilog2i(k), l, ii, jj, mirror = f->mirror;
  size_t Mhh = (size_t)t->Mx*t->Mh, M = (size_t)t->Mx*t->Mz, m;
  double eta = 0.5*t->rho*t->cT, w1 = (k - 0.5)/k;
  if(f->kop_B[ik]) return;
  f->kop_B[ik] = (double _Complex *)calloc(Mhh, sizeof(double _Complex));
  f->kop_C[ik] = (double _Complex *)calloc(Mhh, sizeof(double _Complex));
  for(l = 1; l <= k && l <= t->L; l++){
    const fd3d_cplx *W = t->W_hat + (size_t)(l - 1)*Mhh;
    double w = f->nl ? (k - l + 0.5)/k : 1.0;
    for(m = 0; m < Mhh; m++){ f->kop_B[ik][m] += w*W[m]; f->kop_C[ik][m] += (1.0 - w)*W[m]; }
  }
  for(ii = 0; ii < t->Mx; ii++) for(jj = 0; jj < t->Mh; jj++){
    double _Complex ks = t->Ks_hat[(size_t)ii*t->Mz + jj]*t->dt;
    size_t q = (size_t)ii*t->Mh + jj;
    if(f->nl){ f->kop_B[ik][q] += 0.5*k*ks + (w1 - 1.0)*eta; f->kop_C[ik][q] += 0.5*k*ks + (1.0 - w1)*eta; }
    else f->kop_B[ik][q] += k*ks;
  }
  /* self coefficients (offset 0; with the mirror the top row also sees its image at +1 row) */
  memset(t->work, 0, sizeof(double _Complex)*M);
  for(ii = 0; ii < t->Mx; ii++) for(jj = 0; jj < t->Mh; jj++) t->work[(size_t)ii*t->Mz + jj] = f->kop_B[ik][(size_t)ii*t->Mh + jj];
  from_half(t, t->work);
  fd3d_fft2(t, t->work, +1);
  f->kop_self[ik] = creal(t->work[0])/M;
  f->kop_self_top[ik] = f->kop_self[ik] + (mirror ? creal(t->work[1])/M : 0.0);
}

/* operator (half spectrum) applied to a field: out = op * x (element order) */
static void kop_apply(struct fd3d *f, const double _Complex *opq, const double *x, double *out)
{
  struct fd3d_tab *t = f->t;
  size_t M = (size_t)t->Mx*t->Mz;
  int ii, jj, i, j, row0 = f->mirror ? f->nzp : f->next;
  field_hat(f, x, t->work);
  for(ii = 0; ii < t->Mx; ii++) for(jj = 0; jj < t->Mh; jj++) t->work[(size_t)ii*t->Mz + jj] *= opq[(size_t)ii*t->Mh + jj];
  from_half(t, t->work);
  fd3d_fft2(t, t->work, +1);
  for(i = 0; i < t->nx; i++) for(j = 0; j < f->nzp; j++) out[i*f->nzp + j] = creal(t->work[(size_t)i*t->Mz + row0 + j])/M;
}

/* push k fine slices: constant V (first order) or linear from Vold to V (second order) */
static void push_k_lin(struct fd3d *f, const double *Vold, const double *V, int k)
{
  struct fd3d_tab *t = f->t;
  size_t Mhh = (size_t)t->Mx*t->Mh, m;
  int q;
  double _Complex *a = (double _Complex *)malloc(sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  field_hat(f, Vold, a);
  field_hat(f, V, t->work);
  for(q = 1; q <= k; q++){
    double w = f->nl ? (q - 0.5)/k : 1.0;
    fd3d_cplx *slot = f->hist + (size_t)(f->hk % t->L)*Mhh;
    int ii, jj;
    for(ii = 0; ii < t->Mx; ii++) for(jj = 0; jj < t->Mh; jj++){
      size_t fq = (size_t)ii*t->Mz + jj;
      slot[(size_t)ii*t->Mh + jj] = (fd3d_cplx)((1.0 - w)*a[fq] + w*t->work[fq]);
    }
    f->hk = (f->hk + 1) % t->L;
    if(f->nhist < t->L) f->nhist++;
  }
  (void)m;
  free(a);
}

int fd3d_step_k(struct fd3d *f, const double *load, int k)
{
  struct fd3d_tab *t = f->t;
  int N = f->N, i, it, pass, itot = 0, ik;
  double dt = t->dt*k, *Vold, *Hc, *Vit, *Vn, *conv, *th_p = NULL, *Vs = NULL;
  if(k <= 1 && !f->generic) return fd3d_step(f, load);
  if(k > 64 || (k & (k - 1))) return -1;
  ik = ilog2i(k);
  fd3d_kops(f, k);
  Vold = (double *)malloc(sizeof(double)*N); Hc = (double *)malloc(sizeof(double)*N);
  Vit = (double *)malloc(sizeof(double)*N); Vn = (double *)malloc(sizeof(double)*N); conv = (double *)malloc(sizeof(double)*N);
  memcpy(Vold, f->V, sizeof(double)*N);
  known_k(f, k);
  if(load) for(i = 0; i < N; i++) f->H[i] += load[i];
  memcpy(Hc, f->H, sizeof(double)*N);
  if(f->nl){
    kop_apply(f, f->kop_C[ik], Vold, conv);
    for(i = 0; i < N; i++) Hc[i] += conv[i];
  }
  if(f->law == 1){ th_p = (double *)malloc(sizeof(double)*N); Vs = (double *)malloc(sizeof(double)*N); }
  for(pass = 0; pass < (f->law == 1 ? 2 : 1); pass++){
    memcpy(Vit, pass == 0 ? Vold : Vs, sizeof(double)*N);
    if(f->law == 1)
      for(i = 0; i < N; i++) th_p[i] = theta_adv(f, i, f->theta[i], Vold[i], pass == 0 ? Vold[i] : Vs[i], dt);
    for(it = 0; it < 500; it++){
      double ch = 0.0, vm = 0.0;
      kop_apply(f, f->kop_B[ik], Vit, conv);
      for(i = 0; i < N; i++){
        int row = i % f->nzp;
        double kd = (f->mirror && row == 0) ? f->kop_self_top[ik] : f->kop_self[ik];
        double Hp = Hc[i] + conv[i] - kd*Vit[i];
        if(f->law == 0){
          double d0s = f->delta[i] + (f->nl ? 0.5*dt*Vold[i] : 0.0), dm = fabs(d0s);
          double mu = dm < f->d0 ? f->mu_s - (f->mu_s - f->mu_d)*dm/f->d0 : f->mu_d;
          double slope = dm < f->d0 ? f->sigma[i]*(f->mu_s - f->mu_d)*(f->nl ? 0.5 : 1.0)*dt/f->d0 : 0.0;
          Vn[i] = local_sw(f, i, Hp, kd, f->sigma[i]*mu, slope);
        }else
          Vn[i] = local_rs(f, i, Hp, kd, th_p[i], Vit[i]);
        if(fabs(Vn[i] - Vit[i]) > ch) ch = fabs(Vn[i] - Vit[i]);
        if(fabs(Vn[i]) > vm) vm = fabs(Vn[i]);
      }
      /* damped Jacobi: the long-step operator is not diagonally dominant */
      for(i = 0; i < N; i++) Vit[i] += f->kdamp*(Vn[i] - Vit[i]);
      if(!isfinite(ch)) break;
      if(ch <= 1e-12*(vm > 1e-9 ? vm : 1e-9)) break;
    }
    itot += it + 1;
    {
      /* a long step must not change the peak slip rate much (the step selector aims at
         10 percent); a jump by more than a factor of two means the step was too long or
         the iteration settled on a spurious solution */
      double vo = 0.0, vn = 0.0;
      for(i = 0; i < N; i++){ if(fabs(Vold[i]) > vo) vo = fabs(Vold[i]); if(fabs(Vit[i]) > vn) vn = fabs(Vit[i]); }
      if(vn > 2.0*vo + 1e-12) it = 500;
    }
    if(getenv("FD3D_DEBUG") && (it >= 500 || !isfinite(Vit[0]))){
      double vo = 0.0, vn = 0.0; int io = 0, in_ = 0;
      for(i = 0; i < N; i++){ if(fabs(Vold[i]) > vo){ vo = fabs(Vold[i]); io = i; } if(fabs(Vit[i]) > vn){ vn = fabs(Vit[i]); in_ = i; } }
      fprintf(stderr, "fd3d_step_k: k %d pass %d rejected after %d its: vmax old %.3e (cell %d) new %.3e (cell %d)\n", k, pass, it, vo, io, vn, in_);
    }
    if(it >= 500 || !isfinite(Vit[0])){
      /* not converged: leave the state untouched, the caller retries with fewer substeps */
      free(Vold); free(Hc); free(Vit); free(Vn); free(conv); free(th_p); free(Vs);
      return -2;
    }
    if(f->law == 1 && pass == 0) memcpy(Vs, Vit, sizeof(double)*N);
  }
  memcpy(f->V, Vit, sizeof(double)*N);
  if(f->law == 1) memcpy(f->theta, th_p, sizeof(double)*N);
  kop_apply(f, f->kop_B[ik], f->V, conv);
  for(i = 0; i < N; i++) f->tau[i] = Hc[i] + conv[i];
  push_k_lin(f, Vold, f->V, k);
  for(i = 0; i < N; i++) f->delta[i] += (f->nl ? 0.5*(Vold[i] + f->V[i]) : f->V[i])*dt;
  f->time += dt;
  free(Vold); free(Hc); free(Vit); free(Vn); free(conv); free(th_p); free(Vs);
  return itot;
}


/* ---------------- free-surface layer (stage A, dense first version) ---------------- */

/* frames: fault e1 = x, e2 = z (up), e3 = -y; layer e1 = x, e2 = y, e3 = z */
static void lay_resp(const struct fd3d_tab *t, int src_fault, double sx, double sy, double sz, double a1, double a2,
                     double rx, double ry, double rz, int rcv_fault, int L, double *K, double *ks)
{
  double s[6], T1, T2, v, vprev = 0.0, late[2];
  double x1, x2, x3; int l, q;
  if(src_fault){ x1 = rx - sx; x2 = rz - sz; x3 = -(ry - sy); }
  else{ x1 = rx - sx; x2 = ry - sy; x3 = rz - sz; }
  T1 = 3.0*(sqrt((rx - sx)*(rx - sx) + (ry - sy)*(ry - sy) + (rz - sz)*(rz - sz)) + 4.0*(a1 + a2))/t->cT + 10.0*t->dt;
  T2 = 2.0*T1;
  for(q = 0; q < 2; q++){
    fd3d_rect_ramp_tensor(x1, x2, x3, q ? T2 : T1, a1, a2, t->cL, t->cT, t->rho, src_fault ? 1 : 3, s);
    /* traction: fault receiver -> strike traction e1F.sigma.e3F; layer receiver -> sigma_zz */
    if(src_fault) late[q] = rcv_fault ? s[4] : s[1];       /* fault frame: s31 (= e1 . sigma . e3), s22 = sigma_zz */
    else late[q] = rcv_fault ? -s[3] : s[2];               /* layer frame: -s12 = e1 . sigma . (0,-1,0), s33 */
  }
  *ks = (late[1] - late[0])/(T2 - T1);
  for(l = 1; l <= L; l++){
    fd3d_rect_ramp_tensor(x1, x2, x3, l*t->dt, a1, a2, t->cL, t->cT, t->rho, src_fault ? 1 : 3, s);
    if(src_fault) v = rcv_fault ? s[4] : s[1];
    else v = rcv_fault ? -s[3] : s[2];
    K[l - 1] = v - vprev; vprev = v;
  }
}

/* dense LU without pivoting problems for the small layer operator: partial pivoting */
static int lu_factor(int n, double *A, int *piv)
{
  int i, j, k, p;
  for(k = 0; k < n; k++){
    double m = 0.0; p = k;
    for(i = k; i < n; i++) if(fabs(A[(size_t)i*n + k]) > m){ m = fabs(A[(size_t)i*n + k]); p = i; }
    if(m == 0.0) return -1;
    piv[k] = p;
    if(p != k) for(j = 0; j < n; j++){ double tmp = A[(size_t)k*n + j]; A[(size_t)k*n + j] = A[(size_t)p*n + j]; A[(size_t)p*n + j] = tmp; }
    for(i = k + 1; i < n; i++){
      double fct = A[(size_t)i*n + k]/A[(size_t)k*n + k];
      A[(size_t)i*n + k] = fct;
      for(j = k + 1; j < n; j++) A[(size_t)i*n + j] -= fct*A[(size_t)k*n + j];
    }
  }
  return 0;
}
static void lu_solve(int n, const double *A, const int *piv, double *b)
{
  int i, j;
  for(i = 0; i < n; i++){ if(piv[i] != i){ double tmp = b[i]; b[i] = b[piv[i]]; b[piv[i]] = tmp; } }
  for(i = 0; i < n; i++) for(j = 0; j < i; j++) b[i] -= A[(size_t)i*n + j]*b[j];
  for(i = n - 1; i >= 0; i--){ for(j = i + 1; j < n; j++) b[i] -= A[(size_t)i*n + j]*b[j]; b[i] /= A[(size_t)i*n + i]; }
}

struct fd3d_layer *fd3d_layer_create(struct fd3d *f, double x0f, int nxs, int nys, double hx, double hy,
                                     double xs0, double ys0, int L)
{
  const struct fd3d_tab *t = f->t;
  struct fd3d_layer *ly;
  int N = f->N, ns = nxs*nys, a, b, l, nz = f->nzp;
  double dx = t->dx, *K = (double *)malloc(sizeof(double)*L), ks, *K2 = (double *)malloc(sizeof(double)*L), ks2;
  if(!f->mirror) return NULL;                       /* the layer completes the mirror image */
  ly = (struct fd3d_layer *)calloc(1, sizeof(*ly));
  ly->ns = ns; ly->nxs = nxs; ly->nys = nys; ly->L = L; ly->hx = hx; ly->hy = hy; ly->xs0 = xs0; ly->ys0 = ys0; ly->x0f = x0f;
  ly->eta_p = 0.5*t->rho*t->cL;
  ly->Wss = calloc((size_t)L*ns*ns, sizeof(double)); ly->Kss = calloc((size_t)ns*ns, sizeof(double));
  ly->Wfs = calloc((size_t)L*ns*N, sizeof(double)); ly->Kfs = calloc((size_t)ns*N, sizeof(double));
  ly->Wsf = calloc((size_t)L*N*ns, sizeof(double)); ly->Ksf = calloc((size_t)N*ns, sizeof(double));
  for(a = 0; a < ns; a++){                          /* sources: layer elements */
    double sx = xs0 + (a/nys + 0.5)*hx, sy = ys0 + (a%nys + 0.5)*hy;
    for(b = 0; b < ns; b++){
      double rx = xs0 + (b/nys + 0.5)*hx, ry = ys0 + (b%nys + 0.5)*hy;
      lay_resp(t, 0, sx, sy, 0.0, 0.5*hx, 0.5*hy, rx, ry, 0.0, 0, L, K, &ks);
      ly->Kss[(size_t)b*ns + a] = ks;
      for(l = 0; l < L; l++) ly->Wss[((size_t)l*ns + b)*ns + a] = K[l] - ks*t->dt;
    }
    for(b = 0; b < N; b++){
      double rx = x0f + (b/nz + 0.5)*dx, rz = -(b%nz + 0.5)*dx;
      lay_resp(t, 0, sx, sy, 0.0, 0.5*hx, 0.5*hy, rx, 0.0, rz, 1, L, K, &ks);
      ly->Ksf[(size_t)b*ns + a] = ks;
      for(l = 0; l < L; l++) ly->Wsf[((size_t)l*N + b)*ns + a] = K[l] - ks*t->dt;
    }
  }
  for(a = 0; a < N; a++){                           /* sources: fault elements plus their images */
    double sx = x0f + (a/nz + 0.5)*dx, sz = -(a%nz + 0.5)*dx;
    for(b = 0; b < ns; b++){
      double rx = xs0 + (b/nys + 0.5)*hx, ry = ys0 + (b%nys + 0.5)*hy;
      lay_resp(t, 1, sx, 0.0, sz, 0.5*dx, 0.5*dx, rx, ry, 0.0, 0, L, K, &ks);
      lay_resp(t, 1, sx, 0.0, -sz, 0.5*dx, 0.5*dx, rx, ry, 0.0, 0, L, K2, &ks2);
      ly->Kfs[(size_t)b*N + a] = ks + ks2;
      for(l = 0; l < L; l++) ly->Wfs[((size_t)l*ns + b)*N + a] = K[l] + K2[l] - (ks + ks2)*t->dt;
    }
  }
  free(K); free(K2);
  ly->K1ss = malloc(sizeof(double)*(size_t)ns*ns); ly->K1fs = malloc(sizeof(double)*(size_t)ns*N); ly->K1sf = malloc(sizeof(double)*(size_t)N*ns);
  for(a = 0; a < ns*ns; a++) ly->K1ss[a] = ly->Wss[a] + ly->Kss[a]*t->dt;
  for(a = 0; a < ns*N; a++){ ly->K1fs[a] = ly->Wfs[a] + ly->Kfs[a]*t->dt; ly->K1sf[a] = ly->Wsf[a] + ly->Ksf[a]*t->dt; }
  /* the layer's self operator for the current step: first order K1ss, NL (K1ss + E)/2 - E */
  ly->Bss_lu = malloc(sizeof(double)*(size_t)ns*ns); ly->piv = malloc(sizeof(int)*ns);
  for(a = 0; a < ns*ns; a++) ly->Bss_lu[a] = f->nl ? 0.5*ly->K1ss[a] : ly->K1ss[a];
  if(f->nl) for(a = 0; a < ns; a++) ly->Bss_lu[(size_t)a*ns + a] -= 0.5*ly->eta_p;
  if(lu_factor(ns, ly->Bss_lu, ly->piv) != 0){ fd3d_layer_free(ly); return NULL; }
  ly->V = calloc(ns, sizeof(double)); ly->Vold = calloc(ns, sizeof(double)); ly->delta = calloc(ns, sizeof(double));
  ly->H = calloc(ns, sizeof(double)); ly->tau = calloc(ns, sizeof(double));
  ly->hs = calloc((size_t)L*ns, sizeof(double)); ly->hf = calloc((size_t)L*N, sizeof(double));
  return ly;
}

static void surf_free(struct fd3d_layer *ly);
void fd3d_layer_free(struct fd3d_layer *ly)
{
  if(!ly) return;
  if(ly->surf){ surf_free(ly); return; }
  if(ly->fft){
    free(ly->Wss2); free(ly->Kss2); free(ly->K1ss2); free(ly->Wfs_h); free(ly->Kfs_h); free(ly->K1fs_h);
    free(ly->Wsf_h); free(ly->Ksf_h); free(ly->K1sf_h); free(ly->hs2); free(ly->hsx); free(ly->hfx);
    fd_fft_plan_destroy(ly->pf); fd_fft_plan_destroy(ly->pb);
    fd_fft_plan_destroy(ly->pf2x); fd_fft_plan_destroy(ly->pb2x); fd_fft_plan_destroy(ly->pf2y); fd_fft_plan_destroy(ly->pb2y);
  }
  free(ly->Wss); free(ly->Kss); free(ly->Wfs); free(ly->Kfs); free(ly->Wsf); free(ly->Ksf);
  free(ly->K1ss); free(ly->K1fs); free(ly->K1sf); free(ly->Bss_lu); free(ly->piv);
  free(ly->V); free(ly->Vold); free(ly->delta); free(ly->H); free(ly->tau); free(ly->hs); free(ly->hf);
  free(ly);
}

static void mv(int nr, int nc, const double *A, const double *x, double w, double *y)   /* y += w A x */
{
  int i, j;
  for(i = 0; i < nr; i++){ double s = 0.0; for(j = 0; j < nc; j++) s += A[(size_t)i*nc + j]*x[j]; y[i] += w*s; }
}

int fd3d_step_layer(struct fd3d *f, const double *load)
{
  struct fd3d_layer *ly = f->layer;
  const struct fd3d_tab *t = f->t;
  int N = f->N, ns = ly->ns, i, l, it, itot = 0, lmax = ly->L < ly->nhist + 1 ? ly->L : ly->nhist + 1;
  double dt = t->dt, hw = f->nl ? 0.5 : 1.0;
  double *Hx = calloc(N, sizeof(double)), *Hs = malloc(sizeof(double)*ns), *Vs = malloc(sizeof(double)*ns), *Vp = malloc(sizeof(double)*N);
  double *thbuf = (f->law == 1) ? malloc(sizeof(double)*N) : NULL;
  step_begin(f, load);
  memcpy(ly->Vold, ly->V, sizeof(double)*ns);
  /* known parts: histories (lags >= 2), statics, the V_old half of the current step (NL) */
  for(i = 0; i < ns; i++) ly->H[i] = 0.0;
  mv(ns, ns, ly->Kss, ly->delta, 1.0, ly->H);
  mv(ns, N, ly->Kfs, f->delta, 1.0, ly->H);
  if(!f->static_fn) mv(N, ns, ly->Ksf, ly->delta, 1.0, Hx);   /* with the exact Is the layer's static part is in Is already */
  for(l = 2; l <= lmax; l++){
    int slot = (ly->hk - (l - 1) + 2*ly->L) % ly->L;
    mv(ns, ns, ly->Wss + (size_t)(l - 1)*ns*ns, ly->hs + (size_t)slot*ns, 1.0, ly->H);
    mv(ns, N, ly->Wfs + (size_t)(l - 1)*ns*N, ly->hf + (size_t)slot*N, 1.0, ly->H);
    mv(N, ns, ly->Wsf + (size_t)(l - 1)*N*ns, ly->hs + (size_t)slot*ns, 1.0, Hx);
  }
  if(f->nl){
    mv(ns, ns, ly->K1ss, ly->Vold, 0.5, ly->H);
    for(i = 0; i < ns; i++) ly->H[i] += 0.5*ly->eta_p*ly->Vold[i];
    mv(ns, N, ly->K1fs, f->Vold, 0.5, ly->H);
    mv(N, ns, ly->K1sf, ly->Vold, 0.5, Hx);
  }
  /* iterate fault and layer: the fault sees the layer through the first lag, and vice versa */
  memcpy(Vs, ly->Vold, sizeof(double)*ns);
  f->thnew = thbuf;
  for(it = 0; it < 200; it++){
    double ch = 0.0, vm = 1e-30;
    double *Hf = malloc(sizeof(double)*N);
    memcpy(Hf, Hx, sizeof(double)*N);
    mv(N, ns, ly->K1sf, Vs, hw, Hf);
    memcpy(Vp, f->V, sizeof(double)*N);
    itot += step_solve(f, Hf);
    free(Hf);
    for(i = 0; i < ns; i++) Hs[i] = -ly->H[i];
    mv(ns, N, ly->K1fs, f->V, -hw, Hs);
    lu_solve(ns, ly->Bss_lu, ly->piv, Hs);
    for(i = 0; i < ns; i++){ if(fabs(Hs[i] - Vs[i]) > ch) ch = fabs(Hs[i] - Vs[i]); if(fabs(Hs[i]) > vm) vm = fabs(Hs[i]); }
    for(i = 0; i < N; i++){ if(fabs(f->V[i] - Vp[i]) > ch) ch = fabs(f->V[i] - Vp[i]); if(fabs(f->V[i]) > vm) vm = fabs(f->V[i]); }
    memcpy(Vs, Hs, sizeof(double)*ns);
    if(it > 0 && ch <= 1e-12*vm) break;
  }
  f->thnew = NULL;
  if(thbuf){ memcpy(f->theta, thbuf, sizeof(double)*N); free(thbuf); }
  /* finish: the fault with the layer's contribution in its known part, then the layer */
  memcpy(ly->V, Vs, sizeof(double)*ns);
  mv(N, ns, ly->K1sf, ly->V, hw, Hx);
  for(i = 0; i < N; i++) f->H[i] += Hx[i];
  {
    int slot = ly->hk % ly->L;
    for(i = 0; i < ns; i++){
      double vm_ = f->nl ? 0.5*(ly->Vold[i] + ly->V[i]) : ly->V[i];
      ly->hs[(size_t)slot*ns + i] = vm_; ly->delta[i] += vm_*dt;
    }
    for(i = 0; i < N; i++) ly->hf[(size_t)slot*N + i] = f->nl ? 0.5*(f->Vold[i] + f->V[i]) : f->V[i];
    ly->hk = (ly->hk + 1) % ly->L;
    if(ly->nhist < ly->L) ly->nhist++;
  }
  step_end(f);
  /* layer normal traction at the end of the step (should be zero) */
  for(i = 0; i < ns; i++) ly->tau[i] = 0.0;
  free(Hx); free(Hs); free(Vs); free(Vp);
  return itot;
}


/* ---------------- free-surface layer, FFT-structured ---------------- */
#include <time.h>
static double ly_now(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }
double fd3d_layer_prof[4];                         /* known part, fault solves, layer solves, finish */

/* row FFTs along strike: in[row][x] (real, nrow rows of nin values placed at stride st from 0) ->
   out[k][row] (k < Mh) */
static void lx_fwd(const struct fd3d_layer *ly, int nrow, int nin, int st, const double *in, int in_rowmajor_x_fast,
                   double _Complex *out)
{
  double _Complex *a = malloc(sizeof(double _Complex)*ly->Mx), *b = malloc(sizeof(double _Complex)*ly->Mx);
  int r, x, k;
  for(r = 0; r < nrow; r++){
    memset(a, 0, sizeof(double _Complex)*ly->Mx);
    for(x = 0; x < nin; x++) a[x*st] = in_rowmajor_x_fast ? in[(size_t)r*nin + x] : in[(size_t)x*nrow + r];
    fd_fft_exec(ly->pf, a, b);
    for(k = 0; k < ly->Mh; k++) out[(size_t)k*nrow + r] = b[k];
  }
  free(a); free(b);
}
/* inverse: in[k][row] (half spectrum) -> out real values at x = st*i, i < nout, stored [i][row] (element order) */
static void lx_inv(const struct fd3d_layer *ly, int nrow, const double _Complex *in, int nout, int st, double w, double *out)
{
  double _Complex *a = malloc(sizeof(double _Complex)*ly->Mx), *b = malloc(sizeof(double _Complex)*ly->Mx);
  int r, i, k, M = ly->Mx;
  for(r = 0; r < nrow; r++){
    for(k = 0; k < ly->Mh; k++) a[k] = in[(size_t)k*nrow + r];
    for(k = ly->Mh; k < M; k++) a[k] = conj(a[M - k]);
    fd_fft_exec(ly->pb, a, b);
    for(i = 0; i < nout; i++) out[(size_t)i*nrow + r] += w*creal(b[(size_t)i*st])/M;
  }
  free(a); free(b);
}
/* Y[k][nr] += T[k][nr][nc] X[k][nc] */
static void lx_apply(const struct fd3d_layer *ly, int nr, int nc, const double _Complex *T, const double _Complex *X, double _Complex *Y)
{
  int k, r, c;
  for(k = 0; k < ly->Mh; k++){
    const double _Complex *Tk = T + (size_t)k*nr*nc, *Xk = X + (size_t)k*nc;
    double _Complex *Yk = Y + (size_t)k*nr;
    for(r = 0; r < nr; r++){ double _Complex s = 0.0; for(c = 0; c < nc; c++) s += Tk[(size_t)r*nc + c]*Xk[c]; Yk[r] += s; }
  }
}
/* 2-D layer FFT: x[p][q] (element order) -> half spectrum (Mx2 x Mh2) */
static void l2_fwd(const struct fd3d_layer *ly, const double *x, double _Complex *out)
{
  int Mx = ly->Mx2, My = ly->My2, i, j;
  double _Complex *a = calloc((size_t)Mx*My, sizeof(double _Complex)), *col = malloc(sizeof(double _Complex)*(Mx > My ? Mx : My)),
    *o = malloc(sizeof(double _Complex)*(Mx > My ? Mx : My));
  for(i = 0; i < ly->nxs; i++) for(j = 0; j < ly->nys; j++) a[(size_t)i*My + j] = x[(size_t)i*ly->nys + j];
  for(i = 0; i < Mx; i++){ fd_fft_exec(ly->pf2y, a + (size_t)i*My, o); memcpy(a + (size_t)i*My, o, sizeof(double _Complex)*My); }
  for(j = 0; j < ly->Mh2; j++){
    for(i = 0; i < Mx; i++) col[i] = a[(size_t)i*My + j];
    fd_fft_exec(ly->pf2x, col, o);
    for(i = 0; i < Mx; i++) out[(size_t)i*ly->Mh2 + j] = o[i];
  }
  free(a); free(col); free(o);
}
/* inverse of a half spectrum product into element order (adds w times) */
static void l2_inv(const struct fd3d_layer *ly, const double _Complex *in, double w, double *x)
{
  int Mx = ly->Mx2, My = ly->My2, i, j;
  double _Complex *a = calloc((size_t)Mx*My, sizeof(double _Complex)), *col = malloc(sizeof(double _Complex)*(Mx > My ? Mx : My)),
    *o = malloc(sizeof(double _Complex)*(Mx > My ? Mx : My));
  for(i = 0; i < Mx; i++) for(j = 0; j < ly->Mh2; j++) a[(size_t)i*My + j] = in[(size_t)i*ly->Mh2 + j];
  for(i = 0; i < Mx; i++) for(j = ly->Mh2; j < My; j++) a[(size_t)i*My + j] = conj(a[(size_t)((Mx - i) % Mx)*My + (My - j)]);
  for(j = 0; j < My; j++){
    for(i = 0; i < Mx; i++) col[i] = a[(size_t)i*My + j];
    fd_fft_exec(ly->pb2x, col, o);
    for(i = 0; i < Mx; i++) a[(size_t)i*My + j] = o[i];
  }
  for(i = 0; i < ly->nxs; i++){
    fd_fft_exec(ly->pb2y, a + (size_t)i*My, o);
    for(j = 0; j < ly->nys; j++) x[(size_t)i*ly->nys + j] += w*creal(o[j])/((double)Mx*My);
  }
  free(a); free(col); free(o);
}

struct fd3d_layer *fd3d_layer_create_fft(struct fd3d *f, double x0f, int m, int nxs, int nys, double xs0, int L)
{
  const struct fd3d_tab *t = f->t;
  struct fd3d_layer *ly;
  int nx = t->nx, nz = f->nzp, d, dmin, dmax, q, j, l, dp, dq, k;
  double dx = t->dx, hx = dx/m, *K = malloc(sizeof(double)*L), *K2 = malloc(sizeof(double)*L), ks, ks2;
  double _Complex *buf;
  size_t Mh, s1;
  if(!f->mirror || m < 1) return NULL;
  ly = calloc(1, sizeof(*ly));
  ly->fft = 1; ly->m = m; ly->nxs = nxs; ly->nys = nys; ly->ns = nxs*nys; ly->L = L; ly->hx = hx; ly->hy = hx;
  ly->xs0 = xs0; ly->ys0 = -0.5*nys*hx; ly->x0f = x0f; ly->nx = nx; ly->nz = nz; ly->eta_p = 0.5*t->rho*t->cL;
  ly->c0 = xs0 + 0.5*hx - x0f - 0.5*dx;
  ly->Mx = fd_fft_size_for(nxs + m*(nx - 1) + 1); ly->Mh = ly->Mx/2 + 1; Mh = ly->Mh;
  ly->Mx2 = fd_fft_size_for(2*nxs - 1); ly->My2 = fd_fft_size_for(2*nys - 1); ly->Mh2 = ly->My2/2 + 1;
  ly->pf = fd_fft_plan_create(ly->Mx, FD_FFT_FORWARD); ly->pb = fd_fft_plan_create(ly->Mx, FD_FFT_BACKWARD);
  ly->pf2x = fd_fft_plan_create(ly->Mx2, FD_FFT_FORWARD); ly->pb2x = fd_fft_plan_create(ly->Mx2, FD_FFT_BACKWARD);
  ly->pf2y = fd_fft_plan_create(ly->My2, FD_FFT_FORWARD); ly->pb2y = fd_fft_plan_create(ly->My2, FD_FFT_BACKWARD);
  /* layer-layer: offsets (dp, dq), kernel arrays in FFT order, then 2-D half spectra */
  {
    size_t M2 = (size_t)ly->Mx2*ly->My2, H2 = (size_t)ly->Mx2*ly->Mh2;
    double *A = calloc((size_t)(L + 2)*M2, sizeof(double));       /* [l][M2]: K lags 1..L, then Ks, then K1 */
    double *tmp = malloc(sizeof(double)*ly->ns);
    ly->Wss2 = malloc(sizeof(double _Complex)*(size_t)L*H2); ly->Kss2 = malloc(sizeof(double _Complex)*H2); ly->K1ss2 = malloc(sizeof(double _Complex)*H2);
    for(dp = -(nxs - 1); dp <= nxs - 1; dp++) for(dq = -(nys - 1); dq <= nys - 1; dq++){
      size_t pos = (size_t)((dp + ly->Mx2) % ly->Mx2)*ly->My2 + (dq + ly->My2) % ly->My2;
      lay_resp(t, 0, 0.0, 0.0, 0.0, 0.5*hx, 0.5*hx, dp*hx, dq*hx, 0.0, 0, L, K, &ks);
      for(l = 0; l < L; l++) A[(size_t)l*M2 + pos] = K[l] - ks*t->dt;
      A[(size_t)L*M2 + pos] = ks; A[(size_t)(L + 1)*M2 + pos] = K[0];
      if(dp == 0 && dq == 0) ly->self_ss = K[0];
    }
    /* transform each [M2] array: reuse l2 machinery by building element-order arrays is awkward for
       offsets, so transform directly */
    for(l = 0; l < L + 2; l++){
      double _Complex *a = malloc(sizeof(double _Complex)*M2), *col = malloc(sizeof(double _Complex)*(ly->Mx2 > ly->My2 ? ly->Mx2 : ly->My2)),
        *o = malloc(sizeof(double _Complex)*(ly->Mx2 > ly->My2 ? ly->Mx2 : ly->My2));
      double _Complex *dst = (l < L) ? ly->Wss2 + (size_t)l*H2 : (l == L ? ly->Kss2 : ly->K1ss2);
      int i2, j2;
      for(s1 = 0; s1 < M2; s1++) a[s1] = A[(size_t)l*M2 + s1];
      for(i2 = 0; i2 < ly->Mx2; i2++){ fd_fft_exec(ly->pf2y, a + (size_t)i2*ly->My2, o); memcpy(a + (size_t)i2*ly->My2, o, sizeof(double _Complex)*ly->My2); }
      for(j2 = 0; j2 < ly->Mh2; j2++){
        for(i2 = 0; i2 < ly->Mx2; i2++) col[i2] = a[(size_t)i2*ly->My2 + j2];
        fd_fft_exec(ly->pf2x, col, o);
        for(i2 = 0; i2 < ly->Mx2; i2++) dst[(size_t)i2*ly->Mh2 + j2] = o[i2];
      }
      free(a); free(col); free(o);
    }
    free(A); free(tmp);
  }
  /* fault <-> layer: Toeplitz along strike, offsets d on the lattice */
  dmin = -(nxs - 1) < -m*(nx - 1) ? -(nxs - 1) : -m*(nx - 1);
  dmax = (nxs - 1) > m*(nx - 1) ? (nxs - 1) : m*(nx - 1);
  {
    size_t nfs = (size_t)nys*nz;
    double *Afs = calloc((size_t)(L + 2)*ly->Mx*nfs, sizeof(double));   /* [l][x][q][j] */
    double *Asf = calloc((size_t)(L + 2)*ly->Mx*nfs, sizeof(double));   /* [l][x][j][q] */
    buf = malloc(sizeof(double _Complex)*ly->Mx);
    for(d = dmin; d <= dmax; d++){
      size_t pos = (size_t)((d + ly->Mx) % ly->Mx);
      for(q = 0; q < nys; q++){
        double yq = ly->ys0 + (q + 0.5)*hx;
        for(j = 0; j < nz; j++){
          double zj = -(j + 0.5)*dx;
          /* fault (and image) -> layer: receiver at x offset d hx + c0, d = p - m i in [-m(nx-1), nxs-1];
             layer -> fault: d = m i - p in [-(nxs-1), m(nx-1)]; each range fits the FFT length,
             the union may not (filling both from one loop would alias) */
          if(d >= -m*(nx - 1) && d <= nxs - 1){
          lay_resp(t, 1, 0.0, 0.0, zj, 0.5*dx, 0.5*dx, d*hx + ly->c0, yq, 0.0, 0, L, K, &ks);
          lay_resp(t, 1, 0.0, 0.0, -zj, 0.5*dx, 0.5*dx, d*hx + ly->c0, yq, 0.0, 0, L, K2, &ks2);
          for(l = 0; l < L; l++) Afs[(((size_t)l*ly->Mx + pos)*nys + q)*nz + j] = K[l] + K2[l] - (ks + ks2)*t->dt;
          Afs[(((size_t)L*ly->Mx + pos)*nys + q)*nz + j] = ks + ks2;
          Afs[(((size_t)(L + 1)*ly->Mx + pos)*nys + q)*nz + j] = K[0] + K2[0];
          }
          if(d >= -(nxs - 1) && d <= m*(nx - 1)){
          /* layer -> fault: receiver at x offset d hx - c0 */
          lay_resp(t, 0, 0.0, yq, 0.0, 0.5*hx, 0.5*hx, d*hx - ly->c0, 0.0, zj, 1, L, K, &ks);
          for(l = 0; l < L; l++) Asf[(((size_t)l*ly->Mx + pos)*nz + j)*nys + q] = K[l] - ks*t->dt;
          Asf[(((size_t)L*ly->Mx + pos)*nz + j)*nys + q] = ks;
          Asf[(((size_t)(L + 1)*ly->Mx + pos)*nz + j)*nys + q] = K[0];
          }
        }
      }
    }
    ly->Wfs_h = malloc(sizeof(double _Complex)*(size_t)L*Mh*nfs); ly->Kfs_h = malloc(sizeof(double _Complex)*Mh*nfs); ly->K1fs_h = malloc(sizeof(double _Complex)*Mh*nfs);
    ly->Wsf_h = malloc(sizeof(double _Complex)*(size_t)L*Mh*nfs); ly->Ksf_h = malloc(sizeof(double _Complex)*Mh*nfs); ly->K1sf_h = malloc(sizeof(double _Complex)*Mh*nfs);
    for(l = 0; l < L + 2; l++) for(s1 = 0; s1 < nfs; s1++){
      double _Complex *dfs = (l < L) ? ly->Wfs_h + (size_t)l*Mh*nfs : (l == L ? ly->Kfs_h : ly->K1fs_h);
      double _Complex *dsf = (l < L) ? ly->Wsf_h + (size_t)l*Mh*nfs : (l == L ? ly->Ksf_h : ly->K1sf_h);
      double _Complex *a = calloc(ly->Mx, sizeof(double _Complex));
      int x;
      for(x = 0; x < ly->Mx; x++) a[x] = Afs[((size_t)l*ly->Mx + x)*nfs + s1];
      fd_fft_exec(ly->pf, a, buf);
      for(k = 0; k < (int)Mh; k++) dfs[(size_t)k*nfs + s1] = buf[k];
      for(x = 0; x < ly->Mx; x++) a[x] = Asf[((size_t)l*ly->Mx + x)*nfs + s1];
      fd_fft_exec(ly->pf, a, buf);
      for(k = 0; k < (int)Mh; k++) dsf[(size_t)k*nfs + s1] = buf[k];
      free(a);
    }
    free(Afs); free(Asf); free(buf);
  }
  free(K); free(K2);
  ly->V = calloc(ly->ns, sizeof(double)); ly->Vold = calloc(ly->ns, sizeof(double)); ly->delta = calloc(ly->ns, sizeof(double));
  ly->H = calloc(ly->ns, sizeof(double)); ly->tau = calloc(ly->ns, sizeof(double));
  ly->hs2 = calloc((size_t)L*ly->Mx2*ly->Mh2, sizeof(double _Complex));
  ly->hsx = calloc((size_t)L*Mh*nys, sizeof(double _Complex));
  ly->hfx = calloc((size_t)L*Mh*nz, sizeof(double _Complex));
  return ly;
}

/* element-order fault field (i*nz + j) -> rows [j][x] at stride m */
static void fault_rows_hat(const struct fd3d_layer *ly, const double *V, double _Complex *out)
{
  double *r = malloc(sizeof(double)*(size_t)ly->nz*ly->nx); int i, j;
  for(i = 0; i < ly->nx; i++) for(j = 0; j < ly->nz; j++) r[(size_t)j*ly->nx + i] = V[(size_t)i*ly->nz + j];
  lx_fwd(ly, ly->nz, ly->nx, ly->m, r, 1, out);
  free(r);
}
static void layer_rows_hat(const struct fd3d_layer *ly, const double *V, double _Complex *out)
{
  lx_fwd(ly, ly->nys, ly->nxs, 1, V, 0, out);              /* V in element order p*nys + q: column-major for rows */
}

static int fd3d_step_layer_fft(struct fd3d *f, const double *load)
{
  struct fd3d_layer *ly = f->layer;
  const struct fd3d_tab *t = f->t;
  int N = f->N, ns = ly->ns, nz = ly->nz, nys = ly->nys, i, l, it, inner, itot = 0, lmax = ly->L < ly->nhist + 1 ? ly->L : ly->nhist + 1;
  size_t Mh = ly->Mh, H2 = (size_t)ly->Mx2*ly->Mh2, m;
  double dt = t->dt, hw = f->nl ? 0.5 : 1.0;
  double *Hx = calloc(N, sizeof(double)), *Vs = malloc(sizeof(double)*ns), *Vp = malloc(sizeof(double)*N), *rhs = malloc(sizeof(double)*ns);
  double *tmpf = malloc(sizeof(double)*N), *tmps = malloc(sizeof(double)*ns);
  double _Complex *Xs = malloc(sizeof(double _Complex)*Mh*nys), *Xf = malloc(sizeof(double _Complex)*Mh*nz);
  double _Complex *Ys = calloc(Mh*nys, sizeof(double _Complex)), *Yf = calloc(Mh*nz, sizeof(double _Complex));
  double _Complex *S2 = malloc(sizeof(double _Complex)*H2), *P2 = calloc(H2, sizeof(double _Complex));
  double *thbuf = (f->law == 1) ? malloc(sizeof(double)*N) : NULL;
  double diag = hw*ly->self_ss - (f->nl ? 0.5*ly->eta_p : 0.0);
  double tp0 = ly_now(), tp1;
  step_begin(f, load);
  memcpy(ly->Vold, ly->V, sizeof(double)*ns);
  /* known parts in the spectral domain */
  memset(ly->H, 0, sizeof(double)*ns);
  l2_fwd(ly, ly->delta, S2);
  for(m = 0; m < H2; m++) P2[m] = ly->Kss2[m]*S2[m];
  fault_rows_hat(ly, f->delta, Xf); lx_apply(ly, nys, nz, ly->Kfs_h, Xf, Ys);
  if(!f->static_fn){ layer_rows_hat(ly, ly->delta, Xs); lx_apply(ly, nz, nys, ly->Ksf_h, Xs, Yf); }
  for(l = 2; l <= lmax; l++){
    int slot = (ly->hk - (l - 1) + 2*ly->L) % ly->L;
    const double _Complex *W2 = ly->Wss2 + (size_t)(l - 1)*H2, *h2 = ly->hs2 + (size_t)slot*H2;
    for(m = 0; m < H2; m++) P2[m] += W2[m]*h2[m];
    lx_apply(ly, nys, nz, ly->Wfs_h + (size_t)(l - 1)*Mh*nys*nz, ly->hfx + (size_t)slot*Mh*nz, Ys);
    lx_apply(ly, nz, nys, ly->Wsf_h + (size_t)(l - 1)*Mh*nys*nz, ly->hsx + (size_t)slot*Mh*nys, Yf);
  }
  if(f->nl){
    l2_fwd(ly, ly->Vold, S2);
    for(m = 0; m < H2; m++) P2[m] += 0.5*ly->K1ss2[m]*S2[m];
    for(i = 0; i < ns; i++) ly->H[i] += 0.5*ly->eta_p*ly->Vold[i];
    fault_rows_hat(ly, f->Vold, Xf);
    {
      double _Complex *Yt = calloc(Mh*nys, sizeof(double _Complex)); lx_apply(ly, nys, nz, ly->K1fs_h, Xf, Yt);
      for(m = 0; m < Mh*nys; m++) Ys[m] += 0.5*Yt[m];
      free(Yt);
    }
    layer_rows_hat(ly, ly->Vold, Xs);
    {
      double _Complex *Yt = calloc(Mh*nz, sizeof(double _Complex)); lx_apply(ly, nz, nys, ly->K1sf_h, Xs, Yt);
      for(m = 0; m < Mh*nz; m++) Yf[m] += 0.5*Yt[m];
      free(Yt);
    }
  }
  l2_inv(ly, P2, 1.0, ly->H);
  {
    double *Hrow = calloc((size_t)ns, sizeof(double));
    lx_inv(ly, nys, Ys, ly->nxs, 1, 1.0, Hrow);            /* [p][q] = element order */
    for(i = 0; i < ns; i++) ly->H[i] += Hrow[i];
    free(Hrow);
  }
  {
    double *Hrow = calloc((size_t)N, sizeof(double));       /* [i][j] = element order i*nz + j */
    lx_inv(ly, nz, Yf, ly->nx, ly->m, 1.0, Hrow);
    for(i = 0; i < N; i++) Hx[i] += Hrow[i];
    free(Hrow);
  }
  tp1 = ly_now(); fd3d_layer_prof[0] += tp1 - tp0; tp0 = tp1;
  /* iterate fault and layer */
  memcpy(Vs, ly->Vold, sizeof(double)*ns);
  f->thnew = thbuf;
  double *Hf_prev = malloc(sizeof(double)*N); char *mask = malloc(N);
  for(it = 0; it < 200; it++){
    double ch = 0.0, vm = 1e-30, *Hf = malloc(sizeof(double)*N);
    /* fault sees the layer's current step through the first lag */
    memset(Yf, 0, sizeof(double _Complex)*Mh*nz);
    layer_rows_hat(ly, Vs, Xs); lx_apply(ly, nz, nys, ly->K1sf_h, Xs, Yf);
    memcpy(Hf, Hx, sizeof(double)*N);
    memset(tmpf, 0, sizeof(double)*N); lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, tmpf);
    for(i = 0; i < N; i++) Hf[i] += tmpf[i];
    memcpy(Vp, f->V, sizeof(double)*N);
    tp1 = ly_now(); fd3d_layer_prof[2] += tp1 - tp0; tp0 = tp1;
    itot += step_solve_g(f, Hf, it > 0 ? Vp : NULL, NULL);   /* warm start from the previous outer iterate */
    memcpy(Hf_prev, Hf, sizeof(double)*N);
    tp1 = ly_now(); fd3d_layer_prof[1] += tp1 - tp0; tp0 = tp1;
    free(Hf);
    /* layer: (hw K1ss - nl eta/2) Vs = -(H + hw K1fs Vf), Jacobi with the self coefficient */
    memset(Ys, 0, sizeof(double _Complex)*Mh*nys);
    fault_rows_hat(ly, f->V, Xf); lx_apply(ly, nys, nz, ly->K1fs_h, Xf, Ys);
    for(i = 0; i < ns; i++) rhs[i] = -ly->H[i];
    lx_inv(ly, nys, Ys, ly->nxs, 1, -hw, rhs);
    for(inner = 0; inner < 500; inner++){
      double chi = 0.0, vmi = 1e-30;
      l2_fwd(ly, Vs, S2);
      for(m = 0; m < H2; m++) P2[m] = ly->K1ss2[m]*S2[m];
      memset(tmps, 0, sizeof(double)*ns); l2_inv(ly, P2, hw, tmps);
      if(f->nl) for(i = 0; i < ns; i++) tmps[i] -= 0.5*ly->eta_p*Vs[i];
      for(i = 0; i < ns; i++){
        double dv = (rhs[i] - tmps[i])/diag;
        Vs[i] += dv; if(fabs(dv) > chi) chi = fabs(dv); if(fabs(Vs[i]) > vmi) vmi = fabs(Vs[i]);
      }
      if(chi <= 1e-13*vmi) break;
    }
    ly->n_inner += inner + 1;
    for(i = 0; i < ns; i++){ if(fabs(Vs[i] - ly->V[i]) > ch) ch = fabs(Vs[i] - ly->V[i]); if(fabs(Vs[i]) > vm) vm = fabs(Vs[i]); }
    {
      /* would the fault see a different known part with the new layer rates? */
      double dmax = 0.0, smax = 1.0;
      memset(Yf, 0, sizeof(double _Complex)*Mh*nz);
      for(i = 0; i < ns; i++) tmps[i] = Vs[i] - ly->V[i];
      layer_rows_hat(ly, tmps, Xs); lx_apply(ly, nz, nys, ly->K1sf_h, Xs, Yf);
      memset(tmpf, 0, sizeof(double)*N); lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, tmpf);
      for(i = 0; i < N; i++){ if(fabs(tmpf[i]) > dmax) dmax = fabs(tmpf[i]); if(fabs(f->sigma[i]) > smax) smax = fabs(f->sigma[i]); }
      memcpy(ly->V, Vs, sizeof(double)*ns);
      if(dmax <= 1e-12*smax){ it++; break; }      /* converged: the fault solve already used these rates */
    }
    for(i = 0; i < N; i++){ if(fabs(f->V[i] - Vp[i]) > ch) ch = fabs(f->V[i] - Vp[i]); if(fabs(f->V[i]) > vm) vm = fabs(f->V[i]); }
    if(it > 0 && ch <= 1e-12*vm) break;
  }
  ly->n_outer += it + 1; ly->n_steps++;
  free(Hf_prev); free(mask);
  tp1 = ly_now(); fd3d_layer_prof[2] += tp1 - tp0; tp0 = tp1;
  f->thnew = NULL;
  if(thbuf){ memcpy(f->theta, thbuf, sizeof(double)*N); free(thbuf); }
  /* finish: fault known part with the layer's current step, histories, slips */
  memset(Yf, 0, sizeof(double _Complex)*Mh*nz);
  layer_rows_hat(ly, ly->V, Xs); lx_apply(ly, nz, nys, ly->K1sf_h, Xs, Yf);
  memset(tmpf, 0, sizeof(double)*N); lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, tmpf);
  for(i = 0; i < N; i++) f->H[i] += Hx[i] + tmpf[i];
  {
    int slot = ly->hk % ly->L;
    double *vms = malloc(sizeof(double)*ns), *vmf = malloc(sizeof(double)*N);
    for(i = 0; i < ns; i++){ vms[i] = f->nl ? 0.5*(ly->Vold[i] + ly->V[i]) : ly->V[i]; ly->delta[i] += vms[i]*dt; }
    for(i = 0; i < N; i++) vmf[i] = f->nl ? 0.5*(f->Vold[i] + f->V[i]) : f->V[i];
    l2_fwd(ly, vms, ly->hs2 + (size_t)slot*H2);
    layer_rows_hat(ly, vms, ly->hsx + (size_t)slot*Mh*nys);
    fault_rows_hat(ly, vmf, ly->hfx + (size_t)slot*Mh*nz);
    ly->hk = (ly->hk + 1) % ly->L;
    if(ly->nhist < ly->L) ly->nhist++;
    free(vms); free(vmf);
  }
  step_end(f);
  tp1 = ly_now(); fd3d_layer_prof[3] += tp1 - tp0;
  free(Hx); free(Vs); free(Vp); free(rhs); free(tmpf); free(tmps); free(Xs); free(Xf); free(Ys); free(Yf); free(S2); free(P2);
  return itot;
}

/* ---------------- full virtual surface for dipping faults (stage B) ---------------- */

static void frame_fault(double dip_deg, double F[3][3])
{
  double d = dip_deg*M_PI/180.0;
  F[0][0] = 1.0; F[0][1] = 0.0; F[0][2] = 0.0;
  F[1][0] = 0.0; F[1][1] = cos(d); F[1][2] = sin(d);                 /* up dip */
  F[2][0] = 0.0; F[2][1] = -sin(d); F[2][2] = cos(d);                /* e1 x e2 */
}
/* global stress tensor from local s6 (11 22 33 12 31 32) in frame R (rows = local axes) */
static void to_global(double R[3][3], const double s[6], double G[3][3])
{
  double S[3][3] = {{s[0], s[3], s[4]}, {s[3], s[1], s[5]}, {s[4], s[5], s[2]}};
  int i, j, a, b;
  for(i = 0; i < 3; i++) for(j = 0; j < 3; j++){ double v = 0.0; for(a = 0; a < 3; a++) for(b = 0; b < 3; b++) v += R[a][i]*S[a][b]*R[b][j]; G[i][j] = v; }
}
/* tractions at a receiver for unit rate on a source rectangle, all lags plus the static slope:
   src_fault: fault frame (comp fcomp) else surface frame (comp 1, 2, 3);
   rcv_fault: (shear along s, normal on n) else (sxz, syz, szz); K[r][l], ks[r] */
static void surf_resp(const struct fd3d_tab *t, double F[3][3], int src_fault, int comp, const double c[3], double a1, double a2,
                      const double X[3], int rcv_fault, const double sdir[3], int L, double *K, double *ks)
{
  double Ide[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, (*R)[3] = src_fault ? F : Ide, d[3], x[3], s[6], G[3][3], tr[3], prev[3] = {0, 0, 0}, late[2][3], T1, T2;
  int i, l, q, nr = rcv_fault ? 2 : 3;
  for(i = 0; i < 3; i++) d[i] = X[i] - c[i];
  for(i = 0; i < 3; i++) x[i] = R[i][0]*d[0] + R[i][1]*d[1] + R[i][2]*d[2];
  T1 = 3.0*(sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) + 4.0*(a1 + a2))/t->cT + 10.0*t->dt; T2 = 2.0*T1;
#define SURF_TR(G, tr) do{ if(rcv_fault){ double *n = F[2]; \
      tr[0] = 0.0; tr[1] = 0.0; \
      { int u, v; for(u = 0; u < 3; u++) for(v = 0; v < 3; v++){ tr[0] += sdir[u]*G[u][v]*n[v]; tr[1] += n[u]*G[u][v]*n[v]; } } \
    }else{ tr[0] = G[0][2]; tr[1] = G[1][2]; tr[2] = G[2][2]; } }while(0)
  for(q = 0; q < 2; q++){
    fd3d_rect_ramp_tensor(x[0], x[1], x[2], q ? T2 : T1, a1, a2, t->cL, t->cT, t->rho, comp, s);
    to_global(R, s, G); SURF_TR(G, late[q]);
  }
  for(i = 0; i < nr; i++) ks[i] = (late[1][i] - late[0][i])/(T2 - T1);
  for(l = 1; l <= L; l++){
    fd3d_rect_ramp_tensor(x[0], x[1], x[2], l*t->dt, a1, a2, t->cL, t->cT, t->rho, comp, s);
    to_global(R, s, G); SURF_TR(G, tr);
    for(i = 0; i < nr; i++){ K[(size_t)i*L + l - 1] = tr[i] - prev[i]; prev[i] = tr[i]; }
  }
#undef SURF_TR
}

static void fft_x_inplace(const struct fd3d_layer *ly, const double *in, double _Complex *out, double _Complex *buf)
{
  double _Complex *a = malloc(sizeof(double _Complex)*ly->Mx); int x;
  for(x = 0; x < ly->Mx; x++) a[x] = in[x];
  fd_fft_exec(ly->pf, a, buf);
  for(x = 0; x < ly->Mh; x++) out[x] = buf[x];
  free(a);
}

struct fd3d_layer *fd3d_surf_create_rows(struct fd3d *f, double x0f, double ytop, double ztop, double dip_deg, int fcomp,
                                         int m, int nxs, int nys, double xs0, const double *yc, const double *wy, int L, double fsign);
struct fd3d_layer *fd3d_surf_create(struct fd3d *f, double x0f, double ytop, double ztop, double dip_deg, int fcomp,
                                    int m, int nxs, int nys, double xs0, double ys0, int L, double fsign)
{
  double hx = f->t->dx/m, *yc = malloc(sizeof(double)*nys), *wy = malloc(sizeof(double)*nys);
  struct fd3d_layer *ly; int q;
  for(q = 0; q < nys; q++){ yc[q] = ys0 + (q + 0.5)*hx; wy[q] = hx; }
  ly = fd3d_surf_create_rows(f, x0f, ytop, ztop, dip_deg, fcomp, m, nxs, nys, xs0, yc, wy, L, fsign);
  free(yc); free(wy);
  return ly;
}
/* rows of given centres yc[q] and widths wy[q] (any grading; receivers must not fall within about
   1e-2 of a width from another row's edge line unless exactly on it, see fd_plan.md B5): the
   surface-surface tables are then Toeplitz along strike with dense (row) x (row) blocks */
struct fd3d_layer *fd3d_surf_create_rows(struct fd3d *f, double x0f, double ytop, double ztop, double dip_deg, int fcomp,
                                         int m, int nxs, int nys, double xs0, const double *yc, const double *wy, int L, double fsign)
{
  const struct fd3d_tab *t = f->t;
  struct fd3d_layer *ly;
  int nx = t->nx, nz = f->nzp, d, dmin, dmax, q, j, l, c, r, p, k, N = f->N;
  double dx = t->dx, hx = dx/m, F[3][3], sdir[3], *K = malloc(sizeof(double)*3*L), ks[3];
  size_t Mh, nfs;
  if(f->mirror || m < 1 || (fcomp != 1 && fcomp != 2)) return NULL;
  frame_fault(dip_deg, F);
  for(c = 0; c < 3; c++) sdir[c] = fsign*F[fcomp - 1][c];          /* slip direction of the marcher's slip variable */
  ly = calloc(1, sizeof(*ly));
  ly->fft = 1; ly->surf = 1; ly->fcomp = fcomp; ly->dip = dip_deg; ly->ytop = ytop; ly->ztop = ztop;
  ly->m = m; ly->nxs = nxs; ly->nys = nys; ly->ns = nxs*nys; ly->L = L; ly->hx = hx; ly->hy = wy[0];
  ly->xs0 = xs0; ly->ys0 = yc[0] - 0.5*wy[0]; ly->x0f = x0f; ly->nx = nx; ly->nz = nz; ly->eta_p = 0.5*t->rho*t->cL;
  ly->yc = malloc(sizeof(double)*nys); ly->wy = malloc(sizeof(double)*nys);
  memcpy(ly->yc, yc, sizeof(double)*nys); memcpy(ly->wy, wy, sizeof(double)*nys);
  ly->self_c = calloc((size_t)3*nys, sizeof(double));
  ly->c0 = xs0 + 0.5*hx - x0f - 0.5*dx;
  ly->Mx = fd_fft_size_for((nxs + m*(nx - 1) + 1) > (2*nxs - 1) ? (nxs + m*(nx - 1) + 1) : (2*nxs - 1));   /* fault <-> surface and surface <-> surface offsets */
  ly->Mh = ly->Mx/2 + 1; Mh = ly->Mh;
  ly->Mx2 = fd_fft_size_for(2*nxs - 1); ly->My2 = fd_fft_size_for(2*nys - 1); ly->Mh2 = ly->My2/2 + 1;
  nfs = (size_t)nys*nz;
  ly->pf = fd_fft_plan_create(ly->Mx, FD_FFT_FORWARD); ly->pb = fd_fft_plan_create(ly->Mx, FD_FFT_BACKWARD);
  ly->pf2x = fd_fft_plan_create(ly->Mx2, FD_FFT_FORWARD); ly->pb2x = fd_fft_plan_create(ly->Mx2, FD_FFT_BACKWARD);
  ly->pf2y = fd_fft_plan_create(ly->My2, FD_FFT_FORWARD); ly->pb2y = fd_fft_plan_create(ly->My2, FD_FFT_BACKWARD);
  /* surface-surface: pairs (rcv, src) (0,0) (0,1) (1,0) (1,1) (2,2), Toeplitz along strike (offset
     d = p' - p), dense in (receiver row, source row); row widths may differ */
  {
    static const int pr[5] = {0, 0, 1, 1, 2}, ps[5] = {0, 1, 0, 1, 2};
    size_t nss = (size_t)nys*nys;
    double *A = calloc((size_t)5*(L + 2)*ly->Mx*nss, sizeof(double)), cen[3];
    double _Complex *buf = malloc(sizeof(double _Complex)*ly->Mx), *col = malloc(sizeof(double _Complex)*Mh);
    double *line = malloc(sizeof(double)*ly->Mx);
    int qs, qr, x, uniform = 1;
    for(qs = 1; qs < nys; qs++) if(fabs(wy[qs] - wy[0]) > 1e-9*wy[0] || fabs((yc[qs] - yc[qs - 1]) - wy[0]) > 1e-9*wy[0]) uniform = 0;
    for(qs = 0; qs < nys; qs++){
      cen[0] = 0.0; cen[1] = yc[qs]; cen[2] = 0.0;
      for(qr = 0; qr < nys; qr++){
        /* uniform rows: the response depends on qr - qs only; computed once per offset (qs = 0 or
           qr = 0) and copied */
        int copy_from_qs = -1, copy_from_qr = -1;
        if(uniform && qs > 0 && qr > 0){ copy_from_qs = qs - 1; copy_from_qr = qr - 1; }
        for(d = -(nxs - 1); d <= nxs - 1; d++){
          size_t pos = (size_t)((d + ly->Mx) % ly->Mx);
          double X[3] = {d*hx, yc[qr], 0.0};
          if(copy_from_qs >= 0){
            for(p = 0; p < 5; p++){
              size_t base = ((size_t)p*(L + 2));
              for(l = 0; l < L + 2; l++) A[((base + l)*ly->Mx + pos)*nss + (size_t)qr*nys + qs] = A[((base + l)*ly->Mx + pos)*nss + (size_t)copy_from_qr*nys + copy_from_qs];
            }
            if(d == 0 && qr == qs) for(r = 0; r < 3; r++) ly->self_c[(size_t)r*nys + qs] = ly->self_c[(size_t)r*nys];
            continue;
          }
          for(c = 0; c < 3; c++){
            surf_resp(t, F, 0, c + 1, cen, 0.5*hx, 0.5*wy[qs], X, 0, sdir, L, K, ks);
            for(p = 0; p < 5; p++) if(ps[p] == c){
              size_t base = ((size_t)p*(L + 2));
              r = pr[p];
              for(l = 0; l < L; l++) A[((base + l)*ly->Mx + pos)*nss + (size_t)qr*nys + qs] = K[(size_t)r*L + l] - ks[r]*t->dt;
              A[((base + L)*ly->Mx + pos)*nss + (size_t)qr*nys + qs] = ks[r];
              A[((base + L + 1)*ly->Mx + pos)*nss + (size_t)qr*nys + qs] = K[(size_t)r*L];
              if(d == 0 && qr == qs && (p == 0 || p == 3 || p == 4)) ly->self_c[(size_t)r*nys + qs] = K[(size_t)r*L];
            }
          }
        }
      }
    }
    ly->self_sh = ly->self_c[0]; ly->self_op = ly->self_c[(size_t)2*nys];
    for(p = 0; p < 5; p++){
      ly->T2[p][0] = malloc(sizeof(double _Complex)*(size_t)L*Mh*nss); ly->T2[p][1] = malloc(sizeof(double _Complex)*Mh*nss); ly->T2[p][2] = malloc(sizeof(double _Complex)*Mh*nss);
      for(l = 0; l < L + 2; l++) for(size_t s1 = 0; s1 < nss; s1++){
        double _Complex *dst = l < L ? ly->T2[p][0] + (size_t)l*Mh*nss : ly->T2[p][l - L + 1];
        for(x = 0; x < ly->Mx; x++) line[x] = A[(((size_t)p*(L + 2) + l)*ly->Mx + x)*nss + s1];
        fft_x_inplace(ly, line, col, buf);
        for(k = 0; k < (int)Mh; k++) dst[(size_t)k*nss + s1] = col[k];
      }
    }
    free(A); free(buf); free(col); free(line);
  }
  /* fault <-> surface */
  dmin = -(nxs - 1) < -m*(nx - 1) ? -(nxs - 1) : -m*(nx - 1);
  dmax = (nxs - 1) > m*(nx - 1) ? (nxs - 1) : m*(nx - 1);
  {
    double *Afs = calloc((size_t)3*(L + 2)*ly->Mx*nfs, sizeof(double)), *Asf = calloc((size_t)6*(L + 2)*ly->Mx*nfs, sizeof(double));
    double _Complex *buf = malloc(sizeof(double _Complex)*ly->Mx), *col = malloc(sizeof(double _Complex)*Mh);
    double *line = malloc(sizeof(double)*ly->Mx);
    for(d = dmin; d <= dmax; d++){
      size_t pos = (size_t)((d + ly->Mx) % ly->Mx);
      for(q = 0; q < nys; q++){
        double yq = yc[q];
        for(j = 0; j < nz; j++){
          double yj = ytop - (j + 0.5)*dx*F[1][1], zj = ztop - (j + 0.5)*dx*F[1][2];
          double cf[3] = {0.0, yj, zj}, Xs[3] = {d*hx + ly->c0, yq, 0.0}, cs[3] = {0.0, yq, 0.0}, Xf[3] = {d*hx - ly->c0, yj, zj};
          if(d >= -m*(nx - 1) && d <= nxs - 1){        /* fault -> surface: d = p - m i (ranges as in the layer) */
          surf_resp(t, F, 1, fcomp, cf, 0.5*dx, 0.5*dx, Xs, 0, sdir, L, K, ks);
          if(j == 0 && f->next){                       /* extension rows above the surface, tied to row 0 */
            int e; double *K2 = malloc(sizeof(double)*3*L), ks2[3];
            for(e = 0; e < f->next; e++){
              double ce[3] = {0.0, ytop + (e + 0.5)*dx*F[1][1], ztop + (e + 0.5)*dx*F[1][2]};
              surf_resp(t, F, 1, fcomp, ce, 0.5*dx, 0.5*dx, Xs, 0, sdir, L, K2, ks2);
              for(r = 0; r < 3; r++){ ks[r] += ks2[r]; for(l = 0; l < L; l++) K[(size_t)r*L + l] += K2[(size_t)r*L + l]; }
            }
            free(K2);
          }
          for(r = 0; r < 3; r++){ ks[r] *= fsign; for(l = 0; l < L; l++) K[(size_t)r*L + l] *= fsign; }
          for(r = 0; r < 3; r++){
            for(l = 0; l < L; l++) Afs[(((size_t)r*(L + 2) + l)*ly->Mx + pos)*nfs + (size_t)q*nz + j] = K[(size_t)r*L + l] - ks[r]*t->dt;
            Afs[(((size_t)r*(L + 2) + L)*ly->Mx + pos)*nfs + (size_t)q*nz + j] = ks[r];
            Afs[(((size_t)r*(L + 2) + L + 1)*ly->Mx + pos)*nfs + (size_t)q*nz + j] = K[(size_t)r*L];
          }
          }
          if(d >= -(nxs - 1) && d <= m*(nx - 1)) for(c = 0; c < 3; c++){   /* surface -> fault: d = m i - p */
            surf_resp(t, F, 0, c + 1, cs, 0.5*hx, 0.5*wy[q], Xf, 1, sdir, L, K, ks);
            for(r = 0; r < 2; r++){
              size_t base = (size_t)(r*3 + c)*(L + 2);
              for(l = 0; l < L; l++) Asf[((base + l)*ly->Mx + pos)*nfs + (size_t)j*nys + q] = K[(size_t)r*L + l] - ks[r]*t->dt;
              Asf[((base + L)*ly->Mx + pos)*nfs + (size_t)j*nys + q] = ks[r];
              Asf[((base + L + 1)*ly->Mx + pos)*nfs + (size_t)j*nys + q] = K[(size_t)r*L];
            }
          }
        }
      }
    }
    for(r = 0; r < 3; r++){
      ly->Tfs[r][0] = malloc(sizeof(double _Complex)*(size_t)L*Mh*nfs); ly->Tfs[r][1] = malloc(sizeof(double _Complex)*Mh*nfs); ly->Tfs[r][2] = malloc(sizeof(double _Complex)*Mh*nfs);
    }
    for(r = 0; r < 2; r++) for(c = 0; c < 3; c++){
      ly->Tsf[r][c][0] = malloc(sizeof(double _Complex)*(size_t)L*Mh*nfs); ly->Tsf[r][c][1] = malloc(sizeof(double _Complex)*Mh*nfs); ly->Tsf[r][c][2] = malloc(sizeof(double _Complex)*Mh*nfs);
    }
    for(l = 0; l < L + 2; l++) for(size_t s1 = 0; s1 < nfs; s1++){
      int x;
      for(r = 0; r < 3; r++){
        double _Complex *dst = l < L ? ly->Tfs[r][0] + (size_t)l*Mh*nfs : ly->Tfs[r][l - L + 1];
        for(x = 0; x < ly->Mx; x++) line[x] = Afs[(((size_t)r*(L + 2) + l)*ly->Mx + x)*nfs + s1];
        fft_x_inplace(ly, line, col, buf);
        for(k = 0; k < (int)Mh; k++) dst[(size_t)k*nfs + s1] = col[k];
      }
      for(r = 0; r < 2; r++) for(c = 0; c < 3; c++){
        size_t base = (size_t)(r*3 + c)*(L + 2);
        double _Complex *dst = l < L ? ly->Tsf[r][c][0] + (size_t)l*Mh*nfs : ly->Tsf[r][c][l - L + 1];
        for(x = 0; x < ly->Mx; x++) line[x] = Asf[((base + l)*ly->Mx + x)*nfs + s1];
        fft_x_inplace(ly, line, col, buf);
        for(k = 0; k < (int)Mh; k++) dst[(size_t)k*nfs + s1] = col[k];
      }
    }
    free(Afs); free(Asf); free(buf); free(col); free(line);
  }
  free(K);
  ly->V3 = calloc((size_t)3*ly->ns, sizeof(double)); ly->Vold3 = calloc((size_t)3*ly->ns, sizeof(double));
  ly->delta3 = calloc((size_t)3*ly->ns, sizeof(double)); ly->H3 = calloc((size_t)3*ly->ns, sizeof(double));
  ly->hsx3 = calloc((size_t)3*L*Mh*nys, sizeof(double _Complex));
  ly->hfx = calloc((size_t)L*Mh*nz, sizeof(double _Complex));
  ly->sig0 = calloc(N, sizeof(double)); ly->tn = calloc(N, sizeof(double));
  ly->normal_coupling = 1;
  return ly;
}

static void surf_free(struct fd3d_layer *ly)
{
  int p, r, c, s;
  for(p = 0; p < 5; p++) for(s = 0; s < 3; s++) free(ly->T2[p][s]);
  for(r = 0; r < 3; r++) for(s = 0; s < 3; s++) free(ly->Tfs[r][s]);
  for(r = 0; r < 2; r++) for(c = 0; c < 3; c++) for(s = 0; s < 3; s++) free(ly->Tsf[r][c][s]);
  free(ly->V3); free(ly->Vold3); free(ly->delta3); free(ly->H3); free(ly->hsx3);
  free(ly->yc); free(ly->wy); free(ly->self_c); free(ly->hfx); free(ly->sig0); free(ly->tn);
  fd_fft_plan_destroy(ly->pf); fd_fft_plan_destroy(ly->pb);
  fd_fft_plan_destroy(ly->pf2x); fd_fft_plan_destroy(ly->pb2x); fd_fft_plan_destroy(ly->pf2y); fd_fft_plan_destroy(ly->pb2y);
  free(ly);
}

/* Y_r[k][q'] += w sum_c T[r][c][k][q'][q] X_c[k][q] for the surface-surface tables (row spectra) */
static const int SP_R[5] = {0, 0, 1, 1, 2}, SP_S[5] = {0, 1, 0, 1, 2};
static void surf_apply_rows(const struct fd3d_layer *ly, int slot_kind, int lag, const double _Complex *X3, double w, double _Complex *Y3)
{
  size_t nss = (size_t)ly->nys*ly->nys, Mh = ly->Mh, Q = Mh*ly->nys, k; int p, r, s;
  for(p = 0; p < 5; p++){
    const double _Complex *T = slot_kind == 0 ? ly->T2[p][0] + (size_t)lag*Mh*nss : ly->T2[p][slot_kind];
    const double _Complex *X = X3 + (size_t)SP_S[p]*Q; double _Complex *Y = Y3 + (size_t)SP_R[p]*Q;
    for(k = 0; k < Mh; k++){
      const double _Complex *Tk = T + k*nss, *Xk = X + k*ly->nys; double _Complex *Yk = Y + k*ly->nys;
      for(r = 0; r < ly->nys; r++){ double _Complex a = 0.0; for(s = 0; s < ly->nys; s++) a += Tk[(size_t)r*ly->nys + s]*Xk[s]; Yk[r] += w*a; }
    }
  }
}

static int fd3d_step_surf(struct fd3d *f, const double *load)
{
  struct fd3d_layer *ly = f->layer;
  const struct fd3d_tab *t = f->t;
  int N = f->N, ns = ly->ns, nz = ly->nz, nys = ly->nys, i, l, it, inner, c, r, itot = 0, lmax = ly->L < ly->nhist + 1 ? ly->L : ly->nhist + 1;
  size_t Mh = ly->Mh, nfs = (size_t)nys*nz, m, Q = Mh*nys;
  double dt = t->dt, hw = f->nl ? 0.5 : 1.0, Es[3];
  double *Hx = calloc(N, sizeof(double)), *Hn = calloc(N, sizeof(double)), *V3 = malloc(sizeof(double)*3*ns), *Vp = malloc(sizeof(double)*N);
  double *rhs = malloc(sizeof(double)*3*ns), *tmpf = malloc(sizeof(double)*N), *tmps = malloc(sizeof(double)*3*ns);
  double _Complex *Xs = malloc(sizeof(double _Complex)*Mh*nys), *Xf = malloc(sizeof(double _Complex)*Mh*nz);
  double _Complex *Ys = malloc(sizeof(double _Complex)*3*Mh*nys), *Yf = malloc(sizeof(double _Complex)*2*Mh*nz);
  double _Complex *S3 = malloc(sizeof(double _Complex)*3*Q), *P3 = malloc(sizeof(double _Complex)*3*Q);
  double *thbuf = (f->law == 1) ? malloc(sizeof(double)*N) : NULL, *dend = malloc(sizeof(double)*N), *tnst = malloc(sizeof(double)*N);
  Es[0] = Es[1] = 0.5*t->rho*t->cT; Es[2] = 0.5*t->rho*t->cL;
  step_begin(f, load);
  memcpy(ly->Vold3, ly->V3, sizeof(double)*3*ns);
  /* ---- known parts ---- */
  memset(Ys, 0, sizeof(double _Complex)*3*Mh*nys); memset(Yf, 0, sizeof(double _Complex)*2*Mh*nz);
  for(c = 0; c < 3; c++) layer_rows_hat(ly, ly->delta3 + (size_t)c*ns, S3 + (size_t)c*Q);
  surf_apply_rows(ly, 1, 0, S3, 1.0, Ys);                                         /* statics */
  fault_rows_hat(ly, f->delta, Xf);
  for(r = 0; r < 3; r++) lx_apply(ly, nys, nz, ly->Tfs[r][1], Xf, Ys + (size_t)r*Mh*nys);
  for(c = 0; c < 3; c++){
    layer_rows_hat(ly, ly->delta3 + (size_t)c*ns, Xs);
    if(!f->static_fn) lx_apply(ly, nz, nys, ly->Tsf[0][c][1], Xs, Yf);
    if(!ly->static_n_fn) lx_apply(ly, nz, nys, ly->Tsf[1][c][1], Xs, Yf + Mh*nz);
  }
  for(l = 2; l <= lmax; l++){
    int slot = (ly->hk - (l - 1) + 2*ly->L) % ly->L;
    surf_apply_rows(ly, 0, l - 1, ly->hsx3 + (size_t)slot*3*Q, 1.0, Ys);
    for(r = 0; r < 3; r++) lx_apply(ly, nys, nz, ly->Tfs[r][0] + (size_t)(l - 1)*Mh*nfs, ly->hfx + (size_t)slot*Mh*nz, Ys + (size_t)r*Mh*nys);
    for(c = 0; c < 3; c++) for(r = 0; r < 2; r++)
      lx_apply(ly, nz, nys, ly->Tsf[r][c][0] + (size_t)(l - 1)*Mh*nfs, ly->hsx3 + ((size_t)slot*3 + c)*Mh*nys, Yf + (size_t)r*Mh*nz);
  }
  if(f->nl){
    double _Complex *Yt = malloc(sizeof(double _Complex)*3*Mh*nys);
    for(c = 0; c < 3; c++) layer_rows_hat(ly, ly->Vold3 + (size_t)c*ns, S3 + (size_t)c*Q);
    surf_apply_rows(ly, 2, 0, S3, 0.5, Ys);
    fault_rows_hat(ly, f->Vold, Xf);
    memset(Yt, 0, sizeof(double _Complex)*3*Mh*nys);
    for(r = 0; r < 3; r++) lx_apply(ly, nys, nz, ly->Tfs[r][2], Xf, Yt + (size_t)r*Mh*nys);
    for(m = 0; m < 3*Mh*nys; m++) Ys[m] += 0.5*Yt[m];
    free(Yt);
    Yt = malloc(sizeof(double _Complex)*2*Mh*nz); memset(Yt, 0, sizeof(double _Complex)*2*Mh*nz);
    for(c = 0; c < 3; c++){ layer_rows_hat(ly, ly->Vold3 + (size_t)c*ns, Xs); for(r = 0; r < 2; r++) lx_apply(ly, nz, nys, ly->Tsf[r][c][2], Xs, Yt + (size_t)r*Mh*nz); }
    for(m = 0; m < 2*Mh*nz; m++) Yf[m] += 0.5*Yt[m];
    free(Yt);
  }
  memset(ly->H3, 0, sizeof(double)*3*ns);
  for(c = 0; c < 3; c++){
    lx_inv(ly, nys, Ys + (size_t)c*Mh*nys, ly->nxs, 1, 1.0, ly->H3 + (size_t)c*ns);
    if(f->nl) for(i = 0; i < ns; i++) ly->H3[(size_t)c*ns + i] += 0.5*Es[c]*ly->Vold3[(size_t)c*ns + i];
  }
  lx_inv(ly, nz, Yf, ly->nx, ly->m, 1.0, Hx);
  lx_inv(ly, nz, Yf + Mh*nz, ly->nx, ly->m, 1.0, Hn);
  /* ---- iterate fault and surface ---- */
  memcpy(V3, ly->Vold3, sizeof(double)*3*ns);
  f->thnew = thbuf;
  for(it = 0; it < 200; it++){
    double ch = 0.0, vm = 1e-30, *Hf = malloc(sizeof(double)*N);
    memset(Yf, 0, sizeof(double _Complex)*2*Mh*nz);
    for(c = 0; c < 3; c++){ layer_rows_hat(ly, V3 + (size_t)c*ns, Xs); for(r = 0; r < 2; r++) lx_apply(ly, nz, nys, ly->Tsf[r][c][2], Xs, Yf + (size_t)r*Mh*nz); }
    memcpy(Hf, Hx, sizeof(double)*N); lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, Hf);
    memcpy(ly->tn, Hn, sizeof(double)*N); lx_inv(ly, nz, Yf + Mh*nz, ly->nx, ly->m, hw, ly->tn);
    if(ly->static_n_fn){
      /* exact static normal traction of the slip at the end of the step (current fault rates) */
      for(i = 0; i < N; i++) dend[i] = f->delta[i] + (f->nl ? 0.5*dt*(f->Vold[i] + f->V[i]) : dt*f->V[i]);
      ly->static_n_fn(ly->static_n_ctx, dend, tnst);
      for(i = 0; i < N; i++) ly->tn[i] += tnst[i];
    }
    if(ly->normal_coupling) for(i = 0; i < N; i++){ f->sigma[i] = ly->sig0[i] - ly->tn[i]; if(f->sigma[i] < 0.0) f->sigma[i] = 0.0; }
    memcpy(Vp, f->V, sizeof(double)*N);
    itot += step_solve_g(f, Hf, it > 0 ? Vp : NULL, NULL);
    free(Hf);
    /* surface: zero traction, Jacobi with the self coefficients */
    memset(Ys, 0, sizeof(double _Complex)*3*Mh*nys);
    fault_rows_hat(ly, f->V, Xf);
    for(r = 0; r < 3; r++) lx_apply(ly, nys, nz, ly->Tfs[r][2], Xf, Ys + (size_t)r*Mh*nys);
    for(i = 0; i < 3*ns; i++) rhs[i] = -ly->H3[i];
    for(c = 0; c < 3; c++) lx_inv(ly, nys, Ys + (size_t)c*Mh*nys, ly->nxs, 1, -hw, rhs + (size_t)c*ns);
    for(inner = 0; inner < 500; inner++){
      double chi = 0.0, vmi = 1e-30;
      for(c = 0; c < 3; c++) layer_rows_hat(ly, V3 + (size_t)c*ns, S3 + (size_t)c*Q);
      memset(P3, 0, sizeof(double _Complex)*3*Q);
      surf_apply_rows(ly, 2, 0, S3, 1.0, P3);
      memset(tmps, 0, sizeof(double)*3*ns);
      for(c = 0; c < 3; c++){
        lx_inv(ly, nys, P3 + (size_t)c*Q, ly->nxs, 1, hw, tmps + (size_t)c*ns);
        if(f->nl) for(i = 0; i < ns; i++) tmps[(size_t)c*ns + i] -= 0.5*Es[c]*V3[(size_t)c*ns + i];
      }
      for(c = 0; c < 3; c++) for(i = 0; i < ns; i++){
        size_t ii = (size_t)c*ns + i; double diag = hw*ly->self_c[(size_t)c*nys + i % nys] - (f->nl ? 0.5*Es[c] : 0.0);
        double dv = (rhs[ii] - tmps[ii])/diag;
        V3[ii] += dv; if(fabs(dv) > chi) chi = fabs(dv); if(fabs(V3[ii]) > vmi) vmi = fabs(V3[ii]);
      }
      if(chi <= 1e-13*vmi) break;
    }
    ly->n_inner += inner + 1;
    for(i = 0; i < 3*ns; i++){ if(fabs(V3[i] - ly->V3[i]) > ch) ch = fabs(V3[i] - ly->V3[i]); if(fabs(V3[i]) > vm) vm = fabs(V3[i]); }
    memcpy(ly->V3, V3, sizeof(double)*3*ns);
    for(i = 0; i < N; i++){ if(fabs(f->V[i] - Vp[i]) > ch) ch = fabs(f->V[i] - Vp[i]); if(fabs(f->V[i]) > vm) vm = fabs(f->V[i]); }
    if(it > 0 && ch <= 1e-12*vm) break;
  }
  ly->n_outer += it + 1; ly->n_steps++;
  f->thnew = NULL;
  if(thbuf){ memcpy(f->theta, thbuf, sizeof(double)*N); free(thbuf); }
  /* ---- finish ---- */
  memset(Yf, 0, sizeof(double _Complex)*2*Mh*nz);
  for(c = 0; c < 3; c++){ layer_rows_hat(ly, ly->V3 + (size_t)c*ns, Xs); for(r = 0; r < 2; r++) lx_apply(ly, nz, nys, ly->Tsf[r][c][2], Xs, Yf + (size_t)r*Mh*nz); }
  memset(tmpf, 0, sizeof(double)*N); lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, tmpf);
  for(i = 0; i < N; i++) f->H[i] += Hx[i] + tmpf[i];
  memcpy(ly->tn, Hn, sizeof(double)*N); lx_inv(ly, nz, Yf + Mh*nz, ly->nx, ly->m, hw, ly->tn);
  {
    int slot = ly->hk % ly->L;
    double *vms = malloc(sizeof(double)*3*ns), *vmf = malloc(sizeof(double)*N);
    for(i = 0; i < 3*ns; i++){ vms[i] = f->nl ? 0.5*(ly->Vold3[i] + ly->V3[i]) : ly->V3[i]; ly->delta3[i] += vms[i]*dt; }
    for(i = 0; i < N; i++) vmf[i] = f->nl ? 0.5*(f->Vold[i] + f->V[i]) : f->V[i];
    for(c = 0; c < 3; c++) layer_rows_hat(ly, vms + (size_t)c*ns, ly->hsx3 + ((size_t)slot*3 + c)*Mh*nys);
    fault_rows_hat(ly, vmf, ly->hfx + (size_t)slot*Mh*nz);
    ly->hk = (ly->hk + 1) % ly->L;
    if(ly->nhist < ly->L) ly->nhist++;
    free(vms); free(vmf);
  }
  step_end(f);
  if(ly->static_n_fn){ ly->static_n_fn(ly->static_n_ctx, f->delta, tnst); for(i = 0; i < N; i++) ly->tn[i] += tnst[i]; }
  if(ly->normal_coupling) for(i = 0; i < N; i++){ f->sigma[i] = ly->sig0[i] - ly->tn[i]; if(f->sigma[i] < 0.0) f->sigma[i] = 0.0; }
  free(Hx); free(Hn); free(V3); free(Vp); free(rhs); free(tmpf); free(tmps); free(Xs); free(Xf); free(Ys); free(Yf); free(S3); free(P3); free(dend); free(tnst);
  return itot;
}

/* test hook: one surface/fault kernel (see surf_resp) */
void fd3d_surf_kernel(const struct fd3d_tab *t, double dip_deg, int fcomp, int src_fault, int comp, const double c[3], double a1, double a2,
                      const double X[3], int rcv_fault, int L, double *K, double *ks)
{
  double F[3][3], sdir[3]; int i;
  frame_fault(dip_deg, F);
  for(i = 0; i < 3; i++) sdir[i] = F[fcomp - 1][i];
  surf_resp(t, F, src_fault, comp, c, a1, a2, X, rcv_fault, sdir, L, K, ks);
}
