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
    if(comp == 3) L3_all(x1 - c1, x2 - c2, x3, t, cL, cT, q); else L1_all(x1 - c1, x2 - c2, x3, t, cL, cT, q);
    for(k = 0; k < 6; k++) s[k] += sg*q[k];
  }
  for(k = 0; k < 6; k++) s[k] *= -mu/(2.0*cT);
}
/* full stress tensor (11, 22, 33, 12, 31, 32) at (x1, x2, x3) from unit rate since t = 0 on the
   rectangle |xi1| < a1, |xi2| < a2: comp 1 slip along x1, comp 3 opening. Receivers on the
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
      A[(size_t)i*oz + j] = (fd3d_square_ramp(x1, x2, 0.0, T2, a, cL, cT, rho) - fd3d_square_ramp(x1, x2, 0.0, T1, a, cL, cT, rho))/(T2 - T1);
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
        S1[(size_t)i*oz + j] = fd3d_square_ramp(x1, x2, 0.0, l*dt, a, cL, cT, rho);
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

struct fd3d *fd3d_create2(struct fd3d_tab *t, int nzp, int mirror)
{
  struct fd3d *f = (struct fd3d *)calloc(1, sizeof(*f));
  int N = t->nx*nzp;
  if(mirror && t->nz != 2*nzp){ free(f); return NULL; }
  f->nzp = nzp; f->mirror = mirror;
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
  free(f->Vold);
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
  int row0 = f->mirror ? f->nzp : 0, nzp = f->nzp;
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
  return f->t->k1[1][1] + ((f->mirror && j == 0) ? f->t->k1[1][2] : 0.0);
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
            if(!f->mirror) continue;
            if(di == 0) continue;                     /* own image: in kdiag */
            sj = -1 - sj;                             /* image of row -1 - sj */
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

static int step_solve(struct fd3d *f, const double *Hx)
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
    memcpy(V, pass == 0 ? Vold : Vs, sizeof(double)*N);
    if(f->law == 1)
      for(i = 0; i < N; i++) th_p[i] = theta_adv(f, i, f->theta[i], Vold[i], pass == 0 ? Vold[i] : Vs[i], dt);
    for(it = 0; it < 400; it++){
      double ch = 0.0, vm = 0.0;
      apply_k1(f, V, f->work);
      for(i = 0; i < N; i++){
        double Hp = H[i] + hw*f->work[i];
        double kd = f->nl ? 0.5*kdiag(f, i % f->nzp) - 0.5*eta : kdiag(f, i % f->nzp);
        if(f->law == 0){
          /* slip at the collocation time: delta + dt V (first order) or delta + dt (V_old + V)/2 */
          double d0s = f->delta[i] + (f->nl ? 0.5*dt*Vold[i] : 0.0), dm = fabs(d0s);
          double mu = dm < f->d0 ? f->mu_s - (f->mu_s - f->mu_d)*dm/f->d0 : f->mu_d;
          double slope = dm < f->d0 ? f->sigma[i]*(f->mu_s - f->mu_d)*hw*dt/f->d0 : 0.0;
          Vn[i] = local_sw(f, i, Hp, kd, f->sigma[i]*mu, slope);
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
    double *Vm = (double *)malloc(sizeof(double)*N);
    for(i = 0; i < N; i++){
      f->tau[i] = f->H[i] + 0.5*f->work[i] + (0.5*kdiag(f, i % f->nzp) - 0.5*eta)*f->V[i];
      Vm[i] = 0.5*(Vold[i] + f->V[i]);                    /* midpoint rate of the step */
    }
    fd3d_push(f, Vm);
    for(i = 0; i < N; i++) f->delta[i] += Vm[i]*dt;
    free(Vm);
  }else{
    for(i = 0; i < N; i++) f->tau[i] = f->H[i] + f->work[i] + kdiag(f, i % f->nzp)*f->V[i];
    fd3d_push(f, f->V);
    for(i = 0; i < N; i++) f->delta[i] += f->V[i]*dt;
  }
  f->time += dt;
}

int fd3d_step(struct fd3d *f, const double *load)
{
  int it;
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
  int ii, jj, i, j, row0 = f->mirror ? f->nzp : 0;
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

void fd3d_layer_free(struct fd3d_layer *ly)
{
  if(!ly) return;
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
