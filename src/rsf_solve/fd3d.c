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
  free(t->W_hat); free(t->Ks_hat); free(t->work);
  fd_fft_plan_destroy(t->px); fd_fft_plan_destroy(t->ix); fd_fft_plan_destroy(t->pz); fd_fft_plan_destroy(t->iz);
  free(t);
}

/* ---------------- marcher ---------------- */
struct fd3d *fd3d_create(struct fd3d_tab *t)
{
  struct fd3d *f = (struct fd3d *)calloc(1, sizeof(*f));
  int N = t->nx*t->nz;
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
  free(f->hist); free(f->acc); free(f);
}

static void field_hat(const struct fd3d *f, const double *x, double _Complex *out)
{
  const struct fd3d_tab *t = f->t;
  int i, j;
  memset(out, 0, sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  for(i = 0; i < t->nx; i++) for(j = 0; j < t->nz; j++) out[(size_t)i*t->Mz + j] = x[i*t->nz + j];
  fd3d_fft2(t, out, -1);
}

/* tau0 + Ks delta + history (lags >= 2) */
static void known(struct fd3d *f)
{
  const struct fd3d_tab *t = f->t;
  size_t M = (size_t)t->Mx*t->Mz, m;
  int l, i, j, lmax = t->L < f->nhist + 1 ? t->L : f->nhist + 1;
  field_hat(f, f->delta, t->work);
  for(m = 0; m < M; m++) f->acc[m] = t->Ks_hat[m]*t->work[m];
  {
    size_t Mhh = (size_t)t->Mx*t->Mh; int ii, jj;
    double _Complex *ah = (double _Complex *)calloc(Mhh, sizeof(double _Complex));
    for(l = 2; l <= lmax; l++){
      const fd3d_cplx *W = t->W_hat + (size_t)(l - 1)*Mhh;
      const fd3d_cplx *Vh = f->hist + (size_t)((f->hk - (l - 1) + t->L) % t->L)*Mhh;   /* V^(n+1-l) */
      for(m = 0; m < Mhh; m++) ah[m] += (double _Complex)(W[m]*Vh[m]);
    }
    for(ii = 0; ii < t->Mx; ii++) for(jj = 0; jj < t->Mh; jj++) f->acc[(size_t)ii*t->Mz + jj] += ah[(size_t)ii*t->Mh + jj];
    free(ah);
  }
  from_half(t, f->acc);
  fd3d_fft2(t, f->acc, +1);
  for(i = 0; i < t->nx; i++) for(j = 0; j < t->nz; j++)
    f->H[i*t->nz + j] = f->tau0[i*t->nz + j] + creal(f->acc[(size_t)i*t->Mz + j])/M;
}

/* first-lag traction from neighbours (excluding self) */
static void apply_k1(const struct fd3d *f, const double *V, double *out)
{
  const struct fd3d_tab *t = f->t;
  int i, j, di, dj, nx = t->nx, nz = t->nz;
  for(i = 0; i < nx; i++)
    for(j = 0; j < nz; j++){
      double s = 0.0;
      for(di = -1; di <= 1; di++)
        for(dj = -1; dj <= 1; dj++){
          int si = i - di, sj = j - dj;               /* receiver - source = (di, dj) */
          if((di == 0 && dj == 0) || si < 0 || si >= nx || sj < 0 || sj >= nz) continue;
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

int fd3d_step(struct fd3d *f, const double *load)
{
  struct fd3d_tab *t = f->t;
  int N = f->N, i, it, pass, itot = 0;
  double dt = t->dt, k = t->k1[1][1], *Vn = (double *)malloc(sizeof(double)*N);
  double *th_p = NULL, *Vs = NULL;
  size_t M = (size_t)t->Mx*t->Mz;
  known(f);
  if(load) for(i = 0; i < N; i++) f->H[i] += load[i];
  if(f->law == 1){ th_p = (double *)malloc(sizeof(double)*N); Vs = (double *)malloc(sizeof(double)*N); }
  for(pass = 0; pass < (f->law == 1 ? 2 : 1); pass++){
    double *V = (double *)malloc(sizeof(double)*N);
    memcpy(V, pass == 0 ? f->V : Vs, sizeof(double)*N);
    if(f->law == 1)
      for(i = 0; i < N; i++) th_p[i] = theta_adv(f, i, f->theta[i], f->V[i], pass == 0 ? f->V[i] : Vs[i], dt);
    for(it = 0; it < 400; it++){
      double ch = 0.0, vm = 0.0;
      apply_k1(f, V, f->work);
      for(i = 0; i < N; i++){
        double Hp = f->H[i] + f->work[i];
        if(f->law == 0){
          double dm = fabs(f->delta[i]), mu = dm < f->d0 ? f->mu_s - (f->mu_s - f->mu_d)*dm/f->d0 : f->mu_d;
          double slope = dm < f->d0 ? f->sigma[i]*(f->mu_s - f->mu_d)*dt/f->d0 : 0.0;
          Vn[i] = local_sw(f, i, Hp, k, f->sigma[i]*mu, slope);
        }else
          Vn[i] = local_rs(f, i, Hp, k, th_p[i], V[i]);
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
      if(f->law == 1) memcpy(f->theta, th_p, sizeof(double)*N);
    }
    free(V);
  }
  apply_k1(f, f->V, f->work);
  for(i = 0; i < N; i++) f->tau[i] = f->H[i] + f->work[i] + k*f->V[i];
  /* push V^n into the ring */
  field_hat(f, f->V, t->work);
  to_half(t, t->work, f->hist + (size_t)(f->hk % t->L)*t->Mx*t->Mh);
  f->hk = (f->hk + 1) % t->L;
  if(f->nhist < t->L) f->nhist++;
  for(i = 0; i < N; i++) f->delta[i] += f->V[i]*dt;
  f->time += dt;
  free(Vn); free(th_p); free(Vs);
  return itot;
}
