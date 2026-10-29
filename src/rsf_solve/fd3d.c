/* 3-D fully dynamic extension: kernels, tables, marcher; see fd3d.h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd3d.h"

/* ---------------- Tada (2006) Appendix A ---------------- */
static double H(double x)
{
  return x > 0.0 ? 1.0 : 0.0;
}
static double sgn(double x)
{
  return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
}
static double sq(double a)
{
  return a > 0.0 ? sqrt(a) : 0.0;
}

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
static double g4(double R, double t, double c)
{
  return (t*t*t - R*R/(c*c)*t)*H(t - R/c);
}

double fd3d_L31(double x1, double x2, double x3, double t, double cL, double cT)
{
  double p = cT/cL, chi1 = sqrt(x1*x1 + x3*x3), chi2 = sqrt(x2*x2 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double c3 = cT*cT*cT, chi1_2 = chi1*chi1, chi1_4 = chi1_2*chi1_2, chi1_6 = chi1_4*chi1_2, R5 = R*R*R*R*R;
  double v = H(x1)*H(x2)*H(t - fabs(x3)/cT) - (g3(x1, x2, x3, t, cT) + g3(x2, x1, x3, t, cT))/(2.0*M_PI) - 2.0*c3*x1*(3.0*x3*x3 - x1*x1)/(3.0*M_PI*chi1_6)*(g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT));
  if(x3 != 0.0)
    v += - c3*x1*x2*x3*x3/(M_PI*R5*chi1_2)*(g4(R, t, cL) - g4(R, t, cT)) - 2.0*cT*x1*x3*x3/(M_PI*chi1_4)*(p*p*g2(x1, x2, x3, t, cL) - g2(x1, x2, x3, t, cT));
  v += cT*x2/(2.0*M_PI*chi2*chi2)*g2(x2, x1, x3, t, cT);
  return v;
}

/* full quadrant Green's functions (Tada 2006, A18-A27): L_ij/1 (slip along x1) and
   L_ij/3 (opening), returned as s[6] = (11, 22, 33, 12, 31, 32) */
static double g5(double R, double t, double c)
{
  return (t - R/c)*H(t - R/c);
}
static void L1_all(double x1, double x2, double x3, double t, double cL, double cT, double s[6])
{
  double p = cT/cL, chi1 = sqrt(x1*x1 + x3*x3), chi2 = sqrt(x2*x2 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double c3 = cT*cT*cT, chi1_2 = chi1*chi1, R5 = R*R*R*R*R, pi = M_PI;
  double dg1 = g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT), dg4 = g4(R, t, cL) - g4(R, t, cT);
  double g2L = g2(x1, x2, x3, t, cL), g2T = g2(x1, x2, x3, t, cT), g2sT = g2(x2, x1, x3, t, cT);
  s[0] = -2.0*c3*x3*(3.0*x1*x1 - x3*x3)/(3.0*pi*chi1_2*chi1_2*chi1_2)*dg1 - c3*x1*x1*x2*x3/(pi*R5*chi1_2)*dg4 - cT*x3/(pi*chi1_2)*((1.0 - 2.0*p*p + 2.0*p*p*x1*x1/chi1_2)*g2L + (1.0 - 2.0*x1*x1/chi1_2)*g2T);
  s[3] = c3*x1*x3/(pi*R5)*dg4 - cT*x3/(2.0*pi*chi2*chi2)*g2sT;
  s[1] = c3*x2*x3/(pi*R5)*dg4 - (1.0 - 2.0*p*p)*cT*x3/(pi*chi1_2)*g2L;
  s[4] = fd3d_L31(x1, x2, x3, t, cL, cT);
  s[5] = c3/(pi*R5)*(x3*x3 - R*R/3.0)*dg4 + cT/(pi*R)*(2.0*p*p/3.0*g5(R, t, cL) - g5(R, t, cT)/6.0);
  s[2] = 2.0*c3*x3*(3.0*x1*x1 - x3*x3)/(3.0*pi*chi1_2*chi1_2*chi1_2)*dg1 - c3*x2*x3*x3*x3/(pi*R5*chi1_2)*dg4 - cT*x3/(pi*chi1_2)*((1.0 - 2.0*p*p + 2.0*p*p*x3*x3/chi1_2)*g2L + (1.0 - 2.0*x3*x3/chi1_2)*g2T);
}
static double L11_3(double x1, double x2, double x3, double t, double cL, double cT)
{
  double p = cT/cL, chi1 = sqrt(x1*x1 + x3*x3), chi2 = sqrt(x2*x2 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double c3 = cT*cT*cT, chi1_2 = chi1*chi1, R5 = R*R*R*R*R, pi = M_PI;
  double dg1 = g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT), dg4 = g4(R, t, cL) - g4(R, t, cT);
  return (1.0 - 2.0*p*p)/p*H(x1)*H(x2)*H(t - fabs(x3)/cL) - (1.0 - 2.0*p*p)/(2.0*pi*p)*(g3(x1, x2, x3, t, cL) + g3(x2, x1, x3, t, cL)) - 2.0*c3*x1*(3.0*x3*x3 - x1*x1)/(3.0*pi*chi1_2*chi1_2*chi1_2)*dg1 - c3*x1*x2*x3*x3/(pi*R5*chi1_2)*dg4 - 2.0*cT*x1*x3*x3/(pi*chi1_2*chi1_2)*(p*p*g2(x1, x2, x3, t, cL) - g2(x1, x2, x3, t, cT)) + (1.0 - 2.0*p*p)*cT*x2/(pi*chi2*chi2)*g2(x2, x1, x3, t, cL);
}
static void L3_all(double x1, double x2, double x3, double t, double cL, double cT, double s[6])
{
  double p = cT/cL, chi1 = sqrt(x1*x1 + x3*x3), chi2 = sqrt(x2*x2 + x3*x3), R = sqrt(x1*x1 + x2*x2 + x3*x3);
  double c3 = cT*cT*cT, R5 = R*R*R*R*R, pi = M_PI, dg4 = g4(R, t, cL) - g4(R, t, cT), a[6], b[6];
  s[0] = L11_3(x1, x2, x3, t, cL, cT);
  s[1] = L11_3(x2, x1, x3, t, cL, cT);
  s[3] = c3/(pi*R5)*(x3*x3 - R*R/3.0)*dg4 + cT/(pi*R)*((1.0 - 4.0*p*p/3.0)*g5(R, t, cL) - 2.0/3.0*g5(R, t, cT));
  L1_all(x1, x2, x3, t, cL, cT, a);
  s[4] = a[2];          /* L_31/3 = L_33/1 */
  L1_all(x2, x1, x3, t, cL, cT, b);
  s[5] = b[2];          /* L_32/3 by symmetry */
  s[2] = 1.0/p*H(x1)*H(x2)*H(t - fabs(x3)/cL) - 1.0/(2.0*pi*p)*(g3(x1, x2, x3, t, cL) + g3(x2, x1, x3, t, cL)) + 2.0*c3*x1*(3.0*x3*x3 - x1*x1)/(3.0*pi*pow(chi1, 6))*(g1(x1, x2, x3, t, cL) - g1(x1, x2, x3, t, cT)) + 2.0*c3*x2*(3.0*x3*x3 - x2*x2)/(3.0*pi*pow(chi2, 6))*(g1(x2, x1, x3, t, cL) - g1(x2, x1, x3, t, cT)) + c3*x1*x2*x3*x3/(pi*R5)*(1.0/(chi1*chi1) + 1.0/(chi2*chi2))*dg4 + cT*x1/(pi*chi1*chi1)*((2.0*(1.0 - p*p) + 2.0*p*p*x3*x3/(chi1*chi1))*g2(x1, x2, x3, t, cL) - 2.0*x3*x3/(chi1*chi1)*g2(x1, x2, x3, t, cT)) + cT*x2/(pi*chi2*chi2)*((2.0*(1.0 - p*p) + 2.0*p*p*x3*x3/(chi2*chi2))*g2(x2, x1, x3, t, cL) - 2.0*x3*x3/(chi2*chi2)*g2(x2, x1, x3, t, cT));
}
static void rect_tensor0(double x1, double x2, double x3, double t, double a1, double a2, double cL, double cT, double rho, int comp, double s[6])
{
  double mu = rho*cT*cT, q[6];
  int k, i1, i2;
  double c1;
  double c2;
  double sg;
  double r[6];
  for(k = 0; k < 6; k++)
    s[k] = 0.0;
  for(i1 = 0; i1 < 2; i1++)
    for(i2 = 0; i2 < 2; i2++){
      c1 = i1 ? a1 : -a1;
      c2 = i2 ? a2 : -a2;
      sg = (i1 == i2) ? 1.0 : -1.0;
      if(comp == 3)
        L3_all(x1 - c1, x2 - c2, x3, t, cL, cT, q);
      else if(comp == 2){
        /* slip along x2: L_ij/2(x1, x2) = L_i'j'/1(x2, x1) with the indices 1 <-> 2 swapped */

        L1_all(x2 - c2, x1 - c1, x3, t, cL, cT, r);
        q[0] = r[1];
        q[1] = r[0];
        q[2] = r[2];
        q[3] = r[3];
        q[4] = r[5];
        q[5] = r[4];
      }else
        L1_all(x1 - c1, x2 - c2, x3, t, cL, cT, q);
      for(k = 0; k < 6; k++)
        s[k] += sg*q[k];
    }
  for(k = 0; k < 6; k++)
    s[k] *= -mu/(2.0*cT);
}
/* full stress tensor (11, 22, 33, 12, 31, 32) at (x1, x2, x3) from unit rate since t = 0 on the
   rectangle |xi1| < a1, |xi2| < a2: comp 1 slip along x1, comp 2 slip along x2, comp 3 opening. Receivers on the
   extension of an element edge hit removable 0/0: averaged over a tiny symmetric offset */
void fd3d_rect_ramp_tensor(double x1, double x2, double x3, double t, double a1, double a2,
    double cL, double cT, double rho, int comp, double s[6])
{
  double e = 1e-9*(a1 > a2 ? a1 : a2), lo[6], hi[6];
  int k;
  /* in the source plane the quadrant fields are singular along the shared edge line and cancel
     only between quadrants: a larger offset (1e-3 a) avoids the cancellation, off the plane the
     0/0 is removable and 1e-7 a suffices (same rule as tada_quad.py) */
  double d = (fabs(x3) < e ? 1e-3 : 1e-7)*(a1 > a2 ? a1 : a2);
  if((fabs(fabs(x1) - a1) < e) || (fabs(fabs(x2) - a2) < e)){
    rect_tensor0(x1 - d, x2 - d, x3, t, a1, a2, cL, cT, rho, comp, lo);
    rect_tensor0(x1 + d, x2 + d, x3, t, a1, a2, cL, cT, rho, comp, hi);
    for(k = 0; k < 6; k++)
      s[k] = 0.5*(lo[k] + hi[k]);
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
    for(i = 0; i < Mx; i++)
      col[i] = a[(size_t)i*Mz + j];
    fd_fft_exec(sign < 0 ? t->px : t->ix, col, out);
    for(i = 0; i < Mx; i++)
      a[(size_t)i*Mz + j] = out[i];
  }
  free(col);
  free(out);
}

/* rebuild the full spectrum of a real field from its half (double accumulator) */
static void from_half(const struct fd3d_tab *t, double _Complex *full)
{
  int i, j, Mx = t->Mx, Mz = t->Mz;
  for(i = 0; i < Mx; i++)
    for(j = Mz/2 + 1; j < Mz; j++)
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
  if(comp == 1)
    return fd3d_square_ramp(x1, x2, 0.0, t, a, cL, cT, rho);
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
  if(t)
    t->slipcomp = comp;
  return t;
}
#ifdef FD3D_MPI
static int g_fd3d_use_mpi = 0;
static MPI_Comm g_fd3d_comm;
void fd3d_set_comm(MPI_Comm comm)
{
  g_fd3d_comm = comm;
  g_fd3d_use_mpi = 1;
}
#endif
/* row partition of the half spectra over the ranks (block rows of kx) */
static void tab_partition(struct fd3d_tab *t)
{
  int r, size = 1, rank = 0;
  int a;
  int b;
#ifdef FD3D_MPI
  if(g_fd3d_use_mpi){
    t->comm = g_fd3d_comm;
    MPI_Comm_size(t->comm, &size);
    MPI_Comm_rank(t->comm, &rank);
  } else
    t->comm = MPI_COMM_SELF;
#endif
  t->mpi_rank = rank;
  t->mpi_size = size;
  t->rcounts = (int *)malloc(sizeof(int)*size);
  t->rdispls = (int *)malloc(sizeof(int)*size);
  for(r = 0; r < size; r++){
    a = (int)(((long)t->Mx*r)/size);
    b = (int)(((long)t->Mx*(r + 1))/size);
    t->rcounts[r] = (b - a)*t->Mh;
    t->rdispls[r] = a*t->Mh;
    if(r == rank){
      t->xr0 = a;
      t->nxr = b - a;
    }
  }
}
/* gather the row blocks of a half spectrum held in place (own rows filled) on all ranks */
static void tab_allgather(const struct fd3d_tab *t, double _Complex *half)
{
#ifdef FD3D_MPI
  if(t->mpi_size > 1)
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, half, t->rcounts, t->rdispls, MPI_C_DOUBLE_COMPLEX, t->comm);
#else
  (void)t;
  (void)half;
#endif
}

struct fd3d_tab *fd3d_tab_create(int nx, int nz, double dx, double dt, double cL, double cT, double rho, int lmax)
{
  size_t jm;
  struct fd3d_tab *t = (struct fd3d_tab *)calloc(1, sizeof(*t));
  int ox = 2*nx - 1, oz = 2*nz - 1, i, j, l, Lmax = 0;
  double a = dx/2, *S0, *S1, *A, rmax, T1, T2;
  double x1;
  double x2;
  int P;
  int r;
  int la;
  int lb;
  int nl;
  int q;
  size_t Mh;
  _Complex double *mine;
  double *Ks;
  double x1_2;
  double x2_2;
  int di;
  int dj;
  int ii;
  int jj;
  _Complex double *dst;
  int *sc;
  int *sd;
  int *rc;
  int *rd;
  _Complex double *sbuf;
  _Complex double *rbuf;
  size_t o;
  int qa;
  int qn;
  int ll;
  int ra;
  int rb;
  t->nx = nx;
  t->nz = nz;
  t->dx = dx;
  t->dt = dt;
  t->cL = cL;
  t->cT = cT;
  t->rho = rho;
  t->Mx = fd_fft_size_for(ox);
  t->Mz = fd_fft_size_for(oz);
  t->px = fd_fft_plan_create(t->Mx, FD_FFT_FORWARD);
  t->ix = fd_fft_plan_create(t->Mx, FD_FFT_BACKWARD);
  t->pz = fd_fft_plan_create(t->Mz, FD_FFT_FORWARD);
  t->iz = fd_fft_plan_create(t->Mz, FD_FFT_BACKWARD);
  t->work = (double _Complex *)malloc(sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  for(i = -(nx - 1); i <= nx - 1; i++)
    for(j = -(nz - 1); j <= nz - 1; j++){
      rmax = sqrt((fabs(i*dx) + a)*(fabs(i*dx) + a) + (fabs(j*dx) + a)*(fabs(j*dx) + a));
      l = (int)ceil(rmax/cT/dt) + 1;
      if(l > Lmax)
        Lmax = l;
    }
  t->L = Lmax + 1;
  if(lmax > 0 && t->L > lmax)
    t->L = lmax;
  S0 = (double *)calloc((size_t)ox*oz, sizeof(double));
  S1 = (double *)malloc(sizeof(double)*(size_t)ox*oz);
  A = (double *)malloc(sizeof(double)*(size_t)ox*oz);
  /* static: late-time slope */
  rmax = sqrt(2.0)*(nx > nz ? nx : nz)*dx + dx;
  T1 = 1.5*rmax/cT + 10.0*dt;
  T2 = 2.0*T1;
  for(i = 0; i < ox; i++)
    for(j = 0; j < oz; j++){
      x1 = (i - (nx - 1))*dx;
      x2 = (j - (nz - 1))*dx;
      A[(size_t)i*oz + j] = (tab_kernel(tab_comp, x1, x2, T2, a, cL, cT, rho) - tab_kernel(tab_comp, x1, x2, T1, a, cL, cT, rho))/(T2 - T1);
    }
  t->ks_self = A[(size_t)(nx - 1)*oz + (nz - 1)];
  t->ks_off = (double *)malloc(sizeof(double)*(size_t)ox*oz);
  memcpy(t->ks_off, A, sizeof(double)*(size_t)ox*oz);
  t->Ks_hat = (double _Complex *)malloc(sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  put_offsets(t, A, t->Ks_hat);
  fd3d_fft2(t, t->Ks_hat, -1);
  t->Mh = t->Mz/2 + 1;
  tab_partition(t);
  t->W_hat = (fd3d_cplx *)malloc(sizeof(fd3d_cplx)*(size_t)t->L*t->nxr*t->Mh);
  {
    /* lags in blocks over the ranks: rank r computes lags la .. lb - 1 (W[l] = S(l) - S(l - 1) -
       Ks dt, so it also evaluates S(la - 1)), transforms them and sends every rank its kx rows;
       with one rank this is the serial loop */
    P = t->mpi_size;
    Mh = (size_t)t->Mh;
    /* own lags, full half spectra: nl x Mx x Mh */
    la = 1 + (int)(((long)t->L*t->mpi_rank)/P);
    lb = 1 + (int)(((long)t->L*(t->mpi_rank + 1))/P);
    nl = lb - la;
    mine = (double _Complex *)malloc(sizeof(double _Complex)*((size_t)(nl > 0 ? nl : 1))*t->Mx*Mh);
    if(la > 1)
      for(i = 0; i < ox; i++)
        for(j = 0; j < oz; j++)
          S0[(size_t)i*oz + j] = tab_kernel(tab_comp, (i - (nx - 1))*dx, (j - (nz - 1))*dx, (la - 1)*dt, a, cL, cT, rho);
    for(l = la; l < lb; l++){
      Ks = (double *)malloc(sizeof(double)*(size_t)ox*oz);
      memcpy(Ks, A, sizeof(double)*(size_t)ox*oz);
      for(i = 0; i < ox; i++)
        for(j = 0; j < oz; j++){
          x1_2 = (i - (nx - 1))*dx;
          x2_2 = (j - (nz - 1))*dx;
          S1[(size_t)i*oz + j] = tab_kernel(tab_comp, x1_2, x2_2, l*dt, a, cL, cT, rho);
        }
      if(l == 1){
        for(di = -1; di <= 1; di++)
          for(dj = -1; dj <= 1; dj++)
            t->k1[di + 1][dj + 1] = (nx > 1 || di == 0) && (nz > 1 || dj == 0) ? S1[(size_t)(di + nx - 1)*oz + (dj + nz - 1)] : 0.0;
      }
      for(i = 0; i < ox*oz; i++)
        Ks[i] = S1[i] - S0[i] - Ks[i]*dt;      /* W[l] */
      put_offsets(t, Ks, t->work);
      fd3d_fft2(t, t->work, -1);
      {
        dst = mine + (size_t)(l - la)*t->Mx*Mh;
        for(ii = 0; ii < t->Mx; ii++)
          for(jj = 0; jj < t->Mh; jj++)
            dst[(size_t)ii*Mh + jj] = t->work[(size_t)ii*t->Mz + jj];
      }
      memcpy(S0, S1, sizeof(double)*(size_t)ox*oz);
      free(Ks);
    }
    if(P == 1){
      for(l = 0; l < t->L; l++)
        for(jm = 0; jm < (size_t)t->Mx*Mh; jm++)
          t->W_hat[(size_t)l*t->Mx*Mh+ jm] = (fd3d_cplx)mine[(size_t)l*t->Mx*Mh+ jm];
    }
#ifdef FD3D_MPI
    else{
      /* all-to-all: to rank q, for each own lag, q's rows; from rank s, for each of s's lags, own rows */
      sc = calloc(P, sizeof(int));
      sd = calloc(P, sizeof(int));
      rc = calloc(P, sizeof(int));
      rd = calloc(P, sizeof(int));
      sbuf = malloc(sizeof(double _Complex)*((size_t)(nl > 0 ? nl : 1))*t->Mx*Mh);
      rbuf = malloc(sizeof(double _Complex)*(size_t)t->L*t->nxr*Mh);
      o = 0;
      for(q = 0; q < P; q++){
        qa = t->rdispls[q]/t->Mh;
        qn = t->rcounts[q]/t->Mh;
        sd[q] = (int)o;
        for(ll = 0; ll < nl; ll++){
          memcpy(sbuf + o, mine + ((size_t)ll*t->Mx + qa)*Mh, sizeof(double _Complex)*qn*Mh);
          o += (size_t)qn*Mh;
        }
        sc[q] = (int)o - sd[q];
      }
      o = 0;
      for(r = 0; r < P; r++){
        ra = 1 + (int)(((long)t->L*r)/P);
        rb = 1 + (int)(((long)t->L*(r + 1))/P);
        rd[r] = (int)o;
        rc[r] = (rb - ra)*t->nxr*t->Mh;
        o += (size_t)rc[r];
      }
      MPI_Alltoallv(sbuf, sc, sd, MPI_C_DOUBLE_COMPLEX, rbuf, rc, rd, MPI_C_DOUBLE_COMPLEX, t->comm);
      for(jm = 0; jm < (size_t)t->L*t->nxr*Mh; jm++)
        t->W_hat[jm] = (fd3d_cplx)rbuf[jm];   /* lags in order, own rows */
      /* the first-lag stencil comes from the rank that owns lag 1 */
      MPI_Bcast(&t->k1[0][0], 9, MPI_DOUBLE, 0, t->comm);
      free(sc);
      free(sd);
      free(rc);
      free(rd);
      free(sbuf);
      free(rbuf);
    }
#endif
    free(mine);
  }
  free(S0);
  free(S1);
  free(A);
  return t;
}

void fd3d_tab_free(struct fd3d_tab *t)
{
  if(!t)
    return;
  free(t->W_hat);
  free(t->Ks_hat);
  free(t->work);
  free(t->ks_off);
  free(t->rcounts);
  free(t->rdispls);
  fd_fft_plan_destroy(t->px);
  fd_fft_plan_destroy(t->ix);
  fd_fft_plan_destroy(t->pz);
  fd_fft_plan_destroy(t->iz);
  free(t);
}

/* ---------------- marcher ---------------- */
struct fd3d *fd3d_create(struct fd3d_tab *t)
{
  return fd3d_create2(t, t->nz, 0);
}

struct fd3d *fd3d_create3(struct fd3d_tab *t, int nzp, int mirror, int next);
struct fd3d *fd3d_create2(struct fd3d_tab *t, int nzp, int mirror)
{
  return fd3d_create3(t, nzp, mirror, 0);
}
/* next > 0: the fault plane continues above its top row by next rows whose slip rate copies row 0
   (a surface-breaking fault without an edge in the surface; the table must have nzp + next rows) */
struct fd3d *fd3d_create3(struct fd3d_tab *t, int nzp, int mirror, int next)
{
  struct fd3d *f = (struct fd3d *)calloc(1, sizeof(*f));
  int N = t->nx*nzp;
  int r;
  int size;
  int a;
  int b;
  if(mirror && t->nz != 2*nzp){
    free(f);
    return NULL;
  }
  if(next && (mirror || t->nz != nzp + next)){
    free(f);
    return NULL;
  }
  f->nzp = nzp;
  f->mirror = mirror;
  f->next = next;
  f->nk = getenv("FD3D_NK") ? atoi(getenv("FD3D_NK")) : 0;   /* tests; rsf_solve sets it from -fd_nk */
  {                                                /* column blocks for the local solves (MPI) */
    size = t->mpi_size;
    f->ccounts = (int *)malloc(sizeof(int)*size);
    f->cdispls = (int *)malloc(sizeof(int)*size);
    for(r = 0; r < size; r++){
      a = (int)(((long)t->nx*r)/size);
      b = (int)(((long)t->nx*(r + 1))/size);
      f->ccounts[r] = (b - a)*nzp;
      f->cdispls[r] = a*nzp;
      if(r == t->mpi_rank){
        f->c0 = a;
        f->nc = b - a;
      }
    }
  }
  f->kdamp = getenv("FD3D_KDAMP") ? atof(getenv("FD3D_KDAMP")) : 0.5;
  size_t M = (size_t)t->Mx*t->Mz;
  f->t = t;
  f->N = N;
  f->delta = calloc(N, sizeof(double));
  f->V = calloc(N, sizeof(double));
  f->tau = calloc(N, sizeof(double));
  f->tau0 = calloc(N, sizeof(double));
  f->sigma = calloc(N, sizeof(double));
  f->H = calloc(N, sizeof(double));
  f->work = calloc(N, sizeof(double));
  f->hist = (fd3d_cplx *)calloc((size_t)t->L*t->nxr*t->Mh, sizeof(fd3d_cplx));
  f->acc = (double _Complex *)malloc(sizeof(double _Complex)*M);
  return f;
}

void fd3d_free(struct fd3d *f)
{
  int q;
  if(!f)
    return;
  free(f->delta);
  free(f->V);
  free(f->tau);
  free(f->tau0);
  free(f->sigma);
  free(f->H);
  free(f->work);
  {
    for(q = 0; q < 8; q++){
      free(f->kop_B[q]);
      free(f->kop_C[q]);
    }
  }
  free(f->Vold);
  free(f->dpath);
  free(f->coh);
  free(f->ccounts);
  free(f->cdispls);
  if(f->layer)
    fd3d_layer_free(f->layer);
  free(f->hist);
  free(f->acc);
  free(f);
}

static void field_hat(const struct fd3d *f, const double *x, double _Complex *out)
{
  const struct fd3d_tab *t = f->t;
  int i, j;
  int nzp = f->nzp;
  int r;
  memset(out, 0, sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  for(i = 0; i < t->nx; i++)
    for(j = 0; j < nzp; j++){
      if(f->mirror){
        out[(size_t)i*t->Mz + nzp + j] = x[i*nzp + j];
        out[(size_t)i*t->Mz + nzp - 1 - j] = x[i*nzp + j];      /* image */
      }else if(f->next){
        out[(size_t)i*t->Mz + f->next + j] = x[i*nzp + j];
        if(j == 0){
          for(r = 0; r < f->next; r++)
            out[(size_t)i*t->Mz + r] = x[i*nzp];
        }   /* extension rows copy row 0 */
      }else
        out[(size_t)i*t->Mz + j] = x[i*nzp + j];
    }
  fd3d_fft2(t, out, -1);
}

/* tau0 + Ks delta + history (lags >= 2) */
static void known_k(struct fd3d *f, int k);
static void known(struct fd3d *f)
{
  known_k(f, 1);
}

/* tau0 + Ks delta + history for a step of k substeps: the most recent stored
   slice is k + 1 lags old at the end of the step */
static void known_k(struct fd3d *f, int k)
{
  const struct fd3d_tab *t = f->t;
  size_t M = (size_t)t->Mx*t->Mz, m;
  int l, i, j, lmax = t->L < f->nhist + k ? t->L : f->nhist + k;
  int row0 = f->mirror ? f->nzp : f->next, nzp = f->nzp;
  size_t Mhh;
  size_t Mho;
  int ii;
  int jj;
  _Complex double *ah;
  _Complex double *ao;
  const fd3d_cplx *W;
  const fd3d_cplx *Vh;
  if(f->static_fn)
    memset(f->acc, 0, sizeof(double _Complex)*M);
  else{
    field_hat(f, f->delta, t->work);
    for(m = 0; m < M; m++)
      f->acc[m] = t->Ks_hat[m]*t->work[m];
  }
  {
    /* own kx rows: history sum over the lags; then the row blocks of all ranks */
    Mhh = (size_t)t->Mx*t->Mh;
    Mho = (size_t)t->nxr*t->Mh;
    ah = (double _Complex *)calloc(Mhh, sizeof(double _Complex));
    ao = ah + (size_t)t->xr0*t->Mh;
    for(l = k + 1; l <= lmax; l++){
      W = t->W_hat + (size_t)(l - 1)*Mho;
      Vh = f->hist + (size_t)((f->hk - (l - k) + 2*t->L) % t->L)*Mho;
      for(m = 0; m < Mho; m++)
        ao[m] += (double _Complex)(W[m]*Vh[m]);
    }
    tab_allgather(t, ah);
    for(ii = 0; ii < t->Mx; ii++)
      for(jj = 0; jj < t->Mh; jj++)
        f->acc[(size_t)ii*t->Mz + jj] += ah[(size_t)ii*t->Mh + jj];
    free(ah);
  }
  from_half(t, f->acc);
  fd3d_fft2(t, f->acc, +1);
  for(i = 0; i < t->nx; i++)
    for(j = 0; j < nzp; j++)
      f->H[i*nzp + j] = f->tau0[i*nzp + j] + creal(f->acc[(size_t)i*t->Mz + row0 + j])/M;
  if(f->static_fn){
    f->static_fn(f->static_ctx, f->delta, f->work);
    for(i = 0; i < f->N; i++)
      f->H[i] += f->work[i];
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
  double s;
  int si;
  int sj;
  for(i = 0; i < nx; i++)
    for(j = 0; j < nz; j++){
      s = 0.0;
      for(di = -1; di <= 1; di++)
        for(dj = -1; dj <= 1; dj++){
          si = i - di;
          sj = j - dj;               /* receiver - source = (di, dj) */
          if((di == 0 && dj == 0) || si < 0 || si >= nx || sj >= nz)
            continue;
          if(sj < 0){
            if(!f->mirror && !f->next)
              continue;
            if(di == 0)
              continue;                     /* own image or copy: in kdiag */
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
  if(v < 1e-30)
    v = 1e-30;
  return f->Lc[i]/v + (th - f->Lc[i]/v)*exp(-v*dt/f->Lc[i]);
}

/* local solve at element i: Hp + k V = strength(V) sign(V), scalar */
long fd3d_cnt[3];                                  /* profiling: Jacobi sweeps, local RS solves, Newton iterations */
static double local_sw(const struct fd3d *f, int i, double Hp, double k, double str0, double slope)
{
  double A = -k - slope;
  if(f->locked && f->locked[i])
    return 0.0;
  if(fabs(Hp) <= str0)
    return 0.0;
  return sgn(Hp)*(fabs(Hp) - str0)/A;
}
static double local_rs(const struct fd3d *f, int i, double Hp, double k, double theta, double vguess)
{
  double h = fabs(Hp), A = -k, lo = -80.0, hi, u, un, v, F, g, dg;
  int it;
  if(f->locked && f->locked[i])
    return 0.0;
  if(h <= 0.0)
    return 0.0;
  hi = log(h/A) + 1.0;
  u = log(fabs(vguess) > 1e-30 ? fabs(vguess) : 1e-30);
  if(u > hi - 0.5)
    u = hi - 0.5;
  if(u < lo + 1.0)
    u = lo + 1.0;
  for(it = 0; it < 100; it++){
    v = exp(u);
    F = rs_strength(f, i, v, theta);
    g = h - F - A*v;                                   /* decreasing in v */
    dg = -(rs_dstrength(f, i, v, theta) + A)*v;
    if(g > 0)
      lo = u;
    else
      hi = u;
    un = (dg != 0.0) ? u - g/dg : 0.5*(lo + hi);
    if(!isfinite(un) || un <= lo || un >= hi)
      un = 0.5*(lo + hi);
    if(fabs(un - u) < 1e-13){
      u = un;
      break;
    }
    u = un;
  }
  fd3d_cnt[1]++;
  fd3d_cnt[2] += it + 1;
  return sgn(Hp)*exp(u);
}

/* push one fine-step slice of slip rate into the history ring */
void fd3d_push(struct fd3d *f, const double *V)
{
  struct fd3d_tab *t = f->t;
  int ii;
  int jj;
  fd3d_cplx *dst;
  field_hat(f, V, t->work);
  {                                                /* own kx rows */
    dst = f->hist + (size_t)(f->hk % t->L)*t->nxr*t->Mh;
    for(ii = 0; ii < t->nxr; ii++)
      for(jj = 0; jj < t->Mh; jj++)
        dst[(size_t)ii*t->Mh + jj] = (fd3d_cplx)t->work[(size_t)(t->xr0 + ii)*t->Mz + jj];
  }
  f->hk = (f->hk + 1) % t->L;
  if(f->nhist < t->L)
    f->nhist++;
}

/* empty history, zero slip and time (start of a dynamic phase) */
void fd3d_reset(struct fd3d *f)
{
  struct fd3d_tab *t = f->t;
  struct fd3d_layer *ly = f->layer;
  if(ly && ly->surf){                             /* the surface starts from rest at each phase */
    memset(ly->V3, 0, sizeof(double)*3*ly->ns);
    memset(ly->Vold3, 0, sizeof(double)*3*ly->ns);
    memset(ly->delta3, 0, sizeof(double)*3*ly->ns);
    memset(ly->hsx3, 0, sizeof(double _Complex)*3*(size_t)ly->L*ly->Mh*ly->nys);
    memset(ly->hfx, 0, sizeof(double _Complex)*(size_t)ly->L*ly->Mh*ly->nz);
    memset(ly->tn, 0, sizeof(double)*f->N);
    ly->hk = 0;
    ly->nhist = 0;
  }else if(ly){                                    /* the layer starts from rest at each phase */
    memset(ly->V, 0, sizeof(double)*ly->ns);
    memset(ly->Vold, 0, sizeof(double)*ly->ns);
    memset(ly->delta, 0, sizeof(double)*ly->ns);
    if(ly->fft){
      memset(ly->hs2, 0, sizeof(double _Complex)*(size_t)ly->L*ly->Mx2*ly->Mh2);
      memset(ly->hsx, 0, sizeof(double _Complex)*(size_t)ly->L*ly->Mh*ly->nys);
      memset(ly->hfx, 0, sizeof(double _Complex)*(size_t)ly->L*ly->Mh*ly->nz);
    }else{
      memset(ly->hs, 0, sizeof(double)*(size_t)ly->L*ly->ns);
      memset(ly->hf, 0, sizeof(double)*(size_t)ly->L*f->N);
    }
    ly->hk = 0;
    ly->nhist = 0;
  }
  memset(f->hist, 0, sizeof(fd3d_cplx)*(size_t)t->L*t->nxr*t->Mh);
  f->hk = 0;
  f->nhist = 0;
  f->time = 0.0;
  memset(f->delta, 0, sizeof(double)*f->N);
}

/* one fine step, split in three parts so that the free-surface layer can iterate the
   fault solve with an extra known traction (step_solve's Hx) before finishing */
#include <time.h>
static double fd3d_now(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9*ts.tv_nsec;
}
double fd3d_prof3[4];                              /* known part (history), solve, finish (incl. push), steps */
static void step_begin(struct fd3d *f, const double *load)
{
  struct fd3d_tab *t = f->t;
  int N = f->N, i;
  double eta = 0.5*t->rho*t->cT;
  known(f);
  if(load)
    for(i = 0; i < N; i++)
      f->H[i] += load[i];
  if(!f->Vold)
    f->Vold = (double *)malloc(sizeof(double)*N);
  memcpy(f->Vold, f->V, sizeof(double)*N);
  /* NL: neighbour and self terms are halved (the step's midpoint rate is (V_old + V)/2),
     the V_old half goes into H, the self term becomes kdiag/2 - eta/2 */
  if(f->nl){
    apply_k1(f, f->Vold, f->work);
    for(i = 0; i < N; i++)
      f->H[i] += 0.5*(f->work[i] + (kdiag(f, i % f->nzp) + eta)*f->Vold[i]);
  }
}

static int step_solve_g(struct fd3d *f, const double *Hx, const double *guess, const char *mask);
static int step_solve(struct fd3d *f, const double *Hx)
{
  return step_solve_g(f, Hx, NULL, NULL);
}
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
  double *V;
  char *chg;
  char *act;
  int use_act;
  int nzp;
  int nxx;
  int i_lo;
  int i_hi;
  double ch;
  double vm;
  int ii;
  int jj;
  int di;
  int dj;
  double Hp;
  double kd;
  double dm;
  double sg;
  double coh;
  double mu;
  double slope;
  double vc;
  double sd;
  double vs;
  for(i = 0; i < N; i++)
    H[i] = f->H[i] + (Hx ? Hx[i] : 0.0);
  if(f->vpres){                                   /* prescribed slip rates (tests) */
    memcpy(f->V, f->vpres, sizeof(double)*N);
    free(Vn);
    free(H);
    return 0;
  }
  if(f->law == 1){
    th_p = (double *)malloc(sizeof(double)*N);
    Vs = (double *)malloc(sizeof(double)*N);
  }
  for(pass = 0; pass < (f->law == 1 ? 2 : 1); pass++){
    V = (double *)malloc(sizeof(double)*N);
    memcpy(V, pass == 0 ? (guess ? guess : Vold) : Vs, sizeof(double)*N);
    if(f->law == 1)
      for(i = 0; i < N; i++)
        th_p[i] = theta_adv(f, i, f->theta[i], Vold[i], pass == 0 ? Vold[i] : Vs[i], dt);
    /* active set: after the first sweep, a cell is recomputed only if its own rate or a stencil
       neighbour's changed in the previous sweep (otherwise its known part is unchanged and so is
       its solution); this removes most of the cost of the later sweeps */
    chg = (char *)malloc(N);
    act = (char *)malloc(N);
    use_act = !getenv("FD3D_NO_ACTIVE");
    nzp = f->nzp;
    nxx = t->nx;
    i_lo = f->c0*nzp;
    i_hi = (f->c0 + f->nc)*nzp;      /* own columns (all without MPI) */
    memset(chg, 1, N);
    for(it = 0; it < 400; it++){
      ch = 0.0;
      vm = 0.0;
      fd3d_cnt[0]++;
      apply_k1(f, V, f->work);
      if(use_act && it > 0){
        memset(act, 0, N);
        for(i = 0; i < N; i++)
          if(chg[i]){
            ii = i/nzp;
            jj = i%nzp;
            for(di = -1; di <= 1; di++)
              for(dj = -1; dj <= 1; dj++)
                if(ii + di >= 0 && ii + di < nxx && jj + dj >= 0 && jj + dj < nzp)
                  act[(ii + di)*nzp + jj + dj] = 1;
          }
      }else
        memset(act, 1, N);
      for(i = 0; i < N; i++){
        if((mask && !mask[i]) || !act[i] || i < i_lo || i >= i_hi){
          Vn[i] = V[i];
          continue;
        }
        Hp = H[i] + hw*f->work[i];
        kd = f->nl ? 0.5*kdiag(f, i % f->nzp) - 0.5*eta : kdiag(f, i % f->nzp);
        if(f->law == 0){
          /* slip at the collocation time: delta + dt V (first order) or delta + dt (V_old + V)/2 */
          /* slip distance: net slip, or (path = 1) the accumulated slip path (TPV10: backward motion
             adds to it); strength c + mu max(sigma, 0) (TPV10: tension acts as zero normal stress) */
          dm = f->path ? f->dpath[i] + (f->nl ? 0.5*dt*fabs(Vold[i]) : 0.0) : fabs(f->delta[i] + (f->nl ? 0.5*dt*Vold[i] : 0.0));
          sg = f->sigma[i] > 0.0 ? f->sigma[i] : 0.0;
          coh = f->coh ? f->coh[i] : 0.0;
          mu = dm < f->d0 ? f->mu_s - (f->mu_s - f->mu_d)*dm/f->d0 : f->mu_d;
          slope = dm < f->d0 ? sg*(f->mu_s - f->mu_d)*hw*dt/f->d0 : 0.0;
          Vn[i] = local_sw(f, i, Hp, kd, coh + sg*mu, slope);
          if(slope > 0.0){
            /* the strength falls linearly with |V| only until the slip distance reaches d0 (at
               |V| = vc); beyond, it stays at the dynamic level. If the step weakens faster than
               the radiation term can balance (-kd <= slope; the linearised solve is then ill posed),
               or (sw_saturate = 1) the solution passes vc, take the saturated branch (exact for the
               piecewise-linear strength; the prototype and the earlier tests use the linearised
               solve, so the second case is optional) */
            vc = (f->d0 - dm)/(hw*dt);
            sd = coh + sg*f->mu_d;
            if(-kd - slope <= 0.0 || (f->sw_saturate && fabs(Vn[i]) > vc)){
              vs = (fabs(Hp) > sd) ? (fabs(Hp) - sd)/(-kd) : 0.0;
              if(fabs(Hp) > coh + sg*mu)
                Vn[i] = sgn(Hp)*(vs > vc ? vs : vc);
            }
          }
        }else
          Vn[i] = local_rs(f, i, Hp, kd, th_p[i], V[i]);
        if(fabs(Vn[i] - V[i]) > ch)
          ch = fabs(Vn[i] - V[i]);
      }
#ifdef FD3D_MPI
      if(t->mpi_size > 1){                         /* every rank needs all rates for the next sweep */
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, Vn, f->ccounts, f->cdispls, MPI_DOUBLE, t->comm);
        MPI_Allreduce(MPI_IN_PLACE, &ch, 1, MPI_DOUBLE, MPI_MAX, t->comm);
      }
#endif
      /* changed: by more than 1e-3 of the convergence tolerance (Newton from a nearby guess changes
         the last bits even for unchanged input, which must not keep a cell active) */
      for(i = 0; i < N; i++){
        chg[i] = fabs(Vn[i] - V[i]) > 1e-16*(fabs(V[i]) > 1.0 ? fabs(V[i]) : 1.0);
        if(fabs(Vn[i]) > vm)
          vm = fabs(Vn[i]);
      }
      memcpy(V, Vn, sizeof(double)*N);
      if(ch <= 1e-13*(vm > 1.0 ? vm : 1.0))
        break;     /* Jacobi converged (same rule as the prototype) */
    }
    itot += it + 1;
    if(f->law == 1 && pass == 0)
      memcpy(Vs, V, sizeof(double)*N);
    else{
      memcpy(f->V, V, sizeof(double)*N);
      if(f->law == 1)
        memcpy(f->thnew ? f->thnew : f->theta, th_p, sizeof(double)*N);
    }
    free(V);
    free(chg);
    free(act);
  }
  free(Vn);
  free(th_p);
  free(Vs);
  free(H);
  return itot;
}

static void step_end(struct fd3d *f)
{
  struct fd3d_tab *t = f->t;
  int N = f->N, i;
  double dt = t->dt, eta = 0.5*t->rho*t->cT, *Vold = f->Vold;
  double *Vm;
  apply_k1(f, f->V, f->work);
  if(f->nl){
    Vm = (double *)calloc((size_t)N, sizeof(double));   /* calloc: silences -Wmaybe-uninitialized */
    for(i = 0; i < N; i++){
      f->tau[i] = f->H[i] + 0.5*f->work[i] + (0.5*kdiag(f, i % f->nzp) - 0.5*eta)*f->V[i];
      Vm[i] = 0.5*(Vold[i] + f->V[i]);                    /* midpoint rate of the step */
    }
    fd3d_push(f, Vm);
    for(i = 0; i < N; i++)
      f->delta[i] += Vm[i]*dt;
    if(f->path)
      for(i = 0; i < N; i++)
        f->dpath[i] += 0.5*(fabs(Vold[i]) + fabs(f->V[i]))*dt;
    free(Vm);
  }else{
    for(i = 0; i < N; i++)
      f->tau[i] = f->H[i] + f->work[i] + kdiag(f, i % f->nzp)*f->V[i];
    fd3d_push(f, f->V);
    for(i = 0; i < N; i++)
      f->delta[i] += f->V[i]*dt;
    if(f->path)
      for(i = 0; i < N; i++)
        f->dpath[i] += fabs(f->V[i])*dt;
  }
  f->time += dt;
}

static int fd3d_step_layer_fft(struct fd3d *f, const double *load);
static int fd3d_step_surf(struct fd3d *f, const double *load);
int fd3d_step(struct fd3d *f, const double *load)
{
  int it;
  double t0, t1;
  if(f->layer && f->layer->surf)
    return fd3d_step_surf(f, load);
  if(f->layer && f->layer->fft)
    return fd3d_step_layer_fft(f, load);
  if(f->layer)
    return fd3d_step_layer(f, load);
  t0 = fd3d_now();
  step_begin(f, load);
  t1 = fd3d_now();
  fd3d_prof3[0] += t1 - t0;
  t0 = t1;
  it = step_solve(f, NULL);
  t1 = fd3d_now();
  fd3d_prof3[1] += t1 - t0;
  t0 = t1;
  step_end(f);
  t1 = fd3d_now();
  fd3d_prof3[2] += t1 - t0;
  fd3d_prof3[3] += 1.0;
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

static int ilog2i(int k)
{
  int i = 0;
  while((1 << i) < k)
    i++;
  return i;
}

static void fd3d_kops(struct fd3d *f, int k)
{
  struct fd3d_tab *t = f->t;
  int ik = ilog2i(k), l, ii, jj, mirror = f->mirror;
  size_t Mhh = (size_t)t->Mx*t->Mh, M = (size_t)t->Mx*t->Mz, m;
  double eta = 0.5*t->rho*t->cT, w1 = (k - 0.5)/k;
  const fd3d_cplx *W;
  double w;
  _Complex double ks;
  size_t q;
  if(f->kop_B[ik])
    return;
  f->kop_B[ik] = (double _Complex *)calloc(Mhh, sizeof(double _Complex));
  f->kop_C[ik] = (double _Complex *)calloc(Mhh, sizeof(double _Complex));
  for(l = 1; l <= k && l <= t->L; l++){
    W = t->W_hat + (size_t)(l - 1)*Mhh;
    w = f->nl ? (k - l + 0.5)/k : 1.0;
    for(m = 0; m < Mhh; m++){
      f->kop_B[ik][m] += w*W[m];
      f->kop_C[ik][m] += (1.0 - w)*W[m];
    }
  }
  for(ii = 0; ii < t->Mx; ii++)
    for(jj = 0; jj < t->Mh; jj++){
      ks = t->Ks_hat[(size_t)ii*t->Mz + jj]*t->dt;
      q = (size_t)ii*t->Mh + jj;
      if(f->nl){
        f->kop_B[ik][q] += 0.5*k*ks + (w1 - 1.0)*eta;
        f->kop_C[ik][q] += 0.5*k*ks + (1.0 - w1)*eta;
      } else
        f->kop_B[ik][q] += k*ks;
    }
  /* self coefficients (offset 0; with the mirror the top row also sees its image at +1 row) */
  memset(t->work, 0, sizeof(double _Complex)*M);
  for(ii = 0; ii < t->Mx; ii++)
    for(jj = 0; jj < t->Mh; jj++)
      t->work[(size_t)ii*t->Mz + jj] = f->kop_B[ik][(size_t)ii*t->Mh + jj];
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
  for(ii = 0; ii < t->Mx; ii++)
    for(jj = 0; jj < t->Mh; jj++)
      t->work[(size_t)ii*t->Mz + jj] *= opq[(size_t)ii*t->Mh + jj];
  from_half(t, t->work);
  fd3d_fft2(t, t->work, +1);
  for(i = 0; i < t->nx; i++)
    for(j = 0; j < f->nzp; j++)
      out[i*f->nzp + j] = creal(t->work[(size_t)i*t->Mz + row0 + j])/M;
}

/* push k fine slices: constant V (first order) or linear from Vold to V (second order) */
static void push_k_lin(struct fd3d *f, const double *Vold, const double *V, int k)
{
  struct fd3d_tab *t = f->t;
  size_t Mhh = (size_t)t->Mx*t->Mh, m;
  int q;
  double _Complex *a = (double _Complex *)malloc(sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  double w;
  fd3d_cplx *slot;
  int ii;
  int jj;
  size_t fq;
  field_hat(f, Vold, a);
  field_hat(f, V, t->work);
  for(q = 1; q <= k; q++){
    w = f->nl ? (q - 0.5)/k : 1.0;
    slot = f->hist + (size_t)(f->hk % t->L)*Mhh;

    for(ii = 0; ii < t->Mx; ii++)
      for(jj = 0; jj < t->Mh; jj++){
        fq = (size_t)ii*t->Mz + jj;
        slot[(size_t)ii*t->Mh + jj] = (fd3d_cplx)((1.0 - w)*a[fq] + w*t->work[fq]);
      }
    f->hk = (f->hk + 1) % t->L;
    if(f->nhist < t->L)
      f->nhist++;
  }
  (void)m;
  free(a);
}


/* ---------------- Newton-Krylov solve of a long step (rate and state), fd_plan.md L1 ----------------
   F(V) = Hc + B_k V - S(V; theta_p) = 0 for all cells, S the rate-and-state strength (odd in V),
   B_k the long-step operator (FFT convolution, self term included). Newton with backtracking on
   |F|_2; each Newton system J d = -F, J = B_k - diag(S'), by restarted GMRES, right-preconditioned
   with P = diag(kd - S') (the local Jacobian). Locked cells stay at zero. Returns the number of
   Newton iterations, or -1 if not converged. */
static void nk_resid(struct fd3d *f, int ik, const double *Hc, const double *th, const double *V, double *F, double *conv)
{
  int i, N = f->N;
  kop_apply(f, f->kop_B[ik], V, conv);
  for(i = 0; i < N; i++){
    if(f->locked && f->locked[i]){
      F[i] = 0.0;
      continue;
    }
    F[i] = Hc[i] + conv[i] - sgn(V[i])*rs_strength(f, i, fabs(V[i]), th[i]);
  }
}
static double nk_norm2(int n, const double *x)
{
  double s = 0.0;
  int i;
  for(i = 0; i < n; i++)
    s += x[i]*x[i];
  return sqrt(s);
}
static int nk_solve(struct fd3d *f, int ik, const double *Hc, const double *th, double *V)
{
  int N = f->N, i, j, it, newton, m = 40, restarts = 4;
  double *F = malloc(sizeof(double)*N), *conv = malloc(sizeof(double)*N), *Dp = malloc(sizeof(double)*N), *Pinv = malloc(sizeof(double)*N);
  double *d = malloc(sizeof(double)*N), *Vt = malloc(sizeof(double)*N), *Ft = malloc(sizeof(double)*N), *w = malloc(sizeof(double)*N);
  double *Vk = malloc(sizeof(double)*(size_t)(m + 1)*N), *Hh = calloc((size_t)(m + 1)*m, sizeof(double)), *cs = malloc(sizeof(double)*m), *sn = malloc(sizeof(double)*m), *g = malloc(sizeof(double)*(m + 1)), *y = malloc(sizeof(double)*m);
  double hmax = 0.0, fn, vm;
  int ok = -1;
  double finf;
  double lam;
  double fn_new;
  double dinf;
  int row;
  double kd;
  int r;
  double b0;
  double tol;
  double beta;
  double *vj;
  double *vn;
  double hn;
  int k2;
  double h;
  double *vk2;
  double a1;
  double a2;
  double a1_2;
  double a2_2;
  double rr;
  int jj;
  int k2_2;
  int nj;
  double s;
  double s_2;
  for(i = 0; i < N; i++)
    if(fabs(Hc[i]) > hmax)
      hmax = fabs(Hc[i]);
  nk_resid(f, ik, Hc, th, V, F, conv);
  fn = nk_norm2(N, F);
  for(newton = 0; newton < 40; newton++){
    finf = 0.0;
    dinf = 0.0;
    for(i = 0; i < N; i++)
      if(fabs(F[i]) > finf)
        finf = fabs(F[i]);
    vm = 1e-9;
    for(i = 0; i < N; i++)
      if(fabs(V[i]) > vm)
        vm = fabs(V[i]);
    if(finf <= 1e-11*(hmax > 1.0 ? hmax : 1.0)){
      ok = newton;
      break;
    }
    /* Jacobian diagonal parts */
    for(i = 0; i < N; i++){
      row = i % f->nzp;
      kd = (f->mirror && row == 0) ? f->kop_self_top[ik] : f->kop_self[ik];
      if(f->locked && f->locked[i]){
        Dp[i] = 0.0;
        Pinv[i] = 0.0;
        continue;
      }
      Dp[i] = rs_dstrength(f, i, fabs(V[i]), th[i]);
      Pinv[i] = 1.0/(kd - Dp[i]);
    }
    /* GMRES for J d = -F, right preconditioned, x = P^{-1} z */
    memset(d, 0, sizeof(double)*N);
    {
      b0 = fn;
      tol = 1e-6*b0;
      for(r = 0; r < restarts; r++){
        /* residual r0 = -F - J d */
        kop_apply(f, f->kop_B[ik], d, conv);
        for(i = 0; i < N; i++)
          w[i] = (f->locked && f->locked[i]) ? 0.0 : -F[i] - (conv[i] - Dp[i]*d[i]);
        beta = nk_norm2(N, w);
        if(beta <= tol)
          break;
        for(i = 0; i < N; i++)
          Vk[i] = w[i]/beta;
        memset(g, 0, sizeof(double)*(m + 1));
        g[0] = beta;
        for(j = 0; j < m; j++){
          vj = Vk + (size_t)j*N;
          vn = Vk + (size_t)(j + 1)*N;

          for(i = 0; i < N; i++)
            Ft[i] = Pinv[i]*vj[i];
          kop_apply(f, f->kop_B[ik], Ft, conv);
          for(i = 0; i < N; i++)
            vn[i] = (f->locked && f->locked[i]) ? 0.0 : conv[i] - Dp[i]*Ft[i];
          for(k2 = 0; k2 <= j; k2++){                 /* modified Gram-Schmidt */
            h = 0.0;
            vk2 = Vk + (size_t)k2*N;
            for(i = 0; i < N; i++)
              h += vn[i]*vk2[i];
            Hh[(size_t)k2*m + j] = h;
            for(i = 0; i < N; i++)
              vn[i] -= h*vk2[i];
          }
          hn = nk_norm2(N, vn);
          Hh[(size_t)(j + 1)*m + j] = hn;
          if(hn > 0.0)
            for(i = 0; i < N; i++)
              vn[i] /= hn;
          for(k2 = 0; k2 < j; k2++){                  /* apply the previous rotations */
            a1 = Hh[(size_t)k2*m + j];
            a2 = Hh[(size_t)(k2 + 1)*m + j];
            Hh[(size_t)k2*m + j] = cs[k2]*a1 + sn[k2]*a2;
            Hh[(size_t)(k2 + 1)*m + j] = -sn[k2]*a1 + cs[k2]*a2;
          }
          {
            a1_2 = Hh[(size_t)j*m + j];
            a2_2 = Hh[(size_t)(j + 1)*m + j];
            rr = sqrt(a1_2*a1_2 + a2_2*a2_2);
            cs[j] = rr > 0 ? a1_2/rr : 1.0;
            sn[j] = rr > 0 ? a2_2/rr : 0.0;
            Hh[(size_t)j*m + j] = rr;
            Hh[(size_t)(j + 1)*m + j] = 0.0;
            g[j + 1] = -sn[j]*g[j];
            g[j] = cs[j]*g[j];
          }
          if(fabs(g[j + 1]) <= tol){
            j++;
            break;
          }
        }
        if(j > m)
          j = m;
        {                                              /* back substitution, d += P^{-1} V y */
          nj = j;
          for(jj = nj - 1; jj >= 0; jj--){
            s = g[jj];
            for(k2_2 = jj + 1; k2_2 < nj; k2_2++)
              s -= Hh[(size_t)jj*m + k2_2]*y[k2_2];
            y[jj] = s/Hh[(size_t)jj*m + jj];
          }
          for(i = 0; i < N; i++){
            s_2 = 0.0;
            for(k2_2 = 0; k2_2 < nj; k2_2++)
              s_2 += Vk[(size_t)k2_2*N + i]*y[k2_2];
            d[i] += Pinv[i]*s_2;
          }
        }
      }
    }
    /* backtracking line search on |F|_2 */
    lam = 1.0;
    for(it = 0; it < 12; it++){
      for(i = 0; i < N; i++)
        Vt[i] = V[i] + lam*d[i];
      nk_resid(f, ik, Hc, th, Vt, Ft, conv);
      fn_new = nk_norm2(N, Ft);
      if(isfinite(fn_new) && fn_new <= (1.0 - 1e-4*lam)*fn)
        break;
      lam *= 0.5;
    }
    if(it >= 12)
      break;                                /* no decrease: give up */
    for(i = 0; i < N; i++){
      if(fabs(lam*d[i]) > dinf)
        dinf = fabs(lam*d[i]);
      V[i] = Vt[i];
    }
    memcpy(F, Ft, sizeof(double)*N);
    fn = fn_new;
    if(dinf <= 1e-12*vm){
      ok = newton + 1;
      break;
    }
  }
  free(F);
  free(conv);
  free(Dp);
  free(Pinv);
  free(d);
  free(Vt);
  free(Ft);
  free(w);
  free(Vk);
  free(Hh);
  free(cs);
  free(sn);
  free(g);
  free(y);
  return ok;
}

int fd3d_step_k(struct fd3d *f, const double *load, int k)
{
  struct fd3d_tab *t = f->t;
  int N = f->N, i, it, pass, itot = 0, ik;
  double dt = t->dt*k, *Vold, *Hc, *Vit, *Vn, *conv, *th_p = NULL, *Vs = NULL;
  int nit;
  double ch;
  double vm;
  int row;
  double kd;
  double Hp;
  double d0s;
  double dm;
  double mu;
  double slope;
  double vo;
  double vn;
  double vo_2;
  double vn_2;
  int io;
  int in_;
  if(k <= 1 && !f->generic)
    return fd3d_step(f, load);
  if(k > 64 || (k & (k - 1)))
    return -1;
  if(t->mpi_size > 1)
    return -1;                  /* long steps: one rank only (to be rewritten, fd_plan.md) */
  ik = ilog2i(k);
  fd3d_kops(f, k);
  Vold = (double *)malloc(sizeof(double)*N);
  Hc = (double *)malloc(sizeof(double)*N);
  Vit = (double *)malloc(sizeof(double)*N);
  Vn = (double *)malloc(sizeof(double)*N);
  conv = (double *)malloc(sizeof(double)*N);
  memcpy(Vold, f->V, sizeof(double)*N);
  known_k(f, k);
  if(load)
    for(i = 0; i < N; i++)
      f->H[i] += load[i];
  memcpy(Hc, f->H, sizeof(double)*N);
  if(f->nl){
    kop_apply(f, f->kop_C[ik], Vold, conv);
    for(i = 0; i < N; i++)
      Hc[i] += conv[i];
  }
  if(f->law == 1){
    th_p = (double *)malloc(sizeof(double)*N);
    Vs = (double *)malloc(sizeof(double)*N);
  }
  for(pass = 0; pass < (f->law == 1 ? 2 : 1); pass++){
    memcpy(Vit, pass == 0 ? Vold : Vs, sizeof(double)*N);
    if(f->law == 1)
      for(i = 0; i < N; i++)
        th_p[i] = theta_adv(f, i, f->theta[i], Vold[i], pass == 0 ? Vold[i] : Vs[i], dt);
    if(f->law == 1 && f->nk){
      /* Newton-Krylov instead of the damped Jacobi iteration (fd_plan.md L1) */
      nit = nk_solve(f, ik, Hc, th_p, Vit);
      it = (nit >= 0) ? nit : 500;
    }else
      for(it = 0; it < 500; it++){
        ch = 0.0;
        vm = 0.0;
        kop_apply(f, f->kop_B[ik], Vit, conv);
        for(i = 0; i < N; i++){
          row = i % f->nzp;
          kd = (f->mirror && row == 0) ? f->kop_self_top[ik] : f->kop_self[ik];
          Hp = Hc[i] + conv[i] - kd*Vit[i];
          if(f->law == 0){
            d0s = f->delta[i] + (f->nl ? 0.5*dt*Vold[i] : 0.0);
            dm = fabs(d0s);
            mu = dm < f->d0 ? f->mu_s - (f->mu_s - f->mu_d)*dm/f->d0 : f->mu_d;
            slope = dm < f->d0 ? f->sigma[i]*(f->mu_s - f->mu_d)*(f->nl ? 0.5 : 1.0)*dt/f->d0 : 0.0;
            Vn[i] = local_sw(f, i, Hp, kd, f->sigma[i]*mu, slope);
          }else
            Vn[i] = local_rs(f, i, Hp, kd, th_p[i], Vit[i]);
          if(fabs(Vn[i] - Vit[i]) > ch)
            ch = fabs(Vn[i] - Vit[i]);
          if(fabs(Vn[i]) > vm)
            vm = fabs(Vn[i]);
        }
        /* damped Jacobi: the long-step operator is not diagonally dominant */
        for(i = 0; i < N; i++)
          Vit[i] += f->kdamp*(Vn[i] - Vit[i]);
        if(!isfinite(ch))
          break;
        if(ch <= 1e-12*(vm > 1e-9 ? vm : 1e-9))
          break;
      }
    itot += it + 1;
    {
      /* a long step must not change the peak slip rate much (the step selector aims at
         10 percent); a jump by more than a factor of two means the step was too long or
         the iteration settled on a spurious solution */
      vo = 0.0;
      vn = 0.0;
      for(i = 0; i < N; i++){
        if(fabs(Vold[i]) > vo)
          vo = fabs(Vold[i]);
        if(fabs(Vit[i]) > vn)
          vn = fabs(Vit[i]);
      }
      if(vn > 2.0*vo + 1e-12)
        it = 500;
    }
    if(getenv("FD3D_DEBUG") && (it >= 500 || !isfinite(Vit[0]))){
      vo_2 = 0.0;
      vn_2 = 0.0;
      io = 0;
      in_ = 0;
      for(i = 0; i < N; i++){
        if(fabs(Vold[i]) > vo_2){
          vo_2 = fabs(Vold[i]);
          io = i;
        }
        if(fabs(Vit[i]) > vn_2){
          vn_2 = fabs(Vit[i]);
          in_ = i;
        }
      }
      fprintf(stderr, "fd3d_step_k: k %d pass %d rejected after %d its: vmax old %.3e (cell %d) new %.3e (cell %d)\n", k, pass, it, vo_2, io, vn_2, in_);
    }
    if(it >= 500 || !isfinite(Vit[0])){
      /* not converged: leave the state untouched, the caller retries with fewer substeps */
      free(Vold);
      free(Hc);
      free(Vit);
      free(Vn);
      free(conv);
      free(th_p);
      free(Vs);
      return -2;
    }
    if(f->law == 1 && pass == 0)
      memcpy(Vs, Vit, sizeof(double)*N);
  }
  memcpy(f->V, Vit, sizeof(double)*N);
  if(f->law == 1)
    memcpy(f->theta, th_p, sizeof(double)*N);
  kop_apply(f, f->kop_B[ik], f->V, conv);
  for(i = 0; i < N; i++)
    f->tau[i] = Hc[i] + conv[i];
  push_k_lin(f, Vold, f->V, k);
  for(i = 0; i < N; i++)
    f->delta[i] += (f->nl ? 0.5*(Vold[i] + f->V[i]) : f->V[i])*dt;
  f->time += dt;
  free(Vold);
  free(Hc);
  free(Vit);
  free(Vn);
  free(conv);
  free(th_p);
  free(Vs);
  return itot;
}


/* ---------------- free-surface layer (stage A, dense first version) ---------------- */

/* frames: fault e1 = x, e2 = z (up), e3 = -y; layer e1 = x, e2 = y, e3 = z */
static void lay_resp(const struct fd3d_tab *t, int src_fault, double sx, double sy, double sz, double a1, double a2,
    double rx, double ry, double rz, int rcv_fault, int L, double *K, double *ks)
{
  double s[6], T1, T2, v, vprev = 0.0, late[2];
  double x1, x2, x3;
  int l, q;
  if(src_fault){
    x1 = rx - sx;
    x2 = rz - sz;
    x3 = -(ry - sy);
  } else{
    x1 = rx - sx;
    x2 = ry - sy;
    x3 = rz - sz;
  }
  T1 = 3.0*(sqrt((rx - sx)*(rx - sx) + (ry - sy)*(ry - sy) + (rz - sz)*(rz - sz)) + 4.0*(a1 + a2))/t->cT + 10.0*t->dt;
  T2 = 2.0*T1;
  for(q = 0; q < 2; q++){
    fd3d_rect_ramp_tensor(x1, x2, x3, q ? T2 : T1, a1, a2, t->cL, t->cT, t->rho, src_fault ? 1 : 3, s);
    /* traction: fault receiver -> strike traction e1F.sigma.e3F; layer receiver -> sigma_zz */
    if(src_fault)
      late[q] = rcv_fault ? s[4] : s[1];       /* fault frame: s31 (= e1 . sigma . e3), s22 = sigma_zz */
    else
      late[q] = rcv_fault ? -s[3] : s[2];               /* layer frame: -s12 = e1 . sigma . (0,-1,0), s33 */
  }
  *ks = (late[1] - late[0])/(T2 - T1);
  for(l = 1; l <= L; l++){
    fd3d_rect_ramp_tensor(x1, x2, x3, l*t->dt, a1, a2, t->cL, t->cT, t->rho, src_fault ? 1 : 3, s);
    if(src_fault)
      v = rcv_fault ? s[4] : s[1];
    else
      v = rcv_fault ? -s[3] : s[2];
    K[l - 1] = v - vprev;
    vprev = v;
  }
}

/* dense LU without pivoting problems for the small layer operator: partial pivoting */
static int lu_factor(int n, double *A, int *piv)
{
  int i, j, k, p;
  double m;
  double tmp;
  double fct;
  for(k = 0; k < n; k++){
    m = 0.0;
    p = k;
    for(i = k; i < n; i++)
      if(fabs(A[(size_t)i*n + k]) > m){
        m = fabs(A[(size_t)i*n + k]);
        p = i;
      }
    if(m == 0.0)
      return -1;
    piv[k] = p;
    if(p != k)
      for(j = 0; j < n; j++){
        tmp = A[(size_t)k*n + j];
        A[(size_t)k*n + j] = A[(size_t)p*n + j];
        A[(size_t)p*n + j] = tmp;
      }
    for(i = k + 1; i < n; i++){
      fct = A[(size_t)i*n + k]/A[(size_t)k*n + k];
      A[(size_t)i*n + k] = fct;
      for(j = k + 1; j < n; j++)
        A[(size_t)i*n + j] -= fct*A[(size_t)k*n + j];
    }
  }
  return 0;
}
static void lu_solve(int n, const double *A, const int *piv, double *b)
{
  int i, j;
  double tmp;
  for(i = 0; i < n; i++){
    if(piv[i] != i){
      tmp = b[i];
      b[i] = b[piv[i]];
      b[piv[i]] = tmp;
    }
  }
  for(i = 0; i < n; i++)
    for(j = 0; j < i; j++)
      b[i] -= A[(size_t)i*n + j]*b[j];
  for(i = n - 1; i >= 0; i--){
    for(j = i + 1; j < n; j++)
      b[i] -= A[(size_t)i*n + j]*b[j];
    b[i] /= A[(size_t)i*n + i];
  }
}

struct fd3d_layer *fd3d_layer_create(struct fd3d *f, double x0f, int nxs, int nys, double hx, double hy,
    double xs0, double ys0, int L)
{
  const struct fd3d_tab *t = f->t;
  struct fd3d_layer *ly;
  int N = f->N, ns = nxs*nys, a, b, l, nz = f->nzp;
  double dx = t->dx, *K = (double *)malloc(sizeof(double)*L), ks, *K2 = (double *)malloc(sizeof(double)*L), ks2;
  double sx;
  double sy;
  double rx;
  double ry;
  double rx_2;
  double rz;
  double sx_2;
  double sz;
  double rx_3;
  double ry_2;
  if(!f->mirror)
    return NULL;                       /* the layer completes the mirror image */
  ly = (struct fd3d_layer *)calloc(1, sizeof(*ly));
  ly->ns = ns;
  ly->nxs = nxs;
  ly->nys = nys;
  ly->L = L;
  ly->hx = hx;
  ly->hy = hy;
  ly->xs0 = xs0;
  ly->ys0 = ys0;
  ly->x0f = x0f;
  ly->eta_p = 0.5*t->rho*t->cL;
  ly->Wss = calloc((size_t)L*ns*ns, sizeof(double));
  ly->Kss = calloc((size_t)ns*ns, sizeof(double));
  ly->Wfs = calloc((size_t)L*ns*N, sizeof(double));
  ly->Kfs = calloc((size_t)ns*N, sizeof(double));
  ly->Wsf = calloc((size_t)L*N*ns, sizeof(double));
  ly->Ksf = calloc((size_t)N*ns, sizeof(double));
  for(a = 0; a < ns; a++){                          /* sources: layer elements */
    sx = xs0 + (a/nys + 0.5)*hx;
    sy = ys0 + (a%nys + 0.5)*hy;
    for(b = 0; b < ns; b++){
      rx = xs0 + (b/nys + 0.5)*hx;
      ry = ys0 + (b%nys + 0.5)*hy;
      lay_resp(t, 0, sx, sy, 0.0, 0.5*hx, 0.5*hy, rx, ry, 0.0, 0, L, K, &ks);
      ly->Kss[(size_t)b*ns + a] = ks;
      for(l = 0; l < L; l++)
        ly->Wss[((size_t)l*ns + b)*ns + a] = K[l] - ks*t->dt;
    }
    for(b = 0; b < N; b++){
      rx_2 = x0f + (b/nz + 0.5)*dx;
      rz = -(b%nz + 0.5)*dx;
      lay_resp(t, 0, sx, sy, 0.0, 0.5*hx, 0.5*hy, rx_2, 0.0, rz, 1, L, K, &ks);
      ly->Ksf[(size_t)b*ns + a] = ks;
      for(l = 0; l < L; l++)
        ly->Wsf[((size_t)l*N + b)*ns + a] = K[l] - ks*t->dt;
    }
  }
  for(a = 0; a < N; a++){                           /* sources: fault elements plus their images */
    sx_2 = x0f + (a/nz + 0.5)*dx;
    sz = -(a%nz + 0.5)*dx;
    for(b = 0; b < ns; b++){
      rx_3 = xs0 + (b/nys + 0.5)*hx;
      ry_2 = ys0 + (b%nys + 0.5)*hy;
      lay_resp(t, 1, sx_2, 0.0, sz, 0.5*dx, 0.5*dx, rx_3, ry_2, 0.0, 0, L, K, &ks);
      lay_resp(t, 1, sx_2, 0.0, -sz, 0.5*dx, 0.5*dx, rx_3, ry_2, 0.0, 0, L, K2, &ks2);
      ly->Kfs[(size_t)b*N + a] = ks + ks2;
      for(l = 0; l < L; l++)
        ly->Wfs[((size_t)l*ns + b)*N + a] = K[l] + K2[l] - (ks + ks2)*t->dt;
    }
  }
  free(K);
  free(K2);
  ly->K1ss = malloc(sizeof(double)*(size_t)ns*ns);
  ly->K1fs = malloc(sizeof(double)*(size_t)ns*N);
  ly->K1sf = malloc(sizeof(double)*(size_t)N*ns);
  for(a = 0; a < ns*ns; a++)
    ly->K1ss[a] = ly->Wss[a] + ly->Kss[a]*t->dt;
  for(a = 0; a < ns*N; a++){
    ly->K1fs[a] = ly->Wfs[a] + ly->Kfs[a]*t->dt;
    ly->K1sf[a] = ly->Wsf[a] + ly->Ksf[a]*t->dt;
  }
  /* the layer's self operator for the current step: first order K1ss, NL (K1ss + E)/2 - E */
  ly->Bss_lu = malloc(sizeof(double)*(size_t)ns*ns);
  ly->piv = malloc(sizeof(int)*ns);
  for(a = 0; a < ns*ns; a++)
    ly->Bss_lu[a] = f->nl ? 0.5*ly->K1ss[a] : ly->K1ss[a];
  if(f->nl)
    for(a = 0; a < ns; a++)
      ly->Bss_lu[(size_t)a*ns + a] -= 0.5*ly->eta_p;
  if(lu_factor(ns, ly->Bss_lu, ly->piv) != 0){
    fd3d_layer_free(ly);
    return NULL;
  }
  ly->V = calloc(ns, sizeof(double));
  ly->Vold = calloc(ns, sizeof(double));
  ly->delta = calloc(ns, sizeof(double));
  ly->H = calloc(ns, sizeof(double));
  ly->tau = calloc(ns, sizeof(double));
  ly->hs = calloc((size_t)L*ns, sizeof(double));
  ly->hf = calloc((size_t)L*N, sizeof(double));
  return ly;
}

static void surf_free(struct fd3d_layer *ly);
void fd3d_layer_free(struct fd3d_layer *ly)
{
  if(!ly)
    return;
  if(ly->surf){
    surf_free(ly);
    return;
  }
  if(ly->fft){
    free(ly->Wss2);
    free(ly->Kss2);
    free(ly->K1ss2);
    free(ly->Wfs_h);
    free(ly->Kfs_h);
    free(ly->K1fs_h);
    free(ly->Wsf_h);
    free(ly->Ksf_h);
    free(ly->K1sf_h);
    free(ly->hs2);
    free(ly->hsx);
    free(ly->hfx);
    fd_fft_plan_destroy(ly->pf);
    fd_fft_plan_destroy(ly->pb);
    fd_fft_plan_destroy(ly->pf2x);
    fd_fft_plan_destroy(ly->pb2x);
    fd_fft_plan_destroy(ly->pf2y);
    fd_fft_plan_destroy(ly->pb2y);
  }
  free(ly->Wss);
  free(ly->Kss);
  free(ly->Wfs);
  free(ly->Kfs);
  free(ly->Wsf);
  free(ly->Ksf);
  free(ly->K1ss);
  free(ly->K1fs);
  free(ly->K1sf);
  free(ly->Bss_lu);
  free(ly->piv);
  free(ly->V);
  free(ly->Vold);
  free(ly->delta);
  free(ly->H);
  free(ly->tau);
  free(ly->hs);
  free(ly->hf);
  free(ly);
}

static void mv(int nr, int nc, const double *A, const double *x, double w, double *y)   /* y += w A x */
{
  int i, j;
  double s;
  for(i = 0; i < nr; i++){
    s = 0.0;
    for(j = 0; j < nc; j++)
      s += A[(size_t)i*nc + j]*x[j];
    y[i] += w*s;
  }
}

int fd3d_step_layer(struct fd3d *f, const double *load)
{
  struct fd3d_layer *ly = f->layer;
  const struct fd3d_tab *t = f->t;
  int N = f->N, ns = ly->ns, i, l, it, itot = 0, lmax = ly->L < ly->nhist + 1 ? ly->L : ly->nhist + 1;
  double dt = t->dt, hw = f->nl ? 0.5 : 1.0;
  double *Hx = calloc(N, sizeof(double)), *Hs = malloc(sizeof(double)*ns), *Vs = malloc(sizeof(double)*ns), *Vp = malloc(sizeof(double)*N);
  double *thbuf = (f->law == 1) ? malloc(sizeof(double)*N) : NULL;
  int slot;
  double ch;
  double vm;
  double *Hf;
  int slot_2;
  double vm_;
  step_begin(f, load);
  memcpy(ly->Vold, ly->V, sizeof(double)*ns);
  /* known parts: histories (lags >= 2), statics, the V_old half of the current step (NL) */
  for(i = 0; i < ns; i++)
    ly->H[i] = 0.0;
  mv(ns, ns, ly->Kss, ly->delta, 1.0, ly->H);
  mv(ns, N, ly->Kfs, f->delta, 1.0, ly->H);
  if(!f->static_fn)
    mv(N, ns, ly->Ksf, ly->delta, 1.0, Hx);   /* with the exact Is the layer's static part is in Is already */
  for(l = 2; l <= lmax; l++){
    slot = (ly->hk - (l - 1) + 2*ly->L) % ly->L;
    mv(ns, ns, ly->Wss + (size_t)(l - 1)*ns*ns, ly->hs + (size_t)slot*ns, 1.0, ly->H);
    mv(ns, N, ly->Wfs + (size_t)(l - 1)*ns*N, ly->hf + (size_t)slot*N, 1.0, ly->H);
    mv(N, ns, ly->Wsf + (size_t)(l - 1)*N*ns, ly->hs + (size_t)slot*ns, 1.0, Hx);
  }
  if(f->nl){
    mv(ns, ns, ly->K1ss, ly->Vold, 0.5, ly->H);
    for(i = 0; i < ns; i++)
      ly->H[i] += 0.5*ly->eta_p*ly->Vold[i];
    mv(ns, N, ly->K1fs, f->Vold, 0.5, ly->H);
    mv(N, ns, ly->K1sf, ly->Vold, 0.5, Hx);
  }
  /* iterate fault and layer: the fault sees the layer through the first lag, and vice versa */
  memcpy(Vs, ly->Vold, sizeof(double)*ns);
  f->thnew = thbuf;
  for(it = 0; it < 200; it++){
    ch = 0.0;
    vm = 1e-30;
    Hf = malloc(sizeof(double)*N);
    memcpy(Hf, Hx, sizeof(double)*N);
    mv(N, ns, ly->K1sf, Vs, hw, Hf);
    memcpy(Vp, f->V, sizeof(double)*N);
    itot += step_solve(f, Hf);
    free(Hf);
    for(i = 0; i < ns; i++)
      Hs[i] = -ly->H[i];
    mv(ns, N, ly->K1fs, f->V, -hw, Hs);
    lu_solve(ns, ly->Bss_lu, ly->piv, Hs);
    for(i = 0; i < ns; i++){
      if(fabs(Hs[i] - Vs[i]) > ch)
        ch = fabs(Hs[i] - Vs[i]);
      if(fabs(Hs[i]) > vm)
        vm = fabs(Hs[i]);
    }
    for(i = 0; i < N; i++){
      if(fabs(f->V[i] - Vp[i]) > ch)
        ch = fabs(f->V[i] - Vp[i]);
      if(fabs(f->V[i]) > vm)
        vm = fabs(f->V[i]);
    }
    memcpy(Vs, Hs, sizeof(double)*ns);
    if(it > 0 && ch <= 1e-12*vm)
      break;
  }
  f->thnew = NULL;
  if(thbuf){
    memcpy(f->theta, thbuf, sizeof(double)*N);
    free(thbuf);
  }
  /* finish: the fault with the layer's contribution in its known part, then the layer */
  memcpy(ly->V, Vs, sizeof(double)*ns);
  mv(N, ns, ly->K1sf, ly->V, hw, Hx);
  for(i = 0; i < N; i++)
    f->H[i] += Hx[i];
  {
    slot_2 = ly->hk % ly->L;
    for(i = 0; i < ns; i++){
      vm_ = f->nl ? 0.5*(ly->Vold[i] + ly->V[i]) : ly->V[i];
      ly->hs[(size_t)slot_2*ns + i] = vm_;
      ly->delta[i] += vm_*dt;
    }
    for(i = 0; i < N; i++)
      ly->hf[(size_t)slot_2*N + i] = f->nl ? 0.5*(f->Vold[i] + f->V[i]) : f->V[i];
    ly->hk = (ly->hk + 1) % ly->L;
    if(ly->nhist < ly->L)
      ly->nhist++;
  }
  step_end(f);
  /* layer normal traction at the end of the step (should be zero) */
  for(i = 0; i < ns; i++)
    ly->tau[i] = 0.0;
  free(Hx);
  free(Hs);
  free(Vs);
  free(Vp);
  return itot;
}


/* ---------------- free-surface layer, FFT-structured ---------------- */
#include <time.h>
static double ly_now(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9*ts.tv_nsec;
}
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
    for(x = 0; x < nin; x++)
      a[x*st] = in_rowmajor_x_fast ? in[(size_t)r*nin + x] : in[(size_t)x*nrow + r];
    fd_fft_exec(ly->pf, a, b);
    for(k = 0; k < ly->Mh; k++)
      out[(size_t)k*nrow + r] = b[k];
  }
  free(a);
  free(b);
}
/* inverse: in[k][row] (half spectrum) -> out real values at x = st*i, i < nout, stored [i][row] (element order) */
static void lx_inv(const struct fd3d_layer *ly, int nrow, const double _Complex *in, int nout, int st, double w, double *out)
{
  double _Complex *a = malloc(sizeof(double _Complex)*ly->Mx), *b = malloc(sizeof(double _Complex)*ly->Mx);
  int r, i, k, M = ly->Mx;
  for(r = 0; r < nrow; r++){
    for(k = 0; k < ly->Mh; k++)
      a[k] = in[(size_t)k*nrow + r];
    for(k = ly->Mh; k < M; k++)
      a[k] = conj(a[M - k]);
    fd_fft_exec(ly->pb, a, b);
    for(i = 0; i < nout; i++)
      out[(size_t)i*nrow + r] += w*creal(b[(size_t)i*st])/M;
  }
  free(a);
  free(b);
}
/* Y[k][nr] += T[k][nr][nc] X[k][nc] */
static void lx_apply(const struct fd3d_layer *ly, int nr, int nc, const double _Complex *T, const double _Complex *X, double _Complex *Y)
{
  int k, r, c;
  const _Complex double *Tk;
  const _Complex double *Xk;
  _Complex double *Yk;
  _Complex double s;
  for(k = 0; k < ly->Mh; k++){
    Tk = T + (size_t)k*nr*nc;
    Xk = X + (size_t)k*nc;
    Yk = Y + (size_t)k*nr;
    for(r = 0; r < nr; r++){
      s = 0.0;
      for(c = 0; c < nc; c++)
        s += Tk[(size_t)r*nc + c]*Xk[c];
      Yk[r] += s;
    }
  }
}
/* 2-D layer FFT: x[p][q] (element order) -> half spectrum (Mx2 x Mh2) */
static void l2_fwd(const struct fd3d_layer *ly, const double *x, double _Complex *out)
{
  int Mx = ly->Mx2, My = ly->My2, i, j;
  double _Complex *a = calloc((size_t)Mx*My, sizeof(double _Complex)), *col = malloc(sizeof(double _Complex)*(Mx > My ? Mx : My)),
                  *o = malloc(sizeof(double _Complex)*(Mx > My ? Mx : My));
  for(i = 0; i < ly->nxs; i++)
    for(j = 0; j < ly->nys; j++)
      a[(size_t)i*My + j] = x[(size_t)i*ly->nys + j];
  for(i = 0; i < Mx; i++){
    fd_fft_exec(ly->pf2y, a + (size_t)i*My, o);
    memcpy(a + (size_t)i*My, o, sizeof(double _Complex)*My);
  }
  for(j = 0; j < ly->Mh2; j++){
    for(i = 0; i < Mx; i++)
      col[i] = a[(size_t)i*My + j];
    fd_fft_exec(ly->pf2x, col, o);
    for(i = 0; i < Mx; i++)
      out[(size_t)i*ly->Mh2 + j] = o[i];
  }
  free(a);
  free(col);
  free(o);
}
/* inverse of a half spectrum product into element order (adds w times) */
static void l2_inv(const struct fd3d_layer *ly, const double _Complex *in, double w, double *x)
{
  int Mx = ly->Mx2, My = ly->My2, i, j;
  double _Complex *a = calloc((size_t)Mx*My, sizeof(double _Complex)), *col = malloc(sizeof(double _Complex)*(Mx > My ? Mx : My)),
                  *o = malloc(sizeof(double _Complex)*(Mx > My ? Mx : My));
  for(i = 0; i < Mx; i++)
    for(j = 0; j < ly->Mh2; j++)
      a[(size_t)i*My + j] = in[(size_t)i*ly->Mh2 + j];
  for(i = 0; i < Mx; i++)
    for(j = ly->Mh2; j < My; j++)
      a[(size_t)i*My + j] = conj(a[(size_t)((Mx - i) % Mx)*My + (My - j)]);
  for(j = 0; j < My; j++){
    for(i = 0; i < Mx; i++)
      col[i] = a[(size_t)i*My + j];
    fd_fft_exec(ly->pb2x, col, o);
    for(i = 0; i < Mx; i++)
      a[(size_t)i*My + j] = o[i];
  }
  for(i = 0; i < ly->nxs; i++){
    fd_fft_exec(ly->pb2y, a + (size_t)i*My, o);
    for(j = 0; j < ly->nys; j++)
      x[(size_t)i*ly->nys + j] += w*creal(o[j])/((double)Mx*My);
  }
  free(a);
  free(col);
  free(o);
}

struct fd3d_layer *fd3d_layer_create_fft(struct fd3d *f, double x0f, int m, int nxs, int nys, double xs0, int L)
{
  const struct fd3d_tab *t = f->t;
  struct fd3d_layer *ly;
  int nx = t->nx, nz = f->nzp, d, dmin, dmax, q, j, l, dp, dq, k;
  double dx = t->dx, hx = dx/m, *K = malloc(sizeof(double)*L), *K2 = malloc(sizeof(double)*L), ks, ks2;
  double _Complex *buf;
  size_t Mh, s1;
  size_t M2;
  size_t H2;
  double *A;
  double *tmp;
  size_t pos;
  _Complex double *a;
  _Complex double *col;
  _Complex double *o;
  _Complex double *dst;
  int i2;
  int j2;
  size_t nfs;
  double *Afs;
  double *Asf;
  size_t pos_2;
  double yq;
  double zj;
  _Complex double *dfs;
  _Complex double *dsf;
  _Complex double *a_2;
  int x;
  if(!f->mirror || m < 1)
    return NULL;
  ly = calloc(1, sizeof(*ly));
  ly->fft = 1;
  ly->m = m;
  ly->nxs = nxs;
  ly->nys = nys;
  ly->ns = nxs*nys;
  ly->L = L;
  ly->hx = hx;
  ly->hy = hx;
  ly->xs0 = xs0;
  ly->ys0 = -0.5*nys*hx;
  ly->x0f = x0f;
  ly->nx = nx;
  ly->nz = nz;
  ly->eta_p = 0.5*t->rho*t->cL;
  ly->c0 = xs0 + 0.5*hx - x0f - 0.5*dx;
  ly->Mx = fd_fft_size_for(nxs + m*(nx - 1) + 1);
  ly->Mh = ly->Mx/2 + 1;
  Mh = ly->Mh;
  ly->Mx2 = fd_fft_size_for(2*nxs - 1);
  ly->My2 = fd_fft_size_for(2*nys - 1);
  ly->Mh2 = ly->My2/2 + 1;
  ly->pf = fd_fft_plan_create(ly->Mx, FD_FFT_FORWARD);
  ly->pb = fd_fft_plan_create(ly->Mx, FD_FFT_BACKWARD);
  ly->pf2x = fd_fft_plan_create(ly->Mx2, FD_FFT_FORWARD);
  ly->pb2x = fd_fft_plan_create(ly->Mx2, FD_FFT_BACKWARD);
  ly->pf2y = fd_fft_plan_create(ly->My2, FD_FFT_FORWARD);
  ly->pb2y = fd_fft_plan_create(ly->My2, FD_FFT_BACKWARD);
  /* layer-layer: offsets (dp, dq), kernel arrays in FFT order, then 2-D half spectra */
  {
    M2 = (size_t)ly->Mx2*ly->My2;
    H2 = (size_t)ly->Mx2*ly->Mh2;
    A = calloc((size_t)(L + 2)*M2, sizeof(double));       /* [l][M2]: K lags 1..L, then Ks, then K1 */
    tmp = malloc(sizeof(double)*ly->ns);
    ly->Wss2 = malloc(sizeof(double _Complex)*(size_t)L*H2);
    ly->Kss2 = malloc(sizeof(double _Complex)*H2);
    ly->K1ss2 = malloc(sizeof(double _Complex)*H2);
    for(dp = -(nxs - 1); dp <= nxs - 1; dp++)
      for(dq = -(nys - 1); dq <= nys - 1; dq++){
        pos = (size_t)((dp + ly->Mx2) % ly->Mx2)*ly->My2 + (dq + ly->My2) % ly->My2;
        lay_resp(t, 0, 0.0, 0.0, 0.0, 0.5*hx, 0.5*hx, dp*hx, dq*hx, 0.0, 0, L, K, &ks);
        for(l = 0; l < L; l++)
          A[(size_t)l*M2 + pos] = K[l] - ks*t->dt;
        A[(size_t)L*M2 + pos] = ks;
        A[(size_t)(L + 1)*M2 + pos] = K[0];
        if(dp == 0 && dq == 0)
          ly->self_ss = K[0];
      }
    /* transform each [M2] array: reuse l2 machinery by building element-order arrays is awkward for
       offsets, so transform directly */
    for(l = 0; l < L + 2; l++){
      a = malloc(sizeof(double _Complex)*M2);
      col = malloc(sizeof(double _Complex)*(ly->Mx2 > ly->My2 ? ly->Mx2 : ly->My2));
      o = malloc(sizeof(double _Complex)*(ly->Mx2 > ly->My2 ? ly->Mx2 : ly->My2));
      dst = (l < L) ? ly->Wss2 + (size_t)l*H2 : (l == L ? ly->Kss2 : ly->K1ss2);

      for(s1 = 0; s1 < M2; s1++)
        a[s1] = A[(size_t)l*M2 + s1];
      for(i2 = 0; i2 < ly->Mx2; i2++){
        fd_fft_exec(ly->pf2y, a + (size_t)i2*ly->My2, o);
        memcpy(a + (size_t)i2*ly->My2, o, sizeof(double _Complex)*ly->My2);
      }
      for(j2 = 0; j2 < ly->Mh2; j2++){
        for(i2 = 0; i2 < ly->Mx2; i2++)
          col[i2] = a[(size_t)i2*ly->My2 + j2];
        fd_fft_exec(ly->pf2x, col, o);
        for(i2 = 0; i2 < ly->Mx2; i2++)
          dst[(size_t)i2*ly->Mh2 + j2] = o[i2];
      }
      free(a);
      free(col);
      free(o);
    }
    free(A);
    free(tmp);
  }
  /* fault <-> layer: Toeplitz along strike, offsets d on the lattice */
  dmin = -(nxs - 1) < -m*(nx - 1) ? -(nxs - 1) : -m*(nx - 1);
  dmax = (nxs - 1) > m*(nx - 1) ? (nxs - 1) : m*(nx - 1);
  {
    nfs = (size_t)nys*nz;
    Afs = calloc((size_t)(L + 2)*ly->Mx*nfs, sizeof(double));   /* [l][x][q][j] */
    Asf = calloc((size_t)(L + 2)*ly->Mx*nfs, sizeof(double));   /* [l][x][j][q] */
    buf = malloc(sizeof(double _Complex)*ly->Mx);
    for(d = dmin; d <= dmax; d++){
      pos_2 = (size_t)((d + ly->Mx) % ly->Mx);
      for(q = 0; q < nys; q++){
        yq = ly->ys0 + (q + 0.5)*hx;
        for(j = 0; j < nz; j++){
          zj = -(j + 0.5)*dx;
          /* fault (and image) -> layer: receiver at x offset d hx + c0, d = p - m i in [-m(nx-1), nxs-1];
             layer -> fault: d = m i - p in [-(nxs-1), m(nx-1)]; each range fits the FFT length,
             the union may not (filling both from one loop would alias) */
          if(d >= -m*(nx - 1) && d <= nxs - 1){
            lay_resp(t, 1, 0.0, 0.0, zj, 0.5*dx, 0.5*dx, d*hx + ly->c0, yq, 0.0, 0, L, K, &ks);
            lay_resp(t, 1, 0.0, 0.0, -zj, 0.5*dx, 0.5*dx, d*hx + ly->c0, yq, 0.0, 0, L, K2, &ks2);
            for(l = 0; l < L; l++)
              Afs[(((size_t)l*ly->Mx + pos_2)*nys + q)*nz + j] = K[l] + K2[l] - (ks + ks2)*t->dt;
            Afs[(((size_t)L*ly->Mx + pos_2)*nys + q)*nz + j] = ks + ks2;
            Afs[(((size_t)(L + 1)*ly->Mx + pos_2)*nys + q)*nz + j] = K[0] + K2[0];
          }
          if(d >= -(nxs - 1) && d <= m*(nx - 1)){
            /* layer -> fault: receiver at x offset d hx - c0 */
            lay_resp(t, 0, 0.0, yq, 0.0, 0.5*hx, 0.5*hx, d*hx - ly->c0, 0.0, zj, 1, L, K, &ks);
            for(l = 0; l < L; l++)
              Asf[(((size_t)l*ly->Mx + pos_2)*nz + j)*nys + q] = K[l] - ks*t->dt;
            Asf[(((size_t)L*ly->Mx + pos_2)*nz + j)*nys + q] = ks;
            Asf[(((size_t)(L + 1)*ly->Mx + pos_2)*nz + j)*nys + q] = K[0];
          }
        }
      }
    }
    ly->Wfs_h = malloc(sizeof(double _Complex)*(size_t)L*Mh*nfs);
    ly->Kfs_h = malloc(sizeof(double _Complex)*Mh*nfs);
    ly->K1fs_h = malloc(sizeof(double _Complex)*Mh*nfs);
    ly->Wsf_h = malloc(sizeof(double _Complex)*(size_t)L*Mh*nfs);
    ly->Ksf_h = malloc(sizeof(double _Complex)*Mh*nfs);
    ly->K1sf_h = malloc(sizeof(double _Complex)*Mh*nfs);
    for(l = 0; l < L + 2; l++)
      for(s1 = 0; s1 < nfs; s1++){
        dfs = (l < L) ? ly->Wfs_h + (size_t)l*Mh*nfs : (l == L ? ly->Kfs_h : ly->K1fs_h);
        dsf = (l < L) ? ly->Wsf_h + (size_t)l*Mh*nfs : (l == L ? ly->Ksf_h : ly->K1sf_h);
        a_2 = calloc(ly->Mx, sizeof(double _Complex));

        for(x = 0; x < ly->Mx; x++)
          a_2[x] = Afs[((size_t)l*ly->Mx + x)*nfs + s1];
        fd_fft_exec(ly->pf, a_2, buf);
        for(k = 0; k < (int)Mh; k++)
          dfs[(size_t)k*nfs + s1] = buf[k];
        for(x = 0; x < ly->Mx; x++)
          a_2[x] = Asf[((size_t)l*ly->Mx + x)*nfs + s1];
        fd_fft_exec(ly->pf, a_2, buf);
        for(k = 0; k < (int)Mh; k++)
          dsf[(size_t)k*nfs + s1] = buf[k];
        free(a_2);
      }
    free(Afs);
    free(Asf);
    free(buf);
  }
  free(K);
  free(K2);
  ly->V = calloc(ly->ns, sizeof(double));
  ly->Vold = calloc(ly->ns, sizeof(double));
  ly->delta = calloc(ly->ns, sizeof(double));
  ly->H = calloc(ly->ns, sizeof(double));
  ly->tau = calloc(ly->ns, sizeof(double));
  ly->hs2 = calloc((size_t)L*ly->Mx2*ly->Mh2, sizeof(double _Complex));
  ly->hsx = calloc((size_t)L*Mh*nys, sizeof(double _Complex));
  ly->hfx = calloc((size_t)L*Mh*nz, sizeof(double _Complex));
  return ly;
}

/* element-order fault field (i*nz + j) -> rows [j][x] at stride m */
static void fault_rows_hat(const struct fd3d_layer *ly, const double *V, double _Complex *out)
{
  double *r = malloc(sizeof(double)*(size_t)ly->nz*ly->nx);
  int i, j;
  for(i = 0; i < ly->nx; i++)
    for(j = 0; j < ly->nz; j++)
      r[(size_t)j*ly->nx + i] = V[(size_t)i*ly->nz + j];
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
  int slot;
  const _Complex double *W2;
  const _Complex double *h2;
  _Complex double *Yt;
  _Complex double *Yt_2;
  double *Hrow;
  double *Hrow_2;
  double ch;
  double vm;
  double *Hf;
  double chi;
  double vmi;
  double dv;
  double dmax;
  double smax;
  int slot_2;
  double *vms;
  double *vmf;
  step_begin(f, load);
  memcpy(ly->Vold, ly->V, sizeof(double)*ns);
  /* known parts in the spectral domain */
  memset(ly->H, 0, sizeof(double)*ns);
  l2_fwd(ly, ly->delta, S2);
  for(m = 0; m < H2; m++)
    P2[m] = ly->Kss2[m]*S2[m];
  fault_rows_hat(ly, f->delta, Xf);
  lx_apply(ly, nys, nz, ly->Kfs_h, Xf, Ys);
  if(!f->static_fn){
    layer_rows_hat(ly, ly->delta, Xs);
    lx_apply(ly, nz, nys, ly->Ksf_h, Xs, Yf);
  }
  for(l = 2; l <= lmax; l++){
    slot = (ly->hk - (l - 1) + 2*ly->L) % ly->L;
    W2 = ly->Wss2 + (size_t)(l - 1)*H2;
    h2 = ly->hs2 + (size_t)slot*H2;
    for(m = 0; m < H2; m++)
      P2[m] += W2[m]*h2[m];
    lx_apply(ly, nys, nz, ly->Wfs_h + (size_t)(l - 1)*Mh*nys*nz, ly->hfx + (size_t)slot*Mh*nz, Ys);
    lx_apply(ly, nz, nys, ly->Wsf_h + (size_t)(l - 1)*Mh*nys*nz, ly->hsx + (size_t)slot*Mh*nys, Yf);
  }
  if(f->nl){
    l2_fwd(ly, ly->Vold, S2);
    for(m = 0; m < H2; m++)
      P2[m] += 0.5*ly->K1ss2[m]*S2[m];
    for(i = 0; i < ns; i++)
      ly->H[i] += 0.5*ly->eta_p*ly->Vold[i];
    fault_rows_hat(ly, f->Vold, Xf);
    {
      Yt = calloc(Mh*nys, sizeof(double _Complex));
      lx_apply(ly, nys, nz, ly->K1fs_h, Xf, Yt);
      for(m = 0; m < Mh*nys; m++)
        Ys[m] += 0.5*Yt[m];
      free(Yt);
    }
    layer_rows_hat(ly, ly->Vold, Xs);
    {
      Yt_2 = calloc(Mh*nz, sizeof(double _Complex));
      lx_apply(ly, nz, nys, ly->K1sf_h, Xs, Yt_2);
      for(m = 0; m < Mh*nz; m++)
        Yf[m] += 0.5*Yt_2[m];
      free(Yt_2);
    }
  }
  l2_inv(ly, P2, 1.0, ly->H);
  {
    Hrow = calloc((size_t)ns, sizeof(double));
    lx_inv(ly, nys, Ys, ly->nxs, 1, 1.0, Hrow);            /* [p][q] = element order */
    for(i = 0; i < ns; i++)
      ly->H[i] += Hrow[i];
    free(Hrow);
  }
  {
    Hrow_2 = calloc((size_t)N, sizeof(double));       /* [i][j] = element order i*nz + j */
    lx_inv(ly, nz, Yf, ly->nx, ly->m, 1.0, Hrow_2);
    for(i = 0; i < N; i++)
      Hx[i] += Hrow_2[i];
    free(Hrow_2);
  }
  tp1 = ly_now();
  fd3d_layer_prof[0] += tp1 - tp0;
  tp0 = tp1;
  /* iterate fault and layer */
  memcpy(Vs, ly->Vold, sizeof(double)*ns);
  f->thnew = thbuf;
  double *Hf_prev = malloc(sizeof(double)*N);
  char *mask = malloc(N);
  for(it = 0; it < 200; it++){
    ch = 0.0;
    vm = 1e-30;
    Hf = malloc(sizeof(double)*N);
    /* fault sees the layer's current step through the first lag */
    memset(Yf, 0, sizeof(double _Complex)*Mh*nz);
    layer_rows_hat(ly, Vs, Xs);
    lx_apply(ly, nz, nys, ly->K1sf_h, Xs, Yf);
    memcpy(Hf, Hx, sizeof(double)*N);
    memset(tmpf, 0, sizeof(double)*N);
    lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, tmpf);
    for(i = 0; i < N; i++)
      Hf[i] += tmpf[i];
    memcpy(Vp, f->V, sizeof(double)*N);
    tp1 = ly_now();
    fd3d_layer_prof[2] += tp1 - tp0;
    tp0 = tp1;
    itot += step_solve_g(f, Hf, it > 0 ? Vp : NULL, NULL);   /* warm start from the previous outer iterate */
    memcpy(Hf_prev, Hf, sizeof(double)*N);
    tp1 = ly_now();
    fd3d_layer_prof[1] += tp1 - tp0;
    tp0 = tp1;
    free(Hf);
    /* layer: (hw K1ss - nl eta/2) Vs = -(H + hw K1fs Vf), Jacobi with the self coefficient */
    memset(Ys, 0, sizeof(double _Complex)*Mh*nys);
    fault_rows_hat(ly, f->V, Xf);
    lx_apply(ly, nys, nz, ly->K1fs_h, Xf, Ys);
    for(i = 0; i < ns; i++)
      rhs[i] = -ly->H[i];
    lx_inv(ly, nys, Ys, ly->nxs, 1, -hw, rhs);
    for(inner = 0; inner < 500; inner++){
      chi = 0.0;
      vmi = 1e-30;
      l2_fwd(ly, Vs, S2);
      for(m = 0; m < H2; m++)
        P2[m] = ly->K1ss2[m]*S2[m];
      memset(tmps, 0, sizeof(double)*ns);
      l2_inv(ly, P2, hw, tmps);
      if(f->nl)
        for(i = 0; i < ns; i++)
          tmps[i] -= 0.5*ly->eta_p*Vs[i];
      for(i = 0; i < ns; i++){
        dv = (rhs[i] - tmps[i])/diag;
        Vs[i] += dv;
        if(fabs(dv) > chi)
          chi = fabs(dv);
        if(fabs(Vs[i]) > vmi)
          vmi = fabs(Vs[i]);
      }
      if(chi <= 1e-13*vmi)
        break;
    }
    ly->n_inner += inner + 1;
    for(i = 0; i < ns; i++){
      if(fabs(Vs[i] - ly->V[i]) > ch)
        ch = fabs(Vs[i] - ly->V[i]);
      if(fabs(Vs[i]) > vm)
        vm = fabs(Vs[i]);
    }
    {
      /* would the fault see a different known part with the new layer rates? */
      dmax = 0.0;
      smax = 1.0;
      memset(Yf, 0, sizeof(double _Complex)*Mh*nz);
      for(i = 0; i < ns; i++)
        tmps[i] = Vs[i] - ly->V[i];
      layer_rows_hat(ly, tmps, Xs);
      lx_apply(ly, nz, nys, ly->K1sf_h, Xs, Yf);
      memset(tmpf, 0, sizeof(double)*N);
      lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, tmpf);
      for(i = 0; i < N; i++){
        if(fabs(tmpf[i]) > dmax)
          dmax = fabs(tmpf[i]);
        if(fabs(f->sigma[i]) > smax)
          smax = fabs(f->sigma[i]);
      }
      memcpy(ly->V, Vs, sizeof(double)*ns);
      if(dmax <= 1e-12*smax){
        it++;
        break;
      }      /* converged: the fault solve already used these rates */
    }
    for(i = 0; i < N; i++){
      if(fabs(f->V[i] - Vp[i]) > ch)
        ch = fabs(f->V[i] - Vp[i]);
      if(fabs(f->V[i]) > vm)
        vm = fabs(f->V[i]);
    }
    if(it > 0 && ch <= 1e-12*vm)
      break;
  }
  ly->n_outer += it + 1;
  ly->n_steps++;
  free(Hf_prev);
  free(mask);
  tp1 = ly_now();
  fd3d_layer_prof[2] += tp1 - tp0;
  tp0 = tp1;
  f->thnew = NULL;
  if(thbuf){
    memcpy(f->theta, thbuf, sizeof(double)*N);
    free(thbuf);
  }
  /* finish: fault known part with the layer's current step, histories, slips */
  memset(Yf, 0, sizeof(double _Complex)*Mh*nz);
  layer_rows_hat(ly, ly->V, Xs);
  lx_apply(ly, nz, nys, ly->K1sf_h, Xs, Yf);
  memset(tmpf, 0, sizeof(double)*N);
  lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, tmpf);
  for(i = 0; i < N; i++)
    f->H[i] += Hx[i] + tmpf[i];
  {
    slot_2 = ly->hk % ly->L;
    vms = malloc(sizeof(double)*ns);
    vmf = malloc(sizeof(double)*N);
    for(i = 0; i < ns; i++){
      vms[i] = f->nl ? 0.5*(ly->Vold[i] + ly->V[i]) : ly->V[i];
      ly->delta[i] += vms[i]*dt;
    }
    for(i = 0; i < N; i++)
      vmf[i] = f->nl ? 0.5*(f->Vold[i] + f->V[i]) : f->V[i];
    l2_fwd(ly, vms, ly->hs2 + (size_t)slot_2*H2);
    layer_rows_hat(ly, vms, ly->hsx + (size_t)slot_2*Mh*nys);
    fault_rows_hat(ly, vmf, ly->hfx + (size_t)slot_2*Mh*nz);
    ly->hk = (ly->hk + 1) % ly->L;
    if(ly->nhist < ly->L)
      ly->nhist++;
    free(vms);
    free(vmf);
  }
  step_end(f);
  tp1 = ly_now();
  fd3d_layer_prof[3] += tp1 - tp0;
  free(Hx);
  free(Vs);
  free(Vp);
  free(rhs);
  free(tmpf);
  free(tmps);
  free(Xs);
  free(Xf);
  free(Ys);
  free(Yf);
  free(S2);
  free(P2);
  return itot;
}

/* ---------------- full virtual surface for dipping faults (stage B) ---------------- */

static void frame_fault(double dip_deg, double F[3][3])
{
  double d = dip_deg*M_PI/180.0;
  F[0][0] = 1.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[1][0] = 0.0;
  F[1][1] = cos(d);
  F[1][2] = sin(d);                 /* up dip */
  F[2][0] = 0.0;
  F[2][1] = -sin(d);
  F[2][2] = cos(d);                /* e1 x e2 */
}
/* global stress tensor from local s6 (11 22 33 12 31 32) in frame R (rows = local axes) */
static void to_global(double R[3][3], const double s[6], double G[3][3])
{
  double S[3][3] = {{s[0], s[3], s[4]}, {s[3], s[1], s[5]}, {s[4], s[5], s[2]}};
  int i, j, a, b;
  double v;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++){
      v = 0.0;
      for(a = 0; a < 3; a++)
        for(b = 0; b < 3; b++)
          v += R[a][i]*S[a][b]*R[b][j];
      G[i][j] = v;
    }
}
/* tractions at a receiver for unit rate on a source rectangle, all lags plus the static slope:
   src_fault: fault frame (comp fcomp) else surface frame (comp 1, 2, 3);
   rcv_fault: (shear along s, normal on n) else (sxz, syz, szz); K[r][l], ks[r] */
static void surf_resp(const struct fd3d_tab *t, double F[3][3], int src_fault, int comp, const double c[3], double a1, double a2,
    const double X[3], int rcv_fault, const double sdir[3], int L, double *K, double *ks)
{
  double Ide[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, (*R)[3] = src_fault ? F : Ide, d[3], x[3], s[6], G[3][3], tr[3], prev[3] = {0, 0, 0}, late[2][3], T1, T2;
  int i, l, q, nr = rcv_fault ? 2 : 3;
  for(i = 0; i < 3; i++)
    d[i] = X[i] - c[i];
  for(i = 0; i < 3; i++)
    x[i] = R[i][0]*d[0] + R[i][1]*d[1] + R[i][2]*d[2];
  T1 = 3.0*(sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) + 4.0*(a1 + a2))/t->cT + 10.0*t->dt;
  T2 = 2.0*T1;
#define SURF_TR(G, tr) do{ if(rcv_fault){ double *n = F[2]; \
      tr[0] = 0.0; tr[1] = 0.0; \
      { int u, v; for(u = 0; u < 3; u++) for(v = 0; v < 3; v++){ tr[0] += sdir[u]*G[u][v]*n[v]; tr[1] += n[u]*G[u][v]*n[v]; } } \
    }else{ tr[0] = G[0][2]; tr[1] = G[1][2]; tr[2] = G[2][2]; } }while(0)
  for(q = 0; q < 2; q++){
    fd3d_rect_ramp_tensor(x[0], x[1], x[2], q ? T2 : T1, a1, a2, t->cL, t->cT, t->rho, comp, s);
    to_global(R, s, G);
    SURF_TR(G, late[q]);
  }
  for(i = 0; i < nr; i++)
    ks[i] = (late[1][i] - late[0][i])/(T2 - T1);
  for(l = 1; l <= L; l++){
    fd3d_rect_ramp_tensor(x[0], x[1], x[2], l*t->dt, a1, a2, t->cL, t->cT, t->rho, comp, s);
    to_global(R, s, G);
    SURF_TR(G, tr);
    for(i = 0; i < nr; i++){
      K[(size_t)i*L + l - 1] = tr[i] - prev[i];
      prev[i] = tr[i];
    }
  }
#undef SURF_TR
}

static void fft_x_inplace(const struct fd3d_layer *ly, const double *in, double _Complex *out, double _Complex *buf)
{
  double _Complex *a = malloc(sizeof(double _Complex)*ly->Mx);
  int x;
  for(x = 0; x < ly->Mx; x++)
    a[x] = in[x];
  fd_fft_exec(ly->pf, a, buf);
  for(x = 0; x < ly->Mh; x++)
    out[x] = buf[x];
  free(a);
}

struct fd3d_layer *fd3d_surf_create_rows(struct fd3d *f, double x0f, double ytop, double ztop, double dip_deg, int fcomp,
    int m, int nxs, int nys, double xs0, const double *yc, const double *wy, int L, double fsign);
/* several faults sharing one virtual surface (fd_plan.md G3b): the second and further faults
   need only their fault <-> surface tables; the surface-surface tables are built once */
static int g_surf_skip_T2 = 0;
struct fd3d_layer *fd3d_surf_create(struct fd3d *f, double x0f, double ytop, double ztop, double dip_deg, int fcomp,
    int m, int nxs, int nys, double xs0, double ys0, int L, double fsign)
{
  double hx = f->t->dx/m, *yc = malloc(sizeof(double)*nys), *wy = malloc(sizeof(double)*nys);
  struct fd3d_layer *ly;
  int q;
  for(q = 0; q < nys; q++){
    yc[q] = ys0 + (q + 0.5)*hx;
    wy[q] = hx;
  }
  ly = fd3d_surf_create_rows(f, x0f, ytop, ztop, dip_deg, fcomp, m, nxs, nys, xs0, yc, wy, L, fsign);
  free(yc);
  free(wy);
  return ly;
}
/* rows of given centres yc[q] and widths wy[q] (any grading; receivers must not fall within about
   1e-2 of a width from another row's edge line unless exactly on it, see fd_plan.md B5): the
   surface-surface tables are then Toeplitz along strike with dense (row) x (row) blocks */
struct fd3d_layer *fd3d_surf_create_rows(struct fd3d *f, double x0f, double ytop, double ztop, double dip_deg, int fcomp,
    int m, int nxs, int nys, double xs0, const double *yc, const double *wy, int L, double fsign)
{
  static const int pr[5] = {0, 0, 1, 1, 2}, ps[5] = {0, 1, 0, 1, 2};
  const struct fd3d_tab *t = f->t;
  struct fd3d_layer *ly;
  int nx = t->nx, nz = f->nzp, d, dmin, dmax, q, j, l, c, r, p, k, N = f->N;
  double dx = t->dx, hx = dx/m, F[3][3], sdir[3], *K = malloc(sizeof(double)*3*L), ks[3];
  size_t Mh, nfs;
  size_t nss;
  size_t nkk;
  size_t Lx;
  double *Ap;
  double cen[3];
  _Complex double *buf;
  _Complex double *col;
  double *line;
  int qs;
  int qr;
  int x;
  int uniform;
  size_t s1;
  size_t s0;
  size_t pos;
  double X[3];
  double *Aq;
  _Complex double *dst;
  size_t Lx_2;
  size_t nkk_2;
  double *Afs;
  double *Asf;
  _Complex double *buf_2;
  _Complex double *col_2;
  double yq;
  double yj;
  double zj;
  size_t pos_2;
  double cf[3];
  double Xs[3];
  double cs[3];
  double Xf[3];
  int e;
  double *K2;
  double ks2[3];
  double ce[3];
  double *Ar;
  double *Ar_2;
  size_t sfs;
  size_t ssf;
  _Complex double *dst_2;
  _Complex double *dst_3;
  if(f->mirror || m < 1 || (fcomp != 1 && fcomp != 2))
    return NULL;
  frame_fault(dip_deg, F);
  for(c = 0; c < 3; c++)
    sdir[c] = fsign*F[fcomp - 1][c];          /* slip direction of the marcher's slip variable */
  ly = calloc(1, sizeof(*ly));
  ly->fft = 1;
  ly->surf = 1;
  ly->fcomp = fcomp;
  ly->dip = dip_deg;
  ly->ytop = ytop;
  ly->ztop = ztop;
  ly->m = m;
  ly->nxs = nxs;
  ly->nys = nys;
  ly->ns = nxs*nys;
  ly->L = L;
  ly->hx = hx;
  ly->hy = wy[0];
  ly->xs0 = xs0;
  ly->ys0 = yc[0] - 0.5*wy[0];
  ly->x0f = x0f;
  ly->nx = nx;
  ly->nz = nz;
  ly->eta_p = 0.5*t->rho*t->cL;
  ly->yc = malloc(sizeof(double)*nys);
  ly->wy = malloc(sizeof(double)*nys);
  memcpy(ly->yc, yc, sizeof(double)*nys);
  memcpy(ly->wy, wy, sizeof(double)*nys);
  ly->self_c = calloc((size_t)3*nys, sizeof(double));
  ly->c0 = xs0 + 0.5*hx - x0f - 0.5*dx;
  ly->Mx = fd_fft_size_for((nxs + m*(nx - 1) + 1) > (2*nxs - 1) ? (nxs + m*(nx - 1) + 1) : (2*nxs - 1));   /* fault <-> surface and surface <-> surface offsets */
  ly->Mh = ly->Mx/2 + 1;
  Mh = ly->Mh;
  /* MPI (fd_plan.md S3): the tables are kept for this rank's block of along-strike wavenumbers
     only; the rings and work spectra stay full, partial sums are combined by MPI_Allreduce */
  ly->mpi_size = t->mpi_size;
  ly->k0 = (int)(((long)Mh*t->mpi_rank)/t->mpi_size);
  ly->nk = (int)(((long)Mh*(t->mpi_rank + 1))/t->mpi_size) - ly->k0;
#ifdef FD3D_MPI
  ly->comm = t->comm;
#endif
  ly->Mx2 = fd_fft_size_for(2*nxs - 1);
  ly->My2 = fd_fft_size_for(2*nys - 1);
  ly->Mh2 = ly->My2/2 + 1;
  nfs = (size_t)nys*nz;
  ly->pf = fd_fft_plan_create(ly->Mx, FD_FFT_FORWARD);
  ly->pb = fd_fft_plan_create(ly->Mx, FD_FFT_BACKWARD);
  ly->pf2x = fd_fft_plan_create(ly->Mx2, FD_FFT_FORWARD);
  ly->pb2x = fd_fft_plan_create(ly->Mx2, FD_FFT_BACKWARD);
  ly->pf2y = fd_fft_plan_create(ly->My2, FD_FFT_FORWARD);
  ly->pb2y = fd_fft_plan_create(ly->My2, FD_FFT_BACKWARD);
  /* surface-surface: pairs (rcv, src) (0,0) (0,1) (1,0) (1,1) (2,2), Toeplitz along strike (offset
     d = p' - p), dense in (receiver row, source row); row widths may differ. Built one row pair
     at a time (temporary memory O(L Mx)); each rank keeps its wavenumbers (fd_plan.md S3) */
  if(!g_surf_skip_T2) {
    nss = (size_t)nys*nys;
    nkk = (size_t)ly->nk;
    Lx = (size_t)(L + 2)*ly->Mx;
    Ap = malloc(sizeof(double)*5*Lx);
    buf = malloc(sizeof(double _Complex)*ly->Mx);
    col = malloc(sizeof(double _Complex)*Mh);

    uniform = 1;
    for(p = 0; p < 5; p++){
      ly->T2[p][0] = malloc(sizeof(double _Complex)*(size_t)L*nkk*nss + 1);
      ly->T2[p][1] = malloc(sizeof(double _Complex)*nkk*nss + 1);
      ly->T2[p][2] = malloc(sizeof(double _Complex)*nkk*nss + 1);
    }
    for(qs = 1; qs < nys; qs++)
      if(fabs(wy[qs] - wy[0]) > 1e-9*wy[0] || fabs((yc[qs] - yc[qs - 1]) - wy[0]) > 1e-9*wy[0])
        uniform = 0;
    for(qs = 0; qs < nys; qs++){
      cen[0] = 0.0;
      cen[1] = yc[qs];
      cen[2] = 0.0;
      for(qr = 0; qr < nys; qr++){
        s1 = (size_t)qr*nys + qs;
        if(uniform && qs > 0 && qr > 0){
          /* uniform rows: the response depends on qr - qs only, copy the entry of (qr - 1, qs - 1) */
          s0 = (size_t)(qr - 1)*nys + qs - 1;
          for(p = 0; p < 5; p++){
            for(l = 0; l < L; l++)
              for(k = 0; k < (int)nkk; k++)
                ly->T2[p][0][((size_t)l*nkk + k)*nss + s1] = ly->T2[p][0][((size_t)l*nkk + k)*nss + s0];
            for(k = 0; k < (int)nkk; k++){
              ly->T2[p][1][(size_t)k*nss + s1] = ly->T2[p][1][(size_t)k*nss + s0];
              ly->T2[p][2][(size_t)k*nss + s1] = ly->T2[p][2][(size_t)k*nss + s0];
            }
          }
          if(qr == qs)
            for(r = 0; r < 3; r++)
              ly->self_c[(size_t)r*nys + qs] = ly->self_c[(size_t)r*nys];
          continue;
        }
        memset(Ap, 0, sizeof(double)*5*Lx);
        for(d = -(nxs - 1); d <= nxs - 1; d++){
          pos = (size_t)((d + ly->Mx) % ly->Mx);
          X[0] = d*hx;
          X[1] = yc[qr];
          X[2] = 0.0;
          for(c = 0; c < 3; c++){
            surf_resp(t, F, 0, c + 1, cen, 0.5*hx, 0.5*wy[qs], X, 0, sdir, L, K, ks);
            for(p = 0; p < 5; p++)
              if(ps[p] == c){
                Aq = Ap + (size_t)p*Lx;
                r = pr[p];
                for(l = 0; l < L; l++)
                  Aq[(size_t)l*ly->Mx + pos] = K[(size_t)r*L + l] - ks[r]*t->dt;
                Aq[(size_t)L*ly->Mx + pos] = ks[r];
                Aq[(size_t)(L + 1)*ly->Mx + pos] = K[(size_t)r*L];
                if(d == 0 && qr == qs && (p == 0 || p == 3 || p == 4))
                  ly->self_c[(size_t)r*nys + qs] = K[(size_t)r*L];
              }
          }
        }
        for(p = 0; p < 5; p++)
          for(l = 0; l < L + 2; l++){
            dst = l < L ? ly->T2[p][0] + (size_t)l*nkk*nss : ly->T2[p][l - L + 1];
            line = Ap + (size_t)p*Lx + (size_t)l*ly->Mx;
            fft_x_inplace(ly, line, col, buf);
            for(k = ly->k0; k < ly->k0 + ly->nk; k++)
              dst[(size_t)(k - ly->k0)*nss + s1] = col[k];
          }
      }
    }
    ly->self_sh = ly->self_c[0];
    ly->self_op = ly->self_c[(size_t)2*nys];
    free(Ap);
    free(buf);
    free(col);
    (void)x;
  }
  /* fault <-> surface, one (surface row, fault row) pair at a time */
  dmin = -(nxs - 1) < -m*(nx - 1) ? -(nxs - 1) : -m*(nx - 1);
  dmax = (nxs - 1) > m*(nx - 1) ? (nxs - 1) : m*(nx - 1);
  {
    Lx_2 = (size_t)(L + 2)*ly->Mx;
    nkk_2 = (size_t)ly->nk;
    Afs = malloc(sizeof(double)*3*Lx_2);
    Asf = malloc(sizeof(double)*6*Lx_2);
    buf_2 = malloc(sizeof(double _Complex)*ly->Mx);
    col_2 = malloc(sizeof(double _Complex)*Mh);
    for(r = 0; r < 3; r++){
      ly->Tfs[r][0] = malloc(sizeof(double _Complex)*(size_t)L*nkk_2*nfs + 1);
      ly->Tfs[r][1] = malloc(sizeof(double _Complex)*nkk_2*nfs + 1);
      ly->Tfs[r][2] = malloc(sizeof(double _Complex)*nkk_2*nfs + 1);
    }
    for(r = 0; r < 2; r++)
      for(c = 0; c < 3; c++){
        ly->Tsf[r][c][0] = malloc(sizeof(double _Complex)*(size_t)L*nkk_2*nfs + 1);
        ly->Tsf[r][c][1] = malloc(sizeof(double _Complex)*nkk_2*nfs + 1);
        ly->Tsf[r][c][2] = malloc(sizeof(double _Complex)*nkk_2*nfs + 1);
      }
    for(q = 0; q < nys; q++){
      yq = yc[q];
      for(j = 0; j < nz; j++){
        yj = ytop - (j + 0.5)*dx*F[1][1];
        zj = ztop - (j + 0.5)*dx*F[1][2];
        memset(Afs, 0, sizeof(double)*3*Lx_2);
        memset(Asf, 0, sizeof(double)*6*Lx_2);
        for(d = dmin; d <= dmax; d++){
          pos_2 = (size_t)((d + ly->Mx) % ly->Mx);
          cf[0] = 0.0;
          cf[1] = yj;
          cf[2] = zj;
          Xs[0] = d*hx + ly->c0;
          Xs[1] = yq;
          Xs[2] = 0.0;
          cs[0] = 0.0;
          cs[1] = yq;
          cs[2] = 0.0;
          Xf[0] = d*hx - ly->c0;
          Xf[1] = yj;
          Xf[2] = zj;
          if(d >= -m*(nx - 1) && d <= nxs - 1){        /* fault -> surface: d = p - m i */
            surf_resp(t, F, 1, fcomp, cf, 0.5*dx, 0.5*dx, Xs, 0, sdir, L, K, ks);
            if(j == 0 && f->next){                       /* extension rows above the surface, tied to row 0 */
              K2 = malloc(sizeof(double)*3*L);
              for(e = 0; e < f->next; e++){
                ce[0] = 0.0;
                ce[1] = ytop + (e + 0.5)*dx*F[1][1];
                ce[2] = ztop + (e + 0.5)*dx*F[1][2];
                surf_resp(t, F, 1, fcomp, ce, 0.5*dx, 0.5*dx, Xs, 0, sdir, L, K2, ks2);
                for(r = 0; r < 3; r++){
                  ks[r] += ks2[r];
                  for(l = 0; l < L; l++)
                    K[(size_t)r*L + l] += K2[(size_t)r*L + l];
                }
              }
              free(K2);
            }
            for(r = 0; r < 3; r++){
              ks[r] *= fsign;
              for(l = 0; l < L; l++)
                K[(size_t)r*L + l] *= fsign;
            }
            for(r = 0; r < 3; r++){
              Ar = Afs + (size_t)r*Lx_2;
              for(l = 0; l < L; l++)
                Ar[(size_t)l*ly->Mx + pos_2] = K[(size_t)r*L + l] - ks[r]*t->dt;
              Ar[(size_t)L*ly->Mx + pos_2] = ks[r];
              Ar[(size_t)(L + 1)*ly->Mx + pos_2] = K[(size_t)r*L];
            }
          }
          if(d >= -(nxs - 1) && d <= m*(nx - 1))
            for(c = 0; c < 3; c++){   /* surface -> fault: d = m i - p */
              surf_resp(t, F, 0, c + 1, cs, 0.5*hx, 0.5*wy[q], Xf, 1, sdir, L, K, ks);
              for(r = 0; r < 2; r++){
                Ar_2 = Asf + (size_t)(r*3 + c)*Lx_2;
                for(l = 0; l < L; l++)
                  Ar_2[(size_t)l*ly->Mx + pos_2] = K[(size_t)r*L + l] - ks[r]*t->dt;
                Ar_2[(size_t)L*ly->Mx + pos_2] = ks[r];
                Ar_2[(size_t)(L + 1)*ly->Mx + pos_2] = K[(size_t)r*L];
              }
            }
        }
        for(l = 0; l < L + 2; l++){
          sfs = (size_t)q*nz + j;
          ssf = (size_t)j*nys + q;
          for(r = 0; r < 3; r++){
            dst_2 = l < L ? ly->Tfs[r][0] + (size_t)l*nkk_2*nfs : ly->Tfs[r][l - L + 1];
            fft_x_inplace(ly, Afs + (size_t)r*Lx_2 + (size_t)l*ly->Mx, col_2, buf_2);
            for(k = ly->k0; k < ly->k0 + ly->nk; k++)
              dst_2[(size_t)(k - ly->k0)*nfs + sfs] = col_2[k];
          }
          for(r = 0; r < 2; r++)
            for(c = 0; c < 3; c++){
              dst_3 = l < L ? ly->Tsf[r][c][0] + (size_t)l*nkk_2*nfs : ly->Tsf[r][c][l - L + 1];
              fft_x_inplace(ly, Asf + (size_t)(r*3 + c)*Lx_2 + (size_t)l*ly->Mx, col_2, buf_2);
              for(k = ly->k0; k < ly->k0 + ly->nk; k++)
                dst_3[(size_t)(k - ly->k0)*nfs + ssf] = col_2[k];
            }
        }
      }
    }
    free(Afs);
    free(Asf);
    free(buf_2);
    free(col_2);
  }
  free(K);
  ly->V3 = calloc((size_t)3*ly->ns, sizeof(double));
  ly->Vold3 = calloc((size_t)3*ly->ns, sizeof(double));
  ly->delta3 = calloc((size_t)3*ly->ns, sizeof(double));
  ly->H3 = calloc((size_t)3*ly->ns, sizeof(double));
  ly->hsx3 = calloc((size_t)3*L*Mh*nys, sizeof(double _Complex));
  ly->hfx = calloc((size_t)L*Mh*nz, sizeof(double _Complex));
  ly->sig0 = calloc(N, sizeof(double));
  ly->tn = calloc(N, sizeof(double));
  ly->normal_coupling = 1;
  return ly;
}

static void surf_free(struct fd3d_layer *ly)
{
  int p, r, c, s;
  for(p = 0; p < 5; p++)
    for(s = 0; s < 3; s++)
      free(ly->T2[p][s]);
  for(r = 0; r < 3; r++)
    for(s = 0; s < 3; s++)
      free(ly->Tfs[r][s]);
  for(r = 0; r < 2; r++)
    for(c = 0; c < 3; c++)
      for(s = 0; s < 3; s++)
        free(ly->Tsf[r][c][s]);
  free(ly->V3);
  free(ly->Vold3);
  free(ly->delta3);
  free(ly->H3);
  free(ly->hsx3);
  free(ly->yc);
  free(ly->wy);
  free(ly->self_c);
  free(ly->pre_lu);
  free(ly->pre_piv);
  free(ly->hfx);
  free(ly->sig0);
  free(ly->tn);
  fd_fft_plan_destroy(ly->pf);
  fd_fft_plan_destroy(ly->pb);
  fd_fft_plan_destroy(ly->pf2x);
  fd_fft_plan_destroy(ly->pb2x);
  fd_fft_plan_destroy(ly->pf2y);
  fd_fft_plan_destroy(ly->pb2y);
  free(ly);
}

/* Y[k][nr] += T[k - k0][nr][nc] X[k][nc] over this rank's wavenumbers (surface tables) */
static void lxp_apply(const struct fd3d_layer *ly, int nr, int nc, const double _Complex *T, const double _Complex *X, double _Complex *Y)
{
  int k, r, c;
  const _Complex double *Tk;
  const _Complex double *Xk;
  _Complex double *Yk;
  _Complex double s;
  for(k = ly->k0; k < ly->k0 + ly->nk; k++){
    Tk = T + (size_t)(k - ly->k0)*nr*nc;
    Xk = X + (size_t)k*nc;
    Yk = Y + (size_t)k*nr;
    for(r = 0; r < nr; r++){
      s = 0.0;
      for(c = 0; c < nc; c++)
        s += Tk[(size_t)r*nc + c]*Xk[c];
      Yk[r] += s;
    }
  }
}
/* combine the wavenumber blocks of a spectrum filled by lxp_apply / surf_apply_rows */
static void surf_reduce(const struct fd3d_layer *ly, double _Complex *Y, size_t n)
{
#ifdef FD3D_MPI
  if(ly->mpi_size > 1)
    MPI_Allreduce(MPI_IN_PLACE, Y, (int)n, MPI_C_DOUBLE_COMPLEX, MPI_SUM, ly->comm);
#else
  (void)ly;
  (void)Y;
  (void)n;
#endif
}
/* Y_r[k][q'] += w sum_c T[r][c][k][q'][q] X_c[k][q] for the surface-surface tables (row spectra) */
static const int SP_R[5] = {0, 0, 1, 1, 2}, SP_S[5] = {0, 1, 0, 1, 2};
static void surf_apply_rows(const struct fd3d_layer *ly, int slot_kind, int lag, const double _Complex *X3, double w, double _Complex *Y3)
{
  size_t nss = (size_t)ly->nys*ly->nys, Mh = ly->Mh, Q = Mh*ly->nys, k;
  int p, r, s;
  const _Complex double *T;
  const _Complex double *X;
  _Complex double *Y;
  const _Complex double *Tk;
  const _Complex double *Xk;
  _Complex double *Yk;
  _Complex double a;
  for(p = 0; p < 5; p++){
    T = slot_kind == 0 ? ly->T2[p][0] + (size_t)lag*ly->nk*nss : ly->T2[p][slot_kind];
    X = X3 + (size_t)SP_S[p]*Q;
    Y = Y3 + (size_t)SP_R[p]*Q;
    for(k = (size_t)ly->k0; k < (size_t)(ly->k0 + ly->nk); k++){
      Tk = T + (k - ly->k0)*nss;
      Xk = X + k*ly->nys;
      Yk = Y + k*ly->nys;
      for(r = 0; r < ly->nys; r++){
        a = 0.0;
        for(s = 0; s < ly->nys; s++)
          a += Tk[(size_t)r*ly->nys + s]*Xk[s];
        Yk[r] += w*a;
      }
    }
  }
}


/* ---- the surface's zero-traction solve (fd_plan.md H4) ----
   A v = hw IFFT_x( T2(first lag) FFT_x(v) ) - (nl) Es/2 v = rhs over the surface cells (3
   components). Damped Jacobi needed hundreds of sweeps where the surface edge is near the
   faults (H3). Instead: GMRES, right-preconditioned by the circulant approximation of A (the
   along-strike convolution taken periodic on the FFT length), whose blocks per wavenumber
   (3 nys x 3 nys) are factorised once; FD_SURF_JACOBI=1 restores the old iteration. */
extern void zgetrf_(int *m, int *n, double _Complex *a, int *lda, int *ipiv, int *info);
extern void zgetrs_(char *tr, int *n, int *nrhs, double _Complex *a, int *lda, int *ipiv, double _Complex *b, int *ldb, int *info);
static void surf_A(struct fd3d_layer *ly, int nl, double hw, const double Es[3], const double *v, double *out, double _Complex *S3, double _Complex *P3)
{
  int ns = ly->ns, nys = ly->nys, c, i;
  size_t Q = (size_t)ly->Mh*nys;
  for(c = 0; c < 3; c++)
    layer_rows_hat(ly, v + (size_t)c*ns, S3 + (size_t)c*Q);
  memset(P3, 0, sizeof(double _Complex)*3*Q);
  surf_apply_rows(ly, 2, 0, S3, 1.0, P3);
  surf_reduce(ly, P3, 3*Q);
  memset(out, 0, sizeof(double)*3*ns);
  for(c = 0; c < 3; c++){
    lx_inv(ly, nys, P3 + (size_t)c*Q, ly->nxs, 1, hw, out + (size_t)c*ns);
    if(nl)
      for(i = 0; i < ns; i++)
        out[(size_t)c*ns + i] -= 0.5*Es[c]*v[(size_t)c*ns + i];
  }
}
static void surf_precond_setup(struct fd3d_layer *ly, int nl, double hw, const double Es[3])
{
  static const int pr[5] = {0, 0, 1, 1, 2}, ps[5] = {0, 1, 0, 1, 2};
  int nys = ly->nys, n = 3*nys, k, p, qr, qs, info, c;
  size_t nss = (size_t)nys*nys;
  _Complex double *M;
  const _Complex double *T;
  if(ly->pre_lu)
    return;
  ly->pre_lu = malloc(sizeof(double _Complex)*(size_t)(ly->nk > 0 ? ly->nk : 1)*n*n);
  ly->pre_piv = malloc(sizeof(int)*(size_t)(ly->nk > 0 ? ly->nk : 1)*n);
  for(k = 0; k < ly->nk; k++){
    M = ly->pre_lu + (size_t)k*n*n;
    memset(M, 0, sizeof(double _Complex)*n*n);
    for(p = 0; p < 5; p++){
      T = ly->T2[p][2] + (size_t)k*nss;
      for(qr = 0; qr < nys; qr++)
        for(qs = 0; qs < nys; qs++)
          M[(size_t)(ps[p]*nys + qs)*n + (pr[p]*nys + qr)] += hw*T[(size_t)qr*nys + qs];    /* column major: row pr*nys+qr, column ps*nys+qs */
    }
    if(nl)
      for(c = 0; c < 3; c++)
        for(qr = 0; qr < nys; qr++)
          M[(size_t)(c*nys + qr)*n + (c*nys + qr)] -= 0.5*Es[c];
    zgetrf_(&n, &n, M, &n, ly->pre_piv + (size_t)k*n, &info);
  }
}
static void surf_precond(struct fd3d_layer *ly, const double *r, double *z, double _Complex *S3)
{
  int ns = ly->ns, nys = ly->nys, n = 3*nys, c, k, qr, one = 1, info;
  size_t Q = (size_t)ly->Mh*nys;
  char tr = 'N';
  double _Complex *b = malloc(sizeof(double _Complex)*n);
  int own;
  for(c = 0; c < 3; c++)
    layer_rows_hat(ly, r + (size_t)c*ns, S3 + (size_t)c*Q);
  for(k = 0; k < ly->Mh; k++){
    own = (k >= ly->k0) && (k < ly->k0 + ly->nk);
    if(!own){
      for(c = 0; c < 3; c++)
        for(qr = 0; qr < nys; qr++)
          S3[(size_t)c*Q + (size_t)k*nys + qr] = 0.0;
      continue;
    }
    for(c = 0; c < 3; c++)
      for(qr = 0; qr < nys; qr++)
        b[c*nys + qr] = S3[(size_t)c*Q + (size_t)k*nys + qr];
    zgetrs_(&tr, &n, &one, ly->pre_lu + (size_t)(k - ly->k0)*n*n, &n, ly->pre_piv + (size_t)(k - ly->k0)*n, b, &n, &info);
    for(c = 0; c < 3; c++)
      for(qr = 0; qr < nys; qr++)
        S3[(size_t)c*Q + (size_t)k*nys + qr] = b[c*nys + qr];
  }
  surf_reduce(ly, S3, 3*Q);
  memset(z, 0, sizeof(double)*3*ns);
  for(c = 0; c < 3; c++)
    lx_inv(ly, nys, S3 + (size_t)c*Q, ly->nxs, 1, 1.0, z + (size_t)c*ns);
  free(b);
}
/* solves A V3 = rhs (V3: initial guess in, solution out); returns the number of iterations */
static int surf_solve(struct fd3d_layer *ly, int nl, double hw, const double Es[3], const double *rhs, double *V3,
    double _Complex *S3, double _Complex *P3, double *tmps)
{
  int n = 3*ly->ns, nys = ly->nys, i, it = 0, m = 40, j, k2, restart;
  int inner;
  int c;
  int ns;
  double chi;
  double vmi;
  size_t ii;
  double diag;
  double dv;
  double *r;
  double *w;
  double *z;
  double *Vk;
  double *H;
  double cs[64];
  double sn[64];
  double g[65];
  double y[64];
  double bnorm;
  double beta;
  double *vj;
  double *vn;
  double hn;
  double h;
  double *vk;
  double a1;
  double a2;
  double a1_2;
  double a2_2;
  double rr;
  int jj;
  int nj;
  double s;
  if(getenv("FD_SURF_JACOBI")){
    ns = ly->ns;
    for(inner = 0; inner < 500; inner++){
      chi = 0.0;
      vmi = 1e-30;
      surf_A(ly, nl, hw, Es, V3, tmps, S3, P3);
      for(c = 0; c < 3; c++)
        for(i = 0; i < ns; i++){
          ii = (size_t)c*ns + i;
          diag = hw*ly->self_c[(size_t)c*nys + i % nys] - (nl ? 0.5*Es[c] : 0.0);
          dv = (rhs[ii] - tmps[ii])/diag;
          V3[ii] += dv;
          if(fabs(dv) > chi)
            chi = fabs(dv);
          if(fabs(V3[ii]) > vmi)
            vmi = fabs(V3[ii]);
        }
      if(chi <= 1e-13*vmi)
        break;
    }
    return inner + 1;
  }
  surf_precond_setup(ly, nl, hw, Es);
  {
    r = malloc(sizeof(double)*n);
    w = malloc(sizeof(double)*n);
    z = malloc(sizeof(double)*n);
    Vk = malloc(sizeof(double)*(size_t)(m + 1)*n);
    H = calloc((size_t)(m + 1)*m, sizeof(double));
    bnorm = 0.0;
    for(i = 0; i < n; i++)
      bnorm += rhs[i]*rhs[i];
    bnorm = sqrt(bnorm);
    if(bnorm == 0.0){
      memset(V3, 0, sizeof(double)*n);
      free(r);
      free(w);
      free(z);
      free(Vk);
      free(H);
      return 0;
    }
    for(restart = 0; restart < 10; restart++){
      surf_A(ly, nl, hw, Es, V3, w, S3, P3);
      for(i = 0; i < n; i++)
        r[i] = rhs[i] - w[i];
      beta = 0.0;
      for(i = 0; i < n; i++)
        beta += r[i]*r[i];
      beta = sqrt(beta);
      if(beta <= 1e-13*bnorm)
        break;
      for(i = 0; i < n; i++)
        Vk[i] = r[i]/beta;
      memset(g, 0, sizeof(g));
      g[0] = beta;
      for(j = 0; j < m; j++){
        vj = Vk + (size_t)j*n;
        vn = Vk + (size_t)(j + 1)*n;
        surf_precond(ly, vj, z, S3);
        surf_A(ly, nl, hw, Es, z, vn, S3, P3);
        it++;
        for(k2 = 0; k2 <= j; k2++){
          h = 0.0;
          vk = Vk + (size_t)k2*n;
          for(i = 0; i < n; i++)
            h += vn[i]*vk[i];
          H[(size_t)k2*m + j] = h;
          for(i = 0; i < n; i++)
            vn[i] -= h*vk[i];
        }
        hn = 0.0;
        for(i = 0; i < n; i++)
          hn += vn[i]*vn[i];
        hn = sqrt(hn);
        H[(size_t)(j + 1)*m + j] = hn;
        if(hn > 0.0)
          for(i = 0; i < n; i++)
            vn[i] /= hn;
        for(k2 = 0; k2 < j; k2++){
          a1 = H[(size_t)k2*m + j];
          a2 = H[(size_t)(k2 + 1)*m + j];
          H[(size_t)k2*m + j] = cs[k2]*a1 + sn[k2]*a2;
          H[(size_t)(k2 + 1)*m + j] = -sn[k2]*a1 + cs[k2]*a2;
        }
        {
          a1_2 = H[(size_t)j*m + j];
          a2_2 = H[(size_t)(j + 1)*m + j];
          rr = sqrt(a1_2*a1_2 + a2_2*a2_2);
          cs[j] = rr > 0 ? a1_2/rr : 1.0;
          sn[j] = rr > 0 ? a2_2/rr : 0.0;
          H[(size_t)j*m + j] = rr;
          H[(size_t)(j + 1)*m + j] = 0.0;
          g[j + 1] = -sn[j]*g[j];
          g[j] = cs[j]*g[j];
        }
        if(fabs(g[j + 1]) <= 1e-13*bnorm){
          j++;
          break;
        }
      }
      if(j > m)
        j = m;
      {
        nj = j;
        for(jj = nj - 1; jj >= 0; jj--){
          s = g[jj];
          for(k2 = jj + 1; k2 < nj; k2++)
            s -= H[(size_t)jj*m + k2]*y[k2];
          y[jj] = s/H[(size_t)jj*m + jj];
        }
        memset(w, 0, sizeof(double)*n);
        for(k2 = 0; k2 < nj; k2++)
          for(i = 0; i < n; i++)
            w[i] += Vk[(size_t)k2*n + i]*y[k2];
        surf_precond(ly, w, z, S3);
        for(i = 0; i < n; i++)
          V3[i] += z[i];
      }
    }
    free(r);
    free(w);
    free(z);
    free(Vk);
    free(H);
  }
  return it;
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
  int slot;
  _Complex double *Yt;
  double ch;
  double vm;
  double *Hf;
  int slot_2;
  double *vms;
  double *vmf;
  Es[0] = Es[1] = 0.5*t->rho*t->cT;
  Es[2] = 0.5*t->rho*t->cL;
  step_begin(f, load);
  memcpy(ly->Vold3, ly->V3, sizeof(double)*3*ns);
  /* ---- known parts ---- */
  memset(Ys, 0, sizeof(double _Complex)*3*Mh*nys);
  memset(Yf, 0, sizeof(double _Complex)*2*Mh*nz);
  for(c = 0; c < 3; c++)
    layer_rows_hat(ly, ly->delta3 + (size_t)c*ns, S3 + (size_t)c*Q);
  surf_apply_rows(ly, 1, 0, S3, 1.0, Ys);                                         /* statics */
  fault_rows_hat(ly, f->delta, Xf);
  for(r = 0; r < 3; r++)
    lxp_apply(ly, nys, nz, ly->Tfs[r][1], Xf, Ys + (size_t)r*Mh*nys);
  for(c = 0; c < 3; c++){
    layer_rows_hat(ly, ly->delta3 + (size_t)c*ns, Xs);
    if(!f->static_fn)
      lxp_apply(ly, nz, nys, ly->Tsf[0][c][1], Xs, Yf);
    if(!ly->static_n_fn)
      lxp_apply(ly, nz, nys, ly->Tsf[1][c][1], Xs, Yf + Mh*nz);
  }
  for(l = 2; l <= lmax; l++){
    slot = (ly->hk - (l - 1) + 2*ly->L) % ly->L;
    surf_apply_rows(ly, 0, l - 1, ly->hsx3 + (size_t)slot*3*Q, 1.0, Ys);
    for(r = 0; r < 3; r++)
      lxp_apply(ly, nys, nz, ly->Tfs[r][0] + (size_t)(l - 1)*ly->nk*nfs, ly->hfx + (size_t)slot*Mh*nz, Ys + (size_t)r*Mh*nys);
    for(c = 0; c < 3; c++)
      for(r = 0; r < 2; r++)
        lxp_apply(ly, nz, nys, ly->Tsf[r][c][0] + (size_t)(l - 1)*ly->nk*nfs, ly->hsx3 + ((size_t)slot*3 + c)*Mh*nys, Yf + (size_t)r*Mh*nz);
  }
  if(f->nl){
    Yt = malloc(sizeof(double _Complex)*3*Mh*nys);
    for(c = 0; c < 3; c++)
      layer_rows_hat(ly, ly->Vold3 + (size_t)c*ns, S3 + (size_t)c*Q);
    surf_apply_rows(ly, 2, 0, S3, 0.5, Ys);
    fault_rows_hat(ly, f->Vold, Xf);
    memset(Yt, 0, sizeof(double _Complex)*3*Mh*nys);
    for(r = 0; r < 3; r++)
      lxp_apply(ly, nys, nz, ly->Tfs[r][2], Xf, Yt + (size_t)r*Mh*nys);
    for(m = 0; m < 3*Mh*nys; m++)
      Ys[m] += 0.5*Yt[m];
    free(Yt);
    Yt = malloc(sizeof(double _Complex)*2*Mh*nz);
    memset(Yt, 0, sizeof(double _Complex)*2*Mh*nz);
    for(c = 0; c < 3; c++){
      layer_rows_hat(ly, ly->Vold3 + (size_t)c*ns, Xs);
      for(r = 0; r < 2; r++)
        lxp_apply(ly, nz, nys, ly->Tsf[r][c][2], Xs, Yt + (size_t)r*Mh*nz);
    }
    for(m = 0; m < 2*Mh*nz; m++)
      Yf[m] += 0.5*Yt[m];
    free(Yt);
  }
  surf_reduce(ly, Ys, 3*Mh*nys);
  surf_reduce(ly, Yf, 2*Mh*nz);
  memset(ly->H3, 0, sizeof(double)*3*ns);
  for(c = 0; c < 3; c++){
    lx_inv(ly, nys, Ys + (size_t)c*Mh*nys, ly->nxs, 1, 1.0, ly->H3 + (size_t)c*ns);
    if(f->nl)
      for(i = 0; i < ns; i++)
        ly->H3[(size_t)c*ns + i] += 0.5*Es[c]*ly->Vold3[(size_t)c*ns + i];
  }
  lx_inv(ly, nz, Yf, ly->nx, ly->m, 1.0, Hx);
  lx_inv(ly, nz, Yf + Mh*nz, ly->nx, ly->m, 1.0, Hn);
  /* ---- iterate fault and surface ---- */
  memcpy(V3, ly->Vold3, sizeof(double)*3*ns);
  f->thnew = thbuf;
  for(it = 0; it < 200; it++){
    ch = 0.0;
    vm = 1e-30;
    Hf = malloc(sizeof(double)*N);
    memset(Yf, 0, sizeof(double _Complex)*2*Mh*nz);
    for(c = 0; c < 3; c++){
      layer_rows_hat(ly, V3 + (size_t)c*ns, Xs);
      for(r = 0; r < 2; r++)
        lxp_apply(ly, nz, nys, ly->Tsf[r][c][2], Xs, Yf + (size_t)r*Mh*nz);
    }
    surf_reduce(ly, Yf, 2*Mh*nz);
    memcpy(Hf, Hx, sizeof(double)*N);
    lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, Hf);
    memcpy(ly->tn, Hn, sizeof(double)*N);
    lx_inv(ly, nz, Yf + Mh*nz, ly->nx, ly->m, hw, ly->tn);
    if(ly->static_n_fn){
      /* exact static normal traction of the slip at the end of the step (current fault rates) */
      for(i = 0; i < N; i++)
        dend[i] = f->delta[i] + (f->nl ? 0.5*dt*(f->Vold[i] + f->V[i]) : dt*f->V[i]);
      ly->static_n_fn(ly->static_n_ctx, dend, tnst);
      for(i = 0; i < N; i++)
        ly->tn[i] += tnst[i];
    }
    if(ly->normal_coupling)
      for(i = 0; i < N; i++){
        f->sigma[i] = ly->sig0[i] - ly->tn[i];
        if(f->sigma[i] < 0.0)
          f->sigma[i] = 0.0;
      }
    memcpy(Vp, f->V, sizeof(double)*N);
    itot += step_solve_g(f, Hf, it > 0 ? Vp : NULL, NULL);
    free(Hf);
    /* surface: zero traction, Jacobi with the self coefficients */
    memset(Ys, 0, sizeof(double _Complex)*3*Mh*nys);
    fault_rows_hat(ly, f->V, Xf);
    for(r = 0; r < 3; r++)
      lxp_apply(ly, nys, nz, ly->Tfs[r][2], Xf, Ys + (size_t)r*Mh*nys);
    surf_reduce(ly, Ys, 3*Mh*nys);
    for(i = 0; i < 3*ns; i++)
      rhs[i] = -ly->H3[i];
    for(c = 0; c < 3; c++)
      lx_inv(ly, nys, Ys + (size_t)c*Mh*nys, ly->nxs, 1, -hw, rhs + (size_t)c*ns);
    inner = surf_solve(ly, f->nl, hw, Es, rhs, V3, S3, P3, tmps);
    ly->n_inner += inner;
    for(i = 0; i < 3*ns; i++){
      if(fabs(V3[i] - ly->V3[i]) > ch)
        ch = fabs(V3[i] - ly->V3[i]);
      if(fabs(V3[i]) > vm)
        vm = fabs(V3[i]);
    }
    memcpy(ly->V3, V3, sizeof(double)*3*ns);
    for(i = 0; i < N; i++){
      if(fabs(f->V[i] - Vp[i]) > ch)
        ch = fabs(f->V[i] - Vp[i]);
      if(fabs(f->V[i]) > vm)
        vm = fabs(f->V[i]);
    }
    if(it > 0 && ch <= 1e-12*vm)
      break;
  }
  ly->n_outer += it + 1;
  ly->n_steps++;
  f->thnew = NULL;
  if(thbuf){
    memcpy(f->theta, thbuf, sizeof(double)*N);
    free(thbuf);
  }
  /* ---- finish ---- */
  memset(Yf, 0, sizeof(double _Complex)*2*Mh*nz);
  for(c = 0; c < 3; c++){
    layer_rows_hat(ly, ly->V3 + (size_t)c*ns, Xs);
    for(r = 0; r < 2; r++)
      lxp_apply(ly, nz, nys, ly->Tsf[r][c][2], Xs, Yf + (size_t)r*Mh*nz);
  }
  surf_reduce(ly, Yf, 2*Mh*nz);
  memset(tmpf, 0, sizeof(double)*N);
  lx_inv(ly, nz, Yf, ly->nx, ly->m, hw, tmpf);
  for(i = 0; i < N; i++)
    f->H[i] += Hx[i] + tmpf[i];
  memcpy(ly->tn, Hn, sizeof(double)*N);
  lx_inv(ly, nz, Yf + Mh*nz, ly->nx, ly->m, hw, ly->tn);
  {
    slot_2 = ly->hk % ly->L;
    vms = malloc(sizeof(double)*3*ns);
    vmf = malloc(sizeof(double)*N);
    for(i = 0; i < 3*ns; i++){
      vms[i] = f->nl ? 0.5*(ly->Vold3[i] + ly->V3[i]) : ly->V3[i];
      ly->delta3[i] += vms[i]*dt;
    }
    for(i = 0; i < N; i++)
      vmf[i] = f->nl ? 0.5*(f->Vold[i] + f->V[i]) : f->V[i];
    for(c = 0; c < 3; c++)
      layer_rows_hat(ly, vms + (size_t)c*ns, ly->hsx3 + ((size_t)slot_2*3 + c)*Mh*nys);
    fault_rows_hat(ly, vmf, ly->hfx + (size_t)slot_2*Mh*nz);
    ly->hk = (ly->hk + 1) % ly->L;
    if(ly->nhist < ly->L)
      ly->nhist++;
    free(vms);
    free(vmf);
  }
  step_end(f);
  if(ly->static_n_fn){
    ly->static_n_fn(ly->static_n_ctx, f->delta, tnst);
    for(i = 0; i < N; i++)
      ly->tn[i] += tnst[i];
  }
  if(ly->normal_coupling)
    for(i = 0; i < N; i++){
      f->sigma[i] = ly->sig0[i] - ly->tn[i];
      if(f->sigma[i] < 0.0)
        f->sigma[i] = 0.0;
    }
  free(Hx);
  free(Hn);
  free(V3);
  free(Vp);
  free(rhs);
  free(tmpf);
  free(tmps);
  free(Xs);
  free(Xf);
  free(Ys);
  free(Yf);
  free(S3);
  free(P3);
  free(dend);
  free(tnst);
  return itot;
}

/* test hook: one surface/fault kernel (see surf_resp) */
void fd3d_surf_kernel(const struct fd3d_tab *t, double dip_deg, int fcomp, int src_fault, int comp, const double c[3], double a1, double a2,
    const double X[3], int rcv_fault, int L, double *K, double *ks)
{
  double F[3][3], sdir[3];
  int i;
  frame_fault(dip_deg, F);
  for(i = 0; i < 3; i++)
    sdir[i] = F[fcomp - 1][i];
  surf_resp(t, F, src_fault, comp, c, a1, a2, X, rcv_fault, sdir, L, K, ks);
}

/* ---------------- two parallel planes (stepovers), fd_plan.md P1 ---------------- */

/* cross tables from plane q to plane p (normal offset x3 of the receiver plane from the
   source plane): shear traction along the slip and normal traction (tension positive), split
   by kx rows and built over lags in blocks as the self tables */
static void cross_build(struct fd3d_tab *t, double x3, int comp, fd3d_cplx **Wsh, fd3d_cplx **Wno,
    double _Complex **Ksh, double _Complex **Kno, double *k1max, double **ks_off_out)
{
  size_t jm;
  int nx = t->nx, nz = t->nz, ox = 2*nx - 1, oz = 2*nz - 1, i, j, l, r, P = t->mpi_size, which;
  double a = 0.5*t->dx, dt = t->dt, T1, T2, rmax, s[6];
  size_t Mh = (size_t)t->Mh, M2 = (size_t)ox*oz;
  double *A[2], *S0[2], *S1[2];
  int la = 1 + (int)(((long)t->L*t->mpi_rank)/P), lb = 1 + (int)(((long)t->L*(t->mpi_rank + 1))/P), nl = lb - la;
  double v1[2];
  double v2[2];
  _Complex double * mine[2];
  double v[2];
  double *Wl;
  int ii;
  int jj;
  _Complex double *dst;
  fd3d_cplx *Wout;
  int *sc;
  int *sd;
  int *rc;
  int *rd;
  int q;
  _Complex double *sbuf;
  _Complex double *rbuf;
  size_t o;
  int qa;
  int qn;
  int ll;
  int ra;
  int rb;
#define CROSS_K(x1, x2, tt, out2) do{ fd3d_rect_ramp_tensor((x1), (x2), x3, (tt), a, a, t->cL, t->cT, t->rho, comp, s); \
    (out2)[0] = (comp == 1) ? s[4] : s[5]; (out2)[1] = s[2]; }while(0)
  for(which = 0; which < 2; which++){
    A[which] = malloc(sizeof(double)*M2);
    S0[which] = calloc(M2, sizeof(double));
    S1[which] = malloc(sizeof(double)*M2);
  }
  rmax = sqrt(2.0)*(nx > nz ? nx : nz)*t->dx + t->dx + fabs(x3);
  T1 = 1.5*rmax/t->cT + 10.0*dt;
  T2 = 2.0*T1;
  for(i = 0; i < ox; i++)
    for(j = 0; j < oz; j++){
      CROSS_K((i - (nx - 1))*t->dx, (j - (nz - 1))*t->dx, T1, v1);
      CROSS_K((i - (nx - 1))*t->dx, (j - (nz - 1))*t->dx, T2, v2);
      for(which = 0; which < 2; which++)
        A[which][(size_t)i*oz + j] = (v2[which] - v1[which])/(T2 - T1);
    }
  *Ksh = malloc(sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  *Kno = malloc(sizeof(double _Complex)*(size_t)t->Mx*t->Mz);
  if(ks_off_out){
    ks_off_out[0] = malloc(sizeof(double)*M2);
    memcpy(ks_off_out[0], A[0], sizeof(double)*M2);
    ks_off_out[1] = malloc(sizeof(double)*M2);
    memcpy(ks_off_out[1], A[1], sizeof(double)*M2);
  }
  put_offsets(t, A[0], *Ksh);
  fd3d_fft2(t, *Ksh, -1);
  put_offsets(t, A[1], *Kno);
  fd3d_fft2(t, *Kno, -1);
  *Wsh = malloc(sizeof(fd3d_cplx)*(size_t)t->L*t->nxr*Mh);
  *Wno = malloc(sizeof(fd3d_cplx)*(size_t)t->L*t->nxr*Mh);
  *k1max = 0.0;
  {
    for(which = 0; which < 2; which++)
      mine[which] = malloc(sizeof(double _Complex)*((size_t)(nl > 0 ? nl : 1))*t->Mx*Mh);
    if(la > 1)
      for(i = 0; i < ox; i++)
        for(j = 0; j < oz; j++){
          CROSS_K((i - (nx - 1))*t->dx, (j - (nz - 1))*t->dx, (la - 1)*dt, v);
          S0[0][(size_t)i*oz + j] = v[0];
          S0[1][(size_t)i*oz + j] = v[1];
        }
    for(l = la; l < lb; l++){
      for(i = 0; i < ox; i++)
        for(j = 0; j < oz; j++){
          double v[2];
          CROSS_K((i - (nx - 1))*t->dx, (j - (nz - 1))*t->dx, l*dt, v);
          S1[0][(size_t)i*oz + j] = v[0];
          S1[1][(size_t)i*oz + j] = v[1];
          if(l == 1 && fabs(v[0]) > *k1max)
            *k1max = fabs(v[0]);
        }
      for(which = 0; which < 2; which++){
        Wl = malloc(sizeof(double)*M2);
        dst = mine[which] + (size_t)(l - la)*t->Mx*Mh;
        for(i = 0; i < (int)M2; i++)
          Wl[i] = S1[which][i] - S0[which][i] - A[which][i]*dt;
        put_offsets(t, Wl, t->work);
        fd3d_fft2(t, t->work, -1);
        for(ii = 0; ii < t->Mx; ii++)
          for(jj = 0; jj < t->Mh; jj++)
            dst[(size_t)ii*Mh + jj] = t->work[(size_t)ii*t->Mz + jj];
        memcpy(S0[which], S1[which], sizeof(double)*M2);
        free(Wl);
      }
    }
    for(which = 0; which < 2; which++){
      Wout = which ? *Wno : *Wsh;
      if(P == 1){
        for(jm = 0; jm < (size_t)t->L*t->Mx*Mh; jm++)
          Wout[jm] = (fd3d_cplx)mine[which][jm];
      }
#ifdef FD3D_MPI
      else{
        sc = calloc(P, sizeof(int));
        sd = calloc(P, sizeof(int));
        rc = calloc(P, sizeof(int));
        rd = calloc(P, sizeof(int));
        sbuf = malloc(sizeof(double _Complex)*((size_t)(nl > 0 ? nl : 1))*t->Mx*Mh);
        rbuf = malloc(sizeof(double _Complex)*(size_t)t->L*t->nxr*Mh);
        o = 0;
        for(q = 0; q < P; q++){
          qa = t->rdispls[q]/t->Mh;
          qn = t->rcounts[q]/t->Mh;
          sd[q] = (int)o;
          for(ll = 0; ll < nl; ll++){
            memcpy(sbuf + o, mine[which] + ((size_t)ll*t->Mx + qa)*Mh, sizeof(double _Complex)*qn*Mh);
            o += (size_t)qn*Mh;
          }
          sc[q] = (int)o - sd[q];
        }
        o = 0;
        for(r = 0; r < P; r++){
          ra = 1 + (int)(((long)t->L*r)/P);
          rb = 1 + (int)(((long)t->L*(r + 1))/P);
          rd[r] = (int)o;
          rc[r] = (rb - ra)*t->nxr*t->Mh;
          o += (size_t)rc[r];
        }
        MPI_Alltoallv(sbuf, sc, sd, MPI_C_DOUBLE_COMPLEX, rbuf, rc, rd, MPI_C_DOUBLE_COMPLEX, t->comm);
        for(jm = 0; jm < (size_t)t->L*t->nxr*Mh; jm++)
          Wout[jm] = (fd3d_cplx)rbuf[jm];
        free(sc);
        free(sd);
        free(rc);
        free(rd);
        free(sbuf);
        free(rbuf);
      }
#endif
      free(mine[which]);
    }
#ifdef FD3D_MPI
    if(P > 1)
      MPI_Allreduce(MPI_IN_PLACE, k1max, 1, MPI_DOUBLE, MPI_MAX, t->comm);
#endif
  }
  for(which = 0; which < 2; which++){
    free(A[which]);
    free(S0[which]);
    free(S1[which]);
  }
#undef CROSS_K
}

struct fd3d_pair *fd3d_pair_create(struct fd3d *f0, struct fd3d *f1, double W, int comp)
{
  struct fd3d_pair *pp;
  struct fd3d_tab *t = f0->t;
  double k1c[2], kself = fabs(t->k1[1][1]);
  /* both planes on the same table; with a mirror image (vertical planes reaching the free
     surface) both carry it, the table then covering the doubled rows */
  if(f1->t != t || f0->mirror != f1->mirror || f0->nzp != f1->nzp || f0->layer || f1->layer)
    return NULL;
  pp = calloc(1, sizeof(*pp));
  pp->f[0] = f0;
  pp->f[1] = f1;
  pp->W = W;
  /* receiver plane 1 at +W from source plane 0, and plane 0 at -W from plane 1 */
  cross_build(t, W, comp, &pp->Wsh[0], &pp->Wno[0], &pp->Ksh[0], &pp->Kno[0], &k1c[0], pp->ks_off[0]);   /* 0 -> 1 */
  cross_build(t, -W, comp, &pp->Wsh[1], &pp->Wno[1], &pp->Ksh[1], &pp->Kno[1], &k1c[1], pp->ks_off[1]);  /* 1 -> 0 */
  pp->k1c = k1c[0] > k1c[1] ? k1c[0] : k1c[1];
  if(pp->k1c > 1e-6*kself){
    fd3d_pair_free(pp);
    return NULL;
  }   /* the planes interact within one step: not supported */
  pp->normal_coupling = 0;
  pp->sig0[0] = calloc(f0->N, sizeof(double));
  pp->sig0[1] = calloc(f1->N, sizeof(double));
  pp->tn[0] = calloc(f0->N, sizeof(double));
  pp->tn[1] = calloc(f1->N, sizeof(double));
  return pp;
}

void fd3d_pair_free(struct fd3d_pair *pp)
{
  int k;
  if(!pp)
    return;
  for(k = 0; k < 2; k++){
    free(pp->Wsh[k]);
    free(pp->Wno[k]);
    free(pp->Ksh[k]);
    free(pp->Kno[k]);
    free(pp->sig0[k]);
    free(pp->tn[k]);
    free(pp->sbuf[k]);
    free(pp->ks_off[k][0]);
    free(pp->ks_off[k][1]);
  }
  free(pp);
}

/* cross contribution to plane p from plane q = 1 - p: lags >= 2 of q's history plus (without
   exact statics) the static part of q's slip; shear and normal traction in element order */
static void cross_known(struct fd3d_pair *pp, int p, double *sh, double *no)
{
  struct fd3d *fq = pp->f[1 - p], *fp = pp->f[p];
  const struct fd3d_tab *t = fp->t;
  int q = 1 - p, which, l, i, j, lmax = t->L < fq->nhist + 1 ? t->L : fq->nhist + 1, row0 = fp->mirror ? fp->nzp : fp->next, nzp = fp->nzp;
  size_t Mhh = (size_t)t->Mx*t->Mh, Mho = (size_t)t->nxr*t->Mh, M = (size_t)t->Mx*t->Mz, m;
  double _Complex *acc = malloc(sizeof(double _Complex)*M), *dh = NULL;
  const fd3d_cplx *Wt;
  const _Complex double *Ks;
  _Complex double *ah;
  _Complex double *ao;
  double *out;
  int ii;
  int jj;
  const fd3d_cplx *Wl;
  const fd3d_cplx *Vh;
  if(!fp->static_fn){
    dh = malloc(sizeof(double _Complex)*M);
    field_hat(fq, fq->delta, dh);
  }
  for(which = 0; which < 2; which++){
    Wt = which ? pp->Wno[q] : pp->Wsh[q];
    Ks = which ? pp->Kno[q] : pp->Ksh[q];
    ah = calloc(Mhh, sizeof(double _Complex));
    ao = ah + (size_t)t->xr0*t->Mh;
    out = which ? no : sh;

    for(l = 2; l <= lmax; l++){
      Wl = Wt + (size_t)(l - 1)*Mho;
      Vh = fq->hist + (size_t)((fq->hk - (l - 1) + 2*t->L) % t->L)*Mho;
      for(m = 0; m < Mho; m++)
        ao[m] += (double _Complex)(Wl[m]*Vh[m]);
    }
    tab_allgather(t, ah);
    if(dh)
      for(m = 0; m < M; m++)
        acc[m] = Ks[m]*dh[m];
    else
      memset(acc, 0, sizeof(double _Complex)*M);
    for(ii = 0; ii < t->Mx; ii++)
      for(jj = 0; jj < t->Mh; jj++)
        acc[(size_t)ii*t->Mz + jj] += ah[(size_t)ii*t->Mh + jj];
    free(ah);
    from_half(t, acc);
    fd3d_fft2(t, acc, +1);
    for(i = 0; i < t->nx; i++)
      for(j = 0; j < nzp; j++)
        out[i*nzp + j] = creal(acc[(size_t)i*t->Mz + row0 + j])/M;
  }
  free(acc);
  free(dh);
}

/* exact statics of both planes: pp->static_fn fills both outputs once per step, each marcher's
   static_fn then returns its part (and the cross tables' static part is skipped) */
static void pair_static_out0(void *ctx, const double *delta, double *out)
{
  struct fd3d_pair *pp = (struct fd3d_pair *)ctx;
  (void)delta;
  memcpy(out, pp->sbuf[0], sizeof(double)*pp->f[0]->N);
}
static void pair_static_out1(void *ctx, const double *delta, double *out)
{
  struct fd3d_pair *pp = (struct fd3d_pair *)ctx;
  (void)delta;
  memcpy(out, pp->sbuf[1], sizeof(double)*pp->f[1]->N);
}
int fd3d_pair_step(struct fd3d_pair *pp, const double *load0, const double *load1)
{
  int p, i, it = 0;
  double *sh[2], *no[2];
  struct fd3d *f;
  if(pp->static_fn){
    if(!pp->sbuf[0]){
      pp->sbuf[0] = malloc(sizeof(double)*pp->f[0]->N);
      pp->sbuf[1] = malloc(sizeof(double)*pp->f[1]->N);
    }
    pp->static_fn(pp->static_ctx, pp->f[0]->delta, pp->f[1]->delta, pp->sbuf[0], pp->sbuf[1]);
    pp->f[0]->static_fn = pair_static_out0;
    pp->f[0]->static_ctx = pp;
    pp->f[1]->static_fn = pair_static_out1;
    pp->f[1]->static_ctx = pp;
  }
  step_begin(pp->f[0], load0);
  step_begin(pp->f[1], load1);
  for(p = 0; p < 2; p++){
    sh[p] = malloc(sizeof(double)*pp->f[p]->N);
    no[p] = malloc(sizeof(double)*pp->f[p]->N);
    cross_known(pp, p, sh[p], no[p]);
  }
  for(p = 0; p < 2; p++){
    f = pp->f[p];
    if(pp->static_n_fn)
      pp->static_n_fn(pp->static_n_ctx, p, pp->f[1 - p]->delta, f->delta, pp->tn[p]);
    else
      memcpy(pp->tn[p], no[p], sizeof(double)*f->N);
    if(pp->static_n_fn)
      for(i = 0; i < f->N; i++)
        pp->tn[p][i] += no[p][i];
    if(pp->normal_coupling)
      for(i = 0; i < f->N; i++){
        f->sigma[i] = pp->sig0[p][i] - pp->tn[p][i];
        if(f->sigma[i] < 0.0)
          f->sigma[i] = 0.0;
      }
    it += step_solve(f, sh[p]);
    for(i = 0; i < f->N; i++)
      f->H[i] += sh[p][i];
  }
  for(p = 0; p < 2; p++){
    step_end(pp->f[p]);
    free(sh[p]);
    free(no[p]);
  }
  return it;
}

/* ---------------- several planar faults of any orientation (fd_plan.md G1) ----------------
   Each plane p is an fd3d marcher on its own table (whole space), cells (i, j), i along e1
   (strike), j down the dip; centre O_p + (i + 1/2) dx e1 - (j + 1/2) dx e2, e2 up dip,
   e3 = e1 x e2. Slip of plane p: motion of its +e3 side relative to the -e3 side along
   d_p = e1 (comp 1) or e2 (comp 2); shear traction d_p . sigma . e3, normal traction (tension
   positive) e3 . sigma . e3. Cross interactions are dense tables W[l][i][j], built from
   fd3d_rect_ramp_tensor in the source cell's frame and rotated. The current step's cross term
   K1 = S(dt) is kept (intersecting faults): solved by a fixed point between the planes around
   the planes' own local solves. */
static void mat_rot(const double E[3][3], const double s6[6], double S[3][3])
{
  /* s6 in the frame with rows E (s11 s22 s33 s12 s13 s23) to global S = E^T s E */
  double s[3][3] = {{s6[0], s6[3], s6[4]}, {s6[3], s6[1], s6[5]}, {s6[4], s6[5], s6[2]}};
  int a, b, c, d;
  double v;
  for(a = 0; a < 3; a++)
    for(b = 0; b < 3; b++){
      v = 0.0;
      for(c = 0; c < 3; c++)
        for(d = 0; d < 3; d++)
          v += E[c][a]*s[c][d]*E[d][b];
      S[a][b] = v;
    }
}
static void plane_centre(const struct fd3d_mplane *g, int i, int j, double X[3])
{
  int k;
  for(k = 0; k < 3; k++)
    X[k] = g->O[k] + (i + 0.5)*g->dx*g->E[0][k] - (j + 0.5)*g->dx*g->E[1][k];
}
/* traction (shear along the receiver's slip, normal) on receiver plane gp at X from a unit-rate
   ramp of the source cell (iq, jq) of plane gq, at time t */
static void cross_resp1(const struct fd3d_mplane *gq, int iq, int jq, const struct fd3d_mplane *gp, const double X[3],
    double t, double cL, double cT, double rho, int image, double *ts, double *tn)
{
  double Y[3], r[3], x[3], s6[6], S[3][3], d[3], n[3], Sn[3];
  int k, a;
  plane_centre(gq, iq, jq, Y);
  if(image)
    Y[2] = -Y[2];          /* free-surface image (vertical strike-slip planes): mirrored centre, same frame and slip, as fd3d's mirror */
  for(k = 0; k < 3; k++)
    r[k] = X[k] - Y[k];
  for(a = 0; a < 3; a++)
    x[a] = r[0]*gq->E[a][0] + r[1]*gq->E[a][1] + r[2]*gq->E[a][2];
  for(a = 0; a < 3; a++)
    if(fabs(x[a]) < 1e-9*gq->dx)
      x[a] = 0.0;   /* rotation round-off: in-plane receivers stay in-plane */
  if(t <= 0.0){
    *ts = *tn = 0.0;
    return;
  }
  fd3d_rect_ramp_tensor(x[0], x[1], x[2], t, 0.5*gq->dx, 0.5*gq->dx, cL, cT, rho, gq->comp, s6);
  mat_rot(gq->E, s6, S);
  for(k = 0; k < 3; k++){
    d[k] = gp->E[gp->comp - 1][k];
    n[k] = gp->E[2][k];
  }
  for(a = 0; a < 3; a++)
    Sn[a] = S[a][0]*n[0] + S[a][1]*n[1] + S[a][2]*n[2];
  *ts = d[0]*Sn[0] + d[1]*Sn[1] + d[2]*Sn[2];
  *tn = n[0]*Sn[0] + n[1]*Sn[1] + n[2]*Sn[2];
}
static void cross_resp(const struct fd3d_mplane *gq, int iq, int jq, const struct fd3d_mplane *gp, const double X[3],
    double t, double cL, double cT, double rho, double *ts, double *tn)
{
  double ts2;
  double tn2;
  cross_resp1(gq, iq, jq, gp, X, t, cL, cT, rho, 0, ts, tn);
  if(gq->mirror){
    cross_resp1(gq, iq, jq, gp, X, t, cL, cT, rho, 1, &ts2, &tn2);
    *ts += ts2;
    *tn += tn2;
  }
}

struct fd3d_multi *fd3d_multi_create(int np, struct fd3d **f, const struct fd3d_mplane *geom)
{
  struct fd3d_multi *m = calloc(1, sizeof(*m));
  int p, q, Lmax = 0;
  double cL = f[0]->t->cL, cT = f[0]->t->cT, rho = f[0]->t->rho, dt = f[0]->t->dt;
  struct fd3d_mpair *c;
  int nr;
  int ns;
  int i;
  int j;
  int l;
  int L;
  double rmax;
  double T1;
  double T2;
  double X[3];
  double Y[3];
  double rr;
  double zmin;
  double Y_2[3];
  int jj;
  double X_2[3];
  size_t e;
  double s0;
  double n0;
  double s1;
  double n1;
  double a1;
  double b1;
  double a2;
  double b2;
  double As;
  double An;
  int iq;
  int jq;
  double kself;
  double km;
  size_t e_2;
  m->np = np;
  for(p = 0; p < np; p++){
    m->f[p] = f[p];
    m->g[p] = geom[p];
    /* a mirrored plane (half space) must be vertical with strike slip and its top row at the surface */
    if(f[p]->mirror && (!geom[p].mirror || geom[p].comp != 1 || fabs(geom[p].E[2][2]) > 1e-9 || fabs(geom[p].O[2]) > 1e-6*geom[p].dx)){
      free(m);
      return NULL;
    }
    if(f[p]->layer || fabs(f[p]->t->dt - dt) > 1e-12*dt || f[p]->t->mpi_size > 1){
      free(m);
      return NULL;
    }
    /* extension rows above the surface (f->next) are not in the cross tables (fd_plan.md G3b) */
  }
  for(p = 0; p < np; p++)
    for(q = 0; q < np; q++){
      c = &m->c[q][p];               /* source q, receiver p */
      nr = f[p]->N;
      ns = f[q]->N;
      rmax = 0.0;
      if(p == q)
        continue;
      for(i = 0; i < nr; i++){
        plane_centre(&m->g[p], i/f[p]->nzp, i%f[p]->nzp, X);
        for(j = 0; j < ns; j += (ns > 64 ? ns/64 : 1)){
          plane_centre(&m->g[q], j/f[q]->nzp, j%f[q]->nzp, Y);
          rr = sqrt((X[0]-Y[0])*(X[0]-Y[0]) + (X[1]-Y[1])*(X[1]-Y[1]) + (X[2]-Y[2])*(X[2]-Y[2]));
          if(rr > rmax)
            rmax = rr;
        }
      }
      rmax += 2.0*m->g[q].dx;
      if(m->g[q].mirror){                              /* images lie up to twice the depth further */
        zmin = 0.0;
        for(jj = 0; jj < f[q]->nzp; jj++){
          plane_centre(&m->g[q], 0, jj, Y_2);
          if(Y_2[2] < zmin)
            zmin = Y_2[2];
        }
        rmax += 2.0*fabs(zmin);
      }
      L = (int)ceil(rmax/cT/dt) + 12;                  /* past the S arrival from the farthest cell: W vanishes in 3-D */
      T1 = 1.5*rmax/cT + 10.0*dt;
      T2 = 2.0*T1;
      c->L = L;
      c->nr = nr;
      c->ns = ns;
      c->Ws = calloc((size_t)L*nr*ns, sizeof(double));
      c->Wn = calloc((size_t)L*nr*ns, sizeof(double));
      c->Ks = malloc(sizeof(double)*nr*ns);
      c->Kn = malloc(sizeof(double)*nr*ns);
      c->K1s = malloc(sizeof(double)*nr*ns);
      c->K1n = malloc(sizeof(double)*nr*ns);
      for(i = 0; i < nr; i++){
        plane_centre(&m->g[p], i/f[p]->nzp, i%f[p]->nzp, X_2);
        for(j = 0; j < ns; j++){
          e = (size_t)i*ns + j;

          iq = j/f[q]->nzp;
          jq = j%f[q]->nzp;
          cross_resp(&m->g[q], iq, jq, &m->g[p], X_2, T1, cL, cT, rho, &a1, &b1);
          cross_resp(&m->g[q], iq, jq, &m->g[p], X_2, T2, cL, cT, rho, &a2, &b2);
          As = (a2 - a1)/(T2 - T1);
          An = (b2 - b1)/(T2 - T1);
          c->Ks[e] = As;
          c->Kn[e] = An;
          s0 = n0 = 0.0;
          for(l = 1; l <= L; l++){
            cross_resp(&m->g[q], iq, jq, &m->g[p], X_2, l*dt, cL, cT, rho, &s1, &n1);
            if(l == 1){
              c->K1s[e] = s1;
              c->K1n[e] = n1;
            }
            c->Ws[(size_t)(l - 1)*nr*ns + e] = s1 - s0 - As*dt;
            c->Wn[(size_t)(l - 1)*nr*ns + e] = n1 - n0 - An*dt;
            s0 = s1;
            n0 = n1;
          }
        }
      }
      {                                                 /* the cross first lag, kept only where it matters */
        kself = fabs(f[p]->t->k1[1][1]);
        km = 0.0;
        for(e_2 = 0; e_2 < (size_t)nr*ns; e_2++)
          if(fabs(c->K1s[e_2]) > km)
            km = fabs(c->K1s[e_2]);
        c->has_k1 = km > 1e-12*kself;
        if(c->has_k1)
          m->implicit = 1;
      }
      if(L > Lmax)
        Lmax = L;
    }
  for(p = 0; p < np; p++)
    if(f[p]->t->L > Lmax)
      Lmax = f[p]->t->L;   /* the planes' own windows count too (history seeding, quiet window) */
  m->Lb = Lmax + 1;
  for(p = 0; p < np; p++){
    m->ring[p] = calloc((size_t)m->Lb*f[p]->N, sizeof(double));
    m->sig0[p] = calloc(f[p]->N, sizeof(double));
    m->tn[p] = calloc(f[p]->N, sizeof(double));
  }
  m->tol = 1e-12;
  m->maxit = 200;
  return m;
}

void fd3d_multi_free(struct fd3d_multi *m)
{
  int p, q;
  struct fd3d_mpair *c;
  if(!m)
    return;
  for(p = 0; p < m->np; p++){
    free(m->ring[p]);
    free(m->sig0[p]);
    free(m->tn[p]);
    free(m->sbuf[p]);
    for(q = 0; q < m->np; q++){
      c = &m->c[q][p];
      free(c->Ws);
      free(c->Wn);
      free(c->Ks);
      free(c->Kn);
      free(c->K1s);
      free(c->K1n);
    }
  }
  free(m);
}

void fd3d_multi_reset(struct fd3d_multi *m)
{
  int p;
  for(p = 0; p < m->np; p++){
    fd3d_reset(m->f[p]);
    memset(m->ring[p], 0, sizeof(double)*(size_t)m->Lb*m->f[p]->N);
  }
  m->hk = 0;
}

/* push a slice of all planes (history seeding from QD, and after each step) */
void fd3d_multi_push(struct fd3d_multi *m, double **V)
{
  int p;
  for(p = 0; p < m->np; p++){
    memcpy(m->ring[p] + (size_t)(m->hk % m->Lb)*m->f[p]->N, V[p], sizeof(double)*m->f[p]->N);
    fd3d_push(m->f[p], V[p]);
  }
  m->hk++;
}

static void multi_sbuf_out(void *ctx, const double *delta, double *out)
{
  struct fd3d_multi_slot *s = (struct fd3d_multi_slot *)ctx;
  (void)delta;
  memcpy(out, s->m->sbuf[s->p], sizeof(double)*s->m->f[s->p]->N);
}

int fd3d_multi_step(struct fd3d_multi *m, double **load)
{
  int np = m->np, p, q, i, j, l, it, itot = 0;
  double *Kx[FD3D_MULTI_MAX], *Kn[FD3D_MULTI_MAX], *Hs[FD3D_MULTI_MAX], *Vprev[FD3D_MULTI_MAX];
  double * dl[8];
  int nr;
  struct fd3d_mpair *c;
  int ns;
  double a;
  double b;
  const double *rs;
  const double *rn;
  const double *v;
  const double *Ws;
  const double *Wn;
  double a_2;
  double b_2;
  double a_3;
  double b_3;
  double * dl_2[8];
  double *o;
  double ch;
  double vm;
  struct fd3d *f;
  int nr_2;
  double hw;
  double *tnc;
  struct fd3d_mpair *c_2;
  double a_4;
  double b_4;
  struct fd3d *f_2;
  double *sl;
  if(m->static_fn){                                   /* exact statics of all planes from all slips */

    for(p = 0; p < np; p++){
      if(!m->sbuf[p])
        m->sbuf[p] = malloc(sizeof(double)*m->f[p]->N);
      dl[p] = m->f[p]->delta;
    }
    m->static_fn(m->static_ctx, dl, m->sbuf);
    for(p = 0; p < np; p++){
      m->slot[p].m = m;
      m->slot[p].p = p;
      m->f[p]->static_fn = multi_sbuf_out;
      m->f[p]->static_ctx = &m->slot[p];
    }
  }
  for(p = 0; p < np; p++)
    step_begin(m->f[p], load ? load[p] : NULL);
  /* known cross part: statics (tables unless exact) and lags >= 2 */
  for(p = 0; p < np; p++){
    nr = m->f[p]->N;
    Kx[p] = calloc(nr, sizeof(double));
    Kn[p] = calloc(nr, sizeof(double));
    Hs[p] = malloc(sizeof(double)*nr);
    Vprev[p] = malloc(sizeof(double)*nr);
    for(q = 0; q < np; q++){
      c = &m->c[q][p];

      if(p == q)
        continue;
      ns = c->ns;
      if(!m->static_fn)
        for(i = 0; i < nr; i++){
          a = 0.0;
          b = 0.0;
          rs = c->Ks + (size_t)i*ns;
          rn = c->Kn + (size_t)i*ns;
          for(j = 0; j < ns; j++){
            a += rs[j]*m->f[q]->delta[j];
            b += rn[j]*m->f[q]->delta[j];
          }
          Kx[p][i] += a;
          Kn[p][i] += b;
        }
      for(l = 2; l <= c->L && l - 1 <= m->hk; l++){
        v = m->ring[q] + (size_t)((m->hk - (l - 1) + (long)m->Lb*1000) % m->Lb)*ns;
        Ws = c->Ws + (size_t)(l - 1)*nr*ns;
        Wn = c->Wn + (size_t)(l - 1)*nr*ns;
        for(i = 0; i < nr; i++){
          a_2 = 0.0;
          b_2 = 0.0;
          for(j = 0; j < ns; j++){
            a_2 += Ws[(size_t)i*ns + j]*v[j];
            b_2 += Wn[(size_t)i*ns + j]*v[j];
          }
          Kx[p][i] += a_2;
          Kn[p][i] += b_2;
        }
      }
      if(m->f[p]->nl && c->has_k1)
        for(i = 0; i < nr; i++){   /* NL: half of the current-step term from V_old */
          a_3 = 0.0;
          b_3 = 0.0;
          for(j = 0; j < ns; j++){
            a_3 += c->K1s[(size_t)i*ns + j]*m->f[q]->Vold[j];
            b_3 += c->K1n[(size_t)i*ns + j]*m->f[q]->Vold[j];
          }
          Kx[p][i] += 0.5*a_3;
          Kn[p][i] += 0.5*b_3;
        }
    }
  }
  if(m->static_n_fn)
    for(p = 0; p < np; p++){       /* exact static normal traction (all planes' slips) */
      o = malloc(sizeof(double)*m->f[p]->N);
      for(q = 0; q < np; q++)
        dl_2[q] = m->f[q]->delta;
      m->static_n_fn(m->static_n_ctx, p, dl_2, o);
      for(i = 0; i < m->f[p]->N; i++)
        Kn[p][i] += o[i];
      free(o);
    }
  /* solve: fixed point between the planes for the current-step cross terms. The local solve
     advances the state variable (rate and state) to the end of the step, so every repeat
     starts again from the state at the beginning of the step. */
  double *th0[FD3D_MULTI_MAX];
  for(p = 0; p < np; p++){
    memcpy(Vprev[p], m->f[p]->Vold, sizeof(double)*m->f[p]->N);
    th0[p] = NULL;
    if(m->f[p]->law == 1 && m->implicit){
      th0[p] = malloc(sizeof(double)*m->f[p]->N);
      memcpy(th0[p], m->f[p]->theta, sizeof(double)*m->f[p]->N);
    }
  }
  for(it = 0; it < (m->implicit ? m->maxit : 1); it++){
    ch = 0.0;
    vm = 1e-9;
    for(p = 0; p < np; p++){
      f = m->f[p];
      nr_2 = f->N;
      hw = f->nl ? 0.5 : 1.0;
      tnc = calloc(nr_2, sizeof(double));
      for(i = 0; i < nr_2; i++)
        Hs[p][i] = Kx[p][i];
      for(q = 0; q < np; q++){
        c_2 = &m->c[q][p];
        if(p == q || !c_2->has_k1)
          continue;
        for(i = 0; i < nr_2; i++){
          a_4 = 0.0;
          b_4 = 0.0;
          for(j = 0; j < c_2->ns; j++){
            a_4 += c_2->K1s[(size_t)i*c_2->ns + j]*m->f[q]->V[j];
            b_4 += c_2->K1n[(size_t)i*c_2->ns + j]*m->f[q]->V[j];
          }
          Hs[p][i] += hw*a_4;
          tnc[i] += hw*b_4;
        }
      }
      for(i = 0; i < nr_2; i++)
        m->tn[p][i] = Kn[p][i] + tnc[i];
      if(m->normal_coupling)
        for(i = 0; i < nr_2; i++){
          f->sigma[i] = m->sig0[p][i] - m->tn[p][i];
          if(f->sigma[i] < 0.0)
            f->sigma[i] = 0.0;
        }
      free(tnc);
      if(it > 0){
        memcpy(Vprev[p], f->V, sizeof(double)*nr_2);
        if(th0[p])
          memcpy(f->theta, th0[p], sizeof(double)*nr_2);
      }
      itot += step_solve_g(f, Hs[p], it > 0 ? f->V : NULL, NULL);
      for(i = 0; i < nr_2; i++){
        if(fabs(f->V[i] - Vprev[p][i]) > ch)
          ch = fabs(f->V[i] - Vprev[p][i]);
        if(fabs(f->V[i]) > vm)
          vm = fabs(f->V[i]);
      }
    }
    if(m->implicit && it > 0 && ch <= m->tol*vm)
      break;
  }
  m->last_iters = it + 1;
  for(p = 0; p < np; p++){
    f_2 = m->f[p];
    sl = m->ring[p] + (size_t)(m->hk % m->Lb)*f_2->N;
    for(i = 0; i < f_2->N; i++){
      f_2->H[i] += Hs[p][i];
      sl[i] = f_2->nl ? 0.5*(f_2->Vold[i] + f_2->V[i]) : f_2->V[i];
    }
    step_end(f_2);
    free(Kx[p]);
    free(Kn[p]);
    free(Hs[p]);
    free(Vprev[p]);
    free(th0[p]);
  }
  m->hk++;
  return itot;
}

/* ---------------- several faults with one virtual surface (fd_plan.md G3b) ----------------
   Faults sharing the strike direction x and one along-strike lattice (same nx, x0f), each with
   its own dip, top line and slip sign; one surface of graded rows. The surface tables and
   state live in the first fault's layer, each further fault has its own fault <-> surface
   tables (surface-surface tables skipped). Fault <-> fault coupling: the dense whole-space
   tables of fd3d_multi (the virtual surface supplies the free surface for them too). */
struct fd3d_msurf *fd3d_msurf_create(int nf, struct fd3d **f, const double *ytop, const double *ztop, const double *dip_deg,
    const double *fsign, int fcomp, double x0f, int m, int nxs, int nys, double xs0,
    const double *yc, const double *wy, int L)
{
  struct fd3d_msurf *s = calloc(1, sizeof(*s));
  struct fd3d_mplane g[FD3D_MULTI_MAX];
  int k, q;
  double F[3][3];
  struct fd3d_layer *ly;
  double cut;
  int nys_2;
  int nz;
  int nkk;
  int Lk;
  int q2;
  int j;
  int kk;
  int l;
  int r;
  int cc;
  int ncut;
  size_t nfs;
  double taper;
  double dd;
  double w;
  size_t efs;
  size_t esf;
  if(nf < 1 || nf > FD3D_MULTI_MAX){
    free(s);
    return NULL;
  }
  s->nf = nf;
  for(k = 0; k < nf; k++){
    s->f[k] = f[k];
    s->fsg[k] = fsign[k];
    if(f[k]->t->nx != f[0]->t->nx){
      free(s);
      return NULL;
    }        /* one along-strike lattice */
    frame_fault(dip_deg[k], F);
    memset(&g[k], 0, sizeof(g[k]));
    g[k].dx = f[k]->t->dx;
    g[k].comp = fcomp;
    for(q = 0; q < 3; q++){
      g[k].E[0][q] = F[0][q];
      g[k].E[1][q] = F[1][q];
      g[k].E[2][q] = F[2][q];
    }
    g[k].O[0] = x0f;
    g[k].O[1] = ytop[k];
    g[k].O[2] = ztop[k];       /* cell (0, 0): top row, first column */
  }
  if(nf > 1){
    s->mx = fd3d_multi_create(nf, f, g);
    if(!s->mx){
      free(s);
      return NULL;
    }
  }
  for(k = 0; k < nf; k++){
    g_surf_skip_T2 = (k > 0);
    s->ly[k] = fd3d_surf_create_rows(f[k], x0f, ytop[k], ztop[k], dip_deg[k], fcomp, m, nxs, nys, xs0, yc, wy, L, fsign[k]);
    g_surf_skip_T2 = 0;
    if(!s->ly[k] || s->ly[k]->Mx != s->ly[0]->Mx){
      fd3d_msurf_free(s);
      return NULL;
    }
    f[k]->layer = s->ly[k];
    if(getenv("FD_MSURF_CUT")){
      /* experiment (fd_plan.md H2, "local exact, far approximate"): fault k does not interact
         with surface rows farther than the cutoff from its trace (horizontal distance): their
         fault-to-surface tables (all parts) and surface-to-fault tables (dynamic parts) are
         zeroed; the statics stay exact through Is */
      ly = s->ly[k];
      cut = atof(getenv("FD_MSURF_CUT"));
      nys_2 = ly->nys;
      nz = ly->nz;
      nkk = ly->nk;
      Lk = ly->L;
      ncut = 0;
      nfs = (size_t)nys_2*nz;
      taper = getenv("FD_MSURF_TAPER") ? atof(getenv("FD_MSURF_TAPER")) : 0.0;
      for(q2 = 0; q2 < nys_2; q2++){
        /* weight 1 within the cutoff, a cosine taper over the next `taper` metres, 0 beyond */
        dd = fabs(ly->yc[q2] - ytop[k]);
        w = (dd <= cut) ? 1.0 : ((taper > 0.0 && dd < cut + taper) ? 0.5*(1.0 + cos(M_PI*(dd - cut)/taper)) : 0.0);
        if(w >= 1.0)
          continue;
        if(w == 0.0)
          ncut++;
        for(j = 0; j < nz; j++)
          for(kk = 0; kk < nkk; kk++){
            efs = (size_t)kk*nfs + (size_t)q2*nz + j;
            esf = (size_t)kk*nfs + (size_t)j*nys_2 + q2;
            for(r = 0; r < 3; r++){
              for(l = 0; l < Lk; l++)
                ly->Tfs[r][0][(size_t)l*nkk*nfs + efs] *= w;
              ly->Tfs[r][1][efs] *= w;
              ly->Tfs[r][2][efs] *= w;
            }
            for(r = 0; r < 2; r++)
              for(cc = 0; cc < 3; cc++){
                for(l = 0; l < Lk; l++)
                  ly->Tsf[r][cc][0][(size_t)l*nkk*nfs + esf] *= w;
                ly->Tsf[r][cc][2][esf] *= w;
              }
          }
      }
      fprintf(stderr, "fd3d_msurf: fault %d sees %d of %d surface rows (cutoff %g m from its trace)\n", k, nys_2 - ncut, nys_2, cut);
    }
  }
  return s;
}

void fd3d_msurf_free(struct fd3d_msurf *s)
{
  if(!s)
    return;
  if(s->mx)
    fd3d_multi_free(s->mx);
  free(s);
}

int fd3d_msurf_step(struct fd3d_msurf *s, double **load)
{
  struct fd3d_layer *L0 = s->ly[0];
  const struct fd3d_tab *t = s->f[0]->t;
  int nf = s->nf, ns = L0->ns, nys = L0->nys, i, l, it, inner, c, r, k, q, j, itot = 0, lmax = L0->L < L0->nhist + 1 ? L0->L : L0->nhist + 1, nl = s->f[0]->nl;
  size_t Mh = L0->Mh, m, Q = Mh*nys;
  double dt = t->dt, hw = nl ? 0.5 : 1.0, Es[3];
  double *V3 = malloc(sizeof(double)*3*ns), *rhs = malloc(sizeof(double)*3*ns), *tmps = malloc(sizeof(double)*3*ns);
  double _Complex *Xs = malloc(sizeof(double _Complex)*Mh*nys), *Ys = malloc(sizeof(double _Complex)*3*Mh*nys);
  double _Complex *S3 = malloc(sizeof(double _Complex)*3*Q), *P3 = malloc(sizeof(double _Complex)*3*Q);
  double *Hx[FD3D_MULTI_MAX], *Hn[FD3D_MULTI_MAX], *Kx[FD3D_MULTI_MAX], *Kn[FD3D_MULTI_MAX], *Vp[FD3D_MULTI_MAX], *thb[FD3D_MULTI_MAX];
  double _Complex *Xf[FD3D_MULTI_MAX], *Yf[FD3D_MULTI_MAX];
  double * dl[8];
  int N;
  int nz;
  struct fd3d_layer *ly;
  int nz_2;
  size_t nfs;
  int slot;
  _Complex double *Yt;
  _Complex double *Yt2;
  struct fd3d_mpair *cp;
  double sg;
  double sn;
  int nr;
  int nsr;
  double a;
  double b;
  const double *v;
  const double *Ws;
  const double *Wn;
  double a_2;
  double b_2;
  double a_3;
  double b_3;
  double ch;
  double vm;
  struct fd3d *f;
  struct fd3d_layer *ly_2;
  int N_2;
  int nz_3;
  double *Hf;
  double *tnc;
  struct fd3d_mpair *cp_2;
  double a_4;
  double b_4;
  double * dend[8];
  double *o;
  int nq;
  struct fd3d_layer *ly_3;
  double *fk;
  int slot_2;
  double A[5][5];
  double rhsg[5];
  double gam[5];
  double sdot;
  double sdot_2;
  double fct;
  int e;
  double sdot_3;
  int slot_3;
  double *vms;
  struct fd3d *f_2;
  struct fd3d_layer *ly_4;
  int N_3;
  int nz_4;
  double *tmpf;
  double *vmf;
  struct fd3d_mpair *cp_3;
  double a_5;
  struct fd3d *f_3;
  struct fd3d_layer *ly_5;
  Es[0] = Es[1] = 0.5*t->rho*t->cT;
  Es[2] = 0.5*t->rho*t->cL;
  /* exact statics over all faults (shear), as fd3d_multi_step */
  if(s->mx && s->mx->static_fn){
    for(k = 0; k < nf; k++){
      if(!s->mx->sbuf[k])
        s->mx->sbuf[k] = malloc(sizeof(double)*s->f[k]->N);
      dl[k] = s->f[k]->delta;
    }
    s->mx->static_fn(s->mx->static_ctx, dl, s->mx->sbuf);
    for(k = 0; k < nf; k++){
      s->mx->slot[k].m = s->mx;
      s->mx->slot[k].p = k;
      s->f[k]->static_fn = multi_sbuf_out;
      s->f[k]->static_ctx = &s->mx->slot[k];
    }
  }
  for(k = 0; k < nf; k++){
    N = s->f[k]->N;
    nz = s->ly[k]->nz;
    step_begin(s->f[k], load ? load[k] : NULL);
    Hx[k] = calloc(N, sizeof(double));
    Hn[k] = calloc(N, sizeof(double));
    Kx[k] = calloc(N, sizeof(double));
    Kn[k] = calloc(N, sizeof(double));
    Vp[k] = malloc(sizeof(double)*N);
    thb[k] = (s->f[k]->law == 1) ? malloc(sizeof(double)*N) : NULL;
    Xf[k] = malloc(sizeof(double _Complex)*Mh*nz);
    Yf[k] = calloc(2*Mh*nz, sizeof(double _Complex));
  }
  memcpy(L0->Vold3, L0->V3, sizeof(double)*3*ns);
  /* ---- known parts: surface <- surface and all faults; fault k <- surface ---- */
  memset(Ys, 0, sizeof(double _Complex)*3*Mh*nys);
  for(c = 0; c < 3; c++)
    layer_rows_hat(L0, L0->delta3 + (size_t)c*ns, S3 + (size_t)c*Q);
  surf_apply_rows(L0, 1, 0, S3, 1.0, Ys);
  for(k = 0; k < nf; k++){
    ly = s->ly[k];
    nz_2 = ly->nz;
    nfs = (size_t)nys*nz_2;
    fault_rows_hat(ly, s->f[k]->delta, Xf[k]);
    for(r = 0; r < 3; r++)
      lxp_apply(ly, nys, nz_2, ly->Tfs[r][1], Xf[k], Ys + (size_t)r*Mh*nys);
    for(c = 0; c < 3; c++){
      layer_rows_hat(L0, L0->delta3 + (size_t)c*ns, Xs);
      if(!s->f[k]->static_fn)
        lxp_apply(ly, nz_2, nys, ly->Tsf[0][c][1], Xs, Yf[k]);
      if(!ly->static_n_fn)
        lxp_apply(ly, nz_2, nys, ly->Tsf[1][c][1], Xs, Yf[k] + Mh*nz_2);
    }
    for(l = 2; l <= lmax; l++){
      slot = (L0->hk - (l - 1) + 2*L0->L) % L0->L;
      if(k == 0)
        surf_apply_rows(L0, 0, l - 1, L0->hsx3 + (size_t)slot*3*Q, 1.0, Ys);
      for(r = 0; r < 3; r++)
        lxp_apply(ly, nys, nz_2, ly->Tfs[r][0] + (size_t)(l - 1)*ly->nk*nfs, ly->hfx + (size_t)slot*Mh*nz_2, Ys + (size_t)r*Mh*nys);
      for(c = 0; c < 3; c++)
        for(r = 0; r < 2; r++)
          lxp_apply(ly, nz_2, nys, ly->Tsf[r][c][0] + (size_t)(l - 1)*ly->nk*nfs, L0->hsx3 + ((size_t)slot*3 + c)*Mh*nys, Yf[k] + (size_t)r*Mh*nz_2);
    }
    if(nl){
      Yt = calloc(3*Mh*nys, sizeof(double _Complex));
      Yt2 = calloc(2*Mh*nz_2, sizeof(double _Complex));
      fault_rows_hat(ly, s->f[k]->Vold, Xf[k]);
      for(r = 0; r < 3; r++)
        lxp_apply(ly, nys, nz_2, ly->Tfs[r][2], Xf[k], Yt + (size_t)r*Mh*nys);
      for(m = 0; m < 3*Mh*nys; m++)
        Ys[m] += 0.5*Yt[m];
      for(c = 0; c < 3; c++){
        layer_rows_hat(L0, L0->Vold3 + (size_t)c*ns, Xs);
        for(r = 0; r < 2; r++)
          lxp_apply(ly, nz_2, nys, ly->Tsf[r][c][2], Xs, Yt2 + (size_t)r*Mh*nz_2);
      }
      for(m = 0; m < 2*Mh*nz_2; m++)
        Yf[k][m] += 0.5*Yt2[m];
      free(Yt);
      free(Yt2);
    }
    surf_reduce(ly, Yf[k], 2*Mh*nz_2);
    lx_inv(ly, nz_2, Yf[k], ly->nx, ly->m, 1.0, Hx[k]);
    lx_inv(ly, nz_2, Yf[k] + Mh*nz_2, ly->nx, ly->m, 1.0, Hn[k]);
  }
  if(nl){
    for(c = 0; c < 3; c++)
      layer_rows_hat(L0, L0->Vold3 + (size_t)c*ns, S3 + (size_t)c*Q);
    surf_apply_rows(L0, 2, 0, S3, 0.5, Ys);
  }
  surf_reduce(L0, Ys, 3*Mh*nys);
  memset(L0->H3, 0, sizeof(double)*3*ns);
  for(c = 0; c < 3; c++){
    lx_inv(L0, nys, Ys + (size_t)c*Mh*nys, L0->nxs, 1, 1.0, L0->H3 + (size_t)c*ns);
    if(nl)
      for(i = 0; i < ns; i++)
        L0->H3[(size_t)c*ns + i] += 0.5*Es[c]*L0->Vold3[(size_t)c*ns + i];
  }
  /* fault <- fault (dense whole-space tables, slip signs fsg): statics (unless exact), lags >= 2, NL half */
  if(s->mx)
    for(k = 0; k < nf; k++)
      for(q = 0; q < nf; q++){
        cp = &s->mx->c[q][k];
        sg = s->fsg[k]*s->fsg[q];
        sn = s->fsg[q];
        nr = s->f[k]->N;
        nsr = cp->ns;
        if(k == q)
          continue;
        if(!s->mx->static_fn)
          for(i = 0; i < nr; i++){
            a = 0.0;
            b = 0.0;
            for(j = 0; j < nsr; j++){
              a += cp->Ks[(size_t)i*nsr + j]*s->f[q]->delta[j];
              b += cp->Kn[(size_t)i*nsr + j]*s->f[q]->delta[j];
            }
            Kx[k][i] += sg*a;
            Kn[k][i] += sn*b;
          }
        for(l = 2; l <= cp->L && l - 1 <= s->mx->hk; l++){
          v = s->mx->ring[q] + (size_t)((s->mx->hk - (l - 1) + (long)s->mx->Lb*1000) % s->mx->Lb)*nsr;
          Ws = cp->Ws + (size_t)(l - 1)*nr*nsr;
          Wn = cp->Wn + (size_t)(l - 1)*nr*nsr;
          for(i = 0; i < nr; i++){
            a_2 = 0.0;
            b_2 = 0.0;
            for(j = 0; j < nsr; j++){
              a_2 += Ws[(size_t)i*nsr + j]*v[j];
              b_2 += Wn[(size_t)i*nsr + j]*v[j];
            }
            Kx[k][i] += sg*a_2;
            Kn[k][i] += sn*b_2;
          }
        }
        if(nl && cp->has_k1)
          for(i = 0; i < nr; i++){
            a_3 = 0.0;
            b_3 = 0.0;
            for(j = 0; j < nsr; j++){
              a_3 += cp->K1s[(size_t)i*nsr + j]*s->f[q]->Vold[j];
              b_3 += cp->K1n[(size_t)i*nsr + j]*s->f[q]->Vold[j];
            }
            Kx[k][i] += 0.5*sg*a_3;
            Kn[k][i] += 0.5*sn*b_3;
          }
      }
  /* ---- iterate faults and surface ---- */
  memcpy(V3, L0->Vold3, sizeof(double)*3*ns);
  for(k = 0; k < nf; k++)
    s->f[k]->thnew = thb[k];
#define AA_M 5
  int aa_n = 0, aa_i;
  double *aa_dF[AA_M], *aa_dG[AA_M], *aa_f = malloc(sizeof(double)*3*L0->ns), *aa_g = malloc(sizeof(double)*3*L0->ns);
  for(aa_i = 0; aa_i < AA_M; aa_i++){
    aa_dF[aa_i] = calloc(3*L0->ns, sizeof(double));
    aa_dG[aa_i] = calloc(3*L0->ns, sizeof(double));
  }
  for(it = 0; it < 200; it++){
    ch = 0.0;
    vm = 1e-30;
    for(k = 0; k < nf; k++){
      f = s->f[k];
      ly_2 = s->ly[k];
      N_2 = f->N;
      nz_3 = ly_2->nz;
      Hf = malloc(sizeof(double)*N_2);
      tnc = calloc(N_2, sizeof(double));
      memset(Yf[k], 0, sizeof(double _Complex)*2*Mh*nz_3);
      for(c = 0; c < 3; c++){
        layer_rows_hat(L0, V3 + (size_t)c*ns, Xs);
        for(r = 0; r < 2; r++)
          lxp_apply(ly_2, nz_3, nys, ly_2->Tsf[r][c][2], Xs, Yf[k] + (size_t)r*Mh*nz_3);
      }
      surf_reduce(ly_2, Yf[k], 2*Mh*nz_3);
      for(i = 0; i < N_2; i++)
        Hf[i] = Hx[k][i] + Kx[k][i];
      lx_inv(ly_2, nz_3, Yf[k], ly_2->nx, ly_2->m, hw, Hf);
      memcpy(ly_2->tn, Hn[k], sizeof(double)*N_2);
      lx_inv(ly_2, nz_3, Yf[k] + Mh*nz_3, ly_2->nx, ly_2->m, hw, ly_2->tn);
      for(i = 0; i < N_2; i++)
        ly_2->tn[i] += Kn[k][i];
      if(s->mx)
        for(q = 0; q < nf; q++){                 /* the current step's fault <- fault term */
          cp_2 = &s->mx->c[q][k];
          if(q == k || !cp_2->has_k1)
            continue;
          for(i = 0; i < N_2; i++){
            a_4 = 0.0;
            b_4 = 0.0;
            for(j = 0; j < cp_2->ns; j++){
              a_4 += cp_2->K1s[(size_t)i*cp_2->ns + j]*s->f[q]->V[j];
              b_4 += cp_2->K1n[(size_t)i*cp_2->ns + j]*s->f[q]->V[j];
            }
            Hf[i] += hw*s->fsg[k]*s->fsg[q]*a_4;
            tnc[i] += hw*s->fsg[q]*b_4;
          }
        }
      for(i = 0; i < N_2; i++)
        ly_2->tn[i] += tnc[i];
      if(s->static_n_fn){                                /* exact static normal traction at the end of the step (all faults) */
        o = malloc(sizeof(double)*N_2);
        for(q = 0; q < nf; q++){
          nq = s->f[q]->N;
          dend[q] = malloc(sizeof(double)*nq);
          for(i = 0; i < nq; i++)
            dend[q][i] = s->f[q]->delta[i] + (nl ? 0.5*dt*(s->f[q]->Vold[i] + s->f[q]->V[i]) : dt*s->f[q]->V[i]);
        }
        s->static_n_fn(s->static_n_ctx, k, dend, o);
        for(i = 0; i < N_2; i++)
          ly_2->tn[i] += o[i];
        for(q = 0; q < nf; q++)
          free(dend[q]);
        free(o);
      }
      if(ly_2->normal_coupling)
        for(i = 0; i < N_2; i++){
          f->sigma[i] = ly_2->sig0[i] - ly_2->tn[i];
          if(f->sigma[i] < 0.0)
            f->sigma[i] = 0.0;
        }
      memcpy(Vp[k], f->V, sizeof(double)*N_2);
      itot += step_solve_g(f, Hf, it > 0 ? Vp[k] : NULL, NULL);
      free(Hf);
      free(tnc);
    }
    /* surface: zero traction, Jacobi with the self coefficients */
    memset(Ys, 0, sizeof(double _Complex)*3*Mh*nys);
    for(k = 0; k < nf; k++){
      ly_3 = s->ly[k];
      fault_rows_hat(ly_3, s->f[k]->V, Xf[k]);
      for(r = 0; r < 3; r++)
        lxp_apply(ly_3, nys, ly_3->nz, ly_3->Tfs[r][2], Xf[k], Ys + (size_t)r*Mh*nys);
    }
    surf_reduce(L0, Ys, 3*Mh*nys);
    for(i = 0; i < 3*ns; i++)
      rhs[i] = -L0->H3[i];
    for(c = 0; c < 3; c++)
      lx_inv(L0, nys, Ys + (size_t)c*Mh*nys, L0->nxs, 1, -hw, rhs + (size_t)c*ns);
    inner = surf_solve(L0, nl, hw, Es, rhs, V3, S3, P3, tmps);
    L0->n_inner += inner;
    if(!getenv("FD_SURF_NOAA")){
      /* Anderson acceleration of the fault-surface fixed point (fd_plan.md H4): x = the surface
         rates fed into the fault solves, g = the surface rates that come back; with the
         differences of the last aa_m residuals f = g - x, the next x is g - dG gamma where gamma
         minimises |f - dF gamma| (the fixed point is unchanged; FD_SURF_NOAA=1 turns it off) */
      int n3 = 3*ns, mm = aa_n < AA_M ? aa_n : AA_M, a, b2;
      fk = malloc(sizeof(double)*n3);
      for(i = 0; i < n3; i++)
        fk[i] = V3[i] - L0->V3[i];
      if(aa_n > 0){
        slot_2 = (aa_n - 1) % AA_M;                      /* differences against the previous iterate */
        for(i = 0; i < n3; i++){
          aa_dF[slot_2][i] = fk[i] - aa_f[i];
          aa_dG[slot_2][i] = V3[i] - aa_g[i];
        }
      }
      memcpy(aa_f, fk, sizeof(double)*n3);
      memcpy(aa_g, V3, sizeof(double)*n3);
      if(mm > 0){
        for(a = 0; a < mm; a++){
          for(b2 = 0; b2 < mm; b2++){
            sdot = 0.0;
            for(i = 0; i < n3; i++)
              sdot += aa_dF[a][i]*aa_dF[b2][i];
            A[a][b2] = sdot;
          }
          {
            sdot_2 = 0.0;
            for(i = 0; i < n3; i++)
              sdot_2 += aa_dF[a][i]*fk[i];
            rhsg[a] = sdot_2;
          }
          A[a][a] *= 1.0 + 1e-10;
        }
        /* small symmetric system by Gaussian elimination */
        for(a = 0; a < mm; a++){
          for(b2 = a + 1; b2 < mm; b2++){
            fct = A[a][a] != 0.0 ? A[b2][a]/A[a][a] : 0.0;
            for(e = a; e < mm; e++)
              A[b2][e] -= fct*A[a][e];
            rhsg[b2] -= fct*rhsg[a];
          }
        }
        for(a = mm - 1; a >= 0; a--){
          sdot_3 = rhsg[a];
          for(b2 = a + 1; b2 < mm; b2++)
            sdot_3 -= A[a][b2]*gam[b2];
          gam[a] = A[a][a] != 0.0 ? sdot_3/A[a][a] : 0.0;
        }
        for(a = 0; a < mm; a++)
          if(isfinite(gam[a]))
            for(i = 0; i < n3; i++)
              V3[i] -= gam[a]*aa_dG[a][i];
      }
      aa_n++;
      free(fk);
    }
    for(i = 0; i < 3*ns; i++){
      if(fabs(V3[i] - L0->V3[i]) > ch)
        ch = fabs(V3[i] - L0->V3[i]);
      if(fabs(V3[i]) > vm)
        vm = fabs(V3[i]);
    }
    memcpy(L0->V3, V3, sizeof(double)*3*ns);
    for(k = 0; k < nf; k++)
      for(i = 0; i < s->f[k]->N; i++){
        if(fabs(s->f[k]->V[i] - Vp[k][i]) > ch)
          ch = fabs(s->f[k]->V[i] - Vp[k][i]);
        if(fabs(s->f[k]->V[i]) > vm)
          vm = fabs(s->f[k]->V[i]);
      }
    if(it > 0 && ch <= 1e-12*vm)
      break;
  }
  for(aa_i = 0; aa_i < AA_M; aa_i++){
    free(aa_dF[aa_i]);
    free(aa_dG[aa_i]);
  }
  free(aa_f);
  free(aa_g);
#undef AA_M
  L0->n_outer += it + 1;
  L0->n_steps++;
  if(getenv("FD_DEBUG_ITERS") && (L0->n_steps % 20 == 0))
    fprintf(stderr, "fd3d_msurf_step %i: outer iterations %d (mean %.1f), inner mean %.1f\n",
        (int)L0->n_steps, it + 1, (double)L0->n_outer/L0->n_steps, (double)L0->n_inner/L0->n_steps);
  s->last_iters = it + 1;
  for(k = 0; k < nf; k++){
    s->f[k]->thnew = NULL;
    if(thb[k]){
      memcpy(s->f[k]->theta, thb[k], sizeof(double)*s->f[k]->N);
      free(thb[k]);
    }
  }
  /* ---- finish: faults' known parts, rings, step_end ---- */
  {
    slot_3 = L0->hk % L0->L;
    vms = malloc(sizeof(double)*3*ns);
    for(i = 0; i < 3*ns; i++){
      vms[i] = nl ? 0.5*(L0->Vold3[i] + L0->V3[i]) : L0->V3[i];
      L0->delta3[i] += vms[i]*dt;
    }
    for(c = 0; c < 3; c++)
      layer_rows_hat(L0, vms + (size_t)c*ns, L0->hsx3 + ((size_t)slot_3*3 + c)*Mh*nys);
    free(vms);
    for(k = 0; k < nf; k++){
      f_2 = s->f[k];
      ly_4 = s->ly[k];
      N_3 = f_2->N;
      nz_4 = ly_4->nz;
      tmpf = calloc(N_3, sizeof(double));
      vmf = malloc(sizeof(double)*N_3);
      memset(Yf[k], 0, sizeof(double _Complex)*2*Mh*nz_4);
      for(c = 0; c < 3; c++){
        layer_rows_hat(L0, L0->V3 + (size_t)c*ns, Xs);
        for(r = 0; r < 2; r++)
          lxp_apply(ly_4, nz_4, nys, ly_4->Tsf[r][c][2], Xs, Yf[k] + (size_t)r*Mh*nz_4);
      }
      surf_reduce(ly_4, Yf[k], 2*Mh*nz_4);
      lx_inv(ly_4, nz_4, Yf[k], ly_4->nx, ly_4->m, hw, tmpf);
      for(i = 0; i < N_3; i++)
        f_2->H[i] += Hx[k][i] + Kx[k][i] + tmpf[i];
      if(s->mx)
        for(q = 0; q < nf; q++){                 /* the current step's fault <- fault term, final rates */
          cp_3 = &s->mx->c[q][k];
          if(q == k || !cp_3->has_k1)
            continue;
          for(i = 0; i < N_3; i++){
            a_5 = 0.0;
            for(j = 0; j < cp_3->ns; j++)
              a_5 += cp_3->K1s[(size_t)i*cp_3->ns + j]*s->f[q]->V[j];
            f_2->H[i] += hw*s->fsg[k]*s->fsg[q]*a_5;
          }
        }
      for(i = 0; i < N_3; i++)
        vmf[i] = nl ? 0.5*(f_2->Vold[i] + f_2->V[i]) : f_2->V[i];
      fault_rows_hat(ly_4, vmf, ly_4->hfx + (size_t)slot_3*Mh*nz_4);
      if(s->mx)
        memcpy(s->mx->ring[k] + (size_t)(s->mx->hk % s->mx->Lb)*N_3, vmf, sizeof(double)*N_3);
      free(tmpf);
      free(vmf);
    }
    for(k = 0; k < nf; k++){
      s->ly[k]->hk = (L0->hk + 1) % L0->L;
      if(s->ly[k]->nhist < L0->L)
        s->ly[k]->nhist = L0->nhist + (L0->nhist < L0->L ? 1 : 0);
    }
    if(s->mx)
      s->mx->hk++;
  }
  for(k = 0; k < nf; k++){
    f_3 = s->f[k];
    ly_5 = s->ly[k];
    step_end(f_3);
    if(ly_5->normal_coupling){
      /* normal traction at the end of the step: tn as last computed in the iteration */
      for(i = 0; i < f_3->N; i++){
        f_3->sigma[i] = ly_5->sig0[i] - ly_5->tn[i];
        if(f_3->sigma[i] < 0.0)
          f_3->sigma[i] = 0.0;
      }
    }
    free(Hx[k]);
    free(Hn[k]);
    free(Kx[k]);
    free(Kn[k]);
    free(Vp[k]);
    free(Xf[k]);
    free(Yf[k]);
  }
  free(V3);
  free(rhs);
  free(tmps);
  free(Xs);
  free(Ys);
  free(S3);
  free(P3);
  return itot;
}

/* ---------------- H6: per-fault surface strips of their own strike (fd_plan.md H6) ---------------- */

/* traction t = G n on a receiver with normal n, projected on nr directions: out[r] = dirs[r] . t */
static void project_traction(double G[3][3], const double n[3], int nr, double dirs[][3], double *out)
{
  double tv[3];
  int u, v, r;
  for(u = 0; u < 3; u++){
    tv[u] = 0.0;
    for(v = 0; v < 3; v++)
      tv[u] += G[u][v]*n[v];
  }
  for(r = 0; r < nr; r++)
    out[r] = dirs[r][0]*tv[0] + dirs[r][1]*tv[1] + dirs[r][2]*tv[2];
}

/* out = R in for a 3 x 3 matrix R */
static void rot_vec(double R[3][3], const double in[3], double out[3])
{
  int u, v;
  for(u = 0; u < 3; u++){
    out[u] = 0.0;
    for(v = 0; v < 3; v++)
      out[u] += R[u][v]*in[v];
  }
}

/* Generic response of a receiver to a source rectangle, both in global coordinates.
   source: centre c, frame R (rows: the rectangle's axes, R[2] its normal; slip comp 1, 2 along
   R[0], R[1], 3 opening), half sizes a1 (along R[0]) and a2 (along R[1]);
   receiver at X with normal n and nr projection directions dirs[r];
   output: K[r L + l - 1] the step increments (traction at l dt minus at (l - 1) dt, unit slip
   rate from t = 0) and ks[r] the static slope. */
void fd3d_elem_resp(const struct fd3d_tab *t, const double c[3], double R[3][3], int comp, double a1, double a2,
    const double X[3], const double n[3], int nr, double dirs[][3], int L, double *K, double *ks)
{
  double d[3], x[3], s[6], G[3][3], late[2][3], prev[3], tr[3];
  double T1, T2, Tq;
  int i, l, q, r;

  for(i = 0; i < 3; i++)
    d[i] = X[i] - c[i];
  for(i = 0; i < 3; i++){
    x[i] = R[i][0]*d[0] + R[i][1]*d[1] + R[i][2]*d[2];
    if(fabs(x[i]) < 1e-9*(a1 + a2))
      x[i] = 0.0;
  }
  /* static slope from two late times, after all waves have passed */
  T1 = 3.0*(sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) + 4.0*(a1 + a2))/t->cT + 10.0*t->dt;
  T2 = 2.0*T1;
  for(q = 0; q < 2; q++){
    Tq = q ? T2 : T1;
    fd3d_rect_ramp_tensor(x[0], x[1], x[2], Tq, a1, a2, t->cL, t->cT, t->rho, comp, s);
    to_global(R, s, G);
    project_traction(G, n, nr, dirs, late[q]);
  }
  for(r = 0; r < nr; r++)
    ks[r] = (late[1][r] - late[0][r])/(T2 - T1);
  /* step increments */
  for(r = 0; r < 3; r++)
    prev[r] = 0.0;
  for(l = 1; l <= L; l++){
    fd3d_rect_ramp_tensor(x[0], x[1], x[2], l*t->dt, a1, a2, t->cL, t->cT, t->rho, comp, s);
    to_global(R, s, G);
    project_traction(G, n, nr, dirs, tr);
    for(r = 0; r < nr; r++){
      K[(size_t)r*L + l - 1] = tr[r] - prev[r];
      prev[r] = tr[r];
    }
  }
}

/* largest |K2 - K1| and |K1| over nr rows of L lags, accumulated into worst and scale */
static void compare_resp(const double *K1, const double *K2, int nr, int L, double *worst, double *scale)
{
  int r, l;
  size_t e;
  for(r = 0; r < nr; r++)
    for(l = 0; l < L; l++){
      e = (size_t)r*L + l;
      *worst = fmax(*worst, fabs(K2[e] - K1[e]));
      *scale = fmax(*scale, fabs(K1[e]));
    }
}

/* frame check (fd3d_test h6frames): the responses of the one-fault surface code (surf_resp, in
   the fault's local frame x along strike, y horizontal towards the up-dip side, z up) against
   fd3d_elem_resp with everything rotated about z by `strike` degrees and shifted, for fault ->
   surface, surface -> fault and surface -> surface pairs; returns the largest relative
   difference. surf_resp applies no slip sign to a fault source (the tables do it afterwards). */
double fd3d_h6_frame_check(const struct fd3d_tab *t, double dip_deg, int fcomp, double fsign, double strike_deg, double hs)
{
  double F[3][3], sdir[3], Rz[3][3], RF[3][3], RS[3][3], dirs3[3][3], dirs2[2][3], nrm[3];
  double Ide[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, sh[3] = {1234.5, -987.6, 0.0}, nz[3] = {0, 0, 1};
  double fc[3], sc[3], sc2[3], gfc[3], gsc[3], gsc2[3], ks1[3], ks2[3];
  double worst, scale, a, th, *K1, *K2;
  int L, ip, r, c, u;

  L = 30;
  a = 0.5*t->dx;
  worst = scale = 0.0;
  K1 = malloc(sizeof(double)*3*L);
  K2 = malloc(sizeof(double)*3*L);
  /* local -> global: rotation about z by th, then the shift; frames: rows rotated */
  th = strike_deg*M_PI/180.0;
  Rz[0][0] = cos(th);
  Rz[0][1] = -sin(th);
  Rz[0][2] = 0;
  Rz[1][0] = sin(th);
  Rz[1][1] = cos(th);
  Rz[1][2] = 0;
  Rz[2][0] = 0;
  Rz[2][1] = 0;
  Rz[2][2] = 1;
  frame_fault(dip_deg, F);
  for(c = 0; c < 3; c++)
    sdir[c] = fsign*F[fcomp - 1][c];
  for(r = 0; r < 3; r++){
    rot_vec(Rz, F[r], RF[r]);
    rot_vec(Rz, Ide[r], RS[r]);
  }
  for(r = 0; r < 3; r++)
    for(u = 0; u < 3; u++)
      dirs3[r][u] = RS[r][u];
  rot_vec(Rz, sdir, dirs2[0]);
  rot_vec(Rz, F[2], dirs2[1]);
  rot_vec(Rz, F[2], nrm);
  for(ip = 0; ip < 3; ip++){
    /* a few source / receiver pairs in the local frame, and their global positions */
    fc[0] = 300.0*ip;
    fc[1] = -400.0 - 200.0*ip;
    fc[2] = -700.0 - 100.0*ip;
    sc[0] = -250.0 + 500.0*ip;
    sc[1] = 600.0 - 300.0*ip;
    sc[2] = 0.0;
    sc2[0] = sc[0] + 1500.0;
    sc2[1] = sc[1] - 750.0;
    sc2[2] = 0.0;
    for(u = 0; u < 3; u++){
      gfc[u] = sh[u];
      gsc[u] = sh[u];
      gsc2[u] = sh[u];
      for(r = 0; r < 3; r++){
        gfc[u] += Rz[u][r]*fc[r];
        gsc[u] += Rz[u][r]*sc[r];
        gsc2[u] += Rz[u][r]*sc2[r];
      }
    }
    /* fault -> surface (traction sxz, syz, szz in the local frame = projections on the rotated axes) */
    surf_resp(t, F, 1, fcomp, fc, a, a, sc, 0, sdir, L, K1, ks1);
    fd3d_elem_resp(t, gfc, RF, fcomp, a, a, gsc, nz, 3, dirs3, L, K2, ks2);
    compare_resp(K1, K2, 3, L, &worst, &scale);
    /* surface -> fault (shear along sdir, normal on F[2]) and surface -> surface, each component */
    for(c = 1; c <= 3; c++){
      surf_resp(t, F, 0, c, sc, 0.5*hs, 0.5*hs, fc, 1, sdir, L, K1, ks1);
      fd3d_elem_resp(t, gsc, RS, c, 0.5*hs, 0.5*hs, gfc, nrm, 2, dirs2, L, K2, ks2);
      compare_resp(K1, K2, 2, L, &worst, &scale);
      surf_resp(t, F, 0, c, sc, 0.5*hs, 0.5*hs, sc2, 0, sdir, L, K1, ks1);
      fd3d_elem_resp(t, gsc, RS, c, 0.5*hs, 0.5*hs, gsc2, nz, 3, dirs3, L, K2, ks2);
      compare_resp(K1, K2, 3, L, &worst, &scale);
    }
  }
  free(K1);
  free(K2);
  return scale > 0.0 ? worst/scale : -1.0;
}

/* strip k's local frame to global: position and vector */
static void strip_l2g(const struct fd3d_strips *s, int k, const double xl[3], double X[3])
{
  double ct, st;
  ct = cos(s->th[k]);
  st = sin(s->th[k]);
  X[0] = s->O[k][0] + xl[0]*ct - xl[1]*st;
  X[1] = s->O[k][1] + xl[0]*st + xl[1]*ct;
  X[2] = s->O[k][2] + xl[2];
}
static void strip_v2g(const struct fd3d_strips *s, int k, const double vl[3], double v[3])
{
  double ct, st;
  ct = cos(s->th[k]);
  st = sin(s->th[k]);
  v[0] = vl[0]*ct - vl[1]*st;
  v[1] = vl[0]*st + vl[1]*ct;
  v[2] = vl[2];
}

/* horizontal distance from p to the segment a-b */
static double seg2_dist(const double p[2], const double a[2], const double b[2])
{
  double ab[2], ap[2], l2, u;
  ab[0] = b[0] - a[0];
  ab[1] = b[1] - a[1];
  ap[0] = p[0] - a[0];
  ap[1] = p[1] - a[1];
  l2 = ab[0]*ab[0] + ab[1]*ab[1];
  u = l2 > 0.0 ? (ap[0]*ab[0] + ap[1]*ab[1])/l2 : 0.0;
  if(u < 0.0)
    u = 0.0;
  if(u > 1.0)
    u = 1.0;
  return hypot(p[0] - a[0] - u*ab[0], p[1] - a[1] - u*ab[1]);
}

static void dtab_alloc(struct fd3d_dtab *d, int nr, int nc, int L)
{
  d->nr = nr;
  d->nc = nc;
  d->L = L;
  d->W = calloc((size_t)L*nr*nc, sizeof(double));
  d->Ks = calloc((size_t)nr*nc, sizeof(double));
  d->K1 = calloc((size_t)nr*nc, sizeof(double));
}
static void dtab_free(struct fd3d_dtab *d)
{
  free(d->W);
  free(d->Ks);
  free(d->K1);
}

/* one receiver row of a dense table from fd3d_elem_resp: lag increments K (L values) and static
   slope ks of that row, times sgn, into column col */
static void dtab_fill_row(struct fd3d_dtab *d, const struct fd3d_tab *t, int col, int row, const double *K, double ks, double sgn)
{
  int l, L;
  size_t e;
  L = d->L;
  e = (size_t)row*d->nc + col;
  d->Ks[e] = sgn*ks;
  d->K1[e] = sgn*K[0];
  for(l = 0; l < L; l++)
    d->W[(size_t)l*d->nr*d->nc + e] = sgn*(K[l] - ks*t->dt);
}

/* fault <-> fault: the several-plane tables (G1) in global frames; cell (0, 0) at the top row,
   first column: local (x0f, ytop, ztop) */
static int strips_fault_tables(struct fd3d_strips *s, struct fd3d **f, const double *x0f, const double *ytop,
    const double *ztop, const double *dip_deg, int fcomp)
{
  struct fd3d_mplane g[FD3D_MULTI_MAX];
  double F[3][3], o[3];
  int k, r;
  for(k = 0; k < s->nf; k++){
    o[0] = x0f[k];
    o[1] = ytop[k];
    o[2] = ztop[k];
    frame_fault(dip_deg[k], F);
    memset(&g[k], 0, sizeof(g[k]));
    g[k].dx = f[k]->t->dx;
    g[k].comp = fcomp;
    for(r = 0; r < 3; r++)
      strip_v2g(s, k, F[r], g[k].E[r]);
    strip_l2g(s, k, o, g[k].O);
  }
  s->mx = fd3d_multi_create(s->nf, f, g);
  return s->mx ? 0 : 1;
}

/* global centres: strip cells (p, q) at (xs0 + (p + 1/2) hx, yc[q], 0); fault cells (i, j) at
   (x0f + (i + 1/2) dx, ytop - (j + 1/2) dx cos(dip), ztop - (j + 1/2) dx sin(dip)); trace ends */
static void strips_geometry(struct fd3d_strips *s, const double *x0f, const double *ytop, const double *ztop,
    const double *dip_deg, const double *xs0)
{
  struct fd3d_layer *ly;
  const struct fd3d_tab *t;
  double F[3][3], xl[3];
  int k, p, q, i, j;
  for(k = 0; k < s->nf; k++){
    ly = s->ly[k];
    t = s->f[k]->t;
    frame_fault(dip_deg[k], F);
    s->xg[k] = malloc(sizeof(double)*3*ly->ns);
    s->fg[k] = malloc(sizeof(double)*3*s->f[k]->N);
    for(p = 0; p < ly->nxs; p++)
      for(q = 0; q < ly->nys; q++){
        xl[0] = xs0[k] + (p + 0.5)*ly->hx;
        xl[1] = ly->yc[q];
        xl[2] = 0.0;
        strip_l2g(s, k, xl, s->xg[k] + 3*((size_t)p*ly->nys + q));
      }
    for(i = 0; i < t->nx; i++)
      for(j = 0; j < s->f[k]->nzp; j++){
        xl[0] = x0f[k] + (i + 0.5)*t->dx;
        xl[1] = ytop[k] - (j + 0.5)*t->dx*F[1][1];
        xl[2] = ztop[k] - (j + 0.5)*t->dx*F[1][2];
        strip_l2g(s, k, xl, s->fg[k] + 3*((size_t)i*s->f[k]->nzp + j));
      }
    /* the trace: the fault's top line in the surface plane */
    xl[0] = x0f[k];
    xl[1] = ytop[k];
    xl[2] = 0.0;
    strip_l2g(s, k, xl, s->ta[k]);
    xl[0] = x0f[k] + t->nx*t->dx;
    strip_l2g(s, k, xl, s->tb[k]);
  }
}

/* masks: a strip cell is active where its own trace is the nearest (ties to the lower index) */
static void strips_masks(struct fd3d_strips *s)
{
  double *X, dk, dq;
  int k, q, i;
  for(k = 0; k < s->nf; k++){
    s->active[k] = malloc(s->ly[k]->ns);
    for(i = 0; i < s->ly[k]->ns; i++){
      X = s->xg[k] + 3*(size_t)i;
      dk = seg2_dist(X, s->ta[k], s->tb[k]);
      s->active[k][i] = 1;
      for(q = 0; q < s->nf; q++){
        if(q == k)
          continue;
        dq = seg2_dist(X, s->ta[q], s->tb[q]);
        if(dq < dk || (dq == dk && q < k))
          s->active[k][i] = 0;
      }
    }
  }
}

/* do two squares overlap? centres d apart (from the first to the second), half sizes
   (hxi, hyi) along axes ax[0], ax[1] and (hxj, hyj) along ax[2], ax[3]; separating axes */
static int squares_overlap(const double d[2], double ax[4][2], double hxi, double hyi, double hxj, double hyj)
{
  double pd, ri, rj;
  int u;
  for(u = 0; u < 4; u++){
    pd = fabs(d[0]*ax[u][0] + d[1]*ax[u][1]);
    ri = hxi*fabs(ax[0][0]*ax[u][0] + ax[0][1]*ax[u][1]) + hyi*fabs(ax[1][0]*ax[u][0] + ax[1][1]*ax[u][1]);
    rj = hxj*fabs(ax[2][0]*ax[u][0] + ax[2][1]*ax[u][1]) + hyj*fabs(ax[3][0]*ax[u][0] + ax[3][1]*ax[u][1]);
    if(pd > ri + rj)
      return 0;
  }
  return 1;
}

/* trim (fd_plan.md H6 step 4): where turned cells of different strips overlap (their centres on
   either side of the partition boundary), the one whose centre lies closer to the boundary is
   masked, so that no part of the surface is represented twice (overlaps made the turned cases
   unstable); the gaps this leaves along an oblique boundary are smaller than a cell. Squares
   shrunk by 1e-6 to ignore touching edges. FD_STRIPS_NOTRIM=1 turns this off. */
static void strips_trim(struct fd3d_strips *s)
{
  struct fd3d_layer *lk, *lq;
  double *marg[FD3D_MULTI_MAX], ax[4][2], d[2], *X, *Xi, *Xj;
  double dk, dq, dmin, ck, sk, cq, sq, hxi, hyi, hxj, hyj;
  int changed, k, q, i, j;
  if(getenv("FD_STRIPS_NOTRIM"))
    return;
  /* margin of each cell: distance to the nearest other trace minus to its own */
  for(k = 0; k < s->nf; k++){
    marg[k] = malloc(sizeof(double)*s->ly[k]->ns);
    for(i = 0; i < s->ly[k]->ns; i++){
      X = s->xg[k] + 3*(size_t)i;
      dk = 0.0;
      dmin = 1e300;
      for(q = 0; q < s->nf; q++){
        dq = seg2_dist(X, s->ta[q], s->tb[q]);
        if(q == k)
          dk = dq;
        else if(dq < dmin)
          dmin = dq;
      }
      marg[k][i] = dmin - dk;
    }
  }
  changed = 1;
  while(changed){
    changed = 0;
    for(k = 0; k < s->nf; k++)
      for(q = k + 1; q < s->nf; q++){
        lk = s->ly[k];
        lq = s->ly[q];
        ck = cos(s->th[k]);
        sk = sin(s->th[k]);
        cq = cos(s->th[q]);
        sq = sin(s->th[q]);
        ax[0][0] = ck;
        ax[0][1] = sk;
        ax[1][0] = -sk;
        ax[1][1] = ck;
        ax[2][0] = cq;
        ax[2][1] = sq;
        ax[3][0] = -sq;
        ax[3][1] = cq;
        for(i = 0; i < lk->ns; i++){
          if(!s->active[k][i])
            continue;
          Xi = s->xg[k] + 3*(size_t)i;
          hxi = 0.5*lk->hx*(1.0 - 1e-6);
          hyi = 0.5*lk->wy[i % lk->nys]*(1.0 - 1e-6);
          for(j = 0; j < lq->ns && s->active[k][i]; j++){
            if(!s->active[q][j])
              continue;
            Xj = s->xg[q] + 3*(size_t)j;
            hxj = 0.5*lq->hx*(1.0 - 1e-6);
            hyj = 0.5*lq->wy[j % lq->nys]*(1.0 - 1e-6);
            d[0] = Xj[0] - Xi[0];
            d[1] = Xj[1] - Xi[1];
            if(hypot(d[0], d[1]) > (hxi + hyi + hxj + hyj))
              continue;
            if(squares_overlap(d, ax, hxi, hyi, hxj, hyj)){
              if(marg[k][i] < marg[q][j])
                s->active[k][i] = 0;
              else
                s->active[q][j] = 0;
              changed = 1;
            }
          }
        }
      }
  }
  for(k = 0; k < s->nf; k++)
    free(marg[k]);
}

/* dense tables between fault q / strip q (sources) and strip k / fault k (receivers), q != k;
   strip components in their own frame (x, y, opening), the fault source with its slip sign */
static void strips_dense_pair(struct fd3d_strips *s, int q, int k, const double *dip_deg, const double *fsign, int fcomp)
{
  struct fd3d_layer *lk, *lq;
  const struct fd3d_tab *t;
  double Fq[3][3], Fk[3][3], Rq[3][3], Rk[3][3], RSq[3][3], dirs_s[3][3], dk2[2][3], e[3], sl[3], sdk[3], nk[3];
  double ez[3] = {0, 0, 1}, ks[3], *K;
  int Nq, Nk, nsk, nsq, L, r, i, j, c;
  lk = s->ly[k];
  lq = s->ly[q];
  t = s->f[0]->t;
  L = s->L;
  Nq = s->f[q]->N;
  Nk = s->f[k]->N;
  nsk = lk->ns;
  nsq = lq->ns;
  K = malloc(sizeof(double)*3*L);
  frame_fault(dip_deg[q], Fq);
  frame_fault(dip_deg[k], Fk);
  for(r = 0; r < 3; r++){
    e[0] = (r == 0);
    e[1] = (r == 1);
    e[2] = (r == 2);
    strip_v2g(s, q, Fq[r], Rq[r]);
    strip_v2g(s, q, e, RSq[r]);
    strip_v2g(s, k, e, Rk[r]);
  }
  for(r = 0; r < 3; r++)
    for(i = 0; i < 3; i++)
      dirs_s[r][i] = Rk[r][i];
  for(i = 0; i < 3; i++)
    sl[i] = fsign[k]*Fk[fcomp - 1][i];
  strip_v2g(s, k, sl, sdk);
  strip_v2g(s, k, Fk[2], nk);
  for(i = 0; i < 3; i++){
    dk2[0][i] = sdk[i];
    dk2[1][i] = nk[i];
  }
  /* fault q -> strip k (source slip: the marcher's variable, fsign along Fq[fcomp-1]) */
  dtab_alloc(&s->fs[q][k], 3*nsk, Nq, L);
  for(j = 0; j < Nq; j++)
    for(i = 0; i < nsk; i++){
      fd3d_elem_resp(t, s->fg[q] + 3*(size_t)j, Rq, fcomp, 0.5*t->dx, 0.5*t->dx, s->xg[k] + 3*(size_t)i, ez, 3, dirs_s, L, K, ks);
      for(r = 0; r < 3; r++)
        dtab_fill_row(&s->fs[q][k], t, j, r*nsk + i, K + (size_t)r*L, ks[r], fsign[q]);
    }
  /* strip q -> fault k (receiver: shear along fault k's slip direction, normal) */
  dtab_alloc(&s->sf[q][k], 2*Nk, 3*nsq, L);
  for(c = 1; c <= 3; c++)
    for(j = 0; j < nsq; j++)
      for(i = 0; i < Nk; i++){
        fd3d_elem_resp(t, s->xg[q] + 3*(size_t)j, RSq, c, 0.5*lq->hx, 0.5*lq->wy[j % lq->nys], s->fg[k] + 3*(size_t)i, nk, 2, dk2, L, K, ks);
        for(r = 0; r < 2; r++)
          dtab_fill_row(&s->sf[q][k], t, (c - 1)*nsq + j, r*Nk + i, K + (size_t)r*L, ks[r], 1.0);
      }
  /* strip q -> strip k */
  dtab_alloc(&s->ss[q][k], 3*nsk, 3*nsq, L);
  for(c = 1; c <= 3; c++)
    for(j = 0; j < nsq; j++)
      for(i = 0; i < nsk; i++){
        fd3d_elem_resp(t, s->xg[q] + 3*(size_t)j, RSq, c, 0.5*lq->hx, 0.5*lq->wy[j % lq->nys], s->xg[k] + 3*(size_t)i, ez, 3, dirs_s, L, K, ks);
        for(r = 0; r < 3; r++)
          dtab_fill_row(&s->ss[q][k], t, (c - 1)*nsq + j, r*nsk + i, K + (size_t)r*L, ks[r], 1.0);
      }
  free(K);
}

struct fd3d_strips *fd3d_strips_create(int nf, struct fd3d **f, const double *th_deg, const double (*O)[3],
    const double *x0f, const double *ytop, const double *ztop, const double *dip_deg,
    const double *fsign, int fcomp, int m, const int *nxs, const int *nys,
    const double *xs0, const double **yc, const double **wy, int L)
{
  struct fd3d_strips *s;
  const struct fd3d_tab *t;
  int k, q, i;

  s = calloc(1, sizeof(*s));
  t = f[0]->t;
  s->nf = nf;
  s->L = L;
  for(k = 0; k < nf; k++){
    s->f[k] = f[k];
    s->th[k] = th_deg[k]*M_PI/180.0;
    s->fsg[k] = fsign[k];
    for(i = 0; i < 3; i++)
      s->O[k][i] = O[k][i];
  }
  /* fault <-> fault (before the layers are attached: fd3d_multi_create needs bare marchers) */
  if(nf > 1 && strips_fault_tables(s, f, x0f, ytop, ztop, dip_deg, fcomp)){
    free(s);
    return NULL;
  }
  /* each fault's strip: the one-fault surface layer in its local frame */
  for(k = 0; k < nf; k++){
    if(f[k]->t->dt != t->dt || f[k]->mirror || f[k]->next){
      fd3d_strips_free(s);
      return NULL;
    }
    s->ly[k] = fd3d_surf_create_rows(f[k], x0f[k], ytop[k], ztop[k], dip_deg[k], fcomp, m, nxs[k], nys[k],
        xs0[k], yc[k], wy[k], L, fsign[k]);
    if(!s->ly[k]){
      fd3d_strips_free(s);
      return NULL;
    }
    f[k]->layer = s->ly[k];
  }
  strips_geometry(s, x0f, ytop, ztop, dip_deg, xs0);
  strips_masks(s);
  strips_trim(s);
  for(q = 0; q < nf; q++)
    for(k = 0; k < nf; k++)
      if(q != k)
        strips_dense_pair(s, q, k, dip_deg, fsign, fcomp);
  /* real-space history rings of the pushed slices */
  s->Lb = L + 1;
  s->hk = 0;
  for(k = 0; k < nf; k++){
    s->sring[k] = calloc((size_t)s->Lb*3*s->ly[k]->ns, sizeof(double));
    s->fring[k] = calloc((size_t)s->Lb*f[k]->N, sizeof(double));
  }
  return s;
}

void fd3d_strips_free(struct fd3d_strips *s)
{
  int k, q;
  if(!s)
    return;
  if(s->mx)
    fd3d_multi_free(s->mx);
  for(k = 0; k < s->nf; k++){
    if(s->ly[k]){
      if(s->f[k])
        s->f[k]->layer = NULL;
      fd3d_layer_free(s->ly[k]);
    }
    free(s->active[k]);
    free(s->xg[k]);
    free(s->fg[k]);
    free(s->sring[k]);
    free(s->fring[k]);
    for(q = 0; q < s->nf; q++){
      dtab_free(&s->fs[k][q]);
      dtab_free(&s->sf[k][q]);
      dtab_free(&s->ss[k][q]);
    }
  }
  free(s->lu);
  free(s->lpiv);
  free(s->aidx);
  free(s);
}

/* step-2 check (fd3d_test h6strips): two faults; (a) the active cells of all strips cover each
   global cell position once; (b) shared strike: the dense strip-to-strip and fault-to-strip
   entries equal surf_resp evaluated in the common local frame; (c) a rigid rotation of the
   whole configuration leaves every dense table unchanged. Returns 0 if all pass. */
int fd3d_h6_strips_check(int verbose)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = 1000.0, dt = 0.5*dx/cT;
  int nx = 4, nz = 2, m = 4, L = 12, k, i, fail = 0, nxs[2] = {24, 24}, nys[2] = {16, 16};
  double th[2] = {0.0, 0.0}, O[2][3] = {{0, 0, 0}, {0, 0, 0}}, x0f[2] = {0.0, 0.0}, ytop[2] = {0.0, -3000.0}, ztop[2] = {0.0, 0.0}, dip[2] = {60.0, 75.0}, fs[2] = {1.0, -1.0};
  double xs0[2] = {-2000.0, -2000.0}, *yc[2], *wy[2], hs = dx/m;
  struct fd3d_tab *tb = fd3d_tab_create_c(nx, nz, dx, dt, cL, cT, rho, 0, 2);
  struct fd3d *f[2], *g[2];
  struct fd3d_strips *s1, *s2;
  double worst_b = 0.0, scale_b = 0.0, worst_c = 0.0, scale_c = 0.0;
  int ns;
  int cnt;
  int dup;
  int a0;
  int a1;
  double d;
  double F1[3][3];
  double F0[3][3];
  double sd0[3];
  double *K;
  double ks[3];
  int j;
  int r;
  int l;
  int cc;
  int ns_2;
  double ours;
  double ours_2;
  double th2[2];
  double O2[2][3];
  struct fd3d_dtab * a[3];
  struct fd3d_dtab * b[3];
  int u;
  size_t e;
  size_t n;
  for(k = 0; k < 2; k++){
    yc[k] = malloc(sizeof(double)*nys[k]);
    wy[k] = malloc(sizeof(double)*nys[k]);
    for(i = 0; i < nys[k]; i++){
      yc[k][i] = -5000.0 + (i + 0.5)*hs + 1000.0*0;
      wy[k][i] = hs;
    }
    f[k] = fd3d_create(tb);
    g[k] = fd3d_create(tb);
  }
  /* both strips on the same lattice y in [-5000, -1000] in the shared frame */
  s1 = fd3d_strips_create(2, f, th, (const double (*)[3])O, x0f, ytop, ztop, dip, fs, 2, m, nxs, nys, xs0, (const double **)yc, (const double **)wy, L);
  if(!s1)
    return 1;
  /* (a) partition: every lattice position active in exactly one strip */
  {
    ns = s1->ly[0]->ns;
    cnt = 0;
    dup = 0;
    for(i = 0; i < ns; i++){
      a0 = s1->active[0][i];
      a1 = s1->active[1][i];
      d = fabs(s1->xg[0][3*i] - s1->xg[1][3*i]) + fabs(s1->xg[0][3*i + 1] - s1->xg[1][3*i + 1]);
      if(d > 1e-6)
        dup++;
      if(a0 + a1 == 1)
        cnt++;
    }
    if(cnt != ns || dup)
      fail = 1;
    if(verbose)
      printf("h6strips (a): %d of %d lattice cells active in exactly one strip, %d position mismatches\n", cnt, ns, dup);
  }
  /* (b) shared strike: entries against surf_resp in the local frame */
  {
    K = malloc(sizeof(double)*3*L);
    ns_2 = s1->ly[0]->ns;
    frame_fault(dip[1], F1);
    frame_fault(dip[0], F0);
    for(i = 0; i < 3; i++)
      sd0[i] = fs[0]*F0[1][i];
    for(j = 0; j < ns_2; j += 37)
      for(i = 3; i < ns_2; i += 41){
        for(cc = 1; cc <= 3; cc++){
          surf_resp(tb, F1, 0, cc, s1->xg[1] + 3*j, 0.5*hs, 0.5*hs, s1->xg[0] + 3*i, 0, sd0, L, K, ks);
          for(r = 0; r < 3; r++)
            for(l = 0; l < L; l++){
              ours = s1->ss[1][0].W[(size_t)l*s1->ss[1][0].nr*s1->ss[1][0].nc + (size_t)(r*ns_2 + i)*s1->ss[1][0].nc + (cc - 1)*ns_2 + j] + s1->ss[1][0].Ks[(size_t)(r*ns_2 + i)*s1->ss[1][0].nc + (cc - 1)*ns_2 + j]*dt;
              worst_b = fmax(worst_b, fabs(ours - K[(size_t)r*L + l]));
              scale_b = fmax(scale_b, fabs(K[(size_t)r*L + l]));
            }
        }
      }
    for(j = 0; j < f[1]->N; j++)
      for(i = 1; i < ns_2; i += 53){
        surf_resp(tb, F1, 1, 2, s1->fg[1] + 3*j, 0.5*dx, 0.5*dx, s1->xg[0] + 3*i, 0, sd0, L, K, ks);
        for(r = 0; r < 3; r++)
          for(l = 0; l < L; l++){
            ours_2 = s1->fs[1][0].W[(size_t)l*s1->fs[1][0].nr*s1->fs[1][0].nc + (size_t)(r*ns_2 + i)*s1->fs[1][0].nc + j] + s1->fs[1][0].Ks[(size_t)(r*ns_2 + i)*s1->fs[1][0].nc + j]*dt;
            worst_b = fmax(worst_b, fabs(ours_2 - fs[1]*K[(size_t)r*L + l]));
            scale_b = fmax(scale_b, fabs(K[(size_t)r*L + l]));
          }
      }
    free(K);
    if(!(worst_b <= 1e-10*scale_b))
      fail = 1;
    if(verbose)
      printf("h6strips (b): dense entries against surf_resp in the shared frame: %.2e of %.2e\n", worst_b, scale_b);
  }
  /* (c) rigid rotation by 37 degrees about z (and a shift) of everything */
  {
    th2[0] = 37.0;
    th2[1] = 37.0;
    O2[0][0] = 1500.0;
    O2[0][1] = -700.0;
    O2[0][2] = 0.0;
    O2[1][0] = 1500.0;
    O2[1][1] = -700.0;
    O2[1][2] = 0.0;
    s2 = fd3d_strips_create(2, g, th2, (const double (*)[3])O2, x0f, ytop, ztop, dip, fs, 2, m, nxs, nys, xs0, (const double **)yc, (const double **)wy, L);
    if(!s2)
      return 1;
    {
      a[0] = &s1->fs[1][0];
      a[1] = &s1->sf[1][0];
      a[2] = &s1->ss[0][1];
      b[0] = &s2->fs[1][0];
      b[1] = &s2->sf[1][0];
      b[2] = &s2->ss[0][1];

      for(u = 0; u < 3; u++){
        n = (size_t)a[u]->L*a[u]->nr*a[u]->nc;
        for(e = 0; e < n; e++){
          worst_c = fmax(worst_c, fabs(a[u]->W[e] - b[u]->W[e]));
          scale_c = fmax(scale_c, fabs(a[u]->W[e]));
        }
        n = (size_t)a[u]->nr*a[u]->nc;
        for(e = 0; e < n; e++){
          worst_c = fmax(worst_c, fabs(a[u]->Ks[e] - b[u]->Ks[e])*dt);
          worst_c = fmax(worst_c, fabs(a[u]->K1[e] - b[u]->K1[e]));
        }
      }
      for(k = 0; k < 2; k++)
        for(i = 0; i < s1->ly[k]->ns; i++)
          if(s1->active[k][i] != s2->active[k][i])
            worst_c = 1e300;
    }
    if(!(worst_c <= 1e-8*scale_c))
      fail = 1;     /* rotation round-off near the trace: 1.3e-9 measured */
    if(verbose)
      printf("h6strips (c): rigid rotation by 37 degrees: dense tables and masks change by %.2e of %.2e\n", worst_c, scale_c);
    fd3d_strips_free(s2);
  }
  fd3d_strips_free(s1);
  for(k = 0; k < 2; k++){
    fd3d_free(f[k]);
    fd3d_free(g[k]);
    free(yc[k]);
    free(wy[k]);
  }
  fd3d_tab_free(tb);
  return fail;
}

/* y[r] += w sum_c M[r][c] x[c] for an nr x nc row-major matrix */
static void dmv(const double *M, int nr, int nc, const double *x, double w, double *y)
{
  int r, c;
  double a;
  const double *row;
  for(r = 0; r < nr; r++){
    a = 0.0;
    row = M + (size_t)r*nc;
    for(c = 0; c < nc; c++)
      a += row[c]*x[c];
    y[r] += w*a;
  }
}
/* joint operator of all strips' zero-traction conditions on the stacked rates v (strip k at
   off[k], 3 ns_k each), masked cells held at zero: out = A v restricted to active rows */
struct strips_ctx {
  struct fd3d_strips *s;
  int nl;
  double hw, Es[3];
  int off[FD3D_MULTI_MAX + 1];
  double _Complex *S3[FD3D_MULTI_MAX], *P3[FD3D_MULTI_MAX];
  double *tmp;
};
static void strips_mask(struct strips_ctx *x, double *v)
{
  int k, c, i;
  int ns;
  for(k = 0; k < x->s->nf; k++){
    ns = x->s->ly[k]->ns;
    for(c = 0; c < 3; c++)
      for(i = 0; i < ns; i++)
        if(!x->s->active[k][i])
          v[x->off[k] + c*ns + i] = 0.0;
  }
}
static void strips_A(struct strips_ctx *x, const double *v, double *out)
{
  struct fd3d_strips *s = x->s;
  int k, j;
  double *vm = malloc(sizeof(double)*x->off[s->nf]);
  memcpy(vm, v, sizeof(double)*x->off[s->nf]);
  strips_mask(x, vm);
  for(k = 0; k < s->nf; k++){
    surf_A(s->ly[k], x->nl, x->hw, x->Es, vm + x->off[k], out + x->off[k], x->S3[k], x->P3[k]);
    for(j = 0; j < s->nf; j++)
      if(j != k)
        dmv(s->ss[j][k].K1, s->ss[j][k].nr, s->ss[j][k].nc, vm + x->off[j], x->hw, out + x->off[k]);
  }
  strips_mask(x, out);
  free(vm);
}
static void strips_P(struct strips_ctx *x, const double *r, double *z)
{
  struct fd3d_strips *s = x->s;
  int k;
  for(k = 0; k < s->nf; k++){
    surf_precond_setup(s->ly[k], x->nl, x->hw, x->Es);
    surf_precond(s->ly[k], r + x->off[k], z + x->off[k], x->S3[k]);
  }
  strips_mask(x, z);
}
/* direct solve of the joint strip system (fd_plan.md H6 step 4): the first-lag operator over
   the active cells does not change in time, so it is assembled once (columns from strips_A on
   unit vectors) and LU-factorised (dgetrf); each solve is then exact. The GMRES path below,
   preconditioned per strip, converged slowly where the partition cuts strips obliquely (355
   iterations per step at 30 degrees against 15 with a shared strike). FD_STRIPS_GMRES=1 keeps
   GMRES. Dense: (3 active cells)^2 doubles; fine for test sizes, not for production strips. */
extern void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern void dgetrs_(char *tr, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
static int strips_direct(struct strips_ctx *x, const double *rhs, double *V)
{
  struct fd3d_strips *s = x->s;
  int n = x->off[s->nf], i, j, k, c, info, one = 1;
  char tr = 'N';
  double *e;
  double *col;
  int na;
  int ns;
  int na_2;
  double *b;
  if(!s->lu){
    e = calloc(n, sizeof(double));
    col = malloc(sizeof(double)*n);
    na = 0;
    s->aidx = malloc(sizeof(int)*n);
    for(k = 0; k < s->nf; k++){
      ns = s->ly[k]->ns;
      for(c = 0; c < 3; c++)
        for(i = 0; i < ns; i++)
          if(s->active[k][i])
            s->aidx[na++] = x->off[k] + c*ns + i;
    }
    s->na = na;
    s->lu = malloc(sizeof(double)*(size_t)na*na);
    s->lpiv = malloc(sizeof(int)*na);
    for(j = 0; j < na; j++){
      e[s->aidx[j]] = 1.0;
      strips_A(x, e, col);
      for(i = 0; i < na; i++)
        s->lu[(size_t)j*na + i] = col[s->aidx[i]];          /* column major */
      e[s->aidx[j]] = 0.0;
    }
    dgetrf_(&na, &na, s->lu, &na, s->lpiv, &info);
    free(e);
    free(col);
    if(info != 0){
      fprintf(stderr, "fd3d_strips: the strip operator is singular (dgetrf %d)\n", info);
    }
  }
  {
    na_2 = s->na;
    b = malloc(sizeof(double)*na_2);
    for(i = 0; i < na_2; i++)
      b[i] = rhs[s->aidx[i]];
    dgetrs_(&tr, &na_2, &one, s->lu, &na_2, s->lpiv, b, &na_2, &info);
    memset(V, 0, sizeof(double)*n);
    for(i = 0; i < na_2; i++)
      V[s->aidx[i]] = b[i];
    free(b);
  }
  return 1;
}
/* GMRES(40), right preconditioned, for the joint strip system; returns iterations */
static int strips_solve(struct strips_ctx *x, const double *rhs, double *V)
{
  int n = x->off[x->s->nf], i, it = 0, m = 40, j, k2, restart;
  double *vj;
  double *vn;
  double hn;
  double h;
  double *vk;
  double a1;
  double a2;
  double a1_2;
  double a2_2;
  double rr;
  int jj;
  int nj;
  double sm;
  if(!getenv("FD_STRIPS_GMRES"))
    return strips_direct(x, rhs, V);
  double *r = malloc(sizeof(double)*n), *w = malloc(sizeof(double)*n), *z = malloc(sizeof(double)*n);
  double *Vk = malloc(sizeof(double)*(size_t)(m + 1)*n), *H = calloc((size_t)(m + 1)*m, sizeof(double));
  double cs[64], sn[64], g[65], y[64], bnorm = 0.0, beta;
  for(i = 0; i < n; i++)
    bnorm += rhs[i]*rhs[i];
  bnorm = sqrt(bnorm);
  strips_mask(x, V);
  if(bnorm > 0.0)
    for(restart = 0; restart < 10; restart++){
      strips_A(x, V, w);
      for(i = 0; i < n; i++)
        r[i] = rhs[i] - w[i];
      strips_mask(x, r);
      beta = 0.0;
      for(i = 0; i < n; i++)
        beta += r[i]*r[i];
      beta = sqrt(beta);
      if(beta <= 1e-13*bnorm)
        break;
      for(i = 0; i < n; i++)
        Vk[i] = r[i]/beta;
      memset(g, 0, sizeof(g));
      g[0] = beta;
      for(j = 0; j < m; j++){
        vj = Vk + (size_t)j*n;
        vn = Vk + (size_t)(j + 1)*n;
        strips_P(x, vj, z);
        strips_A(x, z, vn);
        it++;
        for(k2 = 0; k2 <= j; k2++){
          h = 0.0;
          vk = Vk + (size_t)k2*n;
          for(i = 0; i < n; i++)
            h += vn[i]*vk[i];
          H[(size_t)k2*m + j] = h;
          for(i = 0; i < n; i++)
            vn[i] -= h*vk[i];
        }
        hn = 0.0;
        for(i = 0; i < n; i++)
          hn += vn[i]*vn[i];
        hn = sqrt(hn);
        H[(size_t)(j + 1)*m + j] = hn;
        if(hn > 0.0)
          for(i = 0; i < n; i++)
            vn[i] /= hn;
        for(k2 = 0; k2 < j; k2++){
          a1 = H[(size_t)k2*m + j];
          a2 = H[(size_t)(k2 + 1)*m + j];
          H[(size_t)k2*m + j] = cs[k2]*a1 + sn[k2]*a2;
          H[(size_t)(k2 + 1)*m + j] = -sn[k2]*a1 + cs[k2]*a2;
        }
        {
          a1_2 = H[(size_t)j*m + j];
          a2_2 = H[(size_t)(j + 1)*m + j];
          rr = sqrt(a1_2*a1_2 + a2_2*a2_2);
          cs[j] = rr > 0 ? a1_2/rr : 1.0;
          sn[j] = rr > 0 ? a2_2/rr : 0.0;
          H[(size_t)j*m + j] = rr;
          H[(size_t)(j + 1)*m + j] = 0.0;
          g[j + 1] = -sn[j]*g[j];
          g[j] = cs[j]*g[j];
        }
        if(fabs(g[j + 1]) <= 1e-13*bnorm){
          j++;
          break;
        }
      }
      if(j > m)
        j = m;
      {
        nj = j;
        for(jj = nj - 1; jj >= 0; jj--){
          sm = g[jj];
          for(k2 = jj + 1; k2 < nj; k2++)
            sm -= H[(size_t)jj*m + k2]*y[k2];
          y[jj] = sm/H[(size_t)jj*m + jj];
        }
        memset(w, 0, sizeof(double)*n);
        for(k2 = 0; k2 < nj; k2++)
          for(i = 0; i < n; i++)
            w[i] += Vk[(size_t)k2*n + i]*y[k2];
        strips_P(x, w, z);
        for(i = 0; i < n; i++)
          V[i] += z[i];
      }
    }
  strips_mask(x, V);
  free(r);
  free(w);
  free(z);
  free(Vk);
  free(H);
  return it;
}

int fd3d_strips_step(struct fd3d_strips *s, double **load)
{
  const struct fd3d_tab *t = s->f[0]->t;
  int nf = s->nf, i, l, it, c, r, k, q, j, itot = 0, nl = s->f[0]->nl, nst = 0;
  double dt = t->dt, hw = nl ? 0.5 : 1.0;
  struct strips_ctx X;
  double *Hx[FD3D_MULTI_MAX], *Hn[FD3D_MULTI_MAX], *Kx[FD3D_MULTI_MAX], *Kn[FD3D_MULTI_MAX], *Vp[FD3D_MULTI_MAX], *thb[FD3D_MULTI_MAX], *Hs[FD3D_MULTI_MAX];
  double _Complex *Xf[FD3D_MULTI_MAX], *Yf[FD3D_MULTI_MAX], *Ys[FD3D_MULTI_MAX], *Xs[FD3D_MULTI_MAX];
  double *V, *rhs;
  long lmaxd = s->hk;
  double * dl[8];
  struct fd3d_layer *ly;
  int N;
  int nz;
  int nys;
  size_t Mh;
  size_t Q;
  struct fd3d_layer *ly_2;
  int ns;
  int nys_2;
  int nz_2;
  size_t Mh_2;
  size_t Q_2;
  size_t nfs;
  size_t mm;
  int lmax;
  int slot;
  _Complex double *Yt;
  _Complex double *Yt2;
  struct fd3d_dtab *a;
  struct fd3d_dtab *b;
  struct fd3d_dtab *d;
  int nsk;
  int nsq;
  int Nq;
  int Nk;
  long slot_2;
  const double *vf;
  const double *vs;
  struct fd3d_mpair *cp;
  double sg;
  double sn;
  int nr;
  int nsr;
  double a_2;
  double b_2;
  const double *v;
  const double *Ws;
  const double *Wn;
  double a_3;
  double b_3;
  double a_4;
  double b_4;
  double ch;
  double vm;
  struct fd3d *f;
  struct fd3d_layer *ly_3;
  int N_2;
  int nz_3;
  int nys_3;
  int ns_2;
  size_t Mh_3;
  double *Hf;
  double *tnc;
  double *hsx;
  struct fd3d_dtab *d_2;
  struct fd3d_mpair *cp_2;
  double a_5;
  double b_5;
  double * dend[8];
  double *o;
  int nq;
  struct fd3d_layer *ly_4;
  int nys_4;
  int ns_3;
  size_t Mh_4;
  size_t Q_3;
  double *Vn;
  struct fd3d *f_2;
  struct fd3d_layer *ly_5;
  int N_3;
  int nz_4;
  int nys_5;
  int ns_4;
  int slot_3;
  size_t Mh_5;
  double *tmpf;
  double *vmf;
  double *vms;
  struct fd3d_mpair *cp_3;
  double a_6;
  struct fd3d *f_3;
  struct fd3d_layer *ly_6;
  memset(&X, 0, sizeof(X));
  X.s = s;
  X.nl = nl;
  X.hw = hw;
  X.Es[0] = X.Es[1] = 0.5*t->rho*t->cT;
  X.Es[2] = 0.5*t->rho*t->cL;
  for(k = 0; k < nf; k++){
    X.off[k] = nst;
    nst += 3*s->ly[k]->ns;
  }
  X.off[nf] = nst;
  V = malloc(sizeof(double)*nst);
  rhs = malloc(sizeof(double)*nst);
  if(s->mx && s->mx->static_fn){
    for(k = 0; k < nf; k++){
      if(!s->mx->sbuf[k])
        s->mx->sbuf[k] = malloc(sizeof(double)*s->f[k]->N);
      dl[k] = s->f[k]->delta;
    }
    s->mx->static_fn(s->mx->static_ctx, dl, s->mx->sbuf);
    for(k = 0; k < nf; k++){
      s->mx->slot[k].m = s->mx;
      s->mx->slot[k].p = k;
      s->f[k]->static_fn = multi_sbuf_out;
      s->f[k]->static_ctx = &s->mx->slot[k];
    }
  }
  for(k = 0; k < nf; k++){
    ly = s->ly[k];
    N = s->f[k]->N;
    nz = ly->nz;
    nys = ly->nys;
    Mh = ly->Mh;
    Q = Mh*nys;
    step_begin(s->f[k], load ? load[k] : NULL);
    Hx[k] = calloc(N, sizeof(double));
    Hn[k] = calloc(N, sizeof(double));
    Kx[k] = calloc(N, sizeof(double));
    Kn[k] = calloc(N, sizeof(double));
    Vp[k] = malloc(sizeof(double)*N);
    thb[k] = (s->f[k]->law == 1) ? malloc(sizeof(double)*N) : NULL;
    Xf[k] = malloc(sizeof(double _Complex)*Mh*nz);
    Yf[k] = calloc(2*Mh*nz, sizeof(double _Complex));
    Ys[k] = calloc(3*Q, sizeof(double _Complex));
    Xs[k] = malloc(sizeof(double _Complex)*Q);
    X.S3[k] = malloc(sizeof(double _Complex)*3*Q);
    X.P3[k] = malloc(sizeof(double _Complex)*3*Q);
    Hs[k] = calloc(3*ly->ns, sizeof(double));
    memcpy(ly->Vold3, ly->V3, sizeof(double)*3*ly->ns);
  }
  /* ---- known parts (as fd3d_msurf_step, per strip with its own fault; plus the dense couplings) ---- */
  for(k = 0; k < nf; k++){
    ly_2 = s->ly[k];
    ns = ly_2->ns;
    nys_2 = ly_2->nys;
    nz_2 = ly_2->nz;
    Mh_2 = ly_2->Mh;
    Q_2 = Mh_2*nys_2;
    nfs = (size_t)nys_2*nz_2;
    lmax = ly_2->L < ly_2->nhist + 1 ? ly_2->L : ly_2->nhist + 1;
    for(c = 0; c < 3; c++)
      layer_rows_hat(ly_2, ly_2->delta3 + (size_t)c*ns, X.S3[k] + (size_t)c*Q_2);
    surf_apply_rows(ly_2, 1, 0, X.S3[k], 1.0, Ys[k]);
    fault_rows_hat(ly_2, s->f[k]->delta, Xf[k]);
    for(r = 0; r < 3; r++)
      lxp_apply(ly_2, nys_2, nz_2, ly_2->Tfs[r][1], Xf[k], Ys[k] + (size_t)r*Mh_2*nys_2);
    for(c = 0; c < 3; c++){
      layer_rows_hat(ly_2, ly_2->delta3 + (size_t)c*ns, Xs[k]);
      if(!s->f[k]->static_fn)
        lxp_apply(ly_2, nz_2, nys_2, ly_2->Tsf[0][c][1], Xs[k], Yf[k]);
      if(!s->static_n_fn)
        lxp_apply(ly_2, nz_2, nys_2, ly_2->Tsf[1][c][1], Xs[k], Yf[k] + Mh_2*nz_2);
    }
    for(l = 2; l <= lmax; l++){
      slot = (ly_2->hk - (l - 1) + 2*ly_2->L) % ly_2->L;
      surf_apply_rows(ly_2, 0, l - 1, ly_2->hsx3 + (size_t)slot*3*Q_2, 1.0, Ys[k]);
      for(r = 0; r < 3; r++)
        lxp_apply(ly_2, nys_2, nz_2, ly_2->Tfs[r][0] + (size_t)(l - 1)*ly_2->nk*nfs, ly_2->hfx + (size_t)slot*Mh_2*nz_2, Ys[k] + (size_t)r*Mh_2*nys_2);
      for(c = 0; c < 3; c++)
        for(r = 0; r < 2; r++)
          lxp_apply(ly_2, nz_2, nys_2, ly_2->Tsf[r][c][0] + (size_t)(l - 1)*ly_2->nk*nfs, ly_2->hsx3 + ((size_t)slot*3 + c)*Mh_2*nys_2, Yf[k] + (size_t)r*Mh_2*nz_2);
    }
    if(nl){
      Yt = calloc(3*Q_2, sizeof(double _Complex));
      Yt2 = calloc(2*Mh_2*nz_2, sizeof(double _Complex));
      fault_rows_hat(ly_2, s->f[k]->Vold, Xf[k]);
      for(r = 0; r < 3; r++)
        lxp_apply(ly_2, nys_2, nz_2, ly_2->Tfs[r][2], Xf[k], Yt + (size_t)r*Mh_2*nys_2);
      for(mm = 0; mm < 3*Q_2; mm++)
        Ys[k][mm] += 0.5*Yt[mm];
      for(c = 0; c < 3; c++){
        layer_rows_hat(ly_2, ly_2->Vold3 + (size_t)c*ns, Xs[k]);
        for(r = 0; r < 2; r++)
          lxp_apply(ly_2, nz_2, nys_2, ly_2->Tsf[r][c][2], Xs[k], Yt2 + (size_t)r*Mh_2*nz_2);
      }
      for(mm = 0; mm < 2*Mh_2*nz_2; mm++)
        Yf[k][mm] += 0.5*Yt2[mm];
      for(c = 0; c < 3; c++)
        layer_rows_hat(ly_2, ly_2->Vold3 + (size_t)c*ns, X.S3[k] + (size_t)c*Q_2);
      surf_apply_rows(ly_2, 2, 0, X.S3[k], 0.5, Ys[k]);
      free(Yt);
      free(Yt2);
    }
    surf_reduce(ly_2, Yf[k], 2*Mh_2*nz_2);
    lx_inv(ly_2, nz_2, Yf[k], ly_2->nx, ly_2->m, 1.0, Hx[k]);
    lx_inv(ly_2, nz_2, Yf[k] + Mh_2*nz_2, ly_2->nx, ly_2->m, 1.0, Hn[k]);
    surf_reduce(ly_2, Ys[k], 3*Q_2);
    for(c = 0; c < 3; c++){
      lx_inv(ly_2, nys_2, Ys[k] + (size_t)c*Q_2, ly_2->nxs, 1, 1.0, Hs[k] + (size_t)c*ns);
      if(nl)
        for(i = 0; i < ns; i++)
          Hs[k][(size_t)c*ns + i] += 0.5*X.Es[c]*ly_2->Vold3[(size_t)c*ns + i];
    }
  }
  /* dense: other faults and strips -> strip k; other strips -> fault k (statics, lags >= 2, NL half) */
  for(k = 0; k < nf; k++)
    for(q = 0; q < nf; q++){
      a = &s->fs[q][k];
      b = &s->ss[q][k];
      d = &s->sf[q][k];
      nsk = 3*s->ly[k]->ns;
      nsq = 3*s->ly[q]->ns;
      Nq = s->f[q]->N;
      Nk = s->f[k]->N;
      if(q == k)
        continue;
      dmv(a->Ks, a->nr, a->nc, s->f[q]->delta, 1.0, Hs[k]);
      dmv(b->Ks, b->nr, b->nc, s->ly[q]->delta3, 1.0, Hs[k]);
      if(!s->f[k]->static_fn)
        dmv(d->Ks, Nk, d->nc, s->ly[q]->delta3, 1.0, Kx[k]);
      if(!s->static_n_fn)
        dmv(d->Ks + (size_t)Nk*d->nc, Nk, d->nc, s->ly[q]->delta3, 1.0, Kn[k]);
      for(l = 2; l <= s->L && l - 1 <= lmaxd; l++){
        slot_2 = (s->hk - (l - 1) + (long)s->Lb*1000) % s->Lb;
        vf = s->fring[q] + (size_t)slot_2*Nq;
        vs = s->sring[q] + (size_t)slot_2*nsq;
        dmv(a->W + (size_t)(l - 1)*a->nr*a->nc, a->nr, a->nc, vf, 1.0, Hs[k]);
        dmv(b->W + (size_t)(l - 1)*b->nr*b->nc, b->nr, b->nc, vs, 1.0, Hs[k]);
        dmv(d->W + (size_t)(l - 1)*d->nr*d->nc, Nk, d->nc, vs, 1.0, Kx[k]);
        dmv(d->W + (size_t)(l - 1)*d->nr*d->nc + (size_t)Nk*d->nc, Nk, d->nc, vs, 1.0, Kn[k]);
      }
      if(nl){
        dmv(a->K1, a->nr, a->nc, s->f[q]->Vold, 0.5, Hs[k]);
        dmv(b->K1, b->nr, b->nc, s->ly[q]->Vold3, 0.5, Hs[k]);
        dmv(d->K1, Nk, d->nc, s->ly[q]->Vold3, 0.5, Kx[k]);
        dmv(d->K1 + (size_t)Nk*d->nc, Nk, d->nc, s->ly[q]->Vold3, 0.5, Kn[k]);
      }
      (void)nsk;
    }
  /* fault <- fault (as fd3d_msurf_step) */
  if(s->mx)
    for(k = 0; k < nf; k++)
      for(q = 0; q < nf; q++){
        cp = &s->mx->c[q][k];
        sg = s->fsg[k]*s->fsg[q];
        sn = s->fsg[q];
        nr = s->f[k]->N;
        nsr = cp->ns;
        if(k == q)
          continue;
        if(!s->mx->static_fn)
          for(i = 0; i < nr; i++){
            a_2 = 0.0;
            b_2 = 0.0;
            for(j = 0; j < nsr; j++){
              a_2 += cp->Ks[(size_t)i*nsr + j]*s->f[q]->delta[j];
              b_2 += cp->Kn[(size_t)i*nsr + j]*s->f[q]->delta[j];
            }
            Kx[k][i] += sg*a_2;
            Kn[k][i] += sn*b_2;
          }
        for(l = 2; l <= cp->L && l - 1 <= s->mx->hk; l++){
          v = s->mx->ring[q] + (size_t)((s->mx->hk - (l - 1) + (long)s->mx->Lb*1000) % s->mx->Lb)*nsr;
          Ws = cp->Ws + (size_t)(l - 1)*nr*nsr;
          Wn = cp->Wn + (size_t)(l - 1)*nr*nsr;
          for(i = 0; i < nr; i++){
            a_3 = 0.0;
            b_3 = 0.0;
            for(j = 0; j < nsr; j++){
              a_3 += Ws[(size_t)i*nsr + j]*v[j];
              b_3 += Wn[(size_t)i*nsr + j]*v[j];
            }
            Kx[k][i] += sg*a_3;
            Kn[k][i] += sn*b_3;
          }
        }
        if(nl && cp->has_k1)
          for(i = 0; i < nr; i++){
            a_4 = 0.0;
            b_4 = 0.0;
            for(j = 0; j < nsr; j++){
              a_4 += cp->K1s[(size_t)i*nsr + j]*s->f[q]->Vold[j];
              b_4 += cp->K1n[(size_t)i*nsr + j]*s->f[q]->Vold[j];
            }
            Kx[k][i] += 0.5*sg*a_4;
            Kn[k][i] += 0.5*sn*b_4;
          }
      }
  /* ---- iterate faults and strips ---- */
  for(k = 0; k < nf; k++){
    memcpy(V + X.off[k], s->ly[k]->Vold3, sizeof(double)*3*s->ly[k]->ns);
    s->f[k]->thnew = thb[k];
  }
  for(it = 0; it < 200; it++){
    ch = 0.0;
    vm = 1e-30;
    for(k = 0; k < nf; k++){
      f = s->f[k];
      ly_3 = s->ly[k];
      N_2 = f->N;
      nz_3 = ly_3->nz;
      nys_3 = ly_3->nys;
      ns_2 = ly_3->ns;
      Mh_3 = ly_3->Mh;
      Hf = malloc(sizeof(double)*N_2);
      tnc = calloc(N_2, sizeof(double));
      hsx = calloc(2*N_2, sizeof(double));
      memset(Yf[k], 0, sizeof(double _Complex)*2*Mh_3*nz_3);
      for(c = 0; c < 3; c++){
        layer_rows_hat(ly_3, V + X.off[k] + (size_t)c*ns_2, Xs[k]);
        for(r = 0; r < 2; r++)
          lxp_apply(ly_3, nz_3, nys_3, ly_3->Tsf[r][c][2], Xs[k], Yf[k] + (size_t)r*Mh_3*nz_3);
      }
      surf_reduce(ly_3, Yf[k], 2*Mh_3*nz_3);
      for(i = 0; i < N_2; i++)
        Hf[i] = Hx[k][i] + Kx[k][i];
      lx_inv(ly_3, nz_3, Yf[k], ly_3->nx, ly_3->m, hw, Hf);
      memcpy(ly_3->tn, Hn[k], sizeof(double)*N_2);
      lx_inv(ly_3, nz_3, Yf[k] + Mh_3*nz_3, ly_3->nx, ly_3->m, hw, ly_3->tn);
      for(i = 0; i < N_2; i++)
        ly_3->tn[i] += Kn[k][i];
      for(q = 0; q < nf; q++)
        if(q != k){               /* the current step's other strips -> fault k */
          d_2 = &s->sf[q][k];
          dmv(d_2->K1, N_2, d_2->nc, V + X.off[q], hw, hsx);
          dmv(d_2->K1 + (size_t)N_2*d_2->nc, N_2, d_2->nc, V + X.off[q], hw, hsx + N_2);
        }
      for(i = 0; i < N_2; i++){
        Hf[i] += hsx[i];
        tnc[i] += hsx[N_2 + i];
      }
      if(s->mx)
        for(q = 0; q < nf; q++){
          cp_2 = &s->mx->c[q][k];
          if(q == k || !cp_2->has_k1)
            continue;
          for(i = 0; i < N_2; i++){
            a_5 = 0.0;
            b_5 = 0.0;
            for(j = 0; j < cp_2->ns; j++){
              a_5 += cp_2->K1s[(size_t)i*cp_2->ns + j]*s->f[q]->V[j];
              b_5 += cp_2->K1n[(size_t)i*cp_2->ns + j]*s->f[q]->V[j];
            }
            Hf[i] += hw*s->fsg[k]*s->fsg[q]*a_5;
            tnc[i] += hw*s->fsg[q]*b_5;
          }
        }
      for(i = 0; i < N_2; i++)
        ly_3->tn[i] += tnc[i];
      if(s->static_n_fn){
        o = malloc(sizeof(double)*N_2);
        for(q = 0; q < nf; q++){
          nq = s->f[q]->N;
          dend[q] = malloc(sizeof(double)*nq);
          for(i = 0; i < nq; i++)
            dend[q][i] = s->f[q]->delta[i] + (nl ? 0.5*dt*(s->f[q]->Vold[i] + s->f[q]->V[i]) : dt*s->f[q]->V[i]);
        }
        s->static_n_fn(s->static_n_ctx, k, dend, o);
        for(i = 0; i < N_2; i++)
          ly_3->tn[i] += o[i];
        for(q = 0; q < nf; q++)
          free(dend[q]);
        free(o);
      }
      if(ly_3->normal_coupling)
        for(i = 0; i < N_2; i++){
          f->sigma[i] = ly_3->sig0[i] - ly_3->tn[i];
          if(f->sigma[i] < 0.0)
            f->sigma[i] = 0.0;
        }
      memcpy(Vp[k], f->V, sizeof(double)*N_2);
      itot += step_solve_g(f, Hf, it > 0 ? Vp[k] : NULL, NULL);
      free(Hf);
      free(tnc);
      free(hsx);
    }
    /* strips: zero traction, jointly */
    for(k = 0; k < nf; k++){
      ly_4 = s->ly[k];
      nys_4 = ly_4->nys;
      ns_3 = ly_4->ns;
      Mh_4 = ly_4->Mh;
      Q_3 = Mh_4*nys_4;
      memset(Ys[k], 0, sizeof(double _Complex)*3*Q_3);
      fault_rows_hat(ly_4, s->f[k]->V, Xf[k]);
      for(r = 0; r < 3; r++)
        lxp_apply(ly_4, nys_4, ly_4->nz, ly_4->Tfs[r][2], Xf[k], Ys[k] + (size_t)r*Mh_4*nys_4);
      surf_reduce(ly_4, Ys[k], 3*Q_3);
      for(i = 0; i < 3*ns_3; i++)
        rhs[X.off[k] + i] = -Hs[k][i];
      for(c = 0; c < 3; c++)
        lx_inv(ly_4, nys_4, Ys[k] + (size_t)c*Q_3, ly_4->nxs, 1, -hw, rhs + X.off[k] + (size_t)c*ns_3);
      for(q = 0; q < nf; q++)
        if(q != k)
          dmv(s->fs[q][k].K1, s->fs[q][k].nr, s->fs[q][k].nc, s->f[q]->V, -hw, rhs + X.off[k]);
    }
    {
      Vn = malloc(sizeof(double)*nst);
      memcpy(Vn, V, sizeof(double)*nst);
      s->n_inner += strips_solve(&X, rhs, Vn);
      for(i = 0; i < nst; i++){
        if(fabs(Vn[i] - V[i]) > ch)
          ch = fabs(Vn[i] - V[i]);
        if(fabs(Vn[i]) > vm)
          vm = fabs(Vn[i]);
      }
      memcpy(V, Vn, sizeof(double)*nst);
      free(Vn);
    }
    for(k = 0; k < nf; k++)
      for(i = 0; i < s->f[k]->N; i++){
        if(fabs(s->f[k]->V[i] - Vp[k][i]) > ch)
          ch = fabs(s->f[k]->V[i] - Vp[k][i]);
        if(fabs(s->f[k]->V[i]) > vm)
          vm = fabs(s->f[k]->V[i]);
      }
    if(it > 0 && ch <= 1e-12*vm)
      break;
  }
  s->last_iters = it + 1;
  s->n_outer += it + 1;
  s->n_steps++;
  for(k = 0; k < nf; k++){
    memcpy(s->ly[k]->V3, V + X.off[k], sizeof(double)*3*s->ly[k]->ns);
    s->f[k]->thnew = NULL;
    if(thb[k]){
      memcpy(s->f[k]->theta, thb[k], sizeof(double)*s->f[k]->N);
      free(thb[k]);
    }
  }
  /* ---- finish: faults' known parts, rings, step_end ---- */
  for(k = 0; k < nf; k++){
    f_2 = s->f[k];
    ly_5 = s->ly[k];
    N_3 = f_2->N;
    nz_4 = ly_5->nz;
    nys_5 = ly_5->nys;
    ns_4 = ly_5->ns;
    slot_3 = ly_5->hk % ly_5->L;
    Mh_5 = ly_5->Mh;
    tmpf = calloc(N_3, sizeof(double));
    vmf = malloc(sizeof(double)*N_3);
    vms = malloc(sizeof(double)*3*ns_4);
    for(i = 0; i < 3*ns_4; i++){
      vms[i] = nl ? 0.5*(ly_5->Vold3[i] + ly_5->V3[i]) : ly_5->V3[i];
      ly_5->delta3[i] += vms[i]*dt;
    }
    for(c = 0; c < 3; c++)
      layer_rows_hat(ly_5, vms + (size_t)c*ns_4, ly_5->hsx3 + ((size_t)slot_3*3 + c)*Mh_5*nys_5);
    memcpy(s->sring[k] + (size_t)(s->hk % s->Lb)*3*ns_4, vms, sizeof(double)*3*ns_4);
    memset(Yf[k], 0, sizeof(double _Complex)*2*Mh_5*nz_4);
    for(c = 0; c < 3; c++){
      layer_rows_hat(ly_5, ly_5->V3 + (size_t)c*ns_4, Xs[k]);
      for(r = 0; r < 2; r++)
        lxp_apply(ly_5, nz_4, nys_5, ly_5->Tsf[r][c][2], Xs[k], Yf[k] + (size_t)r*Mh_5*nz_4);
    }
    surf_reduce(ly_5, Yf[k], 2*Mh_5*nz_4);
    lx_inv(ly_5, nz_4, Yf[k], ly_5->nx, ly_5->m, hw, tmpf);
    for(i = 0; i < N_3; i++)
      f_2->H[i] += Hx[k][i] + Kx[k][i] + tmpf[i];
    for(q = 0; q < nf; q++)
      if(q != k)
        dmv(s->sf[q][k].K1, N_3, s->sf[q][k].nc, s->ly[q]->V3, hw, f_2->H);
    if(s->mx)
      for(q = 0; q < nf; q++){
        cp_3 = &s->mx->c[q][k];
        if(q == k || !cp_3->has_k1)
          continue;
        for(i = 0; i < N_3; i++){
          a_6 = 0.0;
          for(j = 0; j < cp_3->ns; j++)
            a_6 += cp_3->K1s[(size_t)i*cp_3->ns + j]*s->f[q]->V[j];
          f_2->H[i] += hw*s->fsg[k]*s->fsg[q]*a_6;
        }
      }
    for(i = 0; i < N_3; i++)
      vmf[i] = nl ? 0.5*(f_2->Vold[i] + f_2->V[i]) : f_2->V[i];
    fault_rows_hat(ly_5, vmf, ly_5->hfx + (size_t)slot_3*Mh_5*nz_4);
    memcpy(s->fring[k] + (size_t)(s->hk % s->Lb)*N_3, vmf, sizeof(double)*N_3);
    if(s->mx)
      memcpy(s->mx->ring[k] + (size_t)(s->mx->hk % s->mx->Lb)*N_3, vmf, sizeof(double)*N_3);
    ly_5->hk = (ly_5->hk + 1) % ly_5->L;
    if(ly_5->nhist < ly_5->L)
      ly_5->nhist++;
    free(tmpf);
    free(vmf);
    free(vms);
  }
  s->hk++;
  if(s->mx)
    s->mx->hk++;
  for(k = 0; k < nf; k++){
    f_3 = s->f[k];
    ly_6 = s->ly[k];
    step_end(f_3);
    if(ly_6->normal_coupling)
      for(i = 0; i < f_3->N; i++){
        f_3->sigma[i] = ly_6->sig0[i] - ly_6->tn[i];
        if(f_3->sigma[i] < 0.0)
          f_3->sigma[i] = 0.0;
      }
    free(Hx[k]);
    free(Hn[k]);
    free(Kx[k]);
    free(Kn[k]);
    free(Vp[k]);
    free(Xf[k]);
    free(Yf[k]);
    free(Ys[k]);
    free(Xs[k]);
    free(X.S3[k]);
    free(X.P3[k]);
    free(Hs[k]);
  }
  free(V);
  free(rhs);
  return itot;
}

/* check (1) of H6 (fd3d_test h6step): two faults sharing a strike, with strips on the common
   lattice (partitioned) against the several-fault surface (G3b) on the same cells: prescribed
   slip rates (first `kin` steps 1 m/s on both faults, then free with slip weakening if sw) and
   the fault tractions and surface rates compared over nsteps; returns the largest relative
   difference of the fault shear tractions */
double fd3d_h6_step_check(int nsteps, int sw, double th_deg, int verbose)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = 1000.0, dt = 0.5*dx/cT, hs = 250.0;
  int nx = 4, nz = 2, m = 4, L = 20, k, i, n, nxs = 24, nys = 16, ny2[2] = {16, 16}, nx2[2] = {24, 24};   /* small: dense strip tables are (3 ns)^2 L */
  double ytop[2] = {0.0, -3000.0}, ztop[2] = {0.0, 0.0}, dip[2] = {60.0, 75.0}, fs[2] = {1.0, -1.0}, x0f = 0.0, xs0 = -2000.0;
  double x0f2[2] = {0.0, 0.0}, xs02[2] = {-2000.0, -2000.0}, th[2] = {th_deg, th_deg}, O[2][3] = {{0, 0, 0}, {0, 0, 0}};
  double *yc = malloc(sizeof(double)*nys), *wy = malloc(sizeof(double)*nys);
  const double *ycp[2] = {yc, yc}, *wyp[2] = {wy, wy};
  struct fd3d_tab *tb = fd3d_tab_create_c(nx, nz, dx, dt, cL, cT, rho, 0, 2);
  struct fd3d *fa[2], *fb[2];
  struct fd3d_msurf *ms;
  struct fd3d_strips *st;
  double worst = 0.0, scale = 0.0, wv = 0.0, sv = 0.0;
  double *vp[2];
  struct fd3d * ff[2];
  int u;
  int kin;
  for(i = 0; i < nys; i++){
    yc[i] = -4500.0 + (i + 0.5)*hs;
    wy[i] = hs;
  }
  for(k = 0; k < 2; k++){
    fa[k] = fd3d_create(tb);
    fb[k] = fd3d_create(tb);
    vp[k] = malloc(sizeof(double)*fa[k]->N);
    if(sw){
      ff[0] = fa[k];
      ff[1] = fb[k];
      for(u = 0; u < 2; u++){
        ff[u]->law = 0;
        ff[u]->mu_s = 0.677;
        ff[u]->mu_d = 0.525;
        ff[u]->d0 = 0.4;
        for(i = 0; i < ff[u]->N; i++){
          ff[u]->sigma[i] = 120e6;
          ff[u]->tau0[i] = (k == 0 && i < 4) ? 82e6 : 70e6;
        }
      }
    }
  }
  ms = fd3d_msurf_create(2, fa, ytop, ztop, dip, fs, 2, x0f, m, nxs, nys, xs0, yc, wy, L);
  st = fd3d_strips_create(2, fb, th, (const double (*)[3])O, x0f2, ytop, ztop, dip, fs, 2, m, nx2, ny2, xs02, ycp, wyp, L);
  if(!ms || !st){
    fprintf(stderr, "h6step: setup failed (%p %p)\n", (void *)ms, (void *)st);
    return 1e300;
  }
  for(n = 0; n < nsteps; n++){
    kin = (n < 5);
    for(k = 0; k < 2; k++){
      for(i = 0; i < fa[k]->N; i++)
        vp[k][i] = kin ? 1.0 : 0.0;
      if(!sw || kin){
        fa[k]->vpres = vp[k];
        fb[k]->vpres = vp[k];
      } else {
        fa[k]->vpres = NULL;
        fb[k]->vpres = NULL;
      }
    }
    fd3d_msurf_step(ms, NULL);
    fd3d_strips_step(st, NULL);
    for(k = 0; k < 2; k++)
      for(i = 0; i < fa[k]->N; i++){
        worst = fmax(worst, fabs(fa[k]->tau[i] - fb[k]->tau[i]));
        scale = fmax(scale, fabs(fa[k]->tau[i] - fa[k]->tau0[i]));
        wv = fmax(wv, fabs(fa[k]->V[i] - fb[k]->V[i]));
        sv = fmax(sv, fabs(fa[k]->V[i]));
      }
  }
  if(verbose)
    printf("h6step: %d steps (%s), strike %g: max |tau(strips) - tau(G3b)| %.3e of %.3e (traction change), max |V diff| %.3e of %.3e; outer iterations %.1f, GMRES %.1f per step\n",
        nsteps, sw ? "slip weakening after 5 prescribed steps" : "prescribed slip", th_deg, worst, scale, wv, sv, (double)st->n_outer/st->n_steps, (double)st->n_inner/st->n_steps);
  fd3d_strips_free(st);
  fd3d_msurf_free(ms);
  for(k = 0; k < 2; k++){
    fa[k]->layer = NULL;
    fb[k]->layer = NULL;
    fd3d_free(fa[k]);
    fd3d_free(fb[k]);
    free(vp[k]);
  }
  fd3d_tab_free(tb);
  free(yc);
  free(wy);
  return scale > 0.0 ? worst/scale : -1.0;
}

/* H6 step 4 (fd3d_test h6strike): two faults whose strikes differ by `ang` degrees, each with
   a strip of cells dx/m aligned with its own trace (rows +-2.5 km, columns 1.5 km beyond the
   trace ends), partitioned by the nearest trace. Fault 0 along x from the origin, dipping 60
   degrees (dip slip); fault 1 starting 3 km away at (0, -3 km), turned by ang, dipping 75
   degrees. Both slip at 1 m/s for 5 steps, then stop; the fault shear and normal tractions of
   nsteps steps are written to tau (2 x 2 N per step). Returns 0 on success. */
int fd3d_h6_strike_run(double ang, int m, int nsteps, double *tau, int verbose)
{
  double cL = 6000, cT = 3464, rho = 2670, dx = 1000.0, dt = 0.5*dx/cT, hs = dx/m;
  int nx = 4, nz = 2, L = 20, k, i, n, nxs[2], nys[2];
  double th[2] = {0.0, ang}, O[2][3] = {{0, 0, 0}, {0, -3000.0, 0}}, x0f[2] = {0.0, 0.0}, ytop[2] = {0.0, 0.0}, ztop[2] = {0.0, 0.0};
  double dip[2] = {60.0, 75.0}, fs[2] = {1.0, 1.0}, xs0[2] = {-1500.0, -1500.0}, *yc[2], *wy[2];
  struct fd3d_tab *tb = fd3d_tab_create_c(nx, nz, dx, dt, cL, cT, rho, 0, 2);
  struct fd3d *f[2];
  struct fd3d_strips *s;
  double *vp = malloc(sizeof(double)*nx*nz);
  int nb;
  int nc;
  int a0;
  int a1;
  {
    /* cell boundaries on the trace and on the fault ends: margins that are multiples of hs */
    nb = (int)ceil(2500.0/hs - 1e-9);
    nc = (int)ceil(1500.0/hs - 1e-9);
    for(k = 0; k < 2; k++){
      xs0[k] = -nc*hs;
      nxs[k] = (int)floor(nx*dx/hs + 0.5) + 2*nc;
      nys[k] = 2*nb;
      yc[k] = malloc(sizeof(double)*nys[k]);
      wy[k] = malloc(sizeof(double)*nys[k]);
      for(i = 0; i < nys[k]; i++){
        yc[k][i] = -nb*hs + (i + 0.5)*hs;
        wy[k][i] = hs;
      }
      f[k] = fd3d_create(tb);
    }
  }
  s = fd3d_strips_create(2, f, th, (const double (*)[3])O, x0f, ytop, ztop, dip, fs, 2, m, nxs, nys, xs0, (const double **)yc, (const double **)wy, L);
  if(!s){
    fprintf(stderr, "h6strike: setup failed\n");
    return 1;
  }
  if(verbose){
    a0 = 0;
    a1 = 0;
    for(i = 0; i < s->ly[0]->ns; i++)
      a0 += s->active[0][i];
    for(i = 0; i < s->ly[1]->ns; i++)
      a1 += s->active[1][i];
    fprintf(stderr, "h6strike: turn %g deg, cells %g m: strips %d x %d, active %d and %d\n", ang, hs, nxs[0], nys[0], a0, a1);
  }
  for(n = 0; n < nsteps; n++){
    for(i = 0; i < nx*nz; i++)
      vp[i] = (n < 5) ? 1.0 : 0.0;
    for(k = 0; k < 2; k++)
      f[k]->vpres = vp;
    fd3d_strips_step(s, NULL);
    for(k = 0; k < 2; k++)
      for(i = 0; i < f[k]->N; i++){
        tau[((size_t)n*2 + k)*2*f[k]->N + i] = f[k]->tau[i];
        tau[((size_t)n*2 + k)*2*f[k]->N + f[k]->N + i] = s->ly[k]->tn[i];
      }
  }
  if(verbose)
    fprintf(stderr, "h6strike: outer %.1f, GMRES %.1f per step\n", (double)s->n_outer/s->n_steps, (double)s->n_inner/s->n_steps);
  fd3d_strips_free(s);
  for(k = 0; k < 2; k++){
    f[k]->layer = NULL;
    f[k]->vpres = NULL;
    fd3d_free(f[k]);
    free(yc[k]);
    free(wy[k]);
  }
  fd3d_tab_free(tb);
  free(vp);
  return 0;
}
