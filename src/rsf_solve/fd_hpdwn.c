/* see fd_hpdwn.h */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "fd_hpdwn.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* complex square root with non-negative real part */
static double complex csqrt_pos(double complex z)
{
  double complex r;
  r = csqrt(z);
  if(creal(r) < 0.0)
    r = -r;
  return r;
}

/* traction (txz, tzz) at z = 0 of the upgoing incident field of unit slip on source j, per unit
   of the common factor exp(i k (x - xs)), wavenumber k, frequency w */
static void incident_traction(double k, double complex w, double alpha, double beta, double rho,
                              double tx, double tz, double a, double zs, double complex *txz, double complex *tzz)
{
  double lam, mu, nx, nz, m[2][2];
  double complex kb2, nu, Dr[2], Ds[2], g, Gip, u[2], acc, exx, ezz, sxz, szz, cc, seg, ph;
  int c, i, p, q;
  lam = rho*(alpha*alpha - 2.0*beta*beta);
  mu = rho*beta*beta;
  nx = tz;                                       /* normal as in halfplane_dwn.py */
  nz = -tx;
  m[0][0] = 2.0*tx*nx*mu;
  m[0][1] = (tx*nz + tz*nx)*mu;
  m[1][0] = m[0][1];
  m[1][1] = 2.0*tz*nz*mu;
  kb2 = (w/beta)*(w/beta);
  *txz = 0.0;
  *tzz = 0.0;
  for(c = 0; c < 2; c++){                        /* 0: P (alpha), 1: S (beta) */
    nu = csqrt_pos(k*k - (c ? kb2 : (w/alpha)*(w/alpha)));
    Dr[0] = I*k;
    Dr[1] = nu;                                   /* upgoing side: Zsign = -1, d/dz -> +nu */
    Ds[0] = -I*k;
    Ds[1] = -nu;
    g = 1.0/(2.0*nu);
    for(i = 0; i < 2; i++){
      acc = 0.0;
      for(p = 0; p < 2; p++){
        if(c)
          Gip = (kb2*(i == p ? 1.0 : 0.0) + Dr[i]*Dr[p])*g/(rho*w*w);
        else
          Gip = -(Dr[i]*Dr[p])*g/(rho*w*w);
        for(q = 0; q < 2; q++)
          acc += m[p][q]*Ds[q]*Gip;
      }
      u[i] = acc;
    }
    exx = Dr[0]*u[0];
    ezz = Dr[1]*u[1];
    szz = lam*exx + (lam + 2.0*mu)*ezz;
    sxz = mu*(Dr[1]*u[0] + Dr[0]*u[1]);
    cc = -I*k*tx - nu*tz;
    if(cabs(cc*a) < 1e-8)
      seg = 2.0*a;
    else
      seg = 2.0*csinh(cc*a)/cc;
    ph = cexp(-nu*zs);
    *txz += sxz*seg*ph;
    *tzz += szz*seg*ph;
  }
}

int fd_hpdwn_reflected(int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                       int nr, const double *rx, const double *rz, int nt, const double *t,
                       double alpha, double beta, double rho, double T, int nw, double wIT, double L, int M,
                       double *oxx, double *oxz, double *ozz)
{
  double lam, mu, wI, dw, k, fac;
  double complex w, na, nb, PxzA, PzzA, PxzB, PzzB, det, txz, tzz, A, B, ea, eb, uxA, uzA, uxB, uzB, ph, Bw, e;
  double complex *acc;
  double complex S[3];
  int jw, m, i, j, l, c;
  size_t np, o;
  lam = rho*(alpha*alpha - 2.0*beta*beta);
  mu = rho*beta*beta;
  wI = wIT/T;
  dw = 2.0*M_PI/T;
  np = (size_t)nr*ns;
  acc = (double complex *)calloc(3*np*nt, sizeof(double complex));
  if(!acc)
    return 1;
  for(jw = 0; jw <= nw/2; jw++){
    w = jw*dw + I*wI;
    Bw = -1.0/(w*w);                               /* the ramp */
    fac = (jw == 0) ? 1.0 : 2.0;
    for(m = -M; m <= M; m++){
      k = 2.0*M_PI*m/L;
      na = csqrt_pos(k*k - (w/alpha)*(w/alpha));
      nb = csqrt_pos(k*k - (w/beta)*(w/beta));
      PxzA = mu*(I*k*(-na) + (-na)*I*k);
      PzzA = lam*(I*k*I*k + na*na) + 2.0*mu*na*na;
      PxzB = mu*(I*k*I*k - nb*nb);
      PzzB = lam*(I*k*nb - nb*I*k) + 2.0*mu*(-nb*I*k);
      det = PxzA*PzzB - PxzB*PzzA;
      for(j = 0; j < ns; j++){
        incident_traction(k, w, alpha, beta, rho, stx[j], stz[j], sa[j], sz[j], &txz, &tzz);
        A = (-txz*PzzB + tzz*PxzB)/det;
        B = (-tzz*PxzA + txz*PzzA)/det;
        uxA = I*k*A;
        uzA = -na*A;
        uxB = nb*B;
        uzB = I*k*B;
        for(i = 0; i < nr; i++){
          ea = cexp(-na*rz[i]);
          eb = cexp(-nb*rz[i]);
          ph = cexp(I*k*(rx[i] - sx[j]))*Bw*fac/L;
          o = ((size_t)i*ns + j)*3;
          /* sxx, sxz, szz at the receiver, then their time samples */
            S[0] = (((lam + 2.0*mu)*(I*k*uxA) + lam*(-na*uzA))*ea + ((lam + 2.0*mu)*(I*k*uxB) + lam*(-nb*uzB))*eb)*ph;
            S[1] = (mu*((-na*uxA + I*k*uzA)*ea + (-nb*uxB + I*k*uzB)*eb))*ph;
            S[2] = ((lam*(I*k*uxA) + (lam + 2.0*mu)*(-na*uzA))*ea + (lam*(I*k*uxB) + (lam + 2.0*mu)*(-nb*uzB))*eb)*ph;
            for(l = 0; l < nt; l++){
              e = cexp(-I*(jw*dw)*t[l]);
              for(c = 0; c < 3; c++)
                acc[(o + c)*nt + l] += S[c]*e;
            }
        }
      }
    }
  }
  for(i = 0; i < nr; i++)
    for(j = 0; j < ns; j++)
      for(l = 0; l < nt; l++){
        o = ((size_t)i*ns + j);
        fac = exp(wI*t[l])*dw/(2.0*M_PI);
        oxx[o*nt + l] = creal(acc[(o*3 + 0)*nt + l])*fac;
        oxz[o*nt + l] = creal(acc[(o*3 + 1)*nt + l])*fac;
        ozz[o*nt + l] = creal(acc[(o*3 + 2)*nt + l])*fac;
      }
  free(acc);
  return 0;
}
