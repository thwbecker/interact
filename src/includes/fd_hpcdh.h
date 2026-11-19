#ifndef FD_HPCDH_H
#define FD_HPCDH_H
/*
  Reflected part of the 2-D in-plane half-plane response by Cagniard-de Hoop (fd_plan.md M11):
  the same quantities as fd_hpdwn_batch_ext (fd_hpdwn.h), the shear and normal tractions (in
  fd_project's convention) on receiver elements from unit slip rate on the source segments, at
  the lags t_l = l dt, l = 0 .. nt - 1, and the static reflected tractions (ks, kn, by Crouch and
  Starfield's half plane minus their full plane, interact's tdd_coeff), evaluated in the time domain: per entry the four rays (P,P) (P,S) (S,P)
  (S,S), each point response (1/pi) Im[F(q) dq/dt] on its Cagniard path, integrated along the
  source segment with the arrival fronts removed by a cosine substitution; beyond 40 times the
  entry's first reflected arrival, R = K t + c - b / t (the rays cancel in their sum and lose
  digits later). Port of fully_dynamic/proto/halfplane_cdh.py. Frame x, z (z depth), as fd_hpdwn.
  Returns 0 on success.
*/
int fd_hpcdh_batch(int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                   int nr, const double *rx, const double *rz, const double *rtx, const double *rtz,
                   int npr, const double *prx, const double *prz,
                   double alpha, double beta, double rho, double dt, int nt, double *os, double *on, double *ks, double *kn);
#endif
