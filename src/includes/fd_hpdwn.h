#ifndef FD_HPDWN_H
#define FD_HPDWN_H
/*
  Reflected part of the 2-D in-plane half-plane response (fd_plan.md M10), by the discrete
  wavenumber method: the stresses at receivers from unit slip rate (a slip ramp since t = 0) on
  straight source segments, minus the full-plane part, in the half plane z > 0 (z depth) with a
  traction-free surface z = 0. Port of fully_dynamic/halfplane_dwn.py (part "refl").
  Sources j: centre (sx[j], sz[j]), unit tangent (stx[j], stz[j]) in the (x, z) frame, half-length
  sa[j]; receivers i: points (rx[i], rz[i]). Output: sxx, sxz, szz (x, z frame) at times t[l],
  arrays [((i * ns) + j) * nt + l].
  Parameters: T time window, nw frequencies (nw / 2 + 1 used), wIT the damping omega_I T, L the
  source period, M the wavenumbers |m| <= M (k = 2 pi m / L).
*/
int fd_hpdwn_reflected(int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                       int nr, const double *rx, const double *rz, int nt, const double *t,
                       double alpha, double beta, double rho, double T, int nw, double wIT, double L, int M,
                       double *oxx, double *oxz, double *ozz);
/*
  Batch version for the tables (fd_plan.md M10): the reflected tractions on receiver elements
  (centres (rx, rz), unit tangents (rtx, rtz) in the (x, z) frame) from unit slip rate on the
  source segments, at the lags t_l = l dt, l = 0 .. nt - 1, as the shear and normal tractions of
  fd_project (y-up frame: tangent (rtx, -rtz), normal its rot90); output os, on
  [((i * ns) + j) * nt + l]. Per complex frequency the sums over the wavenumbers are matrix
  products (BLAS zgemm); one FFT per pair gives the lags. The window T = N dt with N the FFT size
  for 2.5 nt, the damping omega_I T = 18, the source period L = alpha T plus twice the extent,
  wavenumbers to 2 pi / h with h a quarter of the smallest half-length or of the smallest
  source plus receiver depth. ks, kn (may be NULL): the static reflected tractions (rate of the
  ramp response for t to infinity), from the same sums at a small frequency. Returns 0 on success.
*/
int fd_hpdwn_batch(int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                   int nr, const double *rx, const double *rz, const double *rtx, const double *rtz,
                   double alpha, double beta, double rho, double dt, int nt, double *os, double *on, double *ks, double *kn);
/* the same with the discretisation (window, wavenumbers, source period) chosen from the receivers
   prx, prz (npr of them; all receivers of the pair), so that a subset of the receivers (one MPI
   rank's rows) gets exactly the values it would get in the full batch */
int fd_hpdwn_batch_ext(int ns, const double *sx, const double *sz, const double *stx, const double *stz, const double *sa,
                       int nr, const double *rx, const double *rz, const double *rtx, const double *rtz,
                       int npr, const double *prx, const double *prz,
                       double alpha, double beta, double rho, double dt, int nt, double *os, double *on, double *ks, double *kn);
#endif
