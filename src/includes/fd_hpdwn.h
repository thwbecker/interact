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
#endif
