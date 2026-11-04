#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fd_hpdwn.h"
/* test of fd_hpdwn.c (fd_plan.md M10): reflected half-plane stresses for one source segment and
   one receiver, read from stdin as xs zs tx tz a xr zr T nw wIT L M nt t_1 .. t_nt (x, z frame, z
   depth); writes sxx sxz szz per time; fully_dynamic/test_fd_hpdwn.py compares with halfplane_dwn.py */
int main(void)
{
  double sx, sz, tx, tz, a, rx, rz, T, wIT, L, t[64], oxx[64], oxz[64], ozz[64];
  int nw, M, nt, l;
  if(scanf("%lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %d %d", &sx, &sz, &tx, &tz, &a, &rx, &rz, &T, &nw, &wIT, &L, &M, &nt) != 13)
    return 1;
  for(l = 0; l < nt; l++)
    if(scanf("%lf", t + l) != 1)
      return 1;
  fd_hpdwn_reflected(1, &sx, &sz, &tx, &tz, &a, 1, &rx, &rz, nt, t, 6000.0, 3464.0, 2670.0, T, nw, wIT, L, M, oxx, oxz, ozz);
  for(l = 0; l < nt; l++)
    printf("%.12e %.12e %.12e\n", oxx[l], oxz[l], ozz[l]);
  return 0;
}
