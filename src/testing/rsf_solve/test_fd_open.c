/*
  test of the 2-D in-plane opening kernel (fd_kernels_open.c, fd_plan.md M7 step 1): reads
  receiver x y and time t per line from stdin, writes the physical stress (sxx, sxy, syy) of a
  unit opening-rate segment of half-length a = 250 m (alpha 6000, beta 3464 m/s, rho 2670);
  fully_dynamic/test_fd_open.py compares with inplane_open.py
*/
#include <stdio.h>
#include <stdlib.h>
#include "fd_kernels.h"

int main(int argc, char **argv)
{
  double xr;
  double yr;
  double t;
  double s[3];
  double alpha;
  double beta;
  double rho;
  double a;
  alpha = 6000.0;
  beta = 3464.0;
  rho = 2670.0;
  a = 250.0;
  while(scanf("%lf %lf %lf", &xr, &yr, &t) == 3){
    fd_seg_ramp_inplane_open(xr, yr, t, a, alpha, beta, rho, 0.0, s);
    printf("%.17e %.17e %.17e\n", s[0], s[1], s[2]);
  }
  return 0;
}
