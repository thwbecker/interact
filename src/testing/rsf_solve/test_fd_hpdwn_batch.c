#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fd_hpdwn.h"
/* test of fd_hpdwn_batch (fd_plan.md M10): a 30 degree fault of n elements of h, top at depth top, all
   pairs, nt lags of dt = 0.5 h / beta; writes the static reflected tractions (K i j shear normal)
   and the lag series (S i j l shear normal); fully_dynamic/test_fd_hpdwn_batch.py compares them with
   Crouch and Starfield and with halfplane_dwn.py.
   usage: test_fd_hpdwn_batch [n] [h] [top] [nt] */
int main(int argc, char **argv)
{
  int n, nt, i, j, l;
  double h, top, dip, dt, *x, *z, *tx, *tz, *a, *os, *on, *ks, *kn;
  n = argc > 1 ? atoi(argv[1]) : 4;
  h = argc > 2 ? atof(argv[2]) : 250.0;
  top = argc > 3 ? atof(argv[3]) : 500.0;
  nt = argc > 4 ? atoi(argv[4]) : 40;
  dip = 30.0*M_PI/180.0;
  dt = 0.5*h/3464.0;
  x = malloc(sizeof(double)*n); z = malloc(sizeof(double)*n); tx = malloc(sizeof(double)*n); tz = malloc(sizeof(double)*n); a = malloc(sizeof(double)*n);
  for(i = 0; i < n; i++){
    tx[i] = cos(dip); tz[i] = sin(dip); a[i] = 0.5*h;
    x[i] = (i + 0.5)*h*tx[i]; z[i] = top + (i + 0.5)*h*tz[i];
  }
  os = malloc(sizeof(double)*n*n*nt); on = malloc(sizeof(double)*n*n*nt); ks = malloc(sizeof(double)*n*n); kn = malloc(sizeof(double)*n*n);
  if(fd_hpdwn_batch(n, x, z, tx, tz, a, n, x, z, tx, tz, 6000.0, 3464.0, 2670.0, dt, nt, os, on, ks, kn)){
    fprintf(stderr, "failed\n");
    return 1;
  }
  printf("# n %d h %g top %g nt %d dt %.17e\n", n, h, top, nt, dt);
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++){
      printf("K %d %d %.12e %.12e\n", i, j, ks[i*n + j], kn[i*n + j]);
      for(l = 0; l < nt; l++)
        printf("S %d %d %d %.12e %.12e\n", i, j, l, os[(i*n + j)*nt + l], on[(i*n + j)*nt + l]);
    }
  return 0;
}
