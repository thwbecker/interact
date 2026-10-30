/*
  dense (non-parallel) pair tables of the 2-D fully dynamic code (fd_plan.md N1): a source
  segment along x and a receiver segment rotated by an angle, kinked at a common end point;
  writes the static blocks and W[l] for a few lags (receiver x source) for comparison with
  the prototype's inplane.W_table / antiplane (fully_dynamic/test_fd_dense.py).
  usage: test_fd_dense <angle_deg> <mode 0 in-plane | 1 antiplane> <outfile>
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fd_tables.h"

int main(int argc, char **argv)
{
  double ang = argc > 1 ? atof(argv[1])*M_PI/180.0 : 0.5, alpha = 6000.0, beta = 3464.0, rho = 2670.0, dx = 100.0, dt = 0.5*dx/beta;
  int mode = argc > 2 ? atoi(argv[2]) : 0, ns = 6, nr = 5, l, i, j, lags[5] = {1, 2, 5, 20, 60};
  struct fd_segment src, rcv;
  struct fd_pair *p;
  char err[256];
  FILE *fo = fopen(argc > 3 ? argv[3] : "fd_dense.txt", "w");
  double *Ws = malloc(sizeof(double)*ns*nr), *Wn = malloc(sizeof(double)*ns*nr);
  src.X0[0] = 0.5*dx;
  src.X0[1] = 0.0;
  src.s[0] = 1.0;
  src.s[1] = 0.0;
  src.dx = dx;
  src.a = 0.5*dx;
  src.n = ns;
  src.kind = FD_SEG_FAULT;
  /* receiver starts at the source's far end (x = ns dx), rotated by ang */
  rcv.s[0] = cos(ang);
  rcv.s[1] = sin(ang);
  rcv.dx = dx;
  rcv.a = 0.5*dx;
  rcv.n = nr;
  rcv.kind = FD_SEG_FAULT;
  rcv.X0[0] = ns*dx + 0.5*dx*rcv.s[0];
  rcv.X0[1] = 0.5*dx*rcv.s[1];
  p = fd_pair_build(&src, &rcv, mode, dt, 1.0, alpha, beta, rho, 1e-9, 200, 1, 0, err);   /* tight: no tail cut before lag 60 */
  if(!p){
    fprintf(stderr, "build failed: %s\n", err);
    return 1;
  }
  fprintf(fo, "# dense %d L %d\n", p->dense, p->L);
  fd_pair_static_dense(p, Ws, Wn);
  for(i = 0; i < nr; i++)
    for(j = 0; j < ns; j++)
      fprintf(fo, "0 %d %d %.15e %.15e\n", i, j, Ws[i*ns + j], Wn[i*ns + j]);
  for(l = 0; l < 5; l++){
    if(lags[l] > p->L)
      continue;
    fd_pair_W_dense(p, lags[l], Ws, Wn);
    for(i = 0; i < nr; i++)
      for(j = 0; j < ns; j++)
        fprintf(fo, "%d %d %d %.15e %.15e\n", lags[l], i, j, Ws[i*ns + j], Wn[i*ns + j]);
  }
  fclose(fo);
  fd_pair_free(p);
  free(Ws);
  free(Wn);
  printf("test_fd_dense: angle %g deg, mode %d, L %d, written\n", ang*180.0/M_PI, mode, p ? 0 : 0);
  return 0;
}
