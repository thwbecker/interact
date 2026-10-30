/*
  test of the dense tables with virtual free-surface elements (fd_plan.md M7 step 2): a fault
  of nf elements of 500 m dipping 30 degrees with its top 500 m below a surface line of ns
  elements of 250 m (in-plane); prints, for every (source, receiver) pair of the three
  segments (fault, surface slip, surface opening), the static block and the increments of
  lags 1..L (FD_CT_OFF should be set), compared by fully_dynamic/test_fd_surf.py with
  densesys.py of the prototype.
  usage: test_fd_surf nf ns L > tables.txt
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fd_tables.h"

int main(int argc, char **argv)
{
  struct fd_segment seg[3];
  struct fd_pair *p;
  int nf, ns, L, a, b, i, j, l;
  double alpha, beta, rho, dt, dip, *Ws, *Wn, *Ks, *Kn;
  char err[256];
  nf = argc > 1 ? atoi(argv[1]) : 4;
  ns = argc > 2 ? atoi(argv[2]) : 8;
  L = argc > 3 ? atoi(argv[3]) : 30;
  alpha = 6000.0;
  beta = 3464.0;
  rho = 2670.0;
  dt = 0.5*250.0/beta;
  dip = 30.0*M_PI/180.0;
  /* fault: from (0, -500) down-dip towards +x, elements of 500 m */
  seg[0].kind = FD_SEG_FAULT;
  seg[0].n = nf;
  seg[0].dx = 500.0;
  seg[0].a = 250.0;
  seg[0].s[0] = cos(dip);
  seg[0].s[1] = -sin(dip);
  seg[0].X0[0] = 250.0*cos(dip);
  seg[0].X0[1] = -500.0 - 250.0*sin(dip);
  /* surface line y = 0 from x = -1000, elements of 250 m, twice (slip and opening) */
  for(a = 1; a <= 2; a++){
    seg[a].kind = (a == 1) ? FD_SEG_SURF_SLIP : FD_SEG_SURF_OPEN;
    seg[a].n = ns;
    seg[a].dx = 250.0;
    seg[a].a = 125.0;
    seg[a].s[0] = 1.0;
    seg[a].s[1] = 0.0;
    seg[a].X0[0] = -1000.0 + 125.0;
    seg[a].X0[1] = 0.0;
  }
  printf("# dt %.17e\n", dt);
  for(a = 0; a < 3; a++)
    for(b = 0; b < 3; b++){
      p = fd_pair_build_dense(seg + a, seg + b, FD_MODE_INPLANE, dt, 1.0, alpha, beta, rho, 1e-12, L + 1, 1, 0, err);
      if(!p){
        fprintf(stderr, "pair %d %d: %s\n", a, b, err);
        return 1;
      }
      Ws = malloc(sizeof(double)*seg[a].n*seg[b].n);
      Wn = malloc(sizeof(double)*seg[a].n*seg[b].n);
      Ks = malloc(sizeof(double)*seg[a].n*seg[b].n);
      Kn = malloc(sizeof(double)*seg[a].n*seg[b].n);
      fd_pair_static_dense(p, Ks, Kn);
      for(i = 0; i < seg[b].n; i++)
        for(j = 0; j < seg[a].n; j++)
          printf("K %d %d %d %d %.17e %.17e\n", a, b, i, j, Ks[i*seg[a].n + j], Kn[i*seg[a].n + j]);
      for(l = 1; l <= L; l++){
        if(l > p->L){
          for(i = 0; i < seg[b].n*seg[a].n; i++){
            Ws[i] = 0.0;
            Wn[i] = 0.0;
          }
        }else
          fd_pair_W_dense(p, l, Ws, Wn);
        for(i = 0; i < seg[b].n; i++)
          for(j = 0; j < seg[a].n; j++)
            printf("W %d %d %d %d %d %.17e %.17e\n", a, b, l, i, j, Ws[i*seg[a].n + j], Wn[i*seg[a].n + j]);
      }
      free(Ws);
      free(Wn);
      free(Ks);
      free(Kn);
      fd_pair_free(p);
    }
  fd_tables_open_cache_clear();
  return 0;
}
