/* test of fd_detect_segments: two parallel segments in shuffled element
   order, one isolated element, one segment at an angle; exit 0 on PASS */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fd_segments.h"
int main(void)
{
  int n = 0, i, k, fail = 0, order[64], tmp, cnt[4] = {0, 0, 0, 0};
  double cx[64], cy[64], tx[64], ty[64], hl[64], c = cos(0.3), s = sin(0.3);
  char err[256];
  struct fd_seglist *S;
  double bx[64];
  double by[64];
  double btx[64];
  double bty[64];
  double bhl[64];
  int g;
  double ex;
  double ey;
  for(i = 0; i < 20; i++){
    cx[n] = -2000 + 100 * (i + 0.5);
    cy[n] = 0;
    tx[n] = 1;
    ty[n] = 0;
    hl[n] = 50;
    n++;
  }    /* seg A */
  for(i = 0; i < 12; i++){
    cx[n] = 500 + 100 * (i + 0.5);
    cy[n] = 800;
    tx[n] = 1;
    ty[n] = 0;
    hl[n] = 50;
    n++;
  }      /* seg B, parallel, offset */
  cx[n] = 9000;
  cy[n] = 9000;
  tx[n] = 1;
  ty[n] = 0;
  hl[n] = 50;
  n++;                                                     /* isolated */
  for(i = 0; i < 7; i++){
    cx[n] = 100 * (i + 0.5) * c;
    cy[n] = -3000 + 100 * (i + 0.5) * s;
    tx[n] = c;
    ty[n] = s;
    hl[n] = 50;
    n++;
  }  /* seg C, rotated */
  /* shuffle */
  for(i = 0; i < n; i++)
    order[i] = i;
  srand(7);
  for(i = n - 1; i > 0; i--){
    k = rand() % (i + 1);
    tmp = order[i];
    order[i] = order[k];
    order[k] = tmp;
  }
  {
    for(i = 0; i < n; i++){
      bx[i] = cx[order[i]];
      by[i] = cy[order[i]];
      btx[i] = tx[order[i]];
      bty[i] = ty[order[i]];
      bhl[i] = hl[order[i]];
    }
    S = fd_detect_segments(n, bx, by, btx, bty, bhl, 1e-6, err);
    if(!S){
      printf("detection failed: %s\n", err);
      return 1;
    }
    printf("%d elements -> %d segments:", n, S->nseg);
    for(k = 0; k < S->nseg; k++){
      printf(" [n=%d dx=%.0f X0=(%.1f,%.1f) s=(%.3f,%.3f)]", S->seg[k].n, S->seg[k].dx, S->seg[k].X0[0], S->seg[k].X0[1], S->seg[k].s[0], S->seg[k].s[1]);
      if(S->seg[k].n == 20)
        cnt[0]++;
      else if(S->seg[k].n == 12)
        cnt[1]++;
      else if(S->seg[k].n == 1)
        cnt[2]++;
      else if(S->seg[k].n == 7)
        cnt[3]++;
      /* ordering check: element j must sit at X0 + j dx s */
      for(i = 0; i < S->seg[k].n; i++){
        g = S->gidx[k][i];
        ex = S->seg[k].X0[0] + i * S->seg[k].dx * S->seg[k].s[0];
        ey = S->seg[k].X0[1] + i * S->seg[k].dx * S->seg[k].s[1];
        if(hypot(bx[g] - ex, by[g] - ey) > 1e-6)
          fail = 1;
      }
    }
    printf("\n");
  }
  if(!(S->nseg == 4 && cnt[0] == 1 && cnt[1] == 1 && cnt[2] == 1 && cnt[3] == 1))
    fail = 1;
  fd_seglist_free(S);
  printf("test_fd_segments: %s\n", fail ? "FAIL" : "PASS");
  return fail;
}
