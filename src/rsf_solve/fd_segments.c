/* grouping of 2-D elements into Toeplitz-capable segments, see fd_segments.h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fd_segments.h"

struct item {
  int g;
  double along, across;
};
static int cmp_along(const void *a, const void *b)
{
  double d = ((const struct item *)a)->along - ((const struct item *)b)->along;
  return d < 0 ? -1 : (d > 0 ? 1 : 0);
}

struct fd_seglist *fd_detect_segments(int n, const double *cx, const double *cy, const double *tx, const double *ty,
    const double *hl, double rtol, char *errmsg)
{
  struct fd_seglist *S = (struct fd_seglist *)calloc(1, sizeof(*S));
  int *done = (int *)calloc(n, sizeof(int)), i, j, k, cnt, cap = 8;
  struct item *it = (struct item *)malloc(sizeof(struct item) * n);
  double dx, tol, nx, ny;
  int lo;
  int hi;
  if(errmsg)
    errmsg[0] = 0;
  S->seg = (struct fd_segment *)malloc(sizeof(struct fd_segment) * cap);
  S->gidx = (int **)malloc(sizeof(int *) * cap);
  for(i = 0; i < n; i++){
    if(done[i])
      continue;
    dx = 2.0 * hl[i];
    tol = rtol * dx;
    nx = -ty[i];
    ny = tx[i];
    /* candidates: same tangent and half-length, on the same line */
    cnt = 0;
    for(j = 0; j < n; j++){
      if(done[j])
        continue;
      if(fabs(hl[j] - hl[i]) > tol || (tx[j] * tx[i] + ty[j] * ty[i]) < 1.0 - 1e-10)
        continue;
      if(fabs((cx[j] - cx[i]) * nx + (cy[j] - cy[i]) * ny) > tol)
        continue;
      it[cnt].g = j;
      it[cnt].along = (cx[j] - cx[i]) * tx[i] + (cy[j] - cy[i]) * ty[i];
      cnt++;
    }
    qsort(it, cnt, sizeof(struct item), cmp_along);
    /* the run containing i: contiguous spacing dx around it */
    for(k = 0; k < cnt && it[k].g != i; k++){ /* locate i */
    }
    {
      lo = k;
      hi = k;
      while(lo > 0 && fabs(it[lo].along - it[lo - 1].along - dx) < tol)
        lo--;
      while(hi < cnt - 1 && fabs(it[hi + 1].along - it[hi].along - dx) < tol)
        hi++;
      if(S->nseg == cap){
        cap *= 2;
        S->seg = (struct fd_segment *)realloc(S->seg, sizeof(struct fd_segment) * cap);
        S->gidx = (int **)realloc(S->gidx, sizeof(int *) * cap);
      }
      S->seg[S->nseg].n = hi - lo + 1;
      S->seg[S->nseg].kind = FD_SEG_FAULT;
      S->seg[S->nseg].dx = dx;
      S->seg[S->nseg].a = hl[i];
      S->seg[S->nseg].s[0] = tx[i];
      S->seg[S->nseg].s[1] = ty[i];
      S->seg[S->nseg].X0[0] = cx[it[lo].g];
      S->seg[S->nseg].X0[1] = cy[it[lo].g];
      S->gidx[S->nseg] = (int *)malloc(sizeof(int) * (hi - lo + 1));
      for(j = lo; j <= hi; j++){
        S->gidx[S->nseg][j - lo] = it[j].g;
        done[it[j].g] = 1;
      }
      S->nseg++;
    }
  }
  free(done);
  free(it);
  return S;
}

void fd_seglist_free(struct fd_seglist *s)
{
  int k;
  if(!s)
    return;
  for(k = 0; k < s->nseg; k++)
    free(s->gidx[k]);
  free(s->gidx);
  free(s->seg);
  free(s);
}
