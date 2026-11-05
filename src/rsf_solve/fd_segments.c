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
      S->seg[S->nseg].pc = NULL;
      S->seg[S->nseg].pt = NULL;
      S->seg[S->nseg].pa = NULL;
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
  for(k = 0; k < s->nseg; k++){
    free(s->gidx[k]);
    free(s->seg[k].pc);
    free(s->seg[k].pt);
    free(s->seg[k].pa);
  }
  free(s->gidx);
  free(s->seg);
  free(s);
}

void fd_seg_element(const struct fd_segment *g, int i, double X[2], double t[2], double *a)
{
  if(g->pc){
    X[0] = g->pc[2*i];
    X[1] = g->pc[2*i + 1];
    t[0] = g->pt[2*i];
    t[1] = g->pt[2*i + 1];
    *a = g->pa[i];
    return;
  }
  X[0] = g->X0[0] + i*g->dx*g->s[0];
  X[1] = g->X0[1] + i*g->dx*g->s[1];
  t[0] = g->s[0];
  t[1] = g->s[1];
  *a = g->a;
}

/* end points of a fault segment: e = 0 start of element 0, e = 1 end of the last element */
static void seg_end(const struct fd_segment *g, int e, double P[2])
{
  double X[2], t[2], a;
  int i;
  i = e ? g->n - 1 : 0;
  fd_seg_element(g, i, X, t, &a);
  P[0] = X[0] + (e ? a : -a)*t[0];
  P[1] = X[1] + (e ? a : -a)*t[1];
}

/* append the elements of segment g (reversed if rev) to the polyline arrays at position m */
static void poly_append(const struct fd_segment *g, const int *gi, int rev, double *pc, double *pt, double *pa, int *gx, int m)
{
  int i, q;
  double X[2], t[2], a;
  for(q = 0; q < g->n; q++){
    i = rev ? g->n - 1 - q : q;
    fd_seg_element(g, i, X, t, &a);
    pc[2*(m + q)] = X[0];
    pc[2*(m + q) + 1] = X[1];
    pt[2*(m + q)] = t[0];                       /* the element keeps its tangent (slip sign of the model) */
    pt[2*(m + q) + 1] = t[1];
    pa[m + q] = a;
    gx[m + q] = gi[i];
  }
}

int fd_seglist_merge_polylines(struct fd_seglist *S, double tol, int minrun)
{
  int n0, k, j, e, f, cur, nchain, len, m, *used, *order, *revs, done, found, best, bestrev, shortrun;
  double P[2], Q[2], d, h;
  struct fd_segment *ns;
  int **ng;
  n0 = S->nseg;
  used = (int *)calloc(n0, sizeof(int));
  order = (int *)malloc(sizeof(int)*n0);
  revs = (int *)malloc(sizeof(int)*n0);
  ns = (struct fd_segment *)malloc(sizeof(struct fd_segment)*n0);
  ng = (int **)malloc(sizeof(int *)*n0);
  nchain = 0;
  for(k = 0; k < n0; k++){
    if(used[k] || S->seg[k].kind != FD_SEG_FAULT)
      continue;
    /* grow a chain from k forwards and backwards */
    len = 1;
    order[0] = k;
    revs[0] = 0;
    used[k] = 1;
    for(f = 0; f < 2; f++){                       /* f = 0: after the end, 1: before the start */
      done = 0;
      while(!done){
        cur = f ? order[0] : order[len - 1];
        if(f)
          seg_end(S->seg + cur, revs[0] ? 1 : 0, P);          /* the chain's free start */
        else
          seg_end(S->seg + cur, revs[len - 1] ? 0 : 1, P);    /* the chain's free end */
        h = (S->seg[cur].pa ? S->seg[cur].pa[0] : S->seg[cur].a)*2.0;
        found = 0;
        best = -1;
        bestrev = 0;
        for(j = 0; j < n0 && !found; j++){
          if(used[j] || S->seg[j].kind != FD_SEG_FAULT)
            continue;
          for(e = 0; e < 2; e++){
            seg_end(S->seg + j, e, Q);
            d = hypot(P[0] - Q[0], P[1] - Q[1]);
            if(d < tol*h){
              best = j;
              /* appending after the end: the join must be j's start (e = 0) unless reversed;
                 prepending before the start: the join must be j's end (e = 1) unless reversed */
              bestrev = f ? (e == 0) : (e == 1);
              found = 1;
              break;
            }
          }
        }
        if(!found){
          done = 1;
          continue;
        }
        used[best] = 1;
        if(f){
          memmove(order + 1, order, sizeof(int)*len);
          memmove(revs + 1, revs, sizeof(int)*len);
          order[0] = best;
          revs[0] = bestrev;
        }else{
          order[len] = best;
          revs[len] = bestrev;
        }
        len++;
      }
    }
    /* a chain of long straight runs (a bend between straight segments) stays as it is: their
       Toeplitz self tables are cheaper; chains with short runs (curved faults) are merged */
    shortrun = 0;
    for(j = 0; j < len; j++)
      if(S->seg[order[j]].n < minrun)
        shortrun = 1;
    if(len == 1 || !shortrun){
      for(j = 0; j < len; j++){
        ns[nchain] = S->seg[order[j]];
        ng[nchain] = S->gidx[order[j]];
        S->gidx[order[j]] = NULL;
        nchain++;
      }
      continue;
    }
    m = 0;
    for(j = 0; j < len; j++)
      m += S->seg[order[j]].n;
    memset(ns + nchain, 0, sizeof(struct fd_segment));
    ns[nchain].kind = FD_SEG_FAULT;
    ns[nchain].n = m;
    ns[nchain].pc = (double *)malloc(sizeof(double)*2*m);
    ns[nchain].pt = (double *)malloc(sizeof(double)*2*m);
    ns[nchain].pa = (double *)malloc(sizeof(double)*m);
    ng[nchain] = (int *)malloc(sizeof(int)*m);
    m = 0;
    for(j = 0; j < len; j++){
      poly_append(S->seg + order[j], S->gidx[order[j]], revs[j], ns[nchain].pc, ns[nchain].pt, ns[nchain].pa, ng[nchain], m);
      m += S->seg[order[j]].n;
    }
    /* the straight-run fields describe element 0 (for messages; geometry comes from the arrays) */
    ns[nchain].X0[0] = ns[nchain].pc[0];
    ns[nchain].X0[1] = ns[nchain].pc[1];
    ns[nchain].s[0] = ns[nchain].pt[0];
    ns[nchain].s[1] = ns[nchain].pt[1];
    ns[nchain].a = ns[nchain].pa[0];
    ns[nchain].dx = 2.0*ns[nchain].pa[0];
    nchain++;
  }
  /* non-fault segments keep their place after the faults */
  for(k = 0; k < n0; k++)
    if(S->seg[k].kind != FD_SEG_FAULT){
      ns[nchain] = S->seg[k];
      ng[nchain] = S->gidx[k];
      S->gidx[k] = NULL;
      nchain++;
    }
  for(k = 0; k < n0; k++){
    if(S->gidx[k]){                               /* merged: the old pieces go */
      free(S->gidx[k]);
      free(S->seg[k].pc);
      free(S->seg[k].pt);
      free(S->seg[k].pa);
    }
  }
  free(S->seg);
  free(S->gidx);
  S->seg = ns;
  S->gidx = ng;
  S->nseg = nchain;
  free(used);
  free(order);
  free(revs);
  return nchain;
}
