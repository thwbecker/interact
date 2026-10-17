/* detection of straight runs of equal, collinear 2-D elements (fd_segments.c) */
#ifndef FD_SEGMENTS_H
#define FD_SEGMENTS_H
#include "fd_tables.h"
struct fd_seglist {
  int nseg;
  struct fd_segment *seg;	/* geometry of each run */
  int **gidx;			/* gidx[k][j]: global element index of element j of run k */
};
/* n elements with centres (cx, cy), unit tangents (tx, ty) and half-lengths hl;
   elements are grouped into runs with identical tangent and half-length whose
   centres lie on one line at spacing 2 hl, in order of increasing position
   along the tangent. rtol is relative to the element length. Returns NULL
   if an element cannot be placed (errmsg set). */
struct fd_seglist *fd_detect_segments(int n, const double *cx, const double *cy, const double *tx, const double *ty,
                                      const double *hl, double rtol, char *errmsg);
void fd_seglist_free(struct fd_seglist *s);
#endif
