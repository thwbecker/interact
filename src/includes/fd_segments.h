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
/* virtual free surface of a 2-D in-plane half plane (fd_plan.md M7 step 2c): the line y = 0,
   x in [xmin, xmax], appended as runs of equal elements (slip and opening segments, global
   index -1); traces xtr on element boundaries; elements of about hnear within wnear of a trace
   or a shallow fault top xsh, hfar elsewhere. Returns the surface elements per component. */
int fd_seglist_add_surface(struct fd_seglist *S, double xmin, double xmax, int ntr, const double *xtr, int nsh, const double *xsh,
                           double hnear, double wnear, double hfar);
#endif
