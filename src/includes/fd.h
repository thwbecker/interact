/*
  fully dynamic extension of rsf_solve: data hung off rsf->fd
  (src/rsf_solve/rsf_fd.c). NULL unless -fd_mode is given. See
  fully_dynamic/fd_c_design.md for the design and milestones.
*/
#ifndef FD_H
#define FD_H
#include <stdio.h>
#include "fd_fft.h"
#include "fd_tables.h"
#include "fd_segments.h"
#include "fd_march.h"
#include "fd3d.h"

struct fd_vars {
  /* options */
  PetscInt  mode;		/* 0 off (rsf->fd stays NULL), 1 dynamic phases, -1 hand-over exercise only */
  PetscReal vsw;		/* QD to FD switch slip rate [m/s] */
  PetscReal vback;		/* FD to QD return threshold [m/s], held for one window */
  PetscReal tol;		/* window tolerance relative to the static interaction envelope */
  PetscReal et;			/* collocation offset in units of dt_min (1 end of step, >1 Peirce-Siebrits) */
  PetscReal xi;			/* variable step: dt <= xi dc / V_max */
  PetscInt  kmax;		/* largest substep multiple (power of two) */
  PetscReal cfl;		/* beta dt_min / dx */
  PetscInt  lmax;		/* cap on the window in fine-grid lags, 0 = prototype bound */
  PetscReal gate_tol;		/* static gate tolerance */
  PetscReal length_scale;	/* geometry unit in metres */
  /* medium */
  int mode2d;			/* FD_MODE_INPLANE or FD_MODE_ANTIPLANE */
  PetscBool halfplane;
  PetscReal mu,nu,rho,alpha,beta,dt_min;
  PetscInt  max_steps;		/* abort a phase beyond this many steps */
  PetscInt  conv,conv_block;	/* history convolution: 1 partitioned (M6), 0 direct; partition length */
  PetscReal back_eps,back_vmax;	/* quasi-steady hand-back criterion (0: off) */
  PetscInt  coarse;		/* coarse history for aligned long steps (1: on) */
  PetscReal dt_restart;		/* TS step after a phase [s] */
  /* 3-D (fd3d.c): square elements on one plane */
  PetscBool dim3, mirror3;
  PetscInt nl;			/* 1: second-order time marching (Noda et al. 2020) */
  PetscInt kmax3;		/* 3-D: 1 allows steps of several substeps (experimental) */
  PetscInt layer_w, layer_m, layer_margin;	/* 3-D half space: free-surface layer (half-width and margin in fault elements, refinement) */
  PetscInt surface_opt, surf_ext, surf_m;	/* 3-D half space: 1 forces the full virtual surface; its margin in fault elements; refinement */
  PetscBool surf3;			/* the 3-D marcher uses the full virtual surface */
  PetscInt surf_grade;		/* 1 (default): surface rows fine next to the trace, widening outwards */
  PetscInt surf_extend, next3;	/* extension rows above the surface (-1: 1 if the fault reaches the surface, else 0) */
  PetscBool surf_m_default;		/* -fd_surface_m not given: refine to keep cS dt/hs >= 2 */
  struct fd3d_tab *tab3;
  struct fd3d *f3;
  int nx3, nz3;
  struct fd3d *f3b;			/* two parallel planes (stepovers): the second plane's marcher */
  struct fd3d_pair *pair3;		/* and their coupling; NULL for one plane */
  int n3tot;				/* marcher elements over both planes (2 nx nz) */
  Vec v3a, v3b;			/* work vectors for the exact static operator */
  VecScatter sc3, scX, scV;		/* MPI: full copies on every rank of v3b, the state, the slip rates */
  Vec v3all, Xall, Vall;
  /* segments, sources (real plus images), tables, marcher system */
  struct fd_seglist *segs;
  PetscInt nsrc;
  struct fd_segment *srcseg;
  int *src_of,*src_rev;
  struct fd_pair **pairs;	/* pairs[source*nseg + receiver] */
  struct fd_sys *sys;
  PetscInt *g_of_s;		/* global element index of segment-order index */
  int M,L;			/* padded FFT length, longest window */
  double nbytes;
  /* stations written during phases */
  PetscInt nstat,stat[16];
  int sstat[16];
  FILE *fstat[16];
  /* accepted-step record for the history seeding */
  PetscInt nrec,irec,nrec_used;
  PetscReal *rec_t;
  double *rec_v;		/* nrec x N, segment order */
  /* state */
  PetscInt  nphase;		/* dynamic phases run so far */
  PetscInt  nflag;		/* accepted steps above vsw seen by the monitor */
  PetscReal t_last_flag;	/* time of the last one [s] */
  PetscBool restart_done;	/* mode -1: the TSRestartStep exercise has been run */
};
#endif
