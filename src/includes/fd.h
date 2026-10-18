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
  PetscReal dt_restart;		/* TS step after a phase [s] */
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
