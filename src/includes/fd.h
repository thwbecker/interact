/*
  fully dynamic extension of rsf_solve: data hung off rsf->fd
  (src/rsf_solve/rsf_fd.c). NULL unless -fd_mode is given. See
  fully_dynamic/fd_c_design.md for the design and milestones.
*/
#ifndef FD_H
#define FD_H
#include "fd_fft.h"

struct fd_vars {
  /* options */
  PetscInt  mode;		/* 0 off (rsf->fd stays NULL), 1 dynamic phases, -1 hand-over exercise only */
  PetscReal vsw;		/* QD to FD switch slip rate [m/s] */
  PetscReal vback;		/* FD to QD return threshold [m/s], held for one window */
  PetscReal tol;		/* window tolerance relative to the static interaction envelope */
  PetscReal et;			/* collocation offset in units of dt_min (1 end of step, >1 Peirce-Siebrits) */
  PetscReal xi;			/* variable step: dt <= xi dc / V_max */
  PetscInt  kmax;		/* largest substep multiple (power of two) */
  /* state */
  PetscInt  nphase;		/* dynamic phases run so far */
  PetscInt  nflag;		/* accepted steps above vsw seen by the monitor */
  PetscReal t_last_flag;	/* time of the last one [s] */
  PetscBool restart_done;	/* mode -1: the TSRestartStep exercise has been run */
};
#endif
