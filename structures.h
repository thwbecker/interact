/* 
   DEFINITION OF THE MEDIUM AND FAULT STRUCTURES
   and some other stuff, like the two fault group 
   $Id: structures.h,v 1.51 2011/01/09 02:02:43 becker Exp $
*/

/*

  fault group structure holds all non-local information that interact
  has about groups of fault patches
  
*/
struct fgrp{
  my_boolean active; /* if set, at least one of the patches of this group
		   is slipping */
  float moment, sarea, cent[3], 
    slip[6];/* temporary moment release slip,
	       moment, and slipping area on fault
	       groups, slip in [strike,dip,normal, x,y,z]
	       
	    */
  float strike[3], dip[3], normal[3]; /* mean orientation */
};

/*

  geometrical quantities about a fault group only used to determine
  fault[].pos[] in interact and to convert from patch format to group
  format after that, the memory is freed again */

struct geog{
  int nrflt;
  /* center of mass (average location of mid point) and min,max and
    range in strike and dip direction */
  COMP_PRECISION center[3],
    pmin[2],pmax[2],prange[2];
  // average vectors in strike and dip direction
  // (normalized)
  COMP_PRECISION strike_vec[3],dip_vec[3];
};


/* 

   MEDIUM STRUCTURE 
   
   holds all general information that is not fault related

*/
struct med{
  /* nr of patches */
  int nrflt;
  /* nr of groups */
  int nrgrp;
  // fault group structures
  struct fgrp *fault_group;
  /* 
     interaction matrix 
     and switches
  */
  struct icm **i;// is I matrix is in memory, this holds the values
  I_MATRIX_PREC *val; /* if I or A matrices are in sparse storage, 
			 this holds the values */
  I_MATRIX_PREC imax,imean,/* max and mean of interaction 
			      matrix entries, and cutoff 
			      value for sparse storage */
    i_mat_cutoff;
  int nmat;// integer with matrix dimensions for files
  size_t i_matrix_prec_size; // precision of matrix
  unsigned int *is1, *is2;// pointers for sparse storage, if wanted
  my_boolean suppress_interactions;/* no interaction but
				   within the same group
				   (self-effect) and effect
				   of slip on magic group
				   MASTER_FAULT_GROUP
				*/
  // note that internally, numbering starts with 0
  // such that the first inputted patch is number 0 and not 1
#define MASTER_FAULT_GROUP 0
  /* 
     boundary condition solution items, 
     first for the unconstrained (SVD) part and
     then *_con for the non-negative, constrained part
     of the equations
     
  */
  int nreq,nreq_con;/* number of active equations 
		       that have to 
		       be satisfied */
  int naflt,naflt_con; /* number of active faults 
			  <= nreq */
  int *nameaf,*nameaf_con; /* names of active faults in 
			      these eqs */
  my_boolean *sma,*sma_con;/* active slide modes of the
			   active faults in order of their 
			   activation
			*/
  short int solver_mode; // type of solver used for A x = b system
  A_MATRIX_PREC *b,*xsol,*b_con,*xsol_con; /* b (rhs) and x (solution) 
					      vectors of system */

  // switch that activates all patches along a fault
  // group if one patch is active
  my_boolean whole_fault_activations;
  // switch that deactivates the whole fault if one
  // patch goes inactive
  my_boolean whole_fault_deactivations;
  //
  // if single patches slip, there stess drop can trigger other patches
  // if this flag is set, will keep triggering patches slipping until no more
  // activations take place (default)
  my_boolean keep_slipping;
  //
  // try to restart given a events.bin file
  my_boolean attempt_restart;

  // count the number of active groups for whole fault 
  // mode search for patches that should be slipping
  int nr_active_groups;
  /*
    nr of iteration at fixed time, total iterations, and 
    maximum number of solver iterations that were actually
    needed at one point
  */
  int iter, total_iter, max_iter_realized;
  /* 
     how to compute group local locations: 0: old 1: new
  */
  int group_geom_mode;
  /* 
     geometrical boundaries for fault geometry 
  */
  COMP_PRECISION xmin[3],xmax[3];
  /* fault length and width extrema */
  COMP_PRECISION lmean,wmean;
  /* geometrical boundaries for stress output, plotting */
  COMP_PRECISION pxmin[3],pxmax[3];
  my_boolean *ok;/* this is only used to determine proper points 
		 if we are calculating output along a 
		 fault plane for stress field output in single
		 step operational mode */
  /* number of columns and rows and slices 
     in stress and displacement field output matrices */
  int n[3];
  /* stress and displacement field output matrices */
  float *s,*u;
  // stress and displacement input locations and nr of input locations
  float *xoloc;int olocnr;
  /* 
     timing and printing issues 
  */
  COMP_PRECISION time,dt0,dt,stop_time,time_to_failure;
  COMP_PRECISION print_interval,slip_line_dt,
    old_moment_time,slip_line_time;
  int nr_timesteps;
  /*
    material properties that don't change within the medium
  */
  COMP_PRECISION pressure,cohesion,min_stress_drop;
  /* 
     total seismic moment release 
  */
  COMP_PRECISION total_moment;
  /* 
     files for output 
  */
  FILE *(*flt_stress_out),*events_out,
    *cevents_out,*(*slip_line_out),*stress_stat_out;
  /* 
     files for input, I matrix 
  */
  FILE *i_mat_in;
  // some switches
  my_boolean flt_stress_init,events_init,save_amat,slip_line_init,
    print_bulk_fields,use_slip_files,use_sparse_storage,
    use_old_imat,use_old_amat,save_imat,check_for_interaction_feedback,
    read_oloc_from_file,suppress_nan_output,geomview_output,
    twod_approx_is_plane_stress,print_plane_coord,
    variable_time_step,debug,no_interactions;
  /* calculation mode and state switches */
  int op_mode,op_state,stress_state_init;
  /* 
     internal checks, initialization flags
  */
  my_boolean geometry_init,bc_init,pos_init,int_mat_init,
    bulk_field_init,read_int_mat_from_file,
    moment_release_init;
  //
  // maximum number of fault groups that will lead
  // to individual output files
  //
  short int max_nr_flt_files;
  // how many modes are active?
  short int nr_flt_mode;
  
  // matrix filenames
  char mfname[STRLEN],hfname[STRLEN];
  /*

    re-used arrays for spatial correlation calculations
    spcorr_interval is the time interval between updates
  */
  float *cr_w1,*cr_w2,*cr_w3;
  COMP_PRECISION *cr_xr,*cr_r,spcorr_interval;
  int *cr_nr;

  /* for SVD */
  COMP_PRECISION wcutoff;
  /* 
     for PGPLOTTING X WINDOWS OUTPUT 
  */
#ifdef USE_PGPLOT
#define MAX_NR_X_WINDOWS 3
  /* plotting issues */
  my_boolean nr_xwindows,moment_array_init;
  int lw_def,*active_flt_list,nr_active_flt;
  int x_window[MAX_NR_X_WINDOWS];
  float tloc[3];
  /* for plotting of moment release on fault one */
  float *momrel;
  int grp0_n,grp0_m;
  /* timing issues */
  COMP_PRECISION time_tic_interval,l_time_tic_interval,
    x_plot_time,x_plot_interval;
  COMP_PRECISION x_scroll_time,x_scroll_interval,
    x_scroll_inc;
#endif
  COMP_PRECISION nan;		/* remember to initialize  */
  /*  */
#ifdef USE_PETSC
  PetscMPIInt comm_size, comm_rank;
  PetscInt    rs, re, rn;
  Mat         pA;
  Vec         pb;
  PetscInt    *indices;
#else
  unsigned int comm_size,comm_rank;
#endif
};

/*

  FAULT STRUCTURE
  this structure should hold all variables that are
  associated with an individual fault patch


  WARNING:
  don't forget to change fltswap and the like if you
  add components

*/

struct flt{
  COMP_PRECISION u[3];/* displacements on fault in strike, 
			 dip, and tensile direction */
  COMP_PRECISION s[3];/* strike, normal and 
			 dip stress */
  COMP_PRECISION sinc[3]; /* assuming that the background
			     stress increases linearly,
			     these are the increments in
			     strike, dip, and normal 
			     direction,s=a*t+b
			  */
  COMP_PRECISION x[3]; /* 
			  location of center of patch or point
			  source or centroid location for 
			  triangles
		       */
#ifdef ALLOW_NON_3DQUAD_GEOM
  MODE_TYPE type;/* 
		    element types:
		    0: rectangular element
		    1: point-source 
		    2: triangular element 
		 */
  COMP_PRECISION *xt; /* coordinates of the three nodes of
			 the triangular element */
#endif
  float pos[2]; /* position in the group of patches
		   assuming linear fault */
  COMP_PRECISION l, w, area;  /* 
			   fault HALF width and total area for rectangular patches, 
			   length is area/(4w)
			   otherwise w and l will hold the area of the triangle

			   
			*/
  float strike, dip;    /* 
			   fault orientation in degrees, 
			   strike is in degrees clockwise 
			   from North, dip is in degrees, too, the
			   angle from horizontal so that dip=90 
			   means vertical 

			*/
  COMP_PRECISION cos_alpha,/* cos and sin alpha, which is
			      angle counterclockwise from East,
			      ie. 90 - strike */
    sin_alpha;
  COMP_PRECISION normal[3],
    t_strike[3],t_dip[3]; /* normal and tangential 
			     (strike and dip direction) 
			     unit vectors */
  float mu_d,mu_s; /* 
		      dynamic and static friction 
		      coefficients
		   */
  COMP_PRECISION f_initial,taud;// initial f ratio and stress drop
  COMP_PRECISION cf[2];// coulomb correction for normal stress
#ifdef LATENCY
  float last_activation_time;/* last time of slip for 
				latency operation mode */
#endif
  unsigned int group; /* fault patch can have a 
			 group assigned */
  MODE_TYPE mode[3];/* 
		     activational mode in each direction
		  */
  my_boolean active; /* fault is active in the current
		     solver */
};

//  structure needed by plotting routines
struct pa{
  float x[5],y[5],z[5];
};

/* 
   interaction coefficient structure
   
   this is the building block of the interaction matrix
   storage, if kept in memory
*/
struct icm{
  I_MATRIX_PREC s[NRMODE_DEF][3];
};
/*
  
  for sorting integer lists (permutations)

 */
struct slist{
  COMP_PRECISION x;
  int i;
};


// typedefs for segseg

#define	SEGSEG_DIM 2               /* Dimension of points */
typedef	int tPointi[SEGSEG_DIM];   /* Type integer point */
typedef	double tPointd[SEGSEG_DIM];   /* Type double point */


/* the structures for blockinvert type programs */
#include "blockinvert_structures.h"
