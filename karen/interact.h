/*
  HEADER FILE FOR INTERACT

  and various other routines of the interact suite
  
  includes function declarations, structure definitions,
  and constants

  $Id: interact.h,v 2.13 2000/08/30 01:33:52 becker Exp becker $ 
*/
/* 
   system headers 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
/* other headers */
#ifdef USE_DOUBLE_PRECISION
#include "precision_double.h"
#else
#include "precision_single.h"
#endif
#ifdef USE_PGPLOT
#include "cpgplot.h"
#endif
/* for calls to external FORTRAN routines under SGI */
#ifdef SGI_SUBROUTINE_CONVENTION
#define dc3d dc3d_
#define mydc3d mydc3d_
#endif
/* 
   programming conventions and macros
*/
typedef short int BOOLEAN;
#define TRUE 1
#define FALSE 0
static COMP_PRECISION sqrarg;
#define SQUARE(a) ((sqrarg=(COMP_PRECISION)(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static int isqrarg;
#define ISQUARE(a) ((isqrarg=(int)(a)) == 0 ? 0 : isqrarg*isqrarg)
#define MIN(x, y) (((x)<(y))?(x):(y))
#define MAX(x, y) (((x)>(y))?(x):(y))
#define SIGN(x) (((x)>=0)?(1.0):(-1.0))
#define STRLEN 200
/* real constants */
#define DEG2RAD 0.0174532925199432958
#define NAN 1e20
/* 
   FILENAMES AND OUTPUT OPTIONS 
*/
#define GEOMETRY_FILE "geom.in"
#define BC_FILE "bc.in"
#define FAULT_PROP_FILE "fp.in"
#define FAULT_STRESS_INIT_FILE "fsi.in"
#define STRESS_OUT_FILE "stress.out"
#define STRESS_HDR_FILE "stress.hdr"
#define FAULT_STRESS_OUT_FILE "flt_stress.out"
#define DISP_OUT_FILE "disp.out"
#define DISP_HDR_FILE "disp.hdr"
#define FAULT_DATA_PREFIX "flt"
#define EVENT_FILE_BINARY  "events.bin"
#define EVENT_FILE_ASCII  "events.dat"
#define CEVENT_FILE "cevents.dat"
#define INTERACTION_MATRIX_FILE "i.dat"
#define A_MATRIX_FILE "a.dat"
#define A_CON_MATRIX_FILE "ac.dat"
#define DEBUG_A_MATRIX_ASCII_OUT "a.asc"
#define DEBUG_B_VECTOR_ASCII_OUT "b.asc"
#define DEBUG_X_VECTOR_ASCII_OUT "x.asc"
#define DEBUG_A_CON_MATRIX_ASCII_OUT "ac.asc"
#define DEBUG_B_CON_VECTOR_ASCII_OUT "bc.asc"
#define DEBUG_X_CON_VECTOR_ASCII_OUT "xc.asc"
#define SLIP_LINE_PREFIX "slipline"
#define STRESS_RELATION_FILE "smlin.in"
/* 
   set this to the default 
   maximum number of fault patches
   that will result in separate output files 
*/
#define MAX_NR_FLT_FILES_DEF 50
/*
  default operational mode switches
*/
#define READ_FAULT_PROPERTIES_DEF TRUE
#define READ_INITIAL_FAULT_STRESS_DEF TRUE
#define SUPPRESS_INTERACTIONS_DEF FALSE
#define WHOLE_FAULT_MODE_DEF FALSE
#define NO_OPENING_MODES_DEF FALSE
#define READ_STRESS_RELATION_FACTORS_DEF TRUE

/* 
   indices used in programs 
*/
#define STRIKE 0
#define DIP 1
#define NORMAL 2
#define X 0
#define Y 1
#define Z 2
/* operational modes for interact */
#define ONE_STEP 1
#define SIMULATE_LOADING 2
#define SIMULATE_LOADING_AND_PLOT 3
/* 
   output modes for divide_patch as used by randomflt  
   and makefault
*/
#define PATCH_OUT_MODE 0
#define GEOMVIEW_MODE 1
#define PSXYZ_MODE 2
#define BC_OUT_MODE 3
/* activation modes for the faults */
#define INACTIVE 0
#define STRIKE_SLIP 10
#define STRIKE_SLIP_LEFTLATERAL 11
#define STRIKE_SLIP_RIGHTLATERAL 12
#define DIP_SLIP 20
#define DIP_SLIP_UPWARD 21
#define DIP_SLIP_DOWNWARD 22
#define NORMAL_SLIP 30
#define NORMAL_SLIP_OUTWARD 31
#define NORMAL_SLIP_INWARD 32
#define MAXSDIR_SLIP 40
/* 
   the second and third mode will lead to a correction of the target 
   stress due to normal stress changes during slip events
*/
#define ACTIVATED 1
#define ACTIVATED_AND_POSITIVE_NORMAL_CORRECTION 2
#define ACTIVATED_AND_NEGATIVE_NORMAL_CORRECTION 3


/* operational modes for obtaining the interaction
   coefficients */
#define I_MAT_IN_MEMORY 0
#define I_MAT_ON_FILE 1
#define CALC_I_COEFF_NOW 2
// max size of I matrix in MB
#define IMAT_SIZE_LIM 300.0
/* 
   minimum range for xmin/xmax as used to determine the 
   size of the X plotting window
*/
#define MIN_GEOM_RANGE 5.0
/* 
   offset operators for the interaction matrix 
   and 3d fields 
*/

#ifdef NO_OPENING_MODES
// opening modes are suppressed, this assumes that strike and dip are 0 and 1 logically
#define NRMODE 2
#define NRMODE3 6
#else
// all modes are possible
#define NRMODE 3
#define NRMODE3 9
#endif
#define POSI(i, j, k, l) ((NRMODE3*((i)*medium->nrflt+(j)))+(3*(k))+(l))
#define POSS(i, j, k, l) (((k)*(medium->n[Y]*medium->n[X])+((i)*medium->n[Y])+(j))*6+(l))
#define POSU(i, j, k, l) (((k)*(medium->n[Y]*medium->n[X])+((i)*medium->n[Y])+(j))*3+(l))
/* 
   DEFINITION OF THE MEDIUM AND FAULT STRUCTURES
*/

/*
  fault group structure
  holds all non-local information about 
  groups of fault patches

 */
struct fgrp{
  BOOLEAN active; /* if set, at least one of the patches 
		     of this group is slipping 
		  */
  float moment, slip[3];/* temporary moment release 
			   on fault groups */
};



/* 
   medium, holds all general information that is not
   fault related 
*/
struct med{
  /* nr of patches, nrflt times 3 */
  int nrflt,nrflt3;
  /* nr of groups */
  int nrgrp;
  // fault group structures
  struct fgrp *fault_group;
  /* 
     interaction matrix 
     and switches
  */
  I_MATRIX_PREC *i,imax,imean;/* matrix itself, max and mean 
				 of entries */
  size_t i_matrix_prec_size; // precision of matrix
  BOOLEAN suppress_interactions;/* no interaction but
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
  BOOLEAN *sma,*sma_con;/* active slide modes of the
			   active faults in order of their 
			   activation
			*/
  A_MATRIX_PREC *b,*b_con; /* b or x vector of system */

  // switch that activates all patches along a fault
  // group, if one patch is active
  BOOLEAN whole_fault_activations;
  // count the number of active groups for whole fault 
  // mode search for patches that should be slipping
  int nr_active_groups;

  /* 
     geometrical boundaries for fault geometry 
  */
  COMP_PRECISION xmin[3],xmax[3];
  /* fault length and width extrema */
  COMP_PRECISION lmin,lmax,lmean,wmin,wmax,wmean;
  /* geometrical boundaries for stress output, plotting */
  COMP_PRECISION pxmin[3],pxmax[3];
  /* number of columns and rows and slices 
     in stress and displcement field output matrices */
  int n[3];
  /* stress and displacement field output matrices */
  float *s,*u;
  /* timing and printing issues */
  COMP_PRECISION time,dt0,dt,stop_time;
  COMP_PRECISION print_interval,slip_line_dt,
    old_moment_time,slip_line_time;
  /*
    material properties that don't change within the medium
  */
  COMP_PRECISION cohesion;
  /* 
     total seismic moment release 
  */
  float total_moment;
  /* files for output */
  FILE *(*flt_stress_out),*events_out,
    *cevents_out,*(*slip_line_out);
  /* files for input */
  FILE *i_mat_in;
  BOOLEAN flt_stress_init,events_init,slip_line_init,
    print_bulk_fields;
  /* calculation mode and state switches */
  int op_mode,op_state,stress_state_init;
  /* 
     internal checks and initializxation flags
   */
  BOOLEAN geometry_init,bc_init,pos_init,int_mat_init,
    bulk_field_init,read_int_mat_from_file,
    moment_release_init;
  //
  // maximum number of fault groups that will lead
  // to individual output files
  //
  int max_nr_flt_files;

  /* 
     for PGPLOTTING X WINDOWS OUTPUT 
  */
#ifdef USE_PGPLOT
#define MAX_NR_X_WINDOWS 3
  /* plotting issues */
  BOOLEAN nr_xwindows,moment_array_init;
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
  COMP_PRECISION x[3]; /* location of center of patch */

  float pos[2]; /* position in the group of patches
		   assuming linear fault */
  COMP_PRECISION w, l;
  float strike, dip;    /* degrees, dip 90 is 
			   horizontal,strike east  is 0, 
			   counterclockwise */
  COMP_PRECISION cos_alpha,
    sin_alpha;/* cos and sin of 90 - strike */
  COMP_PRECISION normal[3],
    t_strike[3],t_dip[3]; /* normal and tangential 
			     (strike and dip direction) 
			     unit vectors */
  float mu_d,mu_s; /* dynamic and static friction 
			   coefficients
			  
			*/
  int group; /* fault patch can have a 
		group assigned */
  short int mode;/* fault operation mode  */
  BOOLEAN active; /* fault is active in the current
		       solver */
};

/* 
   FUNCTION DECLARATIONS 
*/
/* 
   general algorithmic flow 
*/
void initialize(struct med **, struct flt **,BOOLEAN,
		int, BOOLEAN,BOOLEAN,COMP_PRECISION,
		COMP_PRECISION *,COMP_PRECISION *,
		BOOLEAN);
void init_files(struct med **,struct flt **);
void terminate(struct med *, struct flt *);
void phelp(void);
void init_parameters(char **, int , BOOLEAN *,BOOLEAN *,
		     BOOLEAN *,BOOLEAN *,COMP_PRECISION *,
		     int *, BOOLEAN *);
/* 
   input 
*/
void read_geometry(struct med **, struct flt **,BOOLEAN);
FILE *myopen(char *,char *);
void read_boundary_conditions(struct med *,struct flt *);
void read_one_step_bc(FILE *,struct med *,struct flt *);
I_MATRIX_PREC ic_from_file(int , int , int , int ,
			   struct med *,struct flt *);
void read_stress_fac(BOOLEAN,
		     COMP_PRECISION *,COMP_PRECISION *);
/* 
   output 
*/
void print_stress(struct med *,struct flt *);
void print_displacement(struct med *,struct flt *);
void print_fault_data(int grp, char *,struct med *,
		      struct flt *);
void print_fault_stress_and_slip(struct med *,struct flt *);
void print_slip_line(struct med *,struct flt *);
void print_stress_on_fault(struct med *,struct flt *,int);
void print_interaction_matrix(struct med *,struct flt *);
void print_a_matrix(A_MATRIX_PREC *,int , char *,char *);
void print_equations(int,BOOLEAN *,int *,A_MATRIX_PREC *,
		    int ,char *,struct flt *);
void print_solutions(int , int *,struct flt *,struct med *,
		     char *);
void flush_moment_stack(struct med *);

void flush_slipline(struct med *,struct flt *);
/*
  X window output 
*/
#ifdef USE_PGPLOT
#ifdef SGI_SUBROUTINE_CONVENTION
#define palett palett_
#endif
void init_plot_window(struct med *, struct flt *);
void close_plot_window(struct med *, struct flt *);
void plot_patch(int ,struct med *,struct flt *,int ,int );
void add_to_plotting_list(int ,int **,int *);
void plot_time_label(struct med *, struct flt *,int );
void plot_projected_patch(int ,struct med *, 
			  struct flt *,int , int );
void plot_time_tics(struct med *, struct flt *,float);
void update_plots(struct med *, struct flt *);
void plot_quake(int, struct med *,struct flt *);
void palett(int *, float *, float *);
void plot_moment_array(struct med *);
#endif
/* 
   geometry 
*/
void rotate_vec(COMP_PRECISION *, COMP_PRECISION *,
		COMP_PRECISION , COMP_PRECISION );
void rotate_mat(COMP_PRECISION [3][3], 
		COMP_PRECISION [3][3],COMP_PRECISION , 
		COMP_PRECISION );
COMP_PRECISION resolved_stress(COMP_PRECISION *,
			       COMP_PRECISION [3][3],
			       COMP_PRECISION *);
COMP_PRECISION resolved_vector(COMP_PRECISION *,
			       COMP_PRECISION *);
void calculate_corners(COMP_PRECISION [4][3],struct flt *);
COMP_PRECISION project_vector(COMP_PRECISION *,
			      COMP_PRECISION *);
void resolve_force(COMP_PRECISION *,
			     COMP_PRECISION [3][3],
			     COMP_PRECISION *);
void calculate_position_of_patch(struct med *,
				 struct flt *);
void intersect(COMP_PRECISION *,COMP_PRECISION *,
	       COMP_PRECISION *,COMP_PRECISION *,
	       COMP_PRECISION *, int *);
int tri_tri_intersect(float *,float *,float *,
                      float *,float *,float *);
COMP_PRECISION distance(COMP_PRECISION *,
			COMP_PRECISION *);
BOOLEAN far_enough(struct flt *,struct flt *);
COMP_PRECISION norm(COMP_PRECISION *,int );
COMP_PRECISION dotp(COMP_PRECISION *,COMP_PRECISION *,
		    int );
void calc_three_stress_components(COMP_PRECISION [3][3],
				  COMP_PRECISION *,
				  COMP_PRECISION *,
				  COMP_PRECISION *,
				  COMP_PRECISION *,
				  COMP_PRECISION *,
				  COMP_PRECISION *,
				  COMP_PRECISION *);
void calc_base_vecs(COMP_PRECISION *,COMP_PRECISION *,
		    COMP_PRECISION *,COMP_PRECISION ,
		    COMP_PRECISION ,COMP_PRECISION ,
		    COMP_PRECISION );

/* 
   stress and displacement calculations 
*/
void initialize_stress_state(struct flt *,struct med *,
			     BOOLEAN,
			     COMP_PRECISION *, COMP_PRECISION *);
void calc_fields(struct med *,struct flt *,BOOLEAN,
		 BOOLEAN,COMP_PRECISION *, COMP_PRECISION *);
void background_stress(COMP_PRECISION [3][3], 
		       COMP_PRECISION *, COMP_PRECISION,
		       COMP_PRECISION *, COMP_PRECISION *);
void background_disp(COMP_PRECISION *, COMP_PRECISION *, 
		     struct med *,
		     COMP_PRECISION *, COMP_PRECISION *);
void update_stress_state(struct flt *,struct med *);
/* 
   interaction matrix and Greens function routines 
   Okada subroutines
*/
void calc_interaction_matrix(struct med *,struct flt *);
#ifdef USEMYDC3D
/*
extern void mydc3d(COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   int*);
*/
extern void mydc3d(COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,COMP_PRECISION*,
		   COMP_PRECISION*,int*);

#else
extern void dc3d(COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,COMP_PRECISION*,
		 COMP_PRECISION*,int*);
#endif

void eval_green(COMP_PRECISION *,struct flt *,
		COMP_PRECISION *,COMP_PRECISION *, 
		COMP_PRECISION [3][3],int *,
		struct med *);
COMP_PRECISION interaction_coefficient(int , int , int , 
				       int ,struct med *,
				       struct flt *,int *);
void get_right_slip(COMP_PRECISION *,int , int , 
		    struct flt *);
int select_i_coeff_calc_mode(struct med *);
int imatrix_size(struct med *);
/* 
   matrix solution related 
*/
#include "svd.h"
#include "nnls.h"
void assemble_a_matrix(A_MATRIX_PREC *,int ,BOOLEAN *,
		       int, int *,struct flt *fault,
		       struct med *);
void assemble_a_matrix_1(A_MATRIX_PREC *,int ,BOOLEAN *,
			 int, int *,struct flt *fault,
			 struct med *);
void assemble_a_matrix_2(A_MATRIX_PREC *,int ,BOOLEAN *,
			 int, int *,struct flt *fault,
			 struct med *);
void assemble_a_matrix_3(A_MATRIX_PREC *,int ,BOOLEAN *,
			 int, int *,struct flt *fault,
			 struct med *);
void init_equation_system(struct med *,struct flt *fault);
void solve(struct med *,struct flt *fault);
void add_to_active_fault_list(int ,int **,int *,
			      BOOLEAN **);
void add_to_right_hand_side(COMP_PRECISION ,
			    A_MATRIX_PREC **,int *);
void add_solution(int,BOOLEAN *,A_MATRIX_PREC *,int *,
		  struct med *,struct flt *,BOOLEAN );
void resize_arrays(int **,int ,int , A_MATRIX_PREC **,
		   BOOLEAN **);

/* 
   fault fracture criterion and activations 
*/
BOOLEAN fracture_criterion(int, struct flt *,
			   COMP_PRECISION *,BOOLEAN *,
			   struct med *);
void activate_faults(struct flt *,struct med *);
void quake(BOOLEAN *,COMP_PRECISION *,int ,
	   struct flt *,struct med *,BOOLEAN);
void fault_criterion(int ,struct flt *,
			struct med *);
/* 
   misc 
*/
void fltswap(struct flt *, struct flt *);
#define MEMERROR {fprintf(stderr,"mem error\n");exit(-1);}
COMP_PRECISION mygauss_randnr(COMP_PRECISION,long *);
double gasdev(long *);
double ran1(long *);
COMP_PRECISION myrand(long *);
int myrandi(int ,long *);
COMP_PRECISION myrandnr(COMP_PRECISION, long *);
void sincos_deg(COMP_PRECISION *,COMP_PRECISION *,
		COMP_PRECISION );
char *name_boolean(BOOLEAN );
BOOLEAN toggle(BOOLEAN variable);
/*
  model generation 
*/
void divide_fault_in_patches(int ,struct flt *,int , int,
			     int,COMP_PRECISION,
			     BOOLEAN);
void optimize(struct flt *,struct med *);

