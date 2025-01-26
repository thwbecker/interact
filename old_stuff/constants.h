/* headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* for sgi */
#define dc3d dc3d_
/* programming stuff */
#define COMP_PRECISION double
#define BOOLEAN short int
#define TRUE 1
#define FALSE 0
#define SQUARE(x) ((x)*(x))
#define MIN(x, y) (((x)<(y))?(x):(y))
#define MAX(x, y) (((x)>(y))?(x):(y))

/* filenames and output options */
#define GEOMETRY_FILE "geom.in"
#define BC_FILE "bc.in"
#define STRESS_OUT_FILE "stress.out"
#define FAULT_STRESS_OUT_FILE "flt_stress.out"
#define DISP_OUT_FILE "disp.out"
#define FAULT_DATA_PREFIX "flt"
#define MAX_NR_FLT_FILES 30
#define EVENT_FILE "events.dat"
/* indices used in programs */
#define STRIKE 0
#define DIP 1
#define NORMAL 2
#define X 0
#define Y 1
#define Z 2
/* real constants */
#define DEG2RAD 0.0174532925199433
#define EPS 1.0e-8
#define NAN 1e20
#define STRLEN 200
/* material properties, alpha = (lambda+mu)/(lambda+2 mu) */
#define ALPHA 0.66666666666666667
#define LAMBDA 1.0
#define PRESSURE 0.33e-3
#define DELTA_MU 0.1
/* operational modes */
#define ONE_STEP 1
#define SIMULATE_LOADING 2
/* offset operator for the interaction matrix and 3d fields */
#define POSI(i, j, k, l) ((((i)*(medium->nrflt)+(j))*9) + (k*3) + (l))
#define POSS(i, j, k, l) (((k)*(medium->n[Y]*medium->n[X])+((i)*medium->n[Y])+(j))*6+(l))
#define POSU(i, j, k, l) (((k)*(medium->n[Y]*medium->n[X])+((i)*medium->n[Y])+(j))*3+(l))
/* structures */
struct med{
  /* nr of faults */
  int nrflt;
  /* interaction matrix */
  COMP_PRECISION *i;
  /* geometrical boundaries for stress output */
  COMP_PRECISION xmin[3],xmax[3];
  /* number of ;;columns and rows and slices in stress output */
  int n[3];
  /* stress and displacement field outputs */
  float *s,*u;
  /* timing and printing issues */
  COMP_PRECISION time,dt,stop_time;
  COMP_PRECISION print_interval,print_time;
  /* total seismic moment release */
  COMP_PRECISION moment;
  
  /* files for output */
  FILE *(*flt_stress_out),*events_out;
  BOOLEAN flt_stress_init,events_init;

  /* calculation mode and state */
  int op_mode,op_state,stress_state_init;
};

struct flt{
  COMP_PRECISION u[3];/* displacements on fault in strike, dip, 
			 and tensile direction */
  COMP_PRECISION s[3];/* strike, normal and 
			 dip stress */
  COMP_PRECISION sinc[3]; /* assuming that the background
				stress increases linearly,
				these are the increments in
				strike, dip, and normal 
				direction,
				s=a*t+b
			     */
				
  COMP_PRECISION x[3]; /* location of center of patch */
  COMP_PRECISION w, l; /* dimensions of patch,
			  half width and half length
		       */
  COMP_PRECISION dip,ca,sa;/* dip, cos and sin of strike
			      degrees, dip 90 is horizontal,
			      angle west is 0 */
  COMP_PRECISION normal[3],
    t_strike[3],t_dip[3]; /* normal and tangential 
			     (strike and dip direction) 
			     unit vectors */

  COMP_PRECISION mu_s,delta_mu; /* static friction coefficient 
				   and difference between static
				   and dynamic friction 
				   coefficient
				*/
};

/* function declarations */

void calc_fields(struct med *,struct flt *,BOOLEAN,BOOLEAN);
void calc_interaction_matrix(struct med *,struct flt *);
void read_geometry(struct med **, struct flt **);
FILE *myopen(char *,char *);
extern void dc3d(double*,double*,double*,double*,double*,
		 double*,double*,
		 double*,double*,double*,double*,double*,
		 double*,double*,
		 double*,double*,double*,double*,double*,
		 double*,double*,
		 double*,double*,double*,double*,int*);
void read_boundary_conditions(struct med *,struct flt *);
void eval_green(COMP_PRECISION *,struct flt *,COMP_PRECISION *,
		COMP_PRECISION *, COMP_PRECISION [3][3],int *,
		struct med *);

void print_stress(struct med *,struct flt *);
void print_displacement(struct med *,struct flt *);
#define MEMERROR {fprintf(stderr,"mem error\n");exit(-1);}
void rotate_vec(COMP_PRECISION *, COMP_PRECISION *,
		COMP_PRECISION , COMP_PRECISION );
void rotate_mat(COMP_PRECISION [3][3], COMP_PRECISION [3][3],
		COMP_PRECISION , COMP_PRECISION );
COMP_PRECISION resolved_stress(COMP_PRECISION *,
			       COMP_PRECISION [3][3],
			       COMP_PRECISION *);
COMP_PRECISION resolved_vector(COMP_PRECISION *,COMP_PRECISION *);
void background_stress(COMP_PRECISION [3][3], COMP_PRECISION *, 
		       COMP_PRECISION);
void background_disp(COMP_PRECISION *, COMP_PRECISION *, 
		     struct med *);
void initialize(struct med **, struct flt **);
void update_stress_state(struct flt *,struct med *);
BOOLEAN activate_faults(struct flt *,struct med *);
COMP_PRECISION project_vector(COMP_PRECISION *,COMP_PRECISION *);
void resolve_force(COMP_PRECISION *,
			     COMP_PRECISION [3][3],
			     COMP_PRECISION *);
void print_stress_on_fault(struct med *,struct flt *,int);
void terminate(struct med *, struct flt *);
void initialize_stress_state(struct flt *,struct med *);
void quake(COMP_PRECISION *,int ,
	   struct flt *,struct med *);
