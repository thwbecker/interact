/* 
   FUNCTION DECLARATIONS 
*/
/* 
   general algorithmic flow 
*/
void initialize(struct med **, struct flt **,BOOLEAN,
		int, BOOLEAN,BOOLEAN,COMP_PRECISION,
		COMP_PRECISION *,COMP_PRECISION *,
		BOOLEAN,BOOLEAN,BOOLEAN);
void init_files(struct med **,struct flt **);
void terminate(struct med *, struct flt *);
void phelp(void);
void init_parameters(char **, int , BOOLEAN *,BOOLEAN *,
		     BOOLEAN *,BOOLEAN *,COMP_PRECISION *,
		     int *, BOOLEAN *, BOOLEAN *,BOOLEAN *);
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
char *comment_on_code(short int );
char *comment_on_code_bc(short int , COMP_PRECISION );
/* 
   output 
*/
void fiddle_with_limits_for_plot(struct med *,int *,
				 BOOLEAN *,
				 COMP_PRECISION *,BOOLEAN);
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
void add_b_to_a_vector_3d(COMP_PRECISION *,
			  COMP_PRECISION *);
void scale_vector_3d(COMP_PRECISION *,COMP_PRECISION );
void a_equals_b_vector(COMP_PRECISION *,COMP_PRECISION *,
		       int );
void a_equals_b_vector_3d(COMP_PRECISION *,
			  COMP_PRECISION *);

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
void get_local_x_on_plane(COMP_PRECISION *,
			  COMP_PRECISION *,
			  COMP_PRECISION *,
			  COMP_PRECISION *,
			  COMP_PRECISION *);
void get_fault_plane_basevec(COMP_PRECISION *,
			     COMP_PRECISION *,
			     COMP_PRECISION *,
			     struct flt *,
			     struct med *);
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
COMP_PRECISION norm_3d(COMP_PRECISION *);

COMP_PRECISION dotp(COMP_PRECISION *,COMP_PRECISION *,
		    int );
COMP_PRECISION dotp_3d(COMP_PRECISION *,COMP_PRECISION *);
void normalize(COMP_PRECISION *,int );
void normalize_3d(COMP_PRECISION *);
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
void adjust_time_step(struct flt *, struct med *);
/* 
   interaction matrix and Greens function routines 
   Okada subroutines
*/
void calc_interaction_matrix(struct med *,struct flt *);

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
void deactivate_group(int, int , struct flt *, 
		      struct med *);
void deactivate_patch(int, struct flt *, struct med *);
void single_dir_slip_check(BOOLEAN *,COMP_PRECISION *,
			   int , int, 
			   struct flt *, struct med *);
void two_dir_slip_check(BOOLEAN *,COMP_PRECISION *,
			int *,int , struct flt *, 
			struct med *);
COMP_PRECISION coulomb_stress(COMP_PRECISION ,
			      COMP_PRECISION ,
			      COMP_PRECISION ,
			      COMP_PRECISION );
BOOLEAN in_coulomb_compress_regime(COMP_PRECISION);
/* 
   misc 
*/
void fltswap(struct flt *, struct flt *);
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

// segseg

char	SegSegInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p );
char	ParallelInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p );
BOOLEAN	Between( tPointi a, tPointi b, tPointi c );
void	Assigndi( tPointd p, tPointi a );
int	Collinear( tPointi a, tPointi b, tPointi c );
int     AreaSign( tPointi a, tPointi b, tPointi c );
