/* adjust_time_step.c */
void adjust_time_step(struct flt *, struct med *);
/* blockinvert.c */
/* calc_cart_from_eigen_stress.c */
unsigned char read_vecs(int, double *, double *, double *, double *, double *);
void fehelp(char **);
/* calc_eigen_from_cart_stress.c */
void fehelp(char **);
/* calc_interaction_matrix.c */
void print_help_local(char *);
/* calc_spatial_correlation.c */
void calc_spatial_correlation(struct flt *, int, int, int, int, double **, double **, int **, double, double, double *, float **, float **, float **);
double correlation_coefficient(double *, double *, int);
/* calc_stress.c */
void initialize_stress_state(struct flt *, struct med *, unsigned char, double *, double *);
void update_stress_state(struct flt *, struct med *);
void calc_fields(struct med *, struct flt *, unsigned char, unsigned char, double *, double *);
void background_stress(double [3][3], double *, double, double *, double *, double);
void background_disp(double *, double *, struct med *, double *, double *);
void get_local_x_on_plane(double *, double *, double *, double *, double *);
void get_fault_plane_basevec(double *, double *, double *, struct flt *, struct med *);
void calc_deviatoric_stress(double [3][3], double [3][3], double *);
/* calc_stress_stat.c */
/* check_feedback.c */
/* check_interaction.c */
unsigned char incrementally_check_interaction_coefficients(struct flt *, int, struct med *, unsigned char, int *, double);
unsigned char check_coulomb_stress_feedback(int, int, struct flt *, struct med *, unsigned char, int, unsigned char, int *, double);
/* compare_fault.c */
int compare_fault_length(struct flt *, struct flt *);
int compare_fault_width(struct flt *, struct flt *);
/* coulomb_stress.c */
double coulomb_stress(double, double, double, double);
double cstress_drop(double, double, double, double);
unsigned char in_coulomb_compress_regime(double);
/* create_random_mu_file.c */
/* create_random_stress_file.c */
void get_random_stress(double *, double *, double *, int, long *, int);
/* divide_fault_in_patches.c */
void create_patches(int, struct flt *, struct flt **, int *, int *, unsigned char, unsigned char, double *);
void determine_segments(int *, int *, struct flt *, unsigned char, double *);
void divide_fault_in_patches(int, struct flt *, struct flt **, int *, int *, unsigned char, unsigned char);
void get_flt_location(struct flt *, double *, double *, double *, int, int);
/* eigensystem.c */
void eigensystem3d(double [3][3], double *, double *);
/* eval_2dsegment.c */
void eval_2dsegment_plane_strain(double *, struct flt *, double *, double *, double [3][3], int *);
void eval_2dsegment_plane_stress(double *, struct flt *, double *, double *, double [3][3], int *);
void eval_2dsegment_plane_strain_basic(double *, struct flt *, double *, double *, double [3][3], int *);
void eval_2dsegment_plane_stress_basic(double *, struct flt *, double *, double *, double [3][3], int *);
void get_2dseg_geo(double *, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *);
void get_2dseg_ffac(double *, double *, double *, double *, double *, double *, double, double, double, double, double, double, double, double, double, double, double *);
void get_2dseg_disp(double *, double *, double *, double, double, double, double, double, double);
void get_2dseg_stress(double [3][3], double *, double *, double, double, double, double);
/* eval_geo_green.c */
void eval_geo_green(double, double, double, double, double *, double *, double *, double [3][3], int *);
/* eval_green.c */
void eval_green(double *, struct flt *, double *, double *, double [3][3], int *);
void eval_green_basic(double *, struct flt *, double *, double *, double [3][3], int *);
/* eval_okada.c */
void eval_rectangle(double *, struct flt *, double *, double *, double [3][3], int *);
void eval_rectangle_basic(double *, double, double, double, double, double *, double *, double [3][3], int *);
void eval_point(double *, struct flt *, double *, double *, double [3][3], int *);
void eval_point_short(double *, double *, double, double, double, double, double *, double *, double [3][3], int *);
void set_stress_and_disp_nan(double [3][3], double *);
/* eval_triangle.c */
void eval_triangle(double *, struct flt *, double *, double *, double [3][3], int *);
void get_gauss_points(double *, double *, double *, int);
/* far_enough.c */
unsigned char far_enough(struct flt *, struct flt *, double);
/* fit_plane.c */
void fit_plane(int, double *, double *, double *, double *, double *, double *, double *, double *, float *, float *, double *, unsigned char, unsigned char);
void points2patch(struct flt *, double *);
/* fltcopy.c */
void fltswap(struct flt *, struct flt *);
void fltcp(struct flt *, struct flt *);
/* fracture_criterion.c */
unsigned char fracture_criterion(int, struct flt *, double *, unsigned char *, struct med *);
void two_dir_slip_check(unsigned char *, double *, unsigned char *, int, struct flt *, struct med *);
void deactivate_patch(int, struct flt *, struct med *);
void deactivate_group(int, int, struct flt *, struct med *);
int calc_absolute_shear_stress_and_inc(double *, double *, int, struct flt *);
int calc_absolute_shear_stress(double *, int, struct flt *);
/* generate_slipdia.c */
/* geo_okada.c */
void get_projected_fault_parameters(double [2][2], double, double *, double *, double *, double *, double *);
/* geometry.c */
double resolved_stress(double *, double [3][3], double *);
void calc_three_stress_components(double [3][3], double *, double *, double *, double *, double *, double *, double *);
void calc_lhemi_proj(double, double, double *);
void resolve_force(double *, double [3][3], double *);
void calc_base_vecs(double *, double *, double *, double, double, double, double);
void get_maxsdir_stress_drops2(double *, double, double *);
void get_maxsdir_stress_drops(double *, double, double *, double, double *, double *);
void calculate_corners(double [4][3], struct flt *, double *, double *);
void calculate_bloated_corners(double [4][3], struct flt *, double);
void calculate_quad_corners(double [4][3], struct flt *, double);
void calculate_point_source_corners(double [4][3], struct flt *, double, double *, double *);
void calculate_seg_corners(double [4][3], struct flt *, double);
double quad_area(double *);
double triangle_area(double *);
double triangle_area_gh(double *, double *);
void get_gh_tri_vec(double *, double *, double *);
void get_gh_quad_vec(double *, double *, double *, double *);
unsigned char check_planar(double *);
void get_alpha_dip_tri_gh(double *, double *, double *, double *, double *);
void calc_group_geometry(struct med *, struct flt *, struct geog *);
void vec_to_angles(double *, double *, double *);
void angles_to_vec(double, double, double *);
double vec_to_strike(double *);
double vec_to_dip(double *);
void check_fault_angles(struct flt *);
void check_angles(double *, double *);
void globalx(double *, double, double, double *);
void calc_centroid_tri(double *, double *);
void calc_mean_quad_coord(double *, double *);
void calc_centroid_quad(double *, double *);
int compare_cp(double *, double *);
void calc_vec_stat(double *, int, double *);
unsigned char patch_is_2d(unsigned short int);
void calculate_position_of_patch(struct med *, struct flt *);
/* geoproject.c */
void geoproject(double *, double *, int, double, double, double, double, double, int);
/* help_and_comments.c */
void phelp(void);
char *comment_on_code(short int);
char *comment_on_code_bc(short int, double);
/* init.c */
void check_parameters_and_init(int, char **, struct med **, struct flt **, unsigned char *, double *, double *);
void initialize(struct med **, struct flt **, unsigned char, int, unsigned char, unsigned char, double, double *, double *, unsigned char, unsigned char, unsigned char, double, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, short int, unsigned char, unsigned char, double, unsigned char, unsigned char);
void init_files(struct med **, struct flt **);
void terminate(struct med *, struct flt *);
void init_parameters(char **, int, unsigned char *, unsigned char *, unsigned char *, unsigned char *, double *, int *, unsigned char *, unsigned char *, unsigned char *, double *, unsigned char *, double *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, short int *, unsigned char *, unsigned char *, double *, unsigned char *, unsigned char *);
void advance_argument(int *, int, char **);
char *name_boolean(unsigned char);
unsigned char toggle(unsigned char *);
void read_stress_fac(unsigned char, double *, double *, double);
/* input.c */
int read_moment_file(float **, float **, float *, float *, int *, unsigned char);
int read_patch_event_file(float *, int *, int *, float *, float *, FILE *, struct med *);
int write_patch_event_file(float, int, int, float, float *, FILE *);
/* interact.c */
void calc_interaction_matrix(struct med *, struct flt *, unsigned char);
double interaction_coefficient(int, int, int, int, struct flt *, int *);
void get_right_slip(double *, int);
double ic_from_file(int, int, int, int, struct med *);
double aij_from_file(int, int, int, FILE *);
int select_i_coeff_calc_mode(struct med *);
int imatrix_size(struct med *);
/* llgeo.c */
void rotate_vec(double *, double *, double, double);
void rotate_vec2d(double *, double *, double, double);
void rotate_mat(double [3][3], double [3][3], double [3][3]);
void rotate_mat_z(double [3][3], double [3][3], double, double);
double project_vector(double *, double *);
double distance_3d(double *, double *);
double distance(double *, double *, int);
float distance_float(float *, float *, int);
double distance_squared_3d(double *, double *);
double distance_squared(double *, double *, int);
float distance_squared_float(float *, float *, int);
double norm(double *, int);
double norm_3d(double *);
double dotp(double *, double *, int);
double dotp_3d(double *, double *);
void normalize(double *, int);
void normalize_3d(double *);
void a_equals_b_vector(double *, double *, int);
void a_equals_b_vector_3d(double *, double *);
void swap_ab_vector_3d(double *, double *);
void scale_vector_3d(double *, double);
void scale_vector(double *, double, int);
void add_b_to_a_vector_3d(double *, double *);
void sub_b_from_a_vector_3d(double *, double *);
void c_eq_a_plus_b_3d(double *, double *, double *);
void c_eq_a_minus_b_3d(double *, double *, double *);
void c_eq_a_minus_b(double *, double *, double *, int);
void cross_product(double *, double *, double *);
double find_max_vec(double *, int);
double find_max_abs_vec(double *, int);
int find_lde_max(double *, int, int, double);
/* lusolve.c */
void lu_driver(double *, double *, double *, int, int, struct med *);
/* main.c */
/* makefault.c */
/* matrixio.c */
void print_matrix(double *, int, int, FILE *, unsigned char);
void print_a_matrix(double *, int, int, FILE *, double *, unsigned char);
void read_a_matrix_from_file(double *, int, int, char *, char *);
void print_a_matrix_to_file(double *, int, int, char *, char *);
void print_system(double *, double *, double *, int, int, FILE *);
void print_sym3x3matrix(double [3][3], FILE *);
void print_interaction_matrix(struct med *, struct flt *);
void print_reduced_interaction_matrix(struct med *, struct flt *);
void print_vector(double *, int, FILE *);
void print_b_vector(double *, int, FILE *, double *, unsigned char);
/* mspectral.c */
void find_range(int *, int *, float, float, float, float, float *, float *, float *, int);
/* myopen.c */
FILE *myopen(char *, char *);
/* mysincos.c */
void my_sincos_deg(double *, double *, double);
void my_sincos_deg_ftn(double *, double *, double *);
void my_sincos(double *, double *, double);
/* nnls.c */
void nnls_driver(double *, double *, double *, int, int);
/* optimize.c */
void optimize(struct flt *, struct med *);
double distsq(struct flt *, struct flt *);
double max_dist(struct flt *, struct med *);
double penalty_dist(struct flt *, struct med *, double);
/* output.c */
void print_slip_line(struct med *, struct flt *);
void flush_slipline(struct med *, struct flt *);
void print_fault_stress_and_slip(struct med *, struct flt *, unsigned char);
void print_fault_data(char *, struct med *, struct flt *);
void print_fault_stress_stat(FILE *, int, struct med *, struct flt *);
void print_group_data_geom(char *, struct med *, struct flt *, int, int);
double select_val_for_print(struct flt *, int);
void print_stress(struct med *, struct flt *);
void print_stress_on_fault(struct med *, struct flt *, int);
void print_displacement(struct med *, struct flt *);
void print_equations(int, unsigned char *, int *, double *, int, char *, struct flt *);
void print_solutions(int, int *, struct flt *, struct med *, char *);
void flush_moment_stack(struct med *);
void fiddle_with_limits_for_plot(struct med *, int *, unsigned char *, double *, unsigned char);
/* patch2bc.c */
/* patch2corners.c */
/* patch2geom.c */
/* patch2group.c */
/* patch2xyz.c */
/* patch2xyzvec.c */
/* period.c */
void period(float *, float *, int, float, float, float *, float *, int, int *, int *, float *);
void avevar(float *, int, float *, float *);
void fasper(float *, float *, int, float, float, float *, float *, int, int *, int *, float *);
void spread(float, float *, int, float, int);
void realft(float *, int, int);
void four1(float *, int, int);
/* plotevents.c */
/* plotting.c */
/* points2patch.c */
/* print_patch_geometry.c */
void print_patch_geometry_and_bc(int, struct flt *, int, double, unsigned char, FILE *, unsigned char);
/* project_stress.c */
/* quake.c */
void quake(unsigned char *, double *, int, struct flt *, struct med *, unsigned char, unsigned char);
void add_quake_stress(int, unsigned char *, double *, struct flt *, struct med *);
/* randgen.c */
double myrand(long *);
double myrandnr(double, long *);
int myrandi(int, long *);
double mygauss_randnr(double, long *);
double mygr_randnr(double, double, long *);
double gasdev(long *);
double ran1(long *);
double ran2(long *);
/* randomflt.c */
void assign_random(double *, double, double, double, long *, int);
void check_input_parameters(int, char **, int *, long *, int *, double *, unsigned char *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, unsigned char *, int *, double *);
/* randomize_list.c */
void randomize_list(int **, int, unsigned char);
int slist_sort(struct slist *, struct slist *);
/* randomize_strike.c */
/* read_bin_events.c */
/* read_boundary_conditions.c */
void read_boundary_conditions(struct med *, struct flt *, double *, double *);
void read_one_step_bc(FILE *, struct med *, struct flt *, double *, double *);
unsigned char slip_type_bc(int);
/* read_geometry.c */
void read_geometry(char *, struct med **, struct flt **, unsigned char, unsigned char, unsigned char);
/* restart.c */
void adjust_medium_for_restart(struct med *, struct flt *);
/* rupture.c */
unsigned char activate_faults(struct flt *, struct med *);
void fault_criterion(int, struct flt *, struct med *);
/* segseg.c */
void intersect(double *, double *, double *, double *, double *, int *);
char SegSegInt(tPointi, tPointi, tPointi, tPointi, tPointd);
char ParallelInt(tPointi, tPointi, tPointi, tPointi, tPointd);
void Assigndi(tPointd, tPointi);
unsigned char Between(tPointi, tPointi, tPointi);
int Collinear(tPointi, tPointi, tPointi);
int AreaSign(tPointi, tPointi, tPointi);
unsigned char crosses(double *, double *);
/* solve.c */
void solve(struct med *, struct flt *);
void init_equation_system(struct med *, struct flt *);
void add_to_active_fault_list(int, int **, int *, unsigned char **);
void add_to_right_hand_side(double, double **, double **, int *);
void add_solution(int, unsigned char *, double *, int *, struct med *, struct flt *, unsigned char, unsigned char, double);
void assemble_a_matrix(double *, int, unsigned char *, int, int *, struct flt *, struct med *);
/* solve_mode_dependend.c */
void assemble_a_matrix_4(double *, int, unsigned char *, int, int *, struct flt *, struct med *);
void add_quake_stress_4(unsigned char *, double *, int, struct flt *, struct med *);
unsigned char check_coulomb_stress_feedback_4(int, int, struct flt *, struct med *, unsigned char, unsigned char, int *, double);
/* sort_events.c */
/* sparse.c */
size_t create_crs_sparse_from_memory(int, double *, double, unsigned int **, unsigned int **, double **);
size_t create_crs_sparse_from_file(int, double, unsigned int **, unsigned int **, double **, FILE *);
size_t create_ccs_sparse_from_memory(int, double *, double, unsigned int **, unsigned int **, double **);
size_t create_ccs_sparse_from_file(int, double, unsigned int **, unsigned int **, double **, FILE *);
/* sparse_nr.c */
size_t create_nrs_sparse(struct med *, double, double *, unsigned char);
size_t create_nrs_sparse_from_memory(int, double *, double, unsigned int **, double **);
size_t create_nrs_sparse_from_file(int, double, unsigned int **, double **, FILE *);
double get_nrs_sparse_el(int, int, unsigned int *, double *);
/* sparse_solve.c */
void sparse_driver(unsigned int *, unsigned int *, double *, double *, double *, int, struct med *);
/* svd.c */
void svd_driver(double *, double *, double *, int, int, double, int, double *);
/* test_optimize.c */
/* test_sparse.c */
/* test_stuff.c */
/* tri2patch.c */
/* trigonometry.c */
double dist_on_sphere(double, double, double, double);
double azimuth(double, double, double, double);
void get_gc_pole(double, double, double, double, double *, double *);
void get_point_on_gc(double, double, double, double, double, double *, double *, double *);
void get_point_on_course(double, double, double, double, double *, double *);
void lonlat2xyz(double, double, double *);
void xyz2lonlat(double *, double *, double *);
/* tritri.c */
int coplanar_tri_tri(float [3], float [3], float [3], float [3], float [3], float [3], float [3]);
int tri_tri_intersect(float [3], float [3], float [3], float [3], float [3], float [3]);
