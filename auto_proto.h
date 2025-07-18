/* adjust_time_step.c */
void adjust_time_step(struct flt *, struct med *);
/* block_checkflt.c */
/* block_compute_vel_from_omega.c */
/* block_eval_geokada.c */
void block_eval_geookada(double *, double *, double *, double [3][3], double, double, double, double, double, double, double, double, double, double, double, int *, unsigned char, unsigned char);
void block_save_solution_and_faults(double *, int, int, struct bflt *, double *, struct prj *, FILE *, unsigned char, unsigned char);
void block_load_solution_and_faults(double **, int *, int *, struct bflt **, double **, struct prj **, FILE *, unsigned char *, unsigned char *);
void block_eval_blockvec(double *, double, int, double *, struct prj *, double *);
/* block_evaluate_solution.c */
/* block_levmarq.c */
void run_lm(struct bmd *, long int *, struct prj, double *, unsigned char, unsigned char, double, double *, unsigned char, unsigned char, double, unsigned char, unsigned char, unsigned char, int, int, char **);
/* block_matrix.c */
void assemble_block_fltdep_matrices(double **, double **, double **, double **, double *, double *, double *, int, int, int, int, int, int, int, double *, struct bflt *, struct bck *, unsigned char, unsigned char);
void assemble_block_a(double **, double *, int *, int, int, struct bck *, int, struct bmd *);
void assemble_block_d(double **, struct bflt *, int, int, int, struct bmd *);
void assemble_block_i(double **, struct bflt *, int, int, int);
void assemble_block_g(double **, struct bflt *, int, int);
void assemble_block_f(double **, struct bflt *, int, int, struct bck *, double *, unsigned char);
void assemble_block_k(double **, int, int, int, int, int, int, int, double *, double *, double *, unsigned char, unsigned char, double *, struct bflt *, int, int, int, int, double);
void assemble_stress_matrix(double *, int, struct bflt *, double *, int, int);
void block_assemble_fit_vector(double **, int, int, int, int, int, int, int, int, int, double *, int, double *, double *, unsigned char, double *, double **, double *, double **, int, struct bck *, int, struct bflt **, unsigned char, double **, unsigned char, unsigned char, unsigned char, double *, double *, double *, struct prj, unsigned char, int, int, double, struct bmd *);
void block_assemble_dyda_matrix(double **, double **, int, int, int, int, int, int, int, int, int, double *, double *, double **, double *, double **, double **, int, struct bck *, int, struct bflt **, unsigned char, unsigned char, unsigned char, double *, double *, double *, struct prj, unsigned char, int, int, double, struct bmd *);
double block_chi_square(double *, double *, double *, int, int, double, double *, double *);
double block_data_norm(double *, double *, int, int, double, struct bmd *);
int block_slip_direction(int, struct bflt *);
void assign_additional_sol_values(double *, int, int, int, struct bflt *, unsigned char, unsigned char, long int *, int);
void change_locking_depths(struct bflt **, int, int, int, int, double *, double *, double *, unsigned char, struct prj, double **, double **, double *, struct bmd *);
void calc_fault_sn_rms(double *, int, int, double *, double *);
void init_block_mods(struct bmd **);
void init_blocks(struct bmd *, int);
void convert_cart_sol(double *, double *, struct bmd *);
void check_solution_vector(double *, int, int, unsigned char, unsigned char);
void copy_block(struct bck *, struct bck *);
/* block_output.c */
void block_output(struct bmd *, unsigned char, char **, struct prj, unsigned char, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char);
void sort_eigen(int *, double *);
void print_horizontal_stress(double *, FILE *);
void print_projected_stress(double *, FILE *);
void calculate_slipsol_sigma(struct bmd *, double *, double *, int);
void calc_geo_euler_pole(double *, double *, double *, double *);
void print_simple_vel(double *, double *, double *, double *, int, char *);
/* block_read_bflt.c */
void read_bflt(struct bmd *, double, struct prj, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char);
void generate_new_fault(struct bflt **, int *, int *, double *, int, double *, int, double, struct prj, double *, unsigned char, double, unsigned char, unsigned char, double *, double, int, double *, unsigned char *, int *);
void get_bflt_intcoeff(struct bflt **, int, double *, int, double *, int, double *, unsigned char);
void free_bflt(struct bflt **, int, int);
void assign_bflt_dip_mode(struct bflt *, double);
void flip_block_code(struct bflt *, double *);
void print_fault_geometry(struct bflt *, int, FILE *);
void assign_fault_locking_depth_parameters(struct bflt *, double, struct prj, unsigned char, int);
int new_fault(struct bflt **, int);
void init_bflt(struct bflt *);
/* block_read_euler.c */
void read_constrained_euler_poles(struct bmd *, char **, int);
/* block_read_gps.c */
void read_gps_velocities(struct bmd *, struct prj *, double *, char **, unsigned char, int *, int *);
void project_gps_coordinates(double *, double *, int, int, int *, struct bck *, int, struct prj *, char **);
void block_read_centroids(double **, int *);
void find_spherical_rotation(struct bmd *, int, double *, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, double *, int);
void init_gps_pbase(struct bmd *);
void init_cart_gps_velsig(struct bmd *, unsigned char);
void calculate_c2p_gps(double *, double *, struct bmd *);
void convert_sig_p2c(double *, double *, double *);
void remove_cvel_from_vel(struct bmd *, double *, int);
char bname(int);
unsigned char node_in_geo_list(double *, int, double *);
/* block_solve.c */
void solve_block(double *, double *, double *, int, int, int, int, unsigned char, int, int, int, int, double *, double, double *, double **, double **, double *, double *, unsigned char, int, unsigned char, int, unsigned char);
void evaluate_block_solution(double *, int, int, double *, double *, double *, struct bmd *);
/* block_stress.c */
void block_scale_stresses(double *, double *, double *, int, int, int, double, double *, unsigned char);
void calc_horizontal_stress_vec(double *, double *, int);
void calc_horizontal_stress(double *, double *, double *, double *);
void calc_dir_diff_vec(double *, double *, double *, int, unsigned char);
double calc_dir_diff(double, double, unsigned char);
void cart_mat_from_horsym(double, double, double, double *);
void rescale_observed_stresses(double *, double *, double *, double, double *, unsigned char, struct bmd *, unsigned char, unsigned char);
/* blockinvert.c */
/* calc_cart_from_eigen_stress.c */
unsigned char read_vecs(int, double *, double *, double *, double *, double *);
void ccfes_help(char **);
/* calc_design_matrix.c */
void print_help_local2(char *, int, int);
void calc_design_matrix(struct med *, struct flt *, int, int);
void print_design_matrix(struct med *, struct flt *, int, int, FILE *);
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
/* calc_stress_stat.c */
/* check_feedback.c */
/* check_interaction.c */
unsigned char incrementally_check_interaction_coefficients(struct flt *, int, struct med *, unsigned char, int *, double);
unsigned char check_coulomb_stress_feedback(int, int, struct flt *, struct med *, unsigned char, int, unsigned char, int *, double);
/* compare_fault.c */
int compare_fault_length(const void *, const void *);
int compare_fault_width(const void *, const void *);
/* compress_interaction_matrix.c */
/* coulomb_stress.c */
double coulomb_stress(double, double, double, double);
double cstress_drop(double, double, double, double);
unsigned char in_coulomb_compress_regime(double);
/* create_random_mu_file.c */
/* create_random_stress_file.c */
void get_random_stress(double *, double *, double *, int, long *, int);
/* divide_fault_in_patches.c */
void create_patches(int, struct flt *, struct flt **, int *, int *, unsigned char, unsigned char, double *, double, double, long *);
void determine_segments(int *, int *, struct flt *, unsigned char, double *);
void divide_fault_in_patches(int, struct flt *, struct flt **, int *, int *, unsigned char, unsigned char, double, double, long *, unsigned char);
void get_flt_location(struct flt *, double *, double *, double *, int, int);
void randomize_strike_dip(double, double, struct flt *, long *);
/* eigensystem.c */
void eigensystem3d(double [3][3], double *, double *);
void calc_eigensystem_sym3d(double *, double *, double *, unsigned char);
void calc_eigensystem_sym2d(double *, double *, double *, unsigned char);
void eispack_driver(double *, int, int, double *, double *, double *, double *, int);
/* eval_2dsegment.c */
void eval_2dsegment_plane_strain(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
void eval_2dsegment_plane_stress(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
void eval_2dsegment_plane_strain_basic(double *, struct flt *, double *, double *, double [3][3], int *);
void eval_2dsegment_plane_stress_basic(double *, struct flt *, double *, double *, double [3][3], int *);
void get_2dseg_geo(double *, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *);
void get_2dseg_ffac(double *, double *, double *, double *, double *, double *, double, double, double, double, double, double, double, double, double, double, double *);
void get_2dseg_disp(double *, double *, double *, double, double, double, double, double, double);
void get_2dseg_stress(double [3][3], double *, double *, double, double, double, double);
void eval_2dsegment_plane_strain_tdd(double *, struct flt *, double *, double *, double [3][3], int, int *, unsigned char);
/* eval_green.c */
void eval_green_and_project_stress_to_fault(struct flt *, int, int, double *, double *);
void eval_green(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char, unsigned char);
void eval_triangle_general(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char, unsigned char);
void eval_tri_multi_point(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char, unsigned char, unsigned char);
void eval_green_basic(double *, struct flt *, double *, double *, double [3][3], int *);
/* eval_iquad.c */
void eval_iquad(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
/* eval_okada.c */
void eval_okada(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
void eval_okada_basic(double *, double, double, double, double, double *, double *, double [3][3], int *);
void eval_point(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
void eval_point_short(double *, double *, double, double, double, double, double *, double *, double [3][3], int *, unsigned char);
void set_stress_and_disp_nan(double [3][3], double *, unsigned char);
/* eval_triangle_gauss.c */
void eval_triangle_gauss(double *, struct flt *, double *, double *, double [3][3], int *);
void get_gauss_points(double *, double *, double *, int);
/* eval_triangle_nw.c */
void eval_triangle_nw(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
void get_tri_prop_based_on_gh(struct flt *);
/* eval_triangle_tgf.c */
/* far_enough.c */
unsigned char far_enough(struct flt *, struct flt *, double);
/* fit_mean_stress.c */
void eval_stress_correction(double *, double *, double *, double *, int, double *, double *, double *, double *, double *, double *);
/* fit_plane.c */
void fit_plane(int, double *, double *, double *, double *, double *, double *, double *, double *, float *, float *, double *, unsigned char, unsigned char);
void points2patch(struct flt *, double *, unsigned char);
/* fit_simple_stress_from_cart.c */
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
/* fstress2eig.c */
/* fstress2hor.c */
/* generate_random_2d.c */
/* generate_slipdia.c */
/* geo_okada.c */
/* geometry.c */
double patch_area(struct flt *);
void calc_lhemi_proj(double, double, double *);
void resolve_force(double *, double [3][3], double *);
void calc_quad_base_vecs(double *, double *, double *, double, double, double, double);
void get_maxsdir_stress_drops2(double *, double, double *);
void get_maxsdir_stress_drops(double *, double, double *, double, double *, double *);
void calculate_vertices(double *, struct flt *, double *, double *);
void calculate_bloated_vertices(double *, struct flt *, double);
int nvert_of_patch(struct flt *);
int ncon_of_subpatch(struct flt *, int);
double projected_slip_major_to_minor_patch(struct flt *, int, int, int);
int vtk_type_of_patch(struct flt *, int);
int number_of_subpatches(struct flt *);
int node_number_of_subelement(struct flt *, int, int);
void calculate_quad_vertices(double *, struct flt *, double);
void calculate_tri_vertices(double *, struct flt *, double);
void calculate_iquad_vertices(double *, struct flt *, double);
void calculate_point_source_vertices(double *, struct flt *, double, double *, double *);
void calculate_seg_vertices(double *, struct flt *, double);
double quad_area(double *);
double triangle_area(double *);
double triangle_area_gh(double *, double *);
void get_gh_tri_vec(double *, double *, double *);
void get_gh_quad_vec(double *, double *, double *, double *);
void check_fault_normal_vectors(struct flt *);
unsigned char check_planar(double *);
void calc_group_geometry(struct med *, struct flt *, struct geog *);
void vec_to_angles(double *, double *, double *);
void angles_to_vec(double, double, double *);
double vec_to_strike(double *);
double vec_to_dip(double *);
void check_fault_angles(struct flt *);
void check_angles(double *, double *);
void fix_azimuth(double *);
void globalx(double *, double, double, double *);
void calc_centroid_tri(double *, double *);
void calc_tri_bary_coord(double *, double *, double, double, double);
void calc_mean_quad_coord(double *, double *);
void calc_centroid_quad(double *, double *);
unsigned char patch_is_2d(unsigned char);
void calculate_position_of_patch(struct med *, struct flt *);
void compute_cartesian_slip(double *, double *, struct flt *);
void background_stress(double [3][3], double *, double, double *, double *, double);
void background_disp(double *, double *, struct med *, double *, double *);
void get_local_x_on_plane(double *, double *, double *, double *, double *);
void get_fault_plane_basevec(double *, double *, double *, struct flt *, struct med *);
void calc_deviatoric_stress(double [3][3], double [3][3], double *, double *);
void get_sub_normal_vectors(struct flt *, int, double *, double *, double *, double *);
unsigned char is_triangular(unsigned char);
void calc_global_strike_dip_from_local(struct flt *, double *, double *, double *);
void calc_global_slip_and_traction_from_local(struct flt *, double *, double *, double *, double *, double *, double *, double *, unsigned char);
/* geoproject.c */
void geoproject(double *, double *, int, double, double, double, double, double, double, double, int);
/* get_projected_fault_parameters.c */
void get_projected_fault_parameters(double [2][2], double, double *, double *, double *, double *, double *, double *);
/* help_and_comments.c */
void phelp(void);
char *comment_on_code(short int);
char *comment_on_code_bc(short int, double);
/* init.c */
void check_parameters_and_init(int, char **, struct med **, struct flt **, unsigned char *, double *, double *);
void initialize(struct med **, struct flt **, unsigned char, int, unsigned char, unsigned char, double, double *, double *, unsigned char, unsigned char, unsigned char, double, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, short int, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, double, unsigned char, unsigned char, unsigned char);
void init_files(struct med **, struct flt **);
void init_parameters(char **, int, unsigned char *, unsigned char *, unsigned char *, unsigned char *, double *, int *, unsigned char *, unsigned char *, unsigned char *, double *, unsigned char *, double *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, short int *, unsigned char *, double *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, double *, unsigned char *, unsigned char *, unsigned char *, int);
void advance_argument(int *, int, char **);
char *name_boolean(unsigned char);
unsigned char toggle(unsigned char *);
void read_stress_fac(unsigned char, double *, double *, double, struct med *);
/* input.c */
int read_moment_file(float **, float **, float *, float *, int *, unsigned char);
int read_patch_event_file(float *, int *, int *, float *, float *, FILE *, struct med *);
int write_patch_event_file(float, int, int, float, float *, FILE *);
void read_rsf(char *, struct med *, struct flt *);
/* interact.c */
void calc_interaction_matrix(struct med *, struct flt *, unsigned char);
double interaction_coefficient(int, int, int, int, struct flt *, int *);
void get_right_slip(double *, int, double);
double ic_from_file(int, int, int, int, struct med *);
double aij_from_file(int, int, int, FILE *);
int select_i_coeff_calc_mode(struct med *);
size_t imatrix_size(struct med *);
/* kdtree.c */
/* levmarq_numrec.c */
void mrqmin(double *, double *, int, double *, int *, int, double **, double **, double *, double *, struct bmd *, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, double, double *, unsigned char, struct prj);
void print_lm_progress(double, double, double, double, int, int, double *, double, int, char **, int, int, unsigned char, unsigned char);
void covsrt(double **, int, int *, int);
void gaussj(double **, int, double **, int);
void mrqcof(double *, double *, int, double *, int *, int, double **, double *, double *, struct bmd *, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, double, double *, unsigned char, struct prj);
void nrerror(char *);
double **matrix(long, long, long, long);
void free_matrix(double **, long, long, long, long);
/* llgeo.c */
void rotate_vec(double *, double *, double, double);
void rotate_vec2d(double *, double *, double, double);
void rotate_mat(double [3][3], double [3][3], double [3][3]);
void rotate_mat_z(double [3][3], double [3][3], double, double);
double tensor3d_norm(double [3][3]);
void zero_3x3_matrix(double [3][3]);
void scale_3x3_matrix(double [3][3], double);
void add_ay_to_3x3_matrix(double [3][3], double [3][3], double);
void add_y_to_3x3_matrix(double [3][3], double [3][3]);
double distance_3d(double *, double *);
double distance(double *, double *, int);
float distance_float(float *, float *, int);
void assign_to_vector(float *, float, int);
void assign_to_vector_3d(float *, float);
double distance_squared_3d(double *, double *);
double distance_squared(double *, double *, int);
float distance_squared_float(float *, float *, int);
double norm(double *, int);
double norm_3d(double *);
double dotp(double *, double *, int);
double dotp_3d(double *, double *);
void calc_Ax_ftn(double *, int, int, double *, double *);
void calc_Ax(double *, int, int, double *, double *);
void calc_AB_ftn(double *, int, int, double *, int, double *);
void normalize(double *, int);
void normalize_3d(double *);
void a_equals_b_vector(double *, double *, int);
void a_equals_b_vector_3d(double *, double *);
void swap_ab_vector_3d(double *, double *);
void scale_vector_3d(double *, double);
void scale_vector(double *, double, int);
void add_b_to_a_vector(double *, double *, int);
void add_b_to_a_vector_3d(double *, double *);
void sub_b_from_a_vector(double *, double *, int);
void sub_b_from_a_vector_3d(double *, double *);
void c_eq_a_plus_b_3d(double *, double *, double *);
void c_eq_a_minus_b_3d(double *, double *, double *);
void c_eq_a_minus_b(double *, double *, double *, int);
void cross_product(double *, double *, double *);
double find_max_vec(double *, int);
double find_max_abs_vec(double *, int);
float find_max_abs_vec_float(float *, int);
int find_lde_max(double *, int, int, double);
double square(double);
double tracemat9(double *);
int inv_compare_flt(const void *, const void *);
int compare_flt(const void *, const void *);
double rms(double *, int);
void my_vecalloc(double **, int, char *);
void my_ivecalloc(int **, int, char *);
void my_vecrealloc(double **, int, char *);
void my_ivecrealloc(int **, int, char *);
double mean(double *, int, int);
double mean_abs(double *, int, int);
double wmean(double *, int, int, double *);
double wmean_abs(double *, int, int, double *);
void calc_vec_stat(double *, int, double *);
double stddev(double, double, int);
int countzero_vec(double *, int);
double reformat_small(double);
/* lusolve.c */
void lu_driver(double *, double *, double *, int, int, struct med *);
/* makefault.c */
/* matrixio.c */
void print_matrix_ftrn(double *, int, int, FILE *, unsigned char);
void print_matrix_C(double *, int, int, FILE *, unsigned char);
void print_matrix_scaled_ftrn(double *, int, int, FILE *, unsigned char, double);
void print_matrix_ftrn_file(double *, int, int, char *, unsigned char);
void print_a_matrix(double *, int, int, FILE *, double *, unsigned char);
void read_a_matrix_from_file(double *, int, int, char *, char *);
void print_a_matrix_to_file(double *, int, int, char *, char *);
void print_system(double *, double *, double *, int, int, FILE *);
void print_sym3x3matrix(double [3][3], FILE *);
void print_interaction_matrix(struct med *, struct flt *, unsigned char);
void print_reduced_interaction_matrix(struct med *, struct flt *);
void print_vector(double *, int, FILE *);
void print_vector_file(double *, int, char *, char *);
void print_vector_row(double *, int, FILE *);
void print_b_vector(double *, int, FILE *, double *, unsigned char);
void print_3x3_matrix(double [3][3], FILE *);
/* mspectral.c */
void find_range(int *, int *, float, float, float, float, float *, float *, float *, int);
/* myopen.c */
FILE *myopen(char *, char *);
/* myprojectsimple.c */
void myprojectsimple(double *, double *, double, double, double, int);
double oblique_setup(double, double, double *, double, double, double *, GMT_LONG);
void make_euler_matrix(double *, double *, double);
void matrix_3v(double *, double *, double *);
void matrix_2v(double *, double *, double *);
void sphere_project_setup(double, double, double *, double, double, double *, double, double *, double *, GMT_LONG);
void oblique_transform(double, double, double *, double *, double *, double *);
/* mysincos.c */
void my_sincos_deg(double *, double *, double);
void my_sincos_deg_ftnd_(double *, double *, double *);
void my_sincos_ftn_(double *, double *, double *);
void my_sincos(double *, double *, double);
void my_sincosd(double *, double *, double);
void my_sincos_degd(double *, double *, double);
/* nnls.c */
void nnls_driver_i(double *, double *, double *, int, int);
/* optimize.c */
void optimize(struct flt *, struct med *);
double distsq(struct flt *, struct flt *);
double max_dist(struct flt *, struct med *);
double penalty_dist(struct flt *, struct med *, double);
/* output.c */
int mysystem(const char *);
void print_slip_line(struct med *, struct flt *);
void flush_slipline(struct med *, struct flt *);
void print_fault_stress_and_slip(struct med *, struct flt *, unsigned char);
void print_fault_data(char *, struct med *, struct flt *);
void print_fault_stress_stat(FILE *, int, struct med *, struct flt *);
void print_group_data_geom(char *, struct med *, struct flt *, int, int, double);
double select_val_for_print(struct flt *, int);
void print_stress(struct med *, struct flt *);
void print_stress_on_fault(struct med *, struct flt *, int);
void print_displacement(struct med *, struct flt *);
void print_equations(int, unsigned char *, int *, double *, int, char *, struct flt *);
void print_solutions(int, int *, struct flt *, struct med *, char *);
void flush_moment_stack(struct med *);
void fiddle_with_limits_for_plot(struct med *, int *, unsigned char *, double *, unsigned char);
void time_report(char *, char *, struct med *);
/* patch2area.c */
/* patch2bc.c */
/* patch2dis3d.c */
void get_dis3d_parameters(double, double, double, double, double, double, double *, double *, double *, double *);
/* patch2geom.c */
/* patch2group.c */
/* patch2poly3d.c */
/* patch2vertices.c */
/* patch2vtk.c */
/* patch2xyz.c */
/* patch2xyzvec.c */
/* patchquad2patchtri.c */
/* period.c */
void period(float *, float *, int, float, float, float *, float *, int, int *, int *, float *);
void avevar(float *, int, float *, float *);
void fasper(float *, float *, int, float, float, float *, float *, int, int *, int *, float *);
void spread(float, float *, int, float, int);
void realft(float *, int, int);
void four1(float *, int, int);
/* petsc_interact.c */
/* petsc_simple_solve.c */
/* plotevents.c */
/* plotting.c */
/* points2patch.c */
unsigned char read_points_local(double *, int *, unsigned char, FILE *);
/* print_patch_geometry.c */
int print_patch_geometry_and_bc(int, struct flt *, int, double, unsigned char, FILE *, unsigned char, double *);
/* project_stress.c */
/* quake.c */
void quake(unsigned char *, double *, int, struct flt *, struct med *, unsigned char, unsigned char);
void add_quake_stress(int, unsigned char *, double *, struct flt *, struct med *);
/* randgen.c */
void assign_random(double *, double, double, double, long *, int);
double myrand(long *);
double myrandnr(double, long *);
int myrandi(int, long *);
double mygauss_randnr(double, long *);
double mypower_randnr(double, double, double, int, long *);
double gasdev(long *);
double ran1(long *);
double ran2(long *);
/* randomflt.c */
void check_input_parameters(int, char **, int *, long *, int *, double *, unsigned char *, unsigned char *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, unsigned char *, int *, double *, unsigned char *, double *, double *);
/* randomize_list.c */
void randomize_list(int **, int, unsigned char);
int slist_sort(const void *, const void *);
/* randomize_strike.c */
/* read_bin_events.c */
/* read_bin_events_old.c */
/* read_boundary_conditions.c */
void read_boundary_conditions(struct med *, struct flt *, double *, double *, unsigned char);
void read_one_step_bc(FILE *, struct med *, struct flt *, double *, double *, unsigned char);
unsigned char slip_type_bc(int);
unsigned char bc_is_directionally_unconstrained(unsigned char);
unsigned char bc_contains_dip_motion(unsigned char);
/* read_fltdat.c */
int read_fltdat(char *, struct flt *, struct med *, unsigned char);
/* read_geometry.c */
void read_geometry(char *, struct med **, struct flt **, unsigned char, unsigned char, unsigned char, unsigned char);
/* read_stress_observations.c */
void read_stress_observations(struct bmd *, double *, double, unsigned char, double **, double, unsigned char);
/* regular_interact_main.c */
/* restart.c */
void adjust_medium_for_restart(struct med *, struct flt *);
/* rhs.c */
void init_equation_system(struct med *, struct flt *);
void add_to_active_fault_list(int, int **, int *, unsigned char **);
void add_to_right_hand_side(double, double **, double **, int *);
/* rsf_solve.c */
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
/* simul_gr.c */
/* solve.c */
int solve(struct med *, struct flt *);
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
/* stress_aux.c */
double stress_misfit(double [3][3], double [3][3]);
void stress_vec_from_hstate(double, double, double, double, int, double *, double *);
void expand_stress_matrix6to9(double *);
void convert_6sym_to_9_matrix(double *, double [3][3]);
double resolved_stress(double *, double [3][3], double *);
void calc_three_stress_components(double [3][3], double *, double *, double *, double *, double *, double *, double *);
/* string_compare.c */
unsigned char strings_match(char *, char *);
/* svd.c */
void svd_driver_lapack(double *, double *, double *, int, int, double *, int, double *, double **, int, unsigned char, unsigned char);
void svd_driver_numrec(double *, double *, double *, int, int, double *, int, double *, unsigned char, double **, double **, int, unsigned char, int, unsigned char, unsigned char, unsigned char);
void assemble_cov_svd(double *, int, double *, double *);
void print_singular_values(double *, int, FILE *);
void indexx(int, double *, int *);
void reduce_a_matrix(double **, int, int, int);
/* terminate.c */
int terminate(struct med *, struct flt *);
/* test_optimize.c */
/* test_solvers.c */
double mat_value(int, int, int);
/* test_sparse.c */
/* test_stuff.c */
/* test_triangle_stress.c */
/* tri2patch.c */
/* trigonometry.c */
double dist_on_sphere(double, double, double, double);
double dist_on_sphere_deg(double, double, double, double);
double azimuth(double, double, double, double);
void get_gc_pole(double, double, double, double, double *, double *);
void get_point_on_gc(double, double, double, double, double, double *, double *, double *);
void get_point_on_course(double, double, double, double, double *, double *);
void lonlat2xyz(double, double, double *);
void lonlat2xyz_deg(double, double, double *);
void xyz2lonlat(double *, double *, double *);
void xyz2lonlat_deg(double *, double *, double *);
void pv2cv(double *, double *, double *);
void cv2pv(double *, double *, double *);
void calculate_polar_base(double, double, double *);
/* tritri.c */
int coplanar_tri_tri(float [3], float [3], float [3], float [3], float [3], float [3], float [3]);
int tri_tri_intersect(float [3], float [3], float [3], float [3], float [3], float [3]);
