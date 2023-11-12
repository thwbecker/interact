/* adjust_time_step.c */
void adjust_time_step(struct flt *, struct med *);
/* block_checkflt.c */
/* block_compute_vel_from_omega.c */
/* block_eval_geokada.c */
void block_eval_geookada(float *, float *, float *, float [3][3], float, float, float, float, float, float, float, float, float, float, float, int *, unsigned short, unsigned short);
void block_save_solution_and_faults(float *, int, int, struct bflt *, float *, struct prj *, FILE *, unsigned short, unsigned short);
void block_load_solution_and_faults(float **, int *, int *, struct bflt **, float **, struct prj **, FILE *, unsigned short *, unsigned short *);
void block_eval_blockvec(float *, float, int, float *, struct prj *, float *);
/* block_evaluate_solution.c */
/* blockinvert.c */
/* block_levmarq.c */
void run_lm(struct bmd *, long int *, struct prj, float *, unsigned short, unsigned short, float, float *, unsigned short, unsigned short, float, unsigned short, unsigned short, unsigned short, int, int, char **);
/* block_matrix.c */
void assemble_block_fltdep_matrices(float **, float **, float **, float **, float *, float *, float *, int, int, int, int, int, int, int, float *, struct bflt *, struct bck *, unsigned short, unsigned short);
void assemble_block_a(float **, float *, int *, int, int, struct bck *, int, struct bmd *);
void assemble_block_d(float **, struct bflt *, int, int, int, struct bmd *);
void assemble_block_i(float **, struct bflt *, int, int, int);
void assemble_block_g(float **, struct bflt *, int, int);
void assemble_block_f(float **, struct bflt *, int, int, struct bck *, float *, unsigned short);
void assemble_block_k(float **, int, int, int, int, int, int, int, float *, float *, float *, unsigned short, unsigned short, float *, struct bflt *, int, int, int, int, float);
void assemble_stress_matrix(float *, int, struct bflt *, float *, int, int);
void block_assemble_fit_vector(float **, int, int, int, int, int, int, int, int, int, float *, int, float *, float *, unsigned short, float *, float **, float *, float **, int, struct bck *, int, struct bflt **, unsigned short, float **, unsigned short, unsigned short, unsigned short, float *, float *, float *, struct prj, unsigned short, int, int, float, struct bmd *);
void block_assemble_dyda_matrix(float **, float **, int, int, int, int, int, int, int, int, int, float *, float *, float **, float *, float **, float **, int, struct bck *, int, struct bflt **, unsigned short, unsigned short, unsigned short, float *, float *, float *, struct prj, unsigned short, int, int, float, struct bmd *);
float block_chi_square(float *, float *, float *, int, int, float, float *, float *);
float block_data_norm(float *, float *, int, int, float, struct bmd *);
int block_slip_direction(int, struct bflt *);
void assign_additional_sol_values(float *, int, int, int, struct bflt *, unsigned short, unsigned short, long int *, int);
void change_locking_depths(struct bflt **, int, int, int, int, float *, float *, float *, unsigned short, struct prj, float **, float **, float *, struct bmd *);
void calc_fault_sn_rms(float *, int, int, float *, float *);
void init_block_mods(struct bmd **);
void init_blocks(struct bmd *, int);
void convert_cart_sol(float *, float *, struct bmd *);
void check_solution_vector(float *, int, int, unsigned short, unsigned short);
void copy_block(struct bck *, struct bck *);
/* block_output.c */
void block_output(struct bmd *, unsigned short, char **, struct prj, unsigned short, unsigned short, float, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short);
void sort_eigen(int *, float *);
void print_horizontal_stress(float *, FILE *);
void print_projected_stress(float *, FILE *);
void calculate_slipsol_sigma(struct bmd *, float *, float *, int);
void calc_geo_euler_pole(float *, float *, float *, float *);
void print_simple_vel(float *, float *, float *, float *, int, char *);
/* block_read_bflt.c */
void read_bflt(struct bmd *, float, struct prj, unsigned short, float, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short);
void generate_new_fault(struct bflt **, int *, int *, float *, int, float *, int, float, struct prj, float *, unsigned short, float, unsigned short, unsigned short, float *, float, int, float *, unsigned short *, int *);
void get_bflt_intcoeff(struct bflt **, int, float *, int, float *, int, float *, unsigned short);
void free_bflt(struct bflt **, int, int);
void assign_bflt_dip_mode(struct bflt *, float);
void flip_block_code(struct bflt *, float *);
void print_fault_geometry(struct bflt *, int, FILE *);
void assign_fault_locking_depth_parameters(struct bflt *, float, struct prj, unsigned short, int);
int new_fault(struct bflt **, int);
void init_bflt(struct bflt *);
/* block_read_euler.c */
void read_constrained_euler_poles(struct bmd *, char **, int);
/* block_read_gps.c */
void read_gps_velocities(struct bmd *, struct prj *, float *, char **, unsigned short, int *, int *);
void project_gps_coordinates(float *, float *, int, int, int *, struct bck *, int, struct prj *, char **);
void block_read_centroids(float **, int *);
void find_spherical_rotation(struct bmd *, int, float *, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, float *, int);
void init_gps_pbase(struct bmd *);
void init_cart_gps_velsig(struct bmd *, unsigned short);
void calculate_c2p_gps(float *, float *, struct bmd *);
void convert_sig_p2c(float *, float *, float *);
void remove_cvel_from_vel(struct bmd *, float *, int);
char bname(int);
unsigned short node_in_geo_list(float *, int, float *);
/* block_solve.c */
void solve_block(float *, float *, float *, int, int, int, int, unsigned short, int, int, int, int, float *, float, float *, float **, float **, float *, float *, unsigned short, int, unsigned short, int, unsigned short);
void evaluate_block_solution(float *, int, int, float *, float *, float *, struct bmd *);
/* block_stress.c */
void block_scale_stresses(float *, float *, float *, int, int, int, float, float *, unsigned short);
void calc_horizontal_stress_vec(float *, float *, int);
void calc_horizontal_stress(float *, float *, float *, float *);
void calc_dir_diff_vec(float *, float *, float *, int, unsigned short);
float calc_dir_diff(float, float, unsigned short);
void cart_mat_from_horsym(float, float, float, float *);
void rescale_observed_stresses(float *, float *, float *, float, float *, unsigned short, struct bmd *, unsigned short, unsigned short);
/* calc_cart_from_eigen_stress.c */
unsigned short read_vecs(int, float *, float *, float *, float *, float *);
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
void calc_spatial_correlation(struct flt *, int, int, int, int, float **, float **, int **, float, float, float *, float **, float **, float **);
float correlation_coefficient(float *, float *, int);
/* calc_stress.c */
void initialize_stress_state(struct flt *, struct med *, unsigned short, float *, float *);
void update_stress_state(struct flt *, struct med *);
void calc_fields(struct med *, struct flt *, unsigned short, unsigned short, float *, float *);
/* calc_stress_stat.c */
/* check_feedback.c */
/* check_interaction.c */
unsigned short incrementally_check_interaction_coefficients(struct flt *, int, struct med *, unsigned short, int *, float);
unsigned short check_coulomb_stress_feedback(int, int, struct flt *, struct med *, unsigned short, int, unsigned short, int *, float);
/* compare_fault.c */
int compare_fault_length(const void *, const void *);
int compare_fault_width(const void *, const void *);
/* coulomb_stress.c */
float coulomb_stress(float, float, float, float);
float cstress_drop(float, float, float, float);
unsigned short in_coulomb_compress_regime(float);
/* create_random_mu_file.c */
/* create_random_stress_file.c */
void get_random_stress(float *, float *, float *, int, long *, int);
/* divide_fault_in_patches.c */
void create_patches(int, struct flt *, struct flt **, int *, int *, unsigned short, unsigned short, float *, float, float, long *);
void determine_segments(int *, int *, struct flt *, unsigned short, float *);
void divide_fault_in_patches(int, struct flt *, struct flt **, int *, int *, unsigned short, unsigned short, float, float, long *);
void get_flt_location(struct flt *, float *, float *, float *, int, int);
void randomize_strike_dip(float, float, struct flt *, long *);
/* eigensystem.c */
void eigensystem3d(float [3][3], float *, float *);
void calc_eigensystem_sym3d(float *, float *, float *, unsigned short);
void calc_eigensystem_sym2d(float *, float *, float *, unsigned short);
void eispack_driver(float *, int, int, float *, float *, float *, float *, int);
/* eval_2dsegment.c */
void eval_2dsegment_plane_strain(float *, struct flt *, float *, float *, float [3][3], int *);
void eval_2dsegment_plane_stress(float *, struct flt *, float *, float *, float [3][3], int *);
void eval_2dsegment_plane_strain_basic(float *, struct flt *, float *, float *, float [3][3], int *);
void eval_2dsegment_plane_stress_basic(float *, struct flt *, float *, float *, float [3][3], int *);
void get_2dseg_geo(float *, float, float *, float *, float *, float *, float *, float *, float *, float *, float *, int *);
void get_2dseg_ffac(float *, float *, float *, float *, float *, float *, float, float, float, float, float, float, float, float, float, float, float *);
void get_2dseg_disp(float *, float *, float *, float, float, float, float, float, float);
void get_2dseg_stress(float [3][3], float *, float *, float, float, float, float);
void eval_2dsegment_plane_strain_tdd(float *, struct flt *, float *, float *, float [3][3], int, int *);
/* eval_green.c */
void eval_green(float *, struct flt *, float *, float *, float [3][3], int *);
void eval_green_basic(float *, struct flt *, float *, float *, float [3][3], int *);
/* eval_okada.c */
void eval_rectangle(float *, struct flt *, float *, float *, float [3][3], int *);
void eval_rectangle_basic(float *, float, float, float, float, float *, float *, float [3][3], int *);
void eval_point(float *, struct flt *, float *, float *, float [3][3], int *);
void eval_point_short(float *, float *, float, float, float, float, float *, float *, float [3][3], int *);
void set_stress_and_disp_nan(float [3][3], float *);
/* eval_triangle.c */
void get_gauss_points(float *, float *, float *, int);
/* far_enough.c */
unsigned short far_enough(struct flt *, struct flt *, float);
/* fit_mean_stress.c */
void eval_stress_correction(float *, float *, float *, float *, int, float *, float *, float *, float *, float *, float *);
/* fit_plane.c */
void fit_plane(int, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, unsigned short, unsigned short);
void points2patch(struct flt *, float *, unsigned short);
/* fit_simple_stress_from_cart.c */
/* fltcopy.c */
void fltswap(struct flt *, struct flt *);
void fltcp(struct flt *, struct flt *);
/* fracture_criterion.c */
unsigned short fracture_criterion(int, struct flt *, float *, unsigned short *, struct med *);
void two_dir_slip_check(unsigned short *, float *, unsigned short *, int, struct flt *, struct med *);
void deactivate_patch(int, struct flt *, struct med *);
void deactivate_group(int, int, struct flt *, struct med *);
int calc_absolute_shear_stress_and_inc(float *, float *, int, struct flt *);
int calc_absolute_shear_stress(float *, int, struct flt *);
/* fstress2eig.c */
/* fstress2hor.c */
/* generate_random_2d.c */
/* generate_slipdia.c */
/* geometry.c */
void calc_lhemi_proj(float, float, float *);
void resolve_force(float *, float [3][3], float *);
void calc_base_vecs(float *, float *, float *, float, float, float, float);
void get_maxsdir_stress_drops2(float *, float, float *);
void get_maxsdir_stress_drops(float *, float, float *, float, float *, float *);
void calculate_corners(float [4][3], struct flt *, float *, float *);
void calculate_bloated_corners(float [4][3], struct flt *, float);
void calculate_quad_corners(float [4][3], struct flt *, float);
void calculate_point_source_corners(float [4][3], struct flt *, float, float *, float *);
void calculate_seg_corners(float [4][3], struct flt *, float);
float quad_area(float *);
float triangle_area(float *);
float triangle_area_gh(float *, float *);
void get_gh_tri_vec(float *, float *, float *);
void get_gh_quad_vec(float *, float *, float *, float *);
unsigned short check_planar(float *);
void get_alpha_dip_tri_gh(float *, float *, float *, float *, float *);
void calc_group_geometry(struct med *, struct flt *, struct geog *);
void vec_to_angles(float *, float *, float *);
void angles_to_vec(float, float, float *);
float vec_to_strike(float *);
float vec_to_dip(float *);
void check_fault_angles(struct flt *);
void check_angles(float *, float *);
void fix_azimuth(float *);
void globalx(float *, float, float, float *);
void calc_centroid_tri(float *, float *);
void calc_mean_quad_coord(float *, float *);
void calc_centroid_quad(float *, float *);
void calculate_position_of_patch(struct med *, struct flt *);
void compute_cartesian_slip(float *, float *, struct flt *);
void background_stress(float [3][3], float *, float, float *, float *, float);
void background_disp(float *, float *, struct med *, float *, float *);
void get_local_x_on_plane(float *, float *, float *, float *, float *);
void get_fault_plane_basevec(float *, float *, float *, struct flt *, struct med *);
void calc_deviatoric_stress(float [3][3], float [3][3], float *, float *);
/* geo_okada.c */
/* get_projected_fault_parameters.c */
void get_projected_fault_parameters(float [2][2], float, float *, float *, float *, float *, float *, float *);
/* help_and_comments.c */
void phelp(void);
char *comment_on_code(short int);
char *comment_on_code_bc(short int, float);
/* init.c */
void check_parameters_and_init(int, char **, struct med **, struct flt **, unsigned short *, float *, float *);
void initialize(struct med **, struct flt **, unsigned short, int, unsigned short, unsigned short, float, float *, float *, unsigned short, unsigned short, unsigned short, float, unsigned short, float, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, short int, unsigned short, float, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, float, unsigned short, unsigned short);
void init_files(struct med **, struct flt **);
void init_parameters(char **, int, unsigned short *, unsigned short *, unsigned short *, unsigned short *, float *, int *, unsigned short *, unsigned short *, unsigned short *, float *, unsigned short *, float *, unsigned short *, unsigned short *, unsigned short *, unsigned short *, unsigned short *, unsigned short *, unsigned short *, short int *, unsigned short *, float *, unsigned short *, unsigned short *, unsigned short *, unsigned short *, unsigned short *, float *, unsigned short *, unsigned short *, int);
void advance_argument(int *, int, char **);
char *name_boolean(unsigned short);
unsigned short toggle(unsigned short *);
void read_stress_fac(unsigned short, float *, float *, float, struct med *);
/* input.c */
int read_moment_file(float **, float **, float *, float *, int *, unsigned short);
int read_patch_event_file(float *, int *, int *, float *, float *, FILE *, struct med *);
int write_patch_event_file(float, int, int, float, float *, FILE *);
/* interact.c */
void calc_interaction_matrix(struct med *, struct flt *, unsigned short);
float interaction_coefficient(int, int, int, int, struct flt *, int *);
void get_right_slip(float *, int, float);
float ic_from_file(int, int, int, int, struct med *);
float aij_from_file(int, int, int, FILE *);
int select_i_coeff_calc_mode(struct med *);
size_t imatrix_size(struct med *);
/* levmarq_numrec.c */
void mrqmin(float *, float *, int, float *, int *, int, float **, float **, float *, float *, struct bmd *, unsigned short, float, unsigned short, unsigned short, unsigned short, unsigned short, float, float *, unsigned short, struct prj);
void print_lm_progress(float, float, float, float, int, int, float *, float, int, char **, int, int, unsigned short, unsigned short);
void covsrt(float **, int, int *, int);
void gaussj(float **, int, float **, int);
void mrqcof(float *, float *, int, float *, int *, int, float **, float *, float *, struct bmd *, unsigned short, float, unsigned short, unsigned short, unsigned short, unsigned short, float, float *, unsigned short, struct prj);
void nrerror(char *);
float **matrix(long, long, long, long);
void free_matrix(float **, long, long, long, long);
/* llgeo.c */
void rotate_vec(float *, float *, float, float);
void rotate_vec2d(float *, float *, float, float);
void rotate_mat(float [3][3], float [3][3], float [3][3]);
void rotate_mat_z(float [3][3], float [3][3], float, float);
float project_vector(float *, float *);
float distance_3d(float *, float *);
float distance(float *, float *, int);
float distance_float(float *, float *, int);
void assign_to_vector(float *, float, int);
void assign_to_vector_3d(float *, float);
float distance_squared_3d(float *, float *);
float distance_squared(float *, float *, int);
float distance_squared_float(float *, float *, int);
float norm(float *, int);
float norm_3d(float *);
float dotp(float *, float *, int);
float dotp_3d(float *, float *);
void calc_Ax_ftn(float *, int, int, float *, float *);
void calc_Ax(float *, int, int, float *, float *);
void calc_AB_ftn(float *, int, int, float *, int, float *);
void normalize(float *, int);
void normalize_3d(float *);
void a_equals_b_vector(float *, float *, int);
void a_equals_b_vector_3d(float *, float *);
void swap_ab_vector_3d(float *, float *);
void scale_vector_3d(float *, float);
void scale_vector(float *, float, int);
void add_b_to_a_vector(float *, float *, int);
void add_b_to_a_vector_3d(float *, float *);
void sub_b_from_a_vector(float *, float *, int);
void sub_b_from_a_vector_3d(float *, float *);
void c_eq_a_plus_b_3d(float *, float *, float *);
void c_eq_a_minus_b_3d(float *, float *, float *);
void c_eq_a_minus_b(float *, float *, float *, int);
void cross_product(float *, float *, float *);
float find_max_vec(float *, int);
float find_max_abs_vec(float *, int);
float find_max_abs_vec_float(float *, int);
int find_lde_max(float *, int, int, float);
float square(float);
float tracemat9(float *);
int inv_compare_flt(const void *, const void *);
int compare_flt(const void *, const void *);
float rms(float *, int);
void my_vecalloc(float **, int, char *);
void my_ivecalloc(int **, int, char *);
void my_vecrealloc(float **, int, char *);
void my_ivecrealloc(int **, int, char *);
float mean(float *, int, int);
float mean_abs(float *, int, int);
float wmean(float *, int, int, float *);
float wmean_abs(float *, int, int, float *);
void calc_vec_stat(float *, int, float *);
float stddev(float, float, int);
int countzero_vec(float *, int);
float reformat_small(float);
/* lusolve.c */
void lu_driver(float *, float *, float *, int, int, struct med *);
/* main.c */
/* makefault.c */
/* matrixio.c */
void print_matrix_ftrn(float *, int, int, FILE *, unsigned short);
void print_matrix_C(float *, int, int, FILE *, unsigned short);
void print_matrix_scaled_ftrn(float *, int, int, FILE *, unsigned short, float);
void print_matrix_ftrn_file(float *, int, int, char *, unsigned short);
void print_a_matrix(float *, int, int, FILE *, float *, unsigned short);
void read_a_matrix_from_file(float *, int, int, char *, char *);
void print_a_matrix_to_file(float *, int, int, char *, char *);
void print_system(float *, float *, float *, int, int, FILE *);
void print_sym3x3matrix(float [3][3], FILE *);
void print_interaction_matrix(struct med *, struct flt *);
void print_reduced_interaction_matrix(struct med *, struct flt *);
void print_vector(float *, int, FILE *);
void print_vector_file(float *, int, char *, char *);
void print_vector_row(float *, int, FILE *);
void print_b_vector(float *, int, FILE *, float *, unsigned short);
void print_3x3_matrix(float [3][3], FILE *);
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
void my_sincos_deg(float *, float *, float);
void my_sincos_deg_ftn(double *, double *, double *);
void my_sincos(float *, float *, float);
/* nnls.c */
void nnls_driver_i(float *, float *, float *, int, int);
/* optimize.c */
void optimize(struct flt *, struct med *);
float distsq(struct flt *, struct flt *);
float max_dist(struct flt *, struct med *);
float penalty_dist(struct flt *, struct med *, float);
/* output.c */
void print_slip_line(struct med *, struct flt *);
void flush_slipline(struct med *, struct flt *);
void print_fault_stress_and_slip(struct med *, struct flt *, unsigned short);
void print_fault_data(char *, struct med *, struct flt *);
void print_fault_stress_stat(FILE *, int, struct med *, struct flt *);
void print_group_data_geom(char *, struct med *, struct flt *, int, int, float);
float select_val_for_print(struct flt *, int);
void print_stress(struct med *, struct flt *);
void print_stress_on_fault(struct med *, struct flt *, int);
void print_displacement(struct med *, struct flt *);
void print_equations(int, unsigned short *, int *, float *, int, char *, struct flt *);
void print_solutions(int, int *, struct flt *, struct med *, char *);
void flush_moment_stack(struct med *);
void fiddle_with_limits_for_plot(struct med *, int *, unsigned short *, float *, unsigned short);
/* patch2bc.c */
/* patch2corners.c */
/* patch2dis3d.c */
void get_dis3d_parameters(float, float, float, float, float, float, float *, float *, float *, float *);
/* patch2geom.c */
/* patch2group.c */
/* patch2poly3d.c */
/* patch2vtk.c */
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
void drawset(int *, int, int, int, struct pa *);
/* plotting.c */
void update_plots(struct med *, struct flt *);
void init_plot_window(struct med *, struct flt *);
void close_plot_window(struct med *, struct flt *);
void plot_time_label(struct med *, struct flt *, int);
void plot_patch(int, struct med *, struct flt *, int, int);
void plot_projected_patch(int, struct med *, struct flt *, int, int);
void plot_time_tics(struct med *, struct flt *, float);
void add_to_plotting_list(int, int **, int *);
void plot_quake(int, struct med *, struct flt *);
void plot_moment_array(struct med *, float *, int);
void psticks(float *, float *, int);
/* points2patch.c */
unsigned short read_points_local(float *, int *, unsigned short, FILE *);
/* print_patch_geometry.c */
void print_patch_geometry_and_bc(int, struct flt *, int, float, unsigned short, FILE *, unsigned short, float *);
/* project_stress.c */
/* quake.c */
void quake(unsigned short *, float *, int, struct flt *, struct med *, unsigned short, unsigned short);
void add_quake_stress(int, unsigned short *, float *, struct flt *, struct med *);
/* randgen.c */
void assign_random(float *, float, float, float, long *, int);
float myrand(long *);
float myrandnr(float, long *);
int myrandi(int, long *);
float mygauss_randnr(float, long *);
float mypower_randnr(float, float, float, int, long *);
double gasdev(long *);
double ran1(long *);
double ran2(long *);
/* randomflt.c */
void check_input_parameters(int, char **, int *, long *, int *, float *, unsigned short *, unsigned short *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, int *, unsigned short *, int *, float *, unsigned short *, float *, float *);
/* randomize_list.c */
void randomize_list(int **, int, unsigned short);
int slist_sort(const void *, const void *);
/* randomize_strike.c */
/* read_bin_events.c */
/* read_boundary_conditions.c */
void read_boundary_conditions(struct med *, struct flt *, float *, float *, unsigned short);
void read_one_step_bc(FILE *, struct med *, struct flt *, float *, float *, unsigned short);
unsigned short slip_type_bc(int);
/* read_fltdat.c */
void read_fltdat(char *, struct flt *, struct med *, unsigned short);
/* read_geometry.c */
void read_geometry(char *, struct med **, struct flt **, unsigned short, unsigned short, unsigned short, unsigned short);
/* read_stress_observations.c */
void read_stress_observations(struct bmd *, float *, float, unsigned short, float **, float, unsigned short);
/* restart.c */
void adjust_medium_for_restart(struct med *, struct flt *);
/* rhs.c */
void init_equation_system(struct med *, struct flt *);
void add_to_active_fault_list(int, int **, int *, unsigned short **);
void add_to_right_hand_side(float, float **, float **, int *);
/* rupture.c */
unsigned short activate_faults(struct flt *, struct med *);
void fault_criterion(int, struct flt *, struct med *);
/* segseg.c */
void intersect(float *, float *, float *, float *, float *, int *);
char SegSegInt(tPointi, tPointi, tPointi, tPointi, tPointd);
char ParallelInt(tPointi, tPointi, tPointi, tPointi, tPointd);
void Assigndi(tPointd, tPointi);
unsigned short Between(tPointi, tPointi, tPointi);
int Collinear(tPointi, tPointi, tPointi);
int AreaSign(tPointi, tPointi, tPointi);
unsigned short crosses(double *, double *);
/* simul_gr.c */
/* solve.c */
int solve(struct med *, struct flt *);
void add_solution(int, unsigned short *, float *, int *, struct med *, struct flt *, unsigned short, unsigned short, float);
void assemble_a_matrix(float *, int, unsigned short *, int, int *, struct flt *, struct med *);
/* solve_mode_dependend.c */
void assemble_a_matrix_4(float *, int, unsigned short *, int, int *, struct flt *, struct med *);
void add_quake_stress_4(unsigned short *, float *, int, struct flt *, struct med *);
unsigned short check_coulomb_stress_feedback_4(int, int, struct flt *, struct med *, unsigned short, unsigned short, int *, float);
/* sort_events.c */
/* sparse.c */
size_t create_crs_sparse_from_memory(int, float *, float, unsigned int **, unsigned int **, float **);
size_t create_crs_sparse_from_file(int, float, unsigned int **, unsigned int **, float **, FILE *);
size_t create_ccs_sparse_from_memory(int, float *, float, unsigned int **, unsigned int **, float **);
size_t create_ccs_sparse_from_file(int, float, unsigned int **, unsigned int **, float **, FILE *);
/* sparse_nr.c */
size_t create_nrs_sparse(struct med *, float, float *, unsigned short);
size_t create_nrs_sparse_from_memory(int, float *, float, unsigned int **, float **);
size_t create_nrs_sparse_from_file(int, float, unsigned int **, float **, FILE *);
float get_nrs_sparse_el(int, int, unsigned int *, float *);
/* sparse_solve.c */
void sparse_driver(unsigned int *, unsigned int *, float *, float *, float *, int, struct med *);
/* stress_aux.c */
float stress_misfit(float [3][3], float [3][3]);
void stress_vec_from_hstate(float, float, float, float, int, float *, float *);
void expand_stress_matrix6to9(float *);
void convert_6sym_to_9_matrix(float *, float [3][3]);
float resolved_stress(float *, float [3][3], float *);
void calc_three_stress_components(float [3][3], float *, float *, float *, float *, float *, float *, float *);
/* string_compare.c */
unsigned short strings_match(char *, char *);
/* svd.c */
void svd_driver_lapack(float *, float *, float *, int, int, float *, int, float *, float **, int, unsigned short, unsigned short);
void svd_driver_numrec(float *, float *, float *, int, int, float *, int, float *, unsigned short, float **, float **, int, unsigned short, int, unsigned short, unsigned short, unsigned short);
void assemble_cov_svd(float *, int, float *, float *);
void print_singular_values(float *, int, FILE *);
void indexx(int, float *, int *);
void reduce_a_matrix(float **, int, int, int);
/* terminate.c */
void terminate(struct med *, struct flt *);
/* test_optimize.c */
/* test_solvers.c */
float mat_value(int, int, int);
/* test_sparse.c */
/* test_stuff.c */
/* tri2patch.c */
/* trigonometry.c */
float dist_on_sphere(float, float, float, float);
float dist_on_sphere_deg(float, float, float, float);
float azimuth(float, float, float, float);
void get_gc_pole(float, float, float, float, float *, float *);
void get_point_on_gc(float, float, float, float, float, float *, float *, float *);
void get_point_on_course(float, float, float, float, float *, float *);
void lonlat2xyz(float, float, float *);
void lonlat2xyz_deg(float, float, float *);
void xyz2lonlat(float *, float *, float *);
void xyz2lonlat_deg(float *, float *, float *);
void pv2cv(float *, float *, float *);
void cv2pv(float *, float *, float *);
void calculate_polar_base(float, float, float *);
/* tritri.c */
int coplanar_tri_tri(float [3], float [3], float [3], float [3], float [3], float [3], float [3]);
int tri_tri_intersect(float [3], float [3], float [3], float [3], float [3], float [3]);
