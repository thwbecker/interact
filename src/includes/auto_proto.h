/* src/calc_interaction_matrix.c */
void print_help_local(char *);
/* src/calc_stress.c */
void initialize_stress_state(struct flt *, struct med *, unsigned char, double *, double *);
void update_stress_state(struct flt *, struct med *);
void calc_fields(struct med *, struct flt *, unsigned char, unsigned char, double *, double *);
/* src/compress_interaction_matrix.c */
/* src/coulomb_stress.c */
double coulomb_stress(double, double, double, double);
double cstress_drop(double, double, double, double);
unsigned char in_coulomb_compress_regime(double);
/* src/eval_green.c */
void eval_green_and_project_stress_to_fault(struct flt *, int, int, double *, double *, unsigned char, unsigned char);
void eval_green(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char, unsigned char, unsigned char);
void eval_triangle_general(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char, unsigned char);
void eval_green_at_receiver(struct flt *, int, int, double *, double *, double [3][3], int *, unsigned char, unsigned char, unsigned char);
double resolve_stress_on_fault_using_ctx(double [3][3],  struct interact_ctx *, int );
void eval_green_basic(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
double ckernel_func(int, int, void *);
/* src/get_projected_fault_parameters.c */
void get_projected_fault_parameters(double [2][2], double, double *, double *, double *, double *, double *, double *);
/* src/init.c */
void check_parameters_and_init_interact(int, char **, struct med **, struct flt **, unsigned char *, double *, double *);
void initialize_interact(struct med **, struct flt **, unsigned char, unsigned char,int, unsigned char, unsigned char, double, double *, double *, unsigned char, unsigned char, unsigned char, double, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, short int, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char);
void init_files_interact(struct med **, struct flt **);
void init_parameters_interact(char **, int, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, double *, int *, unsigned char *, unsigned char *, unsigned char *, double *, unsigned char *, double *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, short int *, unsigned char *, double *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, double *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, int);
void advance_argument(int *, int, char **);
char *name_boolean(unsigned char);
unsigned char toggle(unsigned char *);
void read_stress_fac(unsigned char, double *, double *, double, struct med *);
/* src/interact_main.c */
/* src/interact_matrix_assembly.c */
void calc_interaction_matrix(struct med *, struct flt *, unsigned char);
double interaction_coefficient(int, int, int, int, struct flt *, int *, unsigned char);
void get_right_slip(double *, int, double,struct flt *);
double ic_from_file(int, int, int, int, struct med *);
double aij_from_file(int, int, int, FILE *);
int select_i_coeff_calc_mode(struct med *);
size_t imatrix_size(struct med *);
/* src/print_patch_geometry.c */
int print_patch_geometry_and_bc(int, struct flt *, int, double, unsigned char, FILE *, unsigned char, double *);
/* src/read_fltdat.c */
int read_fltdat(char *, struct flt *, struct med *, unsigned char);
/* src/read_geometry.c */
void read_geometry(char *, struct med **, struct flt **, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char);
/* src/read_stress_observations.c */
void read_stress_observations(struct bmd *, double *, double, unsigned char, double **, double, unsigned char);
/* src/interact/adjust_time_step.c */
void adjust_time_step(struct flt *, struct med *);
/* src/interact/calc_design_matrix.c */
void print_help_local2(char *, int, int);
void calc_design_matrix(struct med *, struct flt *, int, int);
void print_design_matrix(struct med *, struct flt *, int, int, FILE *);
/* src/interact/check_feedback.c */
/* src/interact/check_interaction.c */
unsigned char incrementally_check_interaction_coefficients(struct flt *, int, struct med *, unsigned char, int *, double);
unsigned char check_coulomb_stress_feedback(int, int, struct flt *, struct med *, unsigned char, int, unsigned char, int *, double);
/* src/interact/fracture_criterion.c */
unsigned char fracture_criterion(int, struct flt *, double *, unsigned char *, struct med *);
void two_dir_slip_check(unsigned char *, double *, unsigned char *, int, struct flt *, struct med *);
void deactivate_patch(int, struct flt *, struct med *);
void deactivate_group(int, int, struct flt *, struct med *);
int calc_absolute_shear_stress_and_inc(double *, double *, int, struct flt *);
int calc_absolute_shear_stress(double *, int, struct flt *);
/* src/interact/help_and_comments.c */
void phelp(void);
char *comment_on_code(short int);
char *comment_on_code_bc(short int, double);
/* src/interact/input.c */
int read_moment_file(float **, float **, float *, float *, int *, unsigned char);
int read_patch_event_file(float *, int *, int *, float *, float *, FILE *, struct med *);
int write_patch_event_file(float, int, int, float, float *, FILE *);
void read_rsf(char *, struct med *, struct flt *);
/* src/interact/optimize.c */
void optimize(struct flt *, struct med *);
double distsq(struct flt *, struct flt *);
double max_dist(struct flt *, struct med *);
double penalty_dist(struct flt *, struct med *, double);
/* src/interact/output.c */
int mysystem(const char *);
void print_slip_line(struct med *, struct flt *);
void flush_slipline(struct med *, struct flt *);
void print_fault_geometry_and_normals(struct flt *, int, char *);
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
/* src/interact/quake.c */
void quake(unsigned char *, double *, int, struct flt *, struct med *, unsigned char, unsigned char);
void add_quake_stress(int, unsigned char *, double *, struct flt *, struct med *);
/* src/interact/read_bin_events.c */
/* src/interact/read_boundary_conditions.c */
void read_boundary_conditions(struct med *, struct flt *, double *, double *, unsigned char);
void read_one_step_bc(FILE *, struct med *, struct flt *, double *, double *, unsigned char);
unsigned char slip_type_bc(int);
unsigned char bc_is_directionally_unconstrained(unsigned char);
unsigned char bc_contains_dip_motion(unsigned char);
/* src/interact/restart.c */
void adjust_medium_for_restart(struct med *, struct flt *);
/* src/interact/rhs.c */
void init_equation_system(struct med *, struct flt *);
void add_to_active_fault_list(int, int **, int *, unsigned char **);
void add_to_right_hand_side(double, double **, double **, int *);
/* src/interact/rupture.c */
unsigned char activate_faults(struct flt *, struct med *);
void fault_criterion(int, struct flt *, struct med *);
/* src/interact/solve_mode_dependent.c */
void assemble_a_matrix_4(double *, int, unsigned char *, int, int *, struct flt *, struct med *);
void add_quake_stress_4(unsigned char *, double *, int, struct flt *, struct med *);
unsigned char check_coulomb_stress_feedback_4(int, int, struct flt *, struct med *, unsigned char, unsigned char, int *, double);
/* src/la_and_geo/divide_fault_in_patches.c */
void create_patches(int, struct flt *, struct flt **, int *, int *, unsigned char, unsigned char, double *, double, double, long *);
void determine_segments(int *, int *, struct flt *, unsigned char, double *);
void divide_fault_in_patches(int, struct flt *, struct flt **, int *, int *, unsigned char, unsigned char, double, double, long *, unsigned char);
void get_flt_location(struct flt *, double *, double *, double *, int, int);
void randomize_strike_dip(double, double, struct flt *, long *);
/* src/la_and_geo/eigensystem.c */
void eigensystem3d(double [3][3], double *, double *);
void calc_eigensystem_sym3d(double *, double *, double *, unsigned char);
void calc_eigensystem_sym2d(double *, double *, double *, unsigned char);
void eispack_driver(double *, int, int, double *, double *, double *, double *, int);
/* src/la_and_geo/far_enough.c */
unsigned char far_enough(struct flt *, struct flt *, double);
/* src/la_and_geo/fit_plane.c */
void fit_plane(int, double *, double *, double *, double *, double *, double *, double *, double *, float *, float *, double *, unsigned char, unsigned char);
void points2patch(struct flt *, double *, unsigned char);
/* src/la_and_geo/geometry.c */
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
/* src/la_and_geo/kdtree.c */
void kdtr_node_init(kd_node);
void KDTreeCreate(int, KDTree *);
void KDTreeDestroy(KDTree *);
void KDTreeReset(KDTree);
void KDTreeView(KDTree);
void KDTreeSetPoints(KDTree, int);
void KDTreeGetPoints(KDTree, int *, kd_node *);
void KDTreeInsertPoint(KDTree, double []);
void KDTreeSetup(KDTree);
void KDTreeFindNearest(KDTree, double [], kd_node *, double *);
/* src/la_and_geo/levmarq_numrec.c */
void mrqmin(double *, double *, int, double *, int *, int, double **, double **, double *, double *, struct bmd *, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, double, double *, unsigned char, struct prj);
void print_lm_progress(double, double, double, double, int, int, double *, double, int, char **, int, int, unsigned char, unsigned char);
void covsrt(double **, int, int *, int);
void gaussj(double **, int, double **, int);
void mrqcof(double *, double *, int, double *, int *, int, double **, double *, double *, struct bmd *, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, double, double *, unsigned char, struct prj);
void nrerror(char *);
double **matrix(long, long, long, long);
void free_matrix(double **, long, long, long, long);
/* src/la_and_geo/llgeo.c */
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
double mwfromm0(double);
/* src/la_and_geo/lusolve.c */
void lu_driver(double *, double *, double *, int, int, struct med *);
/* src/la_and_geo/matrixio.c */
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
/* src/la_and_geo/myprojectsimple.c */
/* src/la_and_geo/nnls.c */
void nnls_driver_i(double *, double *, double *, int, int);
/* src/la_and_geo/petsc_simple_solve.c */
/* src/la_and_geo/project_stress.c */
/* src/la_and_geo/segseg.c */
void intersect(double *, double *, double *, double *, double *, int *);
char SegSegInt(tPointi, tPointi, tPointi, tPointi, tPointd);
char ParallelInt(tPointi, tPointi, tPointi, tPointi, tPointd);
void Assigndi(tPointd, tPointi);
unsigned char Between(tPointi, tPointi, tPointi);
int Collinear(tPointi, tPointi, tPointi);
int AreaSign(tPointi, tPointi, tPointi);
unsigned char crosses(double *, double *);
/* src/la_and_geo/solve.c */
int solve(struct med *, struct flt *);
void add_solution(int, unsigned char *, double *, int *, struct med *, struct flt *, unsigned char, unsigned char, double);
void assemble_a_matrix(double *, int, unsigned char *, int, int *, struct flt *, struct med *);
int par_assemble_a_matrix(int, unsigned char *, int, int *, struct flt *, struct med *);
/* src/la_and_geo/sparse.c */
size_t create_crs_sparse_from_memory(int, double *, double, unsigned int **, unsigned int **, double **);
size_t create_crs_sparse_from_file(int, double, unsigned int **, unsigned int **, double **, FILE *);
size_t create_ccs_sparse_from_memory(int, double *, double, unsigned int **, unsigned int **, double **);
size_t create_ccs_sparse_from_file(int, double, unsigned int **, unsigned int **, double **, FILE *);
/* src/la_and_geo/sparse_nr.c */
size_t create_nrs_sparse(struct med *, double, double *, unsigned char);
size_t create_nrs_sparse_from_memory(int, double *, double, unsigned int **, double **);
size_t create_nrs_sparse_from_file(int, double, unsigned int **, double **, FILE *);
double get_nrs_sparse_el(int, int, unsigned int *, double *);
/* src/la_and_geo/sparse_solve.c */
void sparse_driver(unsigned int *, unsigned int *, double *, double *, double *, int, struct med *);
/* src/la_and_geo/svd.c */
void svd_driver_lapack(double *, double *, double *, int, int, double *, int, double *, double **, int, unsigned char, unsigned char);
void svd_driver_numrec(double *, double *, double *, int, int, double *, int, double *, unsigned char, double **, double **, int, unsigned char, int, unsigned char, unsigned char, unsigned char);
void assemble_cov_svd(double *, int, double *, double *);
void print_singular_values(double *, int, FILE *);
void indexx(int, double *, int *);
void reduce_a_matrix(double **, int, int, int);
/* src/la_and_geo/trigonometry.c */
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
/* src/la_and_geo/tritri.c */
int coplanar_tri_tri(float [3], float [3], float [3], float [3], float [3], float [3], float [3]);
int tri_tri_intersect(float [3], float [3], float [3], float [3], float [3], float [3]);
/* src/block/block_checkflt.c */
/* src/block/block_compute_vel_from_omega.c */
/* src/block/block_eval_geokada.c */
void block_eval_geookada(double *, double *, double *, double [3][3], double, double, double, double, double, double, double, double, double, double, double, int *, unsigned char, unsigned char);
void block_save_solution_and_faults(double *, int, int, struct bflt *, double *, struct prj *, FILE *, unsigned char, unsigned char);
void block_load_solution_and_faults(double **, int *, int *, struct bflt **, double **, struct prj **, FILE *, unsigned char *, unsigned char *);
void block_eval_blockvec(double *, double, int, double *, struct prj *, double *);
/* src/block/block_evaluate_solution.c */
/* src/block/block_levmarq.c */
void run_lm(struct bmd *, long int *, struct prj, double *, unsigned char, unsigned char, double, double *, unsigned char, unsigned char, double, unsigned char, unsigned char, unsigned char, int, int, char **);
/* src/block/block_matrix.c */
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
/* src/block/block_output.c */
void block_output(struct bmd *, unsigned char, char **, struct prj, unsigned char, unsigned char, double, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char);
void sort_eigen(int *, double *);
void print_horizontal_stress(double *, FILE *);
void print_projected_stress(double *, FILE *);
void calculate_slipsol_sigma(struct bmd *, double *, double *, int);
void calc_geo_euler_pole(double *, double *, double *, double *);
void print_simple_vel(double *, double *, double *, double *, int, char *);
/* src/block/block_read_bflt.c */
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
/* src/block/block_read_euler.c */
void read_constrained_euler_poles(struct bmd *, char **, int);
/* src/block/block_read_gps.c */
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
/* src/block/block_solve.c */
void solve_block(double *, double *, double *, int, int, int, int, unsigned char, int, int, int, int, double *, double, double *, double **, double **, double *, double *, unsigned char, int, unsigned char, int, unsigned char);
void evaluate_block_solution(double *, int, int, double *, double *, double *, struct bmd *);
/* src/block/block_stress.c */
void block_scale_stresses(double *, double *, double *, int, int, int, double, double *, unsigned char);
void calc_horizontal_stress_vec(double *, double *, int);
void calc_horizontal_stress(double *, double *, double *, double *);
void calc_dir_diff_vec(double *, double *, double *, int, unsigned char);
double calc_dir_diff(double, double, unsigned char);
void cart_mat_from_horsym(double, double, double, double *);
void rescale_observed_stresses(double *, double *, double *, double, double *, unsigned char, struct bmd *, unsigned char, unsigned char);
/* src/block/blockinvert.c */
/* src/block/geo_okada.c */
/* src/green/eval_2dsegment.c */
void eval_2dsegment_plane_strain(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
void eval_2dsegment_plane_stress(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
void eval_2dsegment_plane_strain_basic(double *, struct flt *, double *, double *, double [3][3], int *);
void eval_2dsegment_plane_stress_basic(double *, struct flt *, double *, double *, double [3][3], int *);
void get_2dseg_geo(double *, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *);
void get_2dseg_ffac(double *, double *, double *, double *, double *, double *, double, double, double, double, double, double, double, double, double, double, double *);
void get_2dseg_disp(double *, double *, double *, double, double, double, double, double, double);
void get_2dseg_stress(double [3][3], double *, double *, double, double, double, double);
void eval_2dsegment_plane_strain_tdd(double *, struct flt *, double *, double *, double [3][3], int, int *, unsigned char);
/* src/green/eval_iquad.c */
void eval_iquad(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
/* src/green/eval_okada.c */
void eval_okada(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char, unsigned char);
void eval_okada_basic(double *, double, double, double, double, double *, double *, double [3][3], int *, unsigned char);
void eval_point(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char, unsigned char);
void eval_point_short(double *, double *, double, double, double, double, double *, double *, double [3][3], int *, unsigned char, unsigned char);
void set_stress_and_disp_nan(double [3][3], double *, unsigned char);
/* src/green/eval_triangle_gauss.c */
void eval_triangle_gauss(double *, struct flt *, double *, double *, double [3][3], int *);
void get_gauss_points(double *, double *, double *, int);
/* src/green/eval_triangle_nw.c */
void eval_triangle_nw(double *, struct flt *, double *, double *, double [3][3], int *, unsigned char);
void get_tri_prop_based_on_gh(struct flt *);
/* src/green/eval_triangle_tgf.c */
/* src/util/calc_cart_from_eigen_stress.c */
unsigned char read_vecs(int, double *, double *, double *, double *, double *);
void ccfes_help(char **);
/* src/util/calc_eigen_from_cart_stress.c */
void fehelp(char **);
/* src/util/calc_spatial_correlation.c */
void calc_spatial_correlation(struct flt *, int, int, int, int, double **, double **, int **, double, double, double *, float **, float **, float **);
double correlation_coefficient(double *, double *, int);
/* src/util/calc_stress_stat.c */
/* src/util/compare_fault.c */
int compare_fault_length(const void *, const void *);
int compare_fault_width(const void *, const void *);
/* src/util/create_random_mu_file.c */
/* src/util/create_random_stress_file.c */
void get_random_stress(double *, double *, double *, int, long *, int);
/* src/util/fit_mean_stress.c */
void eval_stress_correction(double *, double *, double *, double *, int, double *, double *, double *, double *, double *, double *);
/* src/util/fit_simple_stress_from_cart.c */
/* src/util/fltcopy.c */
void fltswap(struct flt *, struct flt *);
void fltcp(struct flt *, struct flt *);
/* src/util/fstress2eig.c */
/* src/util/fstress2hor.c */
/* src/util/generate_random_2d.c */
/* src/util/generate_slipdia.c */
/* src/util/geoproject.c */
/* src/util/makefault.c */
/* src/util/mspectral.c */
void find_range(int *, int *, float, float, float, float, float *, float *, float *, int);
/* src/util/myopen.c */
FILE *myopen(char *, char *);
/* src/util/mysincos.c */
void my_sincos_deg(double *, double *, double);
void my_sincos_deg_ftnd_(double *, double *, double *);
void my_sincos_ftn_(double *, double *, double *);
void my_sincos(double *, double *, double);
void my_sincosd(double *, double *, double);
void my_sincos_degd(double *, double *, double);
/* src/util/patch2area.c */
/* src/util/patch2bc.c */
/* src/util/patch2dis3d.c */
void get_dis3d_parameters(double, double, double, double, double, double, double *, double *, double *, double *);
/* src/util/patch2geom.c */
/* src/util/patch2group.c */
/* src/util/patch2poly3d.c */
/* src/util/patch2vertices.c */
/* src/util/patch2vtk.c */
/* src/util/patch2xyz.c */
/* src/util/patch2xyzvec.c */
/* src/util/patchquad2patchtri.c */
/* src/util/period.c */
void period(float *, float *, int, float, float, float *, float *, int, int *, int *, float *);
void avevar(float *, int, float *, float *);
void fasper(float *, float *, int, float, float, float *, float *, int, int *, int *, float *);
void spread(float, float *, int, float, int);
void realft(float *, int, int);
void four1(float *, int, int);
/* src/util/plotevents.c */
/* src/util/plotting.c */
/* src/util/points2patch.c */
int read_points_local(double *, int *, unsigned char, FILE *, int);
/* src/util/randgen.c */
void assign_random(double *, double, double, double, long *, int);
double myrand(long *);
double myrandnr(double, long *);
int myrandi(int, long *);
double mygauss_randnr(double, long *);
double mypower_randnr(double, double, double, int, long *);
double gasdev(long *);
double ran1(long *);
double ran2(long *);
/* src/util/randomflt.c */
void check_input_parameters(int, char **, int *, long *, int *, double *, unsigned char *, unsigned char *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, unsigned char *, int *, double *, unsigned char *, double *, double *);
/* src/util/randomize_list.c */
void randomize_list(int **, int, unsigned char);
int slist_sort(const void *, const void *);
/* src/util/randomize_strike.c */
/* src/util/simul_gr.c */
/* src/util/sort_events.c */
/* src/util/stress_aux.c */
double stress_misfit(double [3][3], double [3][3]);
void stress_vec_from_hstate(double, double, double, double, int, double *, double *);
void expand_stress_matrix6to9(double *);
void convert_6sym_to_9_matrix(double *, double [3][3]);
double resolved_stress(double *, double [3][3], double *);
void calc_three_stress_components(double [3][3], double *, double *, double *, double *, double *, double *, double *);
/* src/util/string_compare.c */
unsigned char strings_match(char *, char *);
/* src/util/terminate.c */
int terminate(struct med *, struct flt *);
/* src/util/tri2patch.c */
