/* src/calc_interaction_matrix.c */
void print_help_local(char *);
/* src/calc_stress.c */
void initialize_stress_state(struct flt *, struct med *, unsigned char, float *, float *);
void update_stress_state(struct flt *, struct med *);
void calc_fields(struct med *, struct flt *, unsigned char, unsigned char, float *, float *);
/* src/compress_interaction_matrix.c */
/* src/coulomb_stress.c */
float coulomb_stress(float, float, float, float);
float cstress_drop(float, float, float, float);
unsigned char in_coulomb_compress_regime(float);
/* src/eval_green.c */
void eval_green_and_project_stress_to_fault(struct flt *, int, int, float *, float *, unsigned char, unsigned char);
void eval_green(float *, struct flt *, float *, float *, float [3][3], int *, unsigned char, unsigned char, unsigned char);
void eval_triangle_general(float *, struct flt *, float *, float *, float [3][3], int *, unsigned char, unsigned char);
void eval_green_at_receiver(struct flt *, int, int, float *, float *, float [3][3], int *, unsigned char, unsigned char, unsigned char);
void eval_green_basic(float *, struct flt *, float *, float *, float [3][3], int *, unsigned char);
double ckernel_func(int, int, void *);
/* src/get_projected_fault_parameters.c */
void get_projected_fault_parameters(float [2][2], float, float *, float *, float *, float *, float *, float *);
/* src/init.c */
void check_parameters_and_init_interact(int, char **, struct med **, struct flt **, unsigned char *, float *, float *);
void initialize_interact(struct med **, struct flt **, unsigned char,unsigned char, int, unsigned char, unsigned char, float, float *, float *, unsigned char, unsigned char, unsigned char, float, unsigned char, float, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, short int, unsigned char, float, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, float, unsigned char, unsigned char, unsigned char, unsigned char);
void init_files_interact(struct med **, struct flt **);
void init_parameters_interact(char **, int, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, float *, int *, unsigned char *, unsigned char *, unsigned char *, float *, unsigned char *, float *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, short int *, unsigned char *, float *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, float *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *, int);
float resolve_stress_on_fault_using_ctx(float [3][3],  struct interact_ctx *, int );
void advance_argument(int *, int, char **);
char *name_boolean(unsigned char);
unsigned char toggle(unsigned char *);
void read_stress_fac(unsigned char, float *, float *, float, struct med *);
/* src/interact_main.c */
/* src/interact_matrix_assembly.c */
void calc_interaction_matrix(struct med *, struct flt *, unsigned char);
float interaction_coefficient(int, int, int, int, struct flt *, int *, unsigned char);
void get_right_slip(float *, int, float,struct flt *);
float ic_from_file(int, int, int, int, struct med *);
float aij_from_file(int, int, int, FILE *);
int select_i_coeff_calc_mode(struct med *);
size_t imatrix_size(struct med *);
/* src/print_patch_geometry.c */
int print_patch_geometry_and_bc(int, struct flt *, int, float, unsigned char, FILE *, unsigned char, float *);
/* src/read_fltdat.c */
int read_fltdat(char *, struct flt *, struct med *, unsigned char);
/* src/read_geometry.c */
void read_geometry(char *, struct med **, struct flt **, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char);
/* src/read_stress_observations.c */
void read_stress_observations(struct bmd *, float *, float, unsigned char, float **, float, unsigned char);
/* src/interact/adjust_time_step.c */
void adjust_time_step(struct flt *, struct med *);
/* src/interact/calc_design_matrix.c */
void print_help_local2(char *, int, int);
void calc_design_matrix(struct med *, struct flt *, int, int);
void print_design_matrix(struct med *, struct flt *, int, int, FILE *);
/* src/interact/check_feedback.c */
/* src/interact/check_interaction.c */
unsigned char incrementally_check_interaction_coefficients(struct flt *, int, struct med *, unsigned char, int *, float);
unsigned char check_coulomb_stress_feedback(int, int, struct flt *, struct med *, unsigned char, int, unsigned char, int *, float);
/* src/interact/fracture_criterion.c */
unsigned char fracture_criterion(int, struct flt *, float *, unsigned char *, struct med *);
void two_dir_slip_check(unsigned char *, float *, unsigned char *, int, struct flt *, struct med *);
void deactivate_patch(int, struct flt *, struct med *);
void deactivate_group(int, int, struct flt *, struct med *);
int calc_absolute_shear_stress_and_inc(float *, float *, int, struct flt *);
int calc_absolute_shear_stress(float *, int, struct flt *);
/* src/interact/help_and_comments.c */
void phelp(void);
char *comment_on_code(short int);
char *comment_on_code_bc(short int, float);
/* src/interact/input.c */
int read_moment_file(float **, float **, float *, float *, int *, unsigned char);
int read_patch_event_file(float *, int *, int *, float *, float *, FILE *, struct med *);
int write_patch_event_file(float, int, int, float, float *, FILE *);
void read_rsf(char *, struct med *, struct flt *);
/* src/interact/optimize.c */
void optimize(struct flt *, struct med *);
float distsq(struct flt *, struct flt *);
float max_dist(struct flt *, struct med *);
float penalty_dist(struct flt *, struct med *, float);
/* src/interact/output.c */
int mysystem(const char *);
void print_slip_line(struct med *, struct flt *);
void flush_slipline(struct med *, struct flt *);
void print_fault_geometry_and_normals(struct flt *, int, char *);
void print_fault_stress_and_slip(struct med *, struct flt *, unsigned char);
void print_fault_data(char *, struct med *, struct flt *);
void print_fault_stress_stat(FILE *, int, struct med *, struct flt *);
void print_group_data_geom(char *, struct med *, struct flt *, int, int, float);
float select_val_for_print(struct flt *, int);
void print_stress(struct med *, struct flt *);
void print_stress_on_fault(struct med *, struct flt *, int);
void print_displacement(struct med *, struct flt *);
void print_equations(int, unsigned char *, int *, float *, int, char *, struct flt *);
void print_solutions(int, int *, struct flt *, struct med *, char *);
void flush_moment_stack(struct med *);
void fiddle_with_limits_for_plot(struct med *, int *, unsigned char *, float *, unsigned char);
void time_report(char *, char *, struct med *);
/* src/interact/quake.c */
void quake(unsigned char *, float *, int, struct flt *, struct med *, unsigned char, unsigned char);
void add_quake_stress(int, unsigned char *, float *, struct flt *, struct med *);
/* src/interact/read_bin_events.c */
/* src/interact/read_boundary_conditions.c */
void read_boundary_conditions(struct med *, struct flt *, float *, float *, unsigned char);
void read_one_step_bc(FILE *, struct med *, struct flt *, float *, float *, unsigned char);
unsigned char slip_type_bc(int);
unsigned char bc_is_directionally_unconstrained(unsigned char);
unsigned char bc_contains_dip_motion(unsigned char);
/* src/interact/restart.c */
void adjust_medium_for_restart(struct med *, struct flt *);
/* src/interact/rhs.c */
void init_equation_system(struct med *, struct flt *);
void add_to_active_fault_list(int, int **, int *, unsigned char **);
void add_to_right_hand_side(float, float **, float **, int *);
/* src/interact/rupture.c */
unsigned char activate_faults(struct flt *, struct med *);
void fault_criterion(int, struct flt *, struct med *);
/* src/interact/solve_mode_dependent.c */
void assemble_a_matrix_4(float *, int, unsigned char *, int, int *, struct flt *, struct med *);
void add_quake_stress_4(unsigned char *, float *, int, struct flt *, struct med *);
unsigned char check_coulomb_stress_feedback_4(int, int, struct flt *, struct med *, unsigned char, unsigned char, int *, float);
/* src/la_and_geo/divide_fault_in_patches.c */
void create_patches(int, struct flt *, struct flt **, int *, int *, unsigned char, unsigned char, float *, float, float, long *);
void determine_segments(int *, int *, struct flt *, unsigned char, float *);
void divide_fault_in_patches(int, struct flt *, struct flt **, int *, int *, unsigned char, unsigned char, float, float, long *, unsigned char);
void get_flt_location(struct flt *, float *, float *, float *, int, int);
void randomize_strike_dip(float, float, struct flt *, long *);
/* src/la_and_geo/eigensystem.c */
void eigensystem3d(float [3][3], float *, float *);
void calc_eigensystem_sym3d(float *, float *, float *, unsigned char);
void calc_eigensystem_sym2d(float *, float *, float *, unsigned char);
void eispack_driver(float *, int, int, float *, float *, float *, float *, int);
/* src/la_and_geo/far_enough.c */
unsigned char far_enough(struct flt *, struct flt *, float);
/* src/la_and_geo/fit_plane.c */
void fit_plane(int, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, unsigned char, unsigned char);
void points2patch(struct flt *, float *, unsigned char);
/* src/la_and_geo/geometry.c */
float patch_area(struct flt *);
void calc_lhemi_proj(float, float, float *);
void resolve_force(float *, float [3][3], float *);
void calc_quad_base_vecs(float *, float *, float *, float, float, float, float);
void get_maxsdir_stress_drops2(float *, float, float *);
void get_maxsdir_stress_drops(float *, float, float *, float, float *, float *);
void calculate_vertices(float *, struct flt *, float *, float *);
void calculate_bloated_vertices(float *, struct flt *, float);
int nvert_of_patch(struct flt *);
int ncon_of_subpatch(struct flt *, int);
float projected_slip_major_to_minor_patch(struct flt *, int, int, int);
int vtk_type_of_patch(struct flt *, int);
int number_of_subpatches(struct flt *);
int node_number_of_subelement(struct flt *, int, int);
void calculate_quad_vertices(float *, struct flt *, float);
void calculate_tri_vertices(float *, struct flt *, float);
void calculate_iquad_vertices(float *, struct flt *, float);
void calculate_point_source_vertices(float *, struct flt *, float, float *, float *);
void calculate_seg_vertices(float *, struct flt *, float);
float quad_area(float *);
float triangle_area(float *);
float triangle_area_gh(float *, float *);
void get_gh_tri_vec(float *, float *, float *);
void get_gh_quad_vec(float *, float *, float *, float *);
void check_fault_normal_vectors(struct flt *);
unsigned char check_planar(float *);
void calc_group_geometry(struct med *, struct flt *, struct geog *);
void vec_to_angles(float *, float *, float *);
void angles_to_vec(float, float, float *);
float vec_to_strike(float *);
float vec_to_dip(float *);
void check_fault_angles(struct flt *);
void check_angles(float *, double *);
void fix_azimuth(double *);
void globalx(float *, float, float, float *);
void calc_centroid_tri(float *, float *);
void calc_tri_bary_coord(float *, float *, float, float, float);
void calc_mean_quad_coord(float *, float *);
void calc_centroid_quad(float *, float *);
unsigned char patch_is_2d(unsigned char);
void calculate_position_of_patch(struct med *, struct flt *);
void compute_cartesian_slip(float *, float *, struct flt *);
void background_stress(float [3][3], float *, float, float *, float *, float);
void background_disp(float *, float *, struct med *, float *, float *);
void get_local_x_on_plane(float *, float *, float *, float *, float *);
void get_fault_plane_basevec(float *, float *, float *, struct flt *, struct med *);
void calc_deviatoric_stress(float [3][3], float [3][3], float *, float *);
void get_sub_normal_vectors(struct flt *, int, float *, float *, float *, float *);
unsigned char is_triangular(unsigned char);
void calc_global_strike_dip_from_local(struct flt *, float *, float *, float *);
void calc_global_slip_and_traction_from_local(struct flt *, float *, float *, float *, float *, float *, float *, float *, unsigned char);
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
void mrqmin(float *, float *, int, float *, int *, int, float **, float **, float *, float *, struct bmd *, unsigned char, float, unsigned char, unsigned char, unsigned char, unsigned char, float, float *, unsigned char, struct prj);
void print_lm_progress(float, float, float, float, int, int, float *, float, int, char **, int, int, unsigned char, unsigned char);
void covsrt(float **, int, int *, int);
void gaussj(float **, int, float **, int);
void mrqcof(float *, float *, int, float *, int *, int, float **, float *, float *, struct bmd *, unsigned char, float, unsigned char, unsigned char, unsigned char, unsigned char, float, float *, unsigned char, struct prj);
void nrerror(char *);
float **matrix(long, long, long, long);
void free_matrix(float **, long, long, long, long);
/* src/la_and_geo/llgeo.c */
void rotate_vec(float *, float *, double, double);
void rotate_vec2d(float *, float *, double, double);
void rotate_mat(float [3][3], float [3][3], float [3][3]);
void rotate_mat_z(float [3][3], float [3][3], double, double);
float tensor3d_norm(float [3][3]);
void zero_3x3_matrix(float [3][3]);
void scale_3x3_matrix(float [3][3], float);
void add_ay_to_3x3_matrix(float [3][3], float [3][3], float);
void add_y_to_3x3_matrix(float [3][3], float [3][3]);
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
float mwfromm0(float);
/* src/la_and_geo/lusolve.c */
void lu_driver(float *, float *, float *, int, int, struct med *);
/* src/la_and_geo/matrixio.c */
void print_matrix_ftrn(float *, int, int, FILE *, unsigned char);
void print_matrix_C(float *, int, int, FILE *, unsigned char);
void print_matrix_scaled_ftrn(float *, int, int, FILE *, unsigned char, float);
void print_matrix_ftrn_file(float *, int, int, char *, unsigned char);
void print_a_matrix(float *, int, int, FILE *, float *, unsigned char);
void read_a_matrix_from_file(float *, int, int, char *, char *);
void print_a_matrix_to_file(float *, int, int, char *, char *);
void print_system(float *, float *, float *, int, int, FILE *);
void print_sym3x3matrix(float [3][3], FILE *);
void print_interaction_matrix(struct med *, struct flt *, unsigned char);
void print_reduced_interaction_matrix(struct med *, struct flt *);
void print_vector(float *, int, FILE *);
void print_vector_file(float *, int, char *, char *);
void print_vector_row(float *, int, FILE *);
void print_b_vector(float *, int, FILE *, float *, unsigned char);
void print_3x3_matrix(float [3][3], FILE *);
/* src/la_and_geo/myprojectsimple.c */
/* src/la_and_geo/nnls.c */
void nnls_driver_i(float *, float *, float *, int, int);
/* src/la_and_geo/petsc_simple_solve.c */
/* src/la_and_geo/project_stress.c */
/* src/la_and_geo/segseg.c */
void intersect(float *, float *, float *, float *, float *, int *);
char SegSegInt(tPointi, tPointi, tPointi, tPointi, tPointd);
char ParallelInt(tPointi, tPointi, tPointi, tPointi, tPointd);
void Assigndi(tPointd, tPointi);
unsigned char Between(tPointi, tPointi, tPointi);
int Collinear(tPointi, tPointi, tPointi);
int AreaSign(tPointi, tPointi, tPointi);
unsigned char crosses(double *, double *);
/* src/la_and_geo/solve.c */
int solve(struct med *, struct flt *);
void add_solution(int, unsigned char *, float *, int *, struct med *, struct flt *, unsigned char, unsigned char, float);
void assemble_a_matrix(float *, int, unsigned char *, int, int *, struct flt *, struct med *);
int par_assemble_a_matrix(int, unsigned char *, int, int *, struct flt *, struct med *);
/* src/la_and_geo/sparse.c */
size_t create_crs_sparse_from_memory(int, float *, float, unsigned int **, unsigned int **, float **);
size_t create_crs_sparse_from_file(int, float, unsigned int **, unsigned int **, float **, FILE *);
size_t create_ccs_sparse_from_memory(int, float *, float, unsigned int **, unsigned int **, float **);
size_t create_ccs_sparse_from_file(int, float, unsigned int **, unsigned int **, float **, FILE *);
/* src/la_and_geo/sparse_nr.c */
size_t create_nrs_sparse(struct med *, float, float *, unsigned char);
size_t create_nrs_sparse_from_memory(int, float *, float, unsigned int **, float **);
size_t create_nrs_sparse_from_file(int, float, unsigned int **, float **, FILE *);
float get_nrs_sparse_el(int, int, unsigned int *, float *);
/* src/la_and_geo/sparse_solve.c */
void sparse_driver(unsigned int *, unsigned int *, float *, float *, float *, int, struct med *);
/* src/la_and_geo/svd.c */
void svd_driver_lapack(float *, float *, float *, int, int, float *, int, float *, float **, int, unsigned char, unsigned char);
void svd_driver_numrec(float *, float *, float *, int, int, float *, int, float *, unsigned char, float **, float **, int, unsigned char, int, unsigned char, unsigned char, unsigned char);
void assemble_cov_svd(float *, int, float *, float *);
void print_singular_values(float *, int, FILE *);
void indexx(int, float *, int *);
void reduce_a_matrix(float **, int, int, int);
/* src/la_and_geo/trigonometry.c */
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
/* src/la_and_geo/tritri.c */
int coplanar_tri_tri(float [3], float [3], float [3], float [3], float [3], float [3], float [3]);
int tri_tri_intersect(float [3], float [3], float [3], float [3], float [3], float [3]);
/* src/block/block_checkflt.c */
/* src/block/block_compute_vel_from_omega.c */
/* src/block/block_eval_geokada.c */
void block_eval_geookada(float *, float *, float *, float [3][3], float, float, float, float, float, float, float, float, float, float, float, int *, unsigned char, unsigned char);
void block_save_solution_and_faults(float *, int, int, struct bflt *, float *, struct prj *, FILE *, unsigned char, unsigned char);
void block_load_solution_and_faults(float **, int *, int *, struct bflt **, float **, struct prj **, FILE *, unsigned char *, unsigned char *);
void block_eval_blockvec(float *, float, int, float *, struct prj *, float *);
/* src/block/block_evaluate_solution.c */
/* src/block/block_levmarq.c */
void run_lm(struct bmd *, long int *, struct prj, float *, unsigned char, unsigned char, float, float *, unsigned char, unsigned char, float, unsigned char, unsigned char, unsigned char, int, int, char **);
/* src/block/block_matrix.c */
void assemble_block_fltdep_matrices(float **, float **, float **, float **, float *, float *, float *, int, int, int, int, int, int, int, float *, struct bflt *, struct bck *, unsigned char, unsigned char);
void assemble_block_a(float **, float *, int *, int, int, struct bck *, int, struct bmd *);
void assemble_block_d(float **, struct bflt *, int, int, int, struct bmd *);
void assemble_block_i(float **, struct bflt *, int, int, int);
void assemble_block_g(float **, struct bflt *, int, int);
void assemble_block_f(float **, struct bflt *, int, int, struct bck *, float *, unsigned char);
void assemble_block_k(float **, int, int, int, int, int, int, int, float *, float *, float *, unsigned char, unsigned char, float *, struct bflt *, int, int, int, int, float);
void assemble_stress_matrix(float *, int, struct bflt *, float *, int, int);
void block_assemble_fit_vector(float **, int, int, int, int, int, int, int, int, int, float *, int, float *, float *, unsigned char, float *, float **, float *, float **, int, struct bck *, int, struct bflt **, unsigned char, float **, unsigned char, unsigned char, unsigned char, float *, float *, float *, struct prj, unsigned char, int, int, float, struct bmd *);
void block_assemble_dyda_matrix(float **, float **, int, int, int, int, int, int, int, int, int, float *, float *, float **, float *, float **, float **, int, struct bck *, int, struct bflt **, unsigned char, unsigned char, unsigned char, float *, float *, float *, struct prj, unsigned char, int, int, float, struct bmd *);
float block_chi_square(float *, float *, float *, int, int, float, float *, float *);
float block_data_norm(float *, float *, int, int, float, struct bmd *);
int block_slip_direction(int, struct bflt *);
void assign_additional_sol_values(float *, int, int, int, struct bflt *, unsigned char, unsigned char, long int *, int);
void change_locking_depths(struct bflt **, int, int, int, int, float *, float *, float *, unsigned char, struct prj, float **, float **, float *, struct bmd *);
void calc_fault_sn_rms(float *, int, int, float *, float *);
void init_block_mods(struct bmd **);
void init_blocks(struct bmd *, int);
void convert_cart_sol(float *, float *, struct bmd *);
void check_solution_vector(float *, int, int, unsigned char, unsigned char);
void copy_block(struct bck *, struct bck *);
/* src/block/block_output.c */
void block_output(struct bmd *, unsigned char, char **, struct prj, unsigned char, unsigned char, float, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char);
void sort_eigen(int *, float *);
void print_horizontal_stress(float *, FILE *);
void print_projected_stress(float *, FILE *);
void calculate_slipsol_sigma(struct bmd *, float *, float *, int);
void calc_geo_euler_pole(float *, float *, float *, float *);
void print_simple_vel(float *, float *, float *, float *, int, char *);
/* src/block/block_read_bflt.c */
void read_bflt(struct bmd *, float, struct prj, unsigned char, float, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char);
void generate_new_fault(struct bflt **, int *, int *, float *, int, float *, int, float, struct prj, float *, unsigned char, float, unsigned char, unsigned char, float *, float, int, float *, unsigned char *, int *);
void get_bflt_intcoeff(struct bflt **, int, float *, int, float *, int, float *, unsigned char);
void free_bflt(struct bflt **, int, int);
void assign_bflt_dip_mode(struct bflt *, float);
void flip_block_code(struct bflt *, float *);
void print_fault_geometry(struct bflt *, int, FILE *);
void assign_fault_locking_depth_parameters(struct bflt *, float, struct prj, unsigned char, int);
int new_fault(struct bflt **, int);
void init_bflt(struct bflt *);
/* src/block/block_read_euler.c */
void read_constrained_euler_poles(struct bmd *, char **, int);
/* src/block/block_read_gps.c */
void read_gps_velocities(struct bmd *, struct prj *, float *, char **, unsigned char, int *, int *);
void project_gps_coordinates(float *, float *, int, int, int *, struct bck *, int, struct prj *, char **);
void block_read_centroids(float **, int *);
void find_spherical_rotation(struct bmd *, int, float *, unsigned char, unsigned char, unsigned char, unsigned char, unsigned char, float *, int);
void init_gps_pbase(struct bmd *);
void init_cart_gps_velsig(struct bmd *, unsigned char);
void calculate_c2p_gps(float *, float *, struct bmd *);
void convert_sig_p2c(float *, float *, float *);
void remove_cvel_from_vel(struct bmd *, float *, int);
char bname(int);
unsigned char node_in_geo_list(float *, int, float *);
/* src/block/block_solve.c */
void solve_block(float *, float *, float *, int, int, int, int, unsigned char, int, int, int, int, float *, float, float *, float **, float **, float *, float *, unsigned char, int, unsigned char, int, unsigned char);
void evaluate_block_solution(float *, int, int, float *, float *, float *, struct bmd *);
/* src/block/block_stress.c */
void block_scale_stresses(float *, float *, float *, int, int, int, float, float *, unsigned char);
void calc_horizontal_stress_vec(float *, float *, int);
void calc_horizontal_stress(float *, float *, float *, float *);
void calc_dir_diff_vec(float *, float *, float *, int, unsigned char);
float calc_dir_diff(float, float, unsigned char);
void cart_mat_from_horsym(float, float, float, float *);
void rescale_observed_stresses(float *, float *, float *, float, float *, unsigned char, struct bmd *, unsigned char, unsigned char);
/* src/block/blockinvert.c */
/* src/block/geo_okada.c */
/* src/green/eval_2dsegment.c */
void eval_2dsegment_plane_strain(float *, struct flt *, float *, float *, float [3][3], int *, unsigned char);
void eval_2dsegment_plane_stress(float *, struct flt *, float *, float *, float [3][3], int *, unsigned char);
void eval_2dsegment_plane_strain_basic(float *, struct flt *, float *, float *, float [3][3], int *);
void eval_2dsegment_plane_stress_basic(float *, struct flt *, float *, float *, float [3][3], int *);
void get_2dseg_geo(float *, float, float *, float *, float *, float *, float *, float *, float *, float *, float *, int *);
void get_2dseg_ffac(float *, float *, float *, float *, float *, float *, float, float, float, float, float, float, float, float, float, float, float *);
void get_2dseg_disp(float *, float *, float *, float, float, float, float, float, float);
void get_2dseg_stress(float [3][3], float *, float *, float, float, float, float);
void eval_2dsegment_plane_strain_tdd(float *, struct flt *, float *, float *, float [3][3], int, int *, unsigned char);
/* src/green/eval_iquad.c */
void eval_iquad(float *, struct flt *, float *, float *, float [3][3], int *, unsigned char);
/* src/green/eval_okada.c */
void eval_okada(float *, struct flt *, float *, float *, float [3][3], int *, unsigned char, unsigned char);
void eval_okada_basic(float *, float, float, float, float, float *, float *, float [3][3], int *, unsigned char);
void eval_point(float *, struct flt *, float *, float *, float [3][3], int *, unsigned char, unsigned char);
void eval_point_short(float *, float *, float, float, float, float, float *, float *, float [3][3], int *, unsigned char, unsigned char);
void set_stress_and_disp_nan(float [3][3], float *, unsigned char);
/* src/green/eval_triangle_gauss.c */
void eval_triangle_gauss(float *, struct flt *, float *, float *, float [3][3], int *);
void get_gauss_points(float *, float *, float *, int);
/* src/green/eval_triangle_nw.c */
void eval_triangle_nw(float *, struct flt *, float *, float *, float [3][3], int *, unsigned char);
void get_tri_prop_based_on_gh(struct flt *);
/* src/green/eval_triangle_tgf.c */
/* src/util/calc_cart_from_eigen_stress.c */
unsigned char read_vecs(int, float *, float *, float *, float *, float *);
void ccfes_help(char **);
/* src/util/calc_eigen_from_cart_stress.c */
void fehelp(char **);
/* src/util/calc_spatial_correlation.c */
void calc_spatial_correlation(struct flt *, int, int, int, int, float **, float **, int **, float, float, float *, float **, float **, float **);
float correlation_coefficient(float *, float *, int);
/* src/util/calc_stress_stat.c */
/* src/util/compare_fault.c */
int compare_fault_length(const void *, const void *);
int compare_fault_width(const void *, const void *);
/* src/util/create_random_mu_file.c */
/* src/util/create_random_stress_file.c */
void get_random_stress(float *, float *, float *, int, long *, int);
/* src/util/fit_mean_stress.c */
void eval_stress_correction(float *, float *, float *, float *, int, float *, float *, float *, float *, float *, float *);
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
void my_sincos_deg(float *, float *, float);
void my_sincos_deg_ftnd_(double *, double *, double *);
void my_sincos_ftn_(float *, float *, float *);
void my_sincos(float *, float *, float);
void my_sincosd(double *, double *, double);
void my_sincos_degd(double *, double *, double);
/* src/util/patch2area.c */
/* src/util/patch2bc.c */
/* src/util/patch2dis3d.c */
void get_dis3d_parameters(float, float, float, float, float, float, float *, float *, float *, float *);
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
int read_points_local(float *, int *, unsigned char, FILE *, int);
/* src/util/randgen.c */
void assign_random(float *, float, float, float, long *, int);
float myrand(long *);
float myrandnr(float, long *);
int myrandi(int, long *);
float mygauss_randnr(float, long *);
float mypower_randnr(float, float, float, int, long *);
double gasdev(long *);
double ran1(long *);
double ran2(long *);
/* src/util/randomflt.c */
void check_input_parameters(int, char **, int *, long *, int *, float *, unsigned char *, unsigned char *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, int *, unsigned char *, int *, float *, unsigned char *, float *, float *);
/* src/util/randomize_list.c */
void randomize_list(int **, int, unsigned char);
int slist_sort(const void *, const void *);
/* src/util/randomize_strike.c */
/* src/util/simul_gr.c */
/* src/util/sort_events.c */
/* src/util/stress_aux.c */
float stress_misfit(float [3][3], float [3][3]);
void stress_vec_from_hstate(float, float, float, float, int, float *, float *);
void expand_stress_matrix6to9(float *);
void convert_6sym_to_9_matrix(float *, float [3][3]);
float resolved_stress(float *, float [3][3], float *);
void calc_three_stress_components(float [3][3], float *, float *, float *, float *, float *, float *, float *);
/* src/util/string_compare.c */
unsigned char strings_match(char *, char *);
/* src/util/terminate.c */
int terminate(struct med *, struct flt *);
/* src/util/tri2patch.c */
