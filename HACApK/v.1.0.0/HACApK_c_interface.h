//
// Declare Fortran functions as external C functions
//
extern void init_hacapk_struct(void **obj_ptr, int *n);
extern void deallocate_hacapk_struct(void **obj_ptr);
extern void set_hacapk_struct_coord(void **obj_ptr, double *x, double *y, double *z);
extern void make_hacapk_struct_hmat(void **obj_ptr, double *ztol);
extern void hacapk_mult_Ax(void **obj_ptr, double *x, double *b);
extern void hacapk_solve_Ab(void **obj_ptr, double *b, double *x, double *ztol);
extern void hacapk_assemble_dense_mat(void **obj_ptr, double *Ad, int *n);
