#ifdef USE_DOUBLE_PRECISION
#include "precision_double.h"
#elif defined USE_MIXED_PRECISION /* no good */
#include "precision_mixed.h"
#else
#include "precision_single.h"
#endif


void myprojectsimple(COMP_PRECISION *, COMP_PRECISION *,COMP_PRECISION ,
		     COMP_PRECISION , COMP_PRECISION , int );

int	solve_right_spherical_triangle();
int	sphere_azim_dist();
void	oblique_transform (COMP_PRECISION , COMP_PRECISION , COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *);

void make_euler_matrix (double *, double *, double); /*  */
void	matrix_3v(double *, double *, double *);
void	matrix_2v(double *, double *, double *);
void sphere_project_setup (COMP_PRECISION, COMP_PRECISION , COMP_PRECISION *,
			   COMP_PRECISION , COMP_PRECISION , COMP_PRECISION *,
			   COMP_PRECISION , COMP_PRECISION *, COMP_PRECISION *, GMT_LONG); 


