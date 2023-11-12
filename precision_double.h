#define COMP_PRECISION double 
#define C_PREC real*8

/* for the stress and displacement summation, adding up fault
   contributions? */
/* #define SUM_ARR_PREC_IN_DOUBLE */


/* for the A matrix */
#define A_MATRIX_PREC_IN_DOUBLE

#define I_MATRIX_PREC double



#define EPS_COMP_PREC 7.0e-15
#define EPS_AMAT_PREC 7.0e-15

#define ONE_CP_FORMAT "%lf"
#define ONE_IP_FORMAT "%lf"
#define TWO_CP_FORMAT "%lf %lf"
#define TWO_IP_FORMAT "%lf %lf"
#define THREE_CP_FORMAT "%lf %lf %lf"
#define THREE_CPI_FORMAT "%lf %lf %lf %i"
#define SIX_CP_FORMAT "%lf %lf %lf %lf %lf %lf"
#define NINE_CP_FORMAT "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define TWELVE_CP_FORMAT "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define FIFTEEN_CP_FORMAT "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define PATCH_CP_FORMAT "%lf %lf %lf %f %f %lf %lf %i"
#define RS_CP_FORMAT "%lf %lf %lf %f %lf %lf %i"
#define FIELD_CP_FORMAT "%lf %lf %i %lf %lf %i %lf %lf %i"
#define IIF_CP_FORMAT "%i %i %lf"
#define IF_CP_FORMAT "%i %lf"
#define FLTDAT_FORMAT "%f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf %i %i"
#define FLTDAT_BLK_FORMAT "%*f %*f %*f %*f %*f %lf %lf %lf %*f %*f %*f %*i %*i"
#define EISPACK_RS rs_
