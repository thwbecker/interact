/* 

   programming conventions and macros
   $Id: macros.h,v 1.22 2003/07/21 22:37:02 becker Exp $

*/
#define TRUE 1
#define FALSE 0

#define SQUARE(a) (((double)a) * ((double)a))
#ifndef MIN
 #define MIN(x, y) (((x)<(y))?(x):(y))
#endif
#ifndef MAX
 #define MAX(x, y) (((x)>(y))?(x):(y))
#endif

#define INVERSE_LN_TWO  1.442695040888963407359924681001892137426645954152985934135449407
#define log2(x)         (log(x) * INVERSE_LN_TWO)


#define TOGV(x) (((x==FALSE)?(TRUE):(FALSE)))
#define MEMERROR(x) {fprintf(stderr,"%s: memory allocation error\nexiting now.\n",x);exit(-1);}
#define READ_ERROR(x) {fprintf(stderr,"read error file \"%s\", exiting\n",x);exit(-1);}
/* 
   real constants 
*/
#define STRLEN 200
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define PIHALF 1.5707963267948966192313216916397514420985846996875529104874722961
#define TWOPI 6.2831853071795864769252867665590057683943387987502116419498891846
#define DEG2RAD 0.0174532925199432957692369076848861271344287188854172545609719144
#define RAD2DEG 57.295779513082320876798154814105170332405472466564321549160243861
#define PIOVERONEEIGHTY DEG2RAD
#define ONEEIGHTYOVERPI  RAD2DEG

//
// the medium structure has a NaN field which will be 
// initialized as sqrt(-1). the nan macro is only used if 
// medium is not passed to subroutines
//
// WARNING: THIS (WHICH SHOULD BE BETTER)
// DID NOT WORK ON IRIX WITH 
// -Ofast: #define NaN (1./0.-1./0.)
//
#define NaN (1.0/0.0)


#define ONE_MEGABYTE 1048576.0


// derived functions

#define RAD2DEGF(x) ((x)*RAD2DEG)
#define DEG2RADF(x) ((x)*DEG2RAD)

// numerical recipes
#define NUMREC_SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define NUMREC_ISWAP(a,b) {itemp=(a);(a)=(b);(b)=itemp;}


#define RADIUS_EARTH 6371.0 // radius of earth in km
