/* 

   programming conventions and macros

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

#define log2(x)         (log(x) * INVERSE_LN_TWO)


#define TOGV(x) (((x==FALSE)?(TRUE):(FALSE)))
#define MEMERROR(x) {fprintf(stderr,"%s: memory allocation error\nexiting now.\n",x);exit(-1);}
#define PMEMERROR(x) {if(medium->comm_rank==0)fprintf(stderr,"%s: memory allocation error\nexiting now.\n",x);exit(-1);}
#define READ_ERROR(x) {fprintf(stderr,"read error file \"%s\", exiting\n",x);exit(-1);}
#define PERROR(x) {if(medium->comm_rank==0)fprintf(stderr,"ERROR: %s\n",x);exit(-1);}
#define HEADNODE if(medium->comm_rank==0)


#define RAD2DEGF(x) ((x)*RAD2DEG)
#define DEG2RADF(x) ((x)*DEG2RAD)

// numerical recipes
#define NUMREC_SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define NUMREC_ISWAP(a,b) {itemp=(a);(a)=(b);(b)=itemp;}


