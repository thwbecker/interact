#ifdef SINGLE_PREC
#define COMP_PRECISION float
#define DATA_FORMAT "%f"
#else
#define DATA_FORMAT "%lf"
#define COMP_PRECISION double
#endif
void dwnnls(double *, int *, int *, int *, 
	    int *, int *, double *, 
	    double *, double *, int *,int *, double *);

void wnnls(float *, int *, int *, int *, 
	    int *, int *, float *, 
	    float *, float *, int *,int *, float *);

void dwnnls_(double *, int *, int *, int *, 
	    int *, int *, double *, 
	    double *, double *, int *,int *, double *);

void wnnls_(float *, int *, int *, int *, 
	    int *, int *, float *, 
	    float *, float *, int *,int *, float *);

void nnls_driver(int ,int , int ,COMP_PRECISION *,
		 COMP_PRECISION *, COMP_PRECISION *,
		 COMP_PRECISION *,int *);
void nnls_driver_ftn(int *,int *, int *,
		     COMP_PRECISION *,
		     COMP_PRECISION *, 
		     COMP_PRECISION *,
		     COMP_PRECISION *,
		     int *);
void nnls_driver_ftn__(int *,int *, int *,
		       COMP_PRECISION *,
		       COMP_PRECISION *, 
		       COMP_PRECISION *,
		       COMP_PRECISION *,
		       int *);
