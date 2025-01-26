//
// header for randomflt.c
//



// max number of iterations

#define MAXITER 1000000
#define MAX_SEIS_DEPTH 15e3	/* maximum depth [m] for -gra option */
void check_input_parameters(int , char **, int *,long *, int *,
			    COMP_PRECISION *, my_boolean *,
			    my_boolean *,
			    COMP_PRECISION *, 
			    COMP_PRECISION *, COMP_PRECISION *,
			    COMP_PRECISION *, COMP_PRECISION *,
			    COMP_PRECISION *,COMP_PRECISION *,
			    COMP_PRECISION *, COMP_PRECISION *,
			    COMP_PRECISION *,COMP_PRECISION *,
			    COMP_PRECISION *, COMP_PRECISION *,
			    COMP_PRECISION *, COMP_PRECISION *,
			    int *,my_boolean *,
			    int *,
			    COMP_PRECISION *,my_boolean *,COMP_PRECISION *, COMP_PRECISION *);
