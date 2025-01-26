void nnls_driver(A_MATRIX_PREC *, A_MATRIX_PREC *,int );
#ifdef SGI_SUBROUTINE_CONVENTION
#define law_nnls law_nnls_
#endif


extern void law_nnls(A_MATRIX_PREC *,int *,int *, int *,
		     A_MATRIX_PREC *,A_MATRIX_PREC *,
		     A_MATRIX_PREC *,A_MATRIX_PREC *,
		     A_MATRIX_PREC *,int *, int *);
