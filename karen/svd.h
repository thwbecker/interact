#ifdef SGI_SUBROUTINE_CONVENTION
#define svdcmp svdcmp_
#define svbksb svbksb_
#endif
void svd_driver(A_MATRIX_PREC *,A_MATRIX_PREC *,int );
void svbksb(A_MATRIX_PREC *, A_MATRIX_PREC *,
		   A_MATRIX_PREC *,int *,int*,int*,int*,
		   A_MATRIX_PREC *,  A_MATRIX_PREC *);
void svdcmp(A_MATRIX_PREC *,int *,int*,int*,int*,
		   A_MATRIX_PREC *,A_MATRIX_PREC *);
