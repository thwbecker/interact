//
// numerical recipes SVD related routines
//

extern void pythag(COMP_PRECISION *,COMP_PRECISION *);
extern void svbksb(COMP_PRECISION *, COMP_PRECISION *,
		   COMP_PRECISION *,int *,int*,int*,int*,
		   COMP_PRECISION *,  COMP_PRECISION *);
extern void svdcmp(COMP_PRECISION *,int *,int*,int*,int *,
		   COMP_PRECISION *,COMP_PRECISION *);
extern void cginv(COMP_PRECISION *,COMP_PRECISION *,
		  COMP_PRECISION *,COMP_PRECISION *,
		  COMP_PRECISION *,int *,int *);

#ifdef SGI_SUBROUTINE_CONVENTION

#define svdcmp svdcmp_
#define svbksb svbksb_
#define cginv cginv_
#define pythag pythag_


#endif

