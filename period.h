//
// numerical recipes routines for period.c
//
#define SPECA_PREC float
#define SPPF "%f"
#define SPPF2 "%f %f"

void period(SPECA_PREC *,SPECA_PREC *,int ,SPECA_PREC ,SPECA_PREC ,SPECA_PREC *,
	    SPECA_PREC *,int ,int *,int *,SPECA_PREC *);
void avevar(SPECA_PREC *,int ,SPECA_PREC *,SPECA_PREC *);
void four1(SPECA_PREC *,int ,int );
void realft(SPECA_PREC *,int ,int );
void spread(SPECA_PREC ,SPECA_PREC *,int ,SPECA_PREC , int );
void fasper(SPECA_PREC *,SPECA_PREC *,int ,SPECA_PREC ,SPECA_PREC ,SPECA_PREC *,
	    SPECA_PREC *,int ,int *, int *,SPECA_PREC *);



