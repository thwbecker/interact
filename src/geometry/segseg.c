/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 7.  It is not written to be comprehensible without the 
explanation in that book.

Compile:  gcc -o segseg segseg.c

Written by Joseph O'Rourke.
Last modified: November 1997
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1998 by Joseph O'Rourke.  It may be freely 
redistributed in its entirety provided that this copyright notice is 
not removed.
--------------------------------------------------------------------
*/
#include "interact.h"
#include <stdio.h>
#include <math.h>
#include <limits.h>
#define	SEGSEG_X 0
#define	SEGSEG_Y 1


/*-------------------------------------------------------------------*/

//
// segseg driver
//


void intersect(COMP_PRECISION *xa,COMP_PRECISION *xb,
	       COMP_PRECISION *xc,COMP_PRECISION *xd,
	       COMP_PRECISION *xp, int *code_answer)
     
{
  tPointi a,b,c,d;
  tPointd p;
  char code;
  char ans;
  static double intscale=((double)INT_MAX/100.);
  //
  // this is not elegant but should still be faster than
  // the stuff i programmed for floating point number
  //

  a[SEGSEG_X]=(int)(xa[X]*intscale);
  a[SEGSEG_Y]=(int)(xa[Y]*intscale);
  b[SEGSEG_X]=(int)(xb[X]*intscale);
  b[SEGSEG_Y]=(int)(xb[Y]*intscale);
  c[SEGSEG_X]=(int)(xc[X]*intscale);
  c[SEGSEG_Y]=(int)(xc[Y]*intscale);
  d[SEGSEG_X]=(int)(xd[X]*intscale);
  d[SEGSEG_Y]=(int)(xd[Y]*intscale);

  code = SegSegInt(a,b,c,d,p);
  if(code == 'e')
    *code_answer=2;
  else if(code == 'v')
    *code_answer=3;
  else if(code == '1')
    *code_answer=1;
  else 
    *code_answer=0;
  /* 
     printf("segseg: (%i,%i) (%i,%i) (%i,%i) (%i,%i) -> %i (%c)\n",
     a[SEGSEG_X],a[SEGSEG_Y],b[SEGSEG_X],b[SEGSEG_Y],
     c[SEGSEG_X],c[SEGSEG_Y],d[SEGSEG_X],d[SEGSEG_Y],
     *code_answer,code);
     */
}

/*---------------------------------------------------------------------
SegSegInt: Finds the point of intersection p between two closed
segments ab and cd.  Returns p and a char with the following meaning:
   'e': The segments collinearly overlap, sharing a point.
   'v': An endpoint (vertex) of one segment is on the other segment,
        but 'e' doesn't hold.
   '1': The segments inteorsect properly (i.e., they share a point and
        neither 'v' nor 'e' holds).
   '0': The segments do not intersect (i.e., they share no points).
Note that two collinear segments that share just one point, an endpoint
of each, returns 'e' rather than 'v' as one might expect.
---------------------------------------------------------------------*/
char	SegSegInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p )
{
   double  s, t;       /* The two parameters of the parametric eqns. */
   double num, denom;  /* Numerator and denoninator of equations. */
   char code = '?';    /* Return char characterizing intersection. */

   denom = a[SEGSEG_X] * (double)( d[SEGSEG_Y] - c[SEGSEG_Y] ) +
           b[SEGSEG_X] * (double)( c[SEGSEG_Y] - d[SEGSEG_Y] ) +
           d[SEGSEG_X] * (double)( b[SEGSEG_Y] - a[SEGSEG_Y] ) +
           c[SEGSEG_X] * (double)( a[SEGSEG_Y] - b[SEGSEG_Y] );

   /* If denom is zero, then segments are parallel: handle separately. */
   if (denom == 0.0)
      return  ParallelInt(a, b, c, d, p);

   num =    a[SEGSEG_X] * (double)( d[SEGSEG_Y] - c[SEGSEG_Y] ) +
            c[SEGSEG_X] * (double)( a[SEGSEG_Y] - d[SEGSEG_Y] ) +
            d[SEGSEG_X] * (double)( c[SEGSEG_Y] - a[SEGSEG_Y] );
   if ( (num == 0.0) || (num == denom) ) code = 'v';
   s = num / denom;
   //printf("num=%lf, denom=%lf, s=%lf\n", num, denom, s);

   num = -( a[SEGSEG_X] * (double)( c[SEGSEG_Y] - b[SEGSEG_Y] ) +
            b[SEGSEG_X] * (double)( a[SEGSEG_Y] - c[SEGSEG_Y] ) +
            c[SEGSEG_X] * (double)( b[SEGSEG_Y] - a[SEGSEG_Y] ) );
   if ( (num == 0.0) || (num == denom) ) code = 'v';
   t = num / denom;
   //   printf("num=%lf, denom=%lf, t=%lf\n", num, denom, t);

   if      ( (0.0 < s) && (s < 1.0) &&
             (0.0 < t) && (t < 1.0) )
     code = '1';
   else if ( (0.0 > s) || (s > 1.0) ||
             (0.0 > t) || (t > 1.0) )
     code = '0';

   p[SEGSEG_X] = a[SEGSEG_X] + s * ( b[SEGSEG_X] - a[SEGSEG_X] );
   p[SEGSEG_Y] = a[SEGSEG_Y] + s * ( b[SEGSEG_Y] - a[SEGSEG_Y] );

   return code;
}
char	ParallelInt( tPointi a, tPointi b, tPointi c, tPointi d, tPointd p )
{
   if ( !Collinear( a, b, c) )
      return '0';

   if ( Between( a, b, c ) ) {
      Assigndi( p, c );
      return 'e';
   }
   if ( Between( a, b, d ) ) {
      Assigndi( p, d );
      return 'e';
   }
   if ( Between( c, d, a ) ) {
      Assigndi( p, a );
      return 'e';
   }
   if ( Between( c, d, b ) ) {
      Assigndi( p, b );
      return 'e';
   }
   return '0';
}
void	Assigndi( tPointd p, tPointi a )
{
   int i;
   for ( i = 0; i < SEGSEG_DIM; i++ )
      p[i] = a[i];
}
/*---------------------------------------------------------------------
Returns TRUE iff point c lies on the closed segement ab.
Assumes it is already known that abc are collinear.
---------------------------------------------------------------------*/
my_boolean    Between( tPointi a, tPointi b, tPointi c )
{
   tPointi      ba, ca;

   /* If ab not vertical, check betweenness on x; else on y. */
   if ( a[SEGSEG_X] != b[SEGSEG_X] )
      return ((a[SEGSEG_X] <= c[SEGSEG_X]) && (c[SEGSEG_X] <= b[SEGSEG_X])) ||
             ((a[SEGSEG_X] >= c[SEGSEG_X]) && (c[SEGSEG_X] >= b[SEGSEG_X]));
   else
      return ((a[SEGSEG_Y] <= c[SEGSEG_Y]) && (c[SEGSEG_Y] <= b[SEGSEG_Y])) ||
             ((a[SEGSEG_Y] >= c[SEGSEG_Y]) && (c[SEGSEG_Y] >= b[SEGSEG_Y]));
}

int Collinear( tPointi a, tPointi b, tPointi c )
{
   return AreaSign( a, b, c ) == 0;
}
int     AreaSign( tPointi a, tPointi b, tPointi c )
{
    double area2;

    area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
            ( c[0] - a[0] ) * (double)( b[1] - a[1] );

    /* The area should be an integer. */
    if      ( area2 >  0.5 ) return  1;
    else if ( area2 < -0.5 ) return -1;
    else                     return  0;
}

//
// from here old crosses routine
// horrible...

#define CROSS_EPSILON 1.0e-10

my_boolean crosses(double *x, double *y)

{
  double a[4],b[2],c[4],d[2],fac,alpha,beta;
  int i,j,hit;
  my_boolean contact=FALSE;
  
  for(j=0;j<4;j++)a[j]= *(x+j);
  if((a[2] < a[0])||((a[0] == a[2])&&(a[3] < a[1])))
    {fac=a[0];a[0]=a[2];a[2]=fac;fac=a[1];a[1]=a[3];a[3]=fac;}; 
  b[0] = a[2] - a[0];b[1] = a[3] - a[1];
  if ((b[0] != 0)||(b[1] != 0))
    {
      for(j=0;j<4;j++)c[j]= *(y+j);
      if((c[2] < c[0])||((c[0] == c[2])&&(c[3] < c[1])))
	{fac=c[0];c[0]=c[2];c[2]=fac;fac=c[1];c[1]=c[3];c[3]=fac;}; 
      d[0] = c[2] - c[0];d[1] = c[3] - c[1];
      fac = dotp(b,d,2)/(norm(b,2)*norm(d,2));
      if ((fac == 1)||(fac == -1))
	{
	  hit=0;
	  if( (b[0] != 0) && (b[1] != 0) )
	    {
	      alpha = (c[0] + d[0] - a[0])/ b[0];
	      beta =  (c[1] + d[1] - a[1])/ b[1];
	      if(fabs(alpha -beta) < CROSS_EPSILON)
		{if( (alpha >= 0)&&(alpha <= 1))hit++;};
	      alpha = (c[0] - a[0]) / b[0];
	      beta = (c[1] - a[1])/ b[1];
	      if(fabs(alpha-beta) < CROSS_EPSILON)
		{if( (alpha >= 0)&&(alpha <= 1))hit++;};
	    }
	  else
	    {
	      if(b[0] == 0)
		{
		  alpha =  (c[1] + d[1] - a[1])/ b[1];
		  if( (alpha >= 0)&&(alpha <= 1)&&
		     (fabs(c[0]-a[0])<CROSS_EPSILON))
		    hit++;
		  alpha = (c[1] - a[1])/ b[1];
		  if( (alpha >= 0)&&(alpha <= 1))hit++;
		}
	      else
		{
		  alpha = (c[0] + d[0] - a[0])/ b[0];
		  if( (alpha >= 0)&&(alpha <= 1)&&
		     (fabs(c[1]-a[1])<CROSS_EPSILON))
		    hit++;
		  alpha = (c[0] - a[0]) / b[0];
		  if( (alpha >= 0)&&(alpha <= 1))hit++;
		};
	    };
	  
	  if(hit==2)
	    {
	      contact=TRUE;
	    };
	}
      else
	{
	  if((d[0] == 0)&&(d[1] ==0))
	    {
	      /* 	      printf("Fault %3i, hat die Laenge Null.\n", j+1);  */
	      contact=TRUE;
	    }
	  else
	    {
	      alpha=(b[0]*(c[1]-a[1])+b[1]*(a[0]-c[0]))/
		(b[1]*d[0]-b[0]*d[1]);
	      if(b[0]!=0)beta=(c[0]+alpha*d[0]-a[0])/b[0];
	      if(b[1]!=0)beta=(c[1]+alpha*d[1]-a[1])/b[1];
	      if(((alpha <= 1.0+CROSS_EPSILON)&&(alpha >= 0.0-CROSS_EPSILON))&&
		 ((beta <= 1.0+CROSS_EPSILON)&&(beta >= 0.0-CROSS_EPSILON)))
		{contact=TRUE;};
	    };
	};
    }
  else
    {
      contact=TRUE;
    };
  
  return(contact);
}


