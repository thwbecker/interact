#include "interact.h"
/* 


check if two 3-D fault segments are far enough from each other 

leeway is the fraction distance at the edges, 
set to unity for no increased size


$Id: far_enough.c,v 1.6 2005/10/26 01:46:31 becker Exp $

 */
//#define SUPER_DEBUG

my_boolean far_enough(struct flt *fault1,struct flt *fault2,
		   COMP_PRECISION leeway)
{
  my_boolean outcome;
  COMP_PRECISION corner[4][3];
  float a[3],b[3],c[3],d[3],e[3],f[3],g[3],h[3];
  COMP_PRECISION corner2[4][3];
  int i,hit=0;

  calculate_bloated_corners(corner,fault1,leeway);
  calculate_bloated_corners(corner2,fault2,leeway);
  for(i=0;i<3;i++){
    a[i]=corner[0][i];b[i]=corner[1][i];
    c[i]=corner[2][i];d[i]=corner[3][i];
    e[i]=corner2[0][i];f[i]=corner2[1][i];
    g[i]=corner2[2][i];h[i]=corner2[3][i];
  }
  // use four fast triangle tests (slow, but ...)
  hit += tri_tri_intersect(a,b,d,e,f,h);
  hit += tri_tri_intersect(a,b,d,f,g,h);
  hit += tri_tri_intersect(b,c,d,e,f,h);
  hit += tri_tri_intersect(b,c,d,f,g,h);
  if(hit)
    outcome=FALSE;
  else
    outcome=TRUE;
  return(outcome);
}
