/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu


  program to evaluate the displacement and stresses due
  to a triangular element by Gauss integration

  unclear if this ever worked or could work?!

  input is the observational point x, the triangular
  coordinates xt, and the displacement in disp

  disp is in strike, dip, and normal format

  output is u_global (displacements) and sm_global
  (stress matrix)

*/
#include "interact.h"
#include "math.h"

#ifdef ALLOW_NON_3DQUAD_GEOM
void eval_triangle_gauss(COMP_PRECISION *x,struct flt *fault,
			 COMP_PRECISION *disp,COMP_PRECISION *u_global, 
			 COMP_PRECISION sm_global[3][3],int *giret)
{
  //
  // can be 1, 3, 7, or 13, for slip
  //
  // 1 and 7, however, include the centroid which we 
  // use to evaluate the stress. since values will be singular
  // there, should use 3 or 13
  //
#define GAUSS_POINTS 3
  static COMP_PRECISION g[GAUSS_POINTS],
    h[GAUSS_POINTS],w[GAUSS_POINTS];
  static my_boolean init=FALSE;
  COMP_PRECISION u[3],sm[3][3],xf[3];
  int i,iret;
  if(!init){
    get_gauss_points(g,h,w,GAUSS_POINTS);
    init=TRUE;
  }
  /*
    sum up several point sources to obtain the effect
    of a triangle
  */
  u_global[INT_X]=u_global[INT_Y]=u_global[INT_Z]=0.0;
  sm_global[INT_X][INT_X]=sm_global[INT_X][INT_Y]=sm_global[INT_X][INT_Z]=
    sm_global[INT_Y][INT_Y]=sm_global[INT_Y][INT_Z]=sm_global[INT_Z][INT_Z]=0.0;
  *giret=0;
  for(i=0;i<GAUSS_POINTS;i++){
    // assign the coordinate in element
    globalx(fault->xn,g[i],h[i],xf);
    /*
      fprintf(stderr,"eval_triangle: gp %2i: x: %g/%g/%g, g/h/w: %g/%g/%g, cen: %g/%g/%g\n",
      i,xf[INT_X],xf[INT_Y],xf[INT_Z],g[i],h[i],w[i],
      fault->x[INT_X],fault->x[INT_Y],fault->x[INT_Z]);
    */
    // add contribution of point source, are weighted by w[i]
    // at Gauss integration point location
    // fault->w holds the area of the triangle
    eval_point_short(x,xf,fault->w*w[i],
		     fault->sin_alpha,
		     fault->cos_alpha,
		     (COMP_PRECISION)fault->dip,
		     disp,u,sm,&iret);
    if(!iret){// sum up contribution
      u_global[INT_X] += u[INT_X];
      u_global[INT_Y] += u[INT_Y];
      u_global[INT_Z] += u[INT_Z];
      sm_global[INT_X][INT_X] += sm[INT_X][INT_X];
      sm_global[INT_X][INT_Y] += sm[INT_X][INT_Y];
      sm_global[INT_X][INT_Z] += sm[INT_X][INT_Z];
      sm_global[INT_Y][INT_Y] += sm[INT_Y][INT_Y];
      sm_global[INT_Y][INT_Z] += sm[INT_Y][INT_Z];
      sm_global[INT_Z][INT_Z] += sm[INT_Z][INT_Z];
    }else{
      fprintf(stderr,"eval_triangle: gauss point %i is undefined: x: (%g, %g, %g) xg: (%g, %g, %g)\n",
	      i,x[INT_X],x[INT_Y],x[INT_Z],xf[INT_X],xf[INT_Y],xf[INT_Z]);
      fprintf(stderr,"eval_triangle: xt: (%g, %g, %g) (%g, %g, %g) (%g, %g, %g)\n",
	      fault->xn[INT_X  ],fault->xn[INT_Y  ],fault->xn[INT_Z  ],
	      fault->xn[3+INT_X],fault->xn[3+INT_Y],fault->xn[3+INT_Z],
	      fault->xn[6+INT_X],fault->xn[6+INT_Y],fault->xn[6+INT_Z]);
      *giret = iret;
    }
  }
  sm_global[INT_Y][INT_X]=sm_global[INT_X][INT_Y];
  sm_global[INT_Z][INT_X]=sm_global[INT_X][INT_Z];
  sm_global[INT_Z][INT_Y]=sm_global[INT_Y][INT_Z];
}
#endif
/*

  initialize points for Gauss integration over triangular
  element
  
  use 3, 7, or 13 for n
  
*/

#define MAX_GAUSS_POINTS 13

void get_gauss_points(COMP_PRECISION *gg, 
		      COMP_PRECISION *gh, 
		      COMP_PRECISION *gw,
		      int nrgp)
{
  static COMP_PRECISION g[MAX_GAUSS_POINTS],
    h[MAX_GAUSS_POINTS],w[MAX_GAUSS_POINTS];
  int i;
  static my_boolean init=FALSE;
  if(!init){
    if(nrgp > MAX_GAUSS_POINTS){
      fprintf(stderr,
	      "init_gauss_points: too manmy Gauss points requested, %i\n",
	      nrgp);
      exit(-1);
    }
    switch(nrgp){
    case 1:{// at centroid, only one
      g[0]=1.0/3.0;
      h[0]=g[0];
      w[0]=1.0;
      break;
    }
    case 3:{ /* triangle formula second order */
      g[0]=1.0/6.0;
      g[1]=2.0/3.0;
      g[2]=g[0];
      h[0]=g[0];
      h[1]=g[0];
      h[2]=g[1];
      w[0]=(1.0/3.0);
      w[1]=w[0];
      w[2]=w[0];
      break;
    }
    case 7:{
      /* triangle formula fifth order */
      g[0]=0.1012865073235;
      g[1]=0.7974269853531;
      g[2]=g[0];
      g[3]=0.4701420641051;
      g[4]=g[3];
      g[5]=0.0597158717898;
      g[6]=1.0/3.0;
      h[0]=g[0];
      h[1]=g[0];
      h[2]=g[1];
      h[3]=g[5];
      h[4]=g[3];
      h[5]=g[3];
      h[6]=g[6];
      w[0]=0.1259391805448;
      w[1]=w[0];
      w[2]=w[0];
      w[3]=0.1323941527885;
      w[4]=w[3];
      w[5]=w[3];
      w[6]=0.225;
      break;
    }
    case 13:{
      /* triangle formula seventh order */
      g[0]=0.0651301029022;
      g[1]=0.8697397941956;
      g[2]=g[0];
      g[3]=0.3128654960049;
      g[4]=0.6384441885698;
      g[5]=0.0486903154253;
      g[6]=g[4];
      g[7]=g[3];
      g[8]=g[5];
      g[9]= 0.2603459660790;
      g[10]=0.4793080678419;
      g[11]=g[9];
      g[12]=1.0/3.0;
      h[0]=g[0];
      h[1]=g[0];
      h[2]=g[1];
      h[3]=g[5];
      h[4]=g[3];
      h[5]=g[4];
      h[6]=g[5];
      h[7]=g[4];
      h[8]=g[3];
      h[9]=g[9];
      h[10]=g[9];
      h[11]=g[10];
      h[12]=g[12];
      w[0]=0.0533472356088;
      w[1]=w[0];
      w[2]=w[0];
      w[3]=0.0771137608903;
      w[4]=w[3];
      w[5]=w[3];
      w[6]=w[3];
      w[7]=w[3];
      w[8]=w[3];
      w[9]=0.1756152574332;
      w[10]=w[9];
      w[11]=w[9];
      w[12]=-0.1495700444677;
      break;
    }
    default:{
      fprintf(stderr,
	      "init_gauss_points: for triangles, %i is undefined\n",
	      nrgp);
      exit(-1);
      break;
    }}
    init=TRUE;
  }
  for(i=0;i<nrgp;i++){
    gg[i]=g[i];
    gh[i]=h[i];
    gw[i]=w[i];
  }
}









