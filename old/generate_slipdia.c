/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: generate_slipdia.c,v 1.4 2003/02/13 22:45:12 becker Exp $
*/
#include "interact.h"
/*

  given a background stress direction, resolve the strike and dip stresses
  on various fault geometries and calculate the resulting plunge vector 
  
  output is in a x-y z system where x-y are a lower hemisphere projection
  of the strike and dip of the plane and z is from the resulting plunge

*/


int main(int argc, char **argv)
{
  struct flt *fault;
  COMP_PRECISION strike,dip,ddip,dstrike,sm[3][3],xhemi[3],stress[3],
    tstress[3],alpha,sin_dip,cos_dip,sin_alpha,cos_alpha;

  ddip=5.0;dstrike=5.0;
  fault=calloc(1,sizeof(struct flt));
  switch(argc){
  case 1:{
    break;
  }
  case 2:{
    sscanf(argv[1],ONE_CP_FORMAT,&dstrike);
    break;
  }
  case 3:{
    sscanf(argv[1],ONE_CP_FORMAT,&dstrike);
    sscanf(argv[2],ONE_CP_FORMAT,&ddip);
    break;
  }
 default:{
    fprintf(stderr,"%s [dstrike ddip]\n",argv[0]);
    fprintf(stderr,"\treads stresses in format s_xx s_xy s_xz s_yy s_yz s_zz from stdin\n");
    fprintf(stderr,"\tloops through all possible strike/dip combinations\n");
    fprintf(stderr,"\tand outputs the largest sheaer stress drop direction motion in a\n");
    fprintf(stderr,"\tlower hemisphere projecton (x-y) and plunge (z) system\n");
    exit(-1);
  }}
  // read in stress state
  fprintf(stderr,"%s: reading s_ij (6 components)\n",argv[0]);
  fscanf(stdin,SIX_CP_FORMAT,&sm[X][X],&sm[X][Y],&sm[X][Z],
	 &sm[Y][Y],&sm[Y][Z],&sm[Z][Z]);
  sm[Y][X] = sm[X][Y];
  sm[Z][X] = sm[X][Z];
  sm[Z][Y] = sm[Y][Z];
  fprintf(stderr,"%s: sxx: %g sxx: %g sxx: %g sxx: %g sxx: %g sxx: %g\n",
	  argv[0],sm[X][X],sm[X][Y],sm[X][Z],sm[Y][Y],sm[Y][Z],sm[Z][Z]);
  
  
  for(strike = 0.0;strike <= 360.0;strike += dstrike)
    for(dip=0.0;dip<=90.0;dip+=ddip){
      /*
	convert to radians 
      */
      alpha = 90.0 - strike;
      // get the lower hemisphere projection
      calc_lhemi_proj(dip, strike, xhemi);
      /* 
	 calculate fault base vectors in normal, strike, and dip direction 
      */
      my_sincos_deg(&sin_dip,&cos_dip,dip);
      my_sincos_deg(&sin_alpha,&cos_alpha,alpha);
      calc_base_vecs(fault->t_strike,fault->normal,fault->t_dip,
		     sin_alpha,cos_alpha,sin_dip,cos_dip);
      calc_three_stress_components(sm,fault->normal,fault->t_strike,
				   fault->normal,fault->t_dip,
				   (stress+STRIKE),(stress+NORMAL),
				   (stress+DIP));
      // obtain the target stress drop given the regional stresses
      get_maxsdir_stress_drops2(stress,1.0,tstress);
      //
      // output is x and y in the hemisphere projection and 
      // the resulting plunge angle
      fprintf(stdout,"%g %g %g\n",xhemi[X],xhemi[Y],
	      atan2(tstress[DIP],tstress[STRIKE])/PI);

    }
  exit(0);
}
