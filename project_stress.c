/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: project_stress.c,v 2.11 2003/02/13 22:45:12 becker Exp $
*/
#include "interact.h"
/*
  converts stress matrix into components resolved on
  fault 

  call without arguments for short description 
*/


int main(int argc, char **argv)
{
  struct flt *fault;
  my_boolean use_rake;
  COMP_PRECISION strike,dip,rake,sm[3][3],x,y,z,st,sn,sd,
    rake_vec[3],trac[3],cr,sr,alpha,sin_dip,cos_dip,sin_alpha,cos_alpha;
  int i;

  fault=calloc(1,sizeof(struct flt));
  switch(argc){
  case 3:{
    sscanf(argv[1],ONE_CP_FORMAT,&strike);
    sscanf(argv[2],ONE_CP_FORMAT,&dip);
    use_rake=FALSE;
    break;
  }
  case 4:{
    sscanf(argv[1],ONE_CP_FORMAT,&strike);
    sscanf(argv[2],ONE_CP_FORMAT,&dip);
    sscanf(argv[3],ONE_CP_FORMAT,&rake);
    use_rake=TRUE;
    break;
  }
  default:{
    fprintf(stderr,"%s strike dip [rake]\n",argv[0]);
    fprintf(stderr,"\treads stresses in format x y z s_xx s_xy s_xz s_yy s_yz s_zz from stdin\n");
    fprintf(stderr,"\tprojects them onto a fault oriented at angle strike in degrees clockwise from north\n");
    fprintf(stderr,"\tand dip in degress downward from horizontal\n\n");
    fprintf(stderr,"\toutput will be to stdout in format x y z s_strike s_dip s_normal\n\n");
    fprintf(stderr,"\tif rake is set (degrees counterclockwise from 0 (strike) to 90 (dip)) output is\n");
    fprintf(stderr,"\tx y z s_rake s_normal instead\n\n");
    exit(-1);
  }}
  fprintf(stderr,"%s: using strike: %g and dip: %g\n",argv[0],strike,dip);
  /* 
     convert to radians 
  */
  alpha = 90.0 - strike;
  alpha *= DEG2RAD;
  if((dip > 90) || (dip < 0)){
    fprintf(stderr,"%s: dip should be between 0 and 90\n",
	    argv[0]);
    exit(-1);
  }
  dip *= DEG2RAD;
  /* 
     calculate fault base vectors in normal, strike, and dip direction 
  */
  my_sincos(&sin_dip,&cos_dip,dip);
  my_sincos(&sin_alpha,&cos_alpha,alpha);
  calc_quad_base_vecs(fault->t_strike,fault->normal,fault->t_dip,
		      sin_alpha,cos_alpha,sin_dip,cos_dip);
  if(!use_rake){
    /* 
       output of three components,
       s_strike s_dip s_normal
    */
   
    while(fscanf(stdin,NINE_CP_FORMAT,&x,&y,&z,&sm[INT_X][INT_X],&sm[INT_X][INT_Y],&sm[INT_X][INT_Z],
		 &sm[INT_Y][INT_Y],&sm[INT_Y][INT_Z],&sm[INT_Z][INT_Z])==9){
      sm[INT_Y][INT_X] = sm[INT_X][INT_Y];
      sm[INT_Z][INT_X] = sm[INT_X][INT_Z];
      sm[INT_Z][INT_Y] = sm[INT_Y][INT_Z];
      calc_three_stress_components(sm,fault->normal,fault->t_strike,
				   fault->normal,fault->t_dip,&st,&sn,&sd);
#ifdef PRINT_BASE_VECTORS
	print_vector(fault->t_strike,3,stderr);
	print_vector(fault->t_dip,3,stderr);
	print_vector(fault->normal,3,stderr);
#endif
      fprintf(stdout,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
	      x,y,z,st,sd,sn);
    }
  }else{
    /* 
       rake is given, output of s_rake and s_normal 
    */
    my_sincos_deg(&sr,&cr,rake);
    for(i=0;i<3;i++){/* obtain the rake vector from 
			the strike and dip vectors 
		     */
      rake_vec[i] = cr * fault->t_strike[i];
      rake_vec[i]+= sr * fault->t_dip[i];
    } 
    while(fscanf(stdin,NINE_CP_FORMAT,
		 &x,&y,&z,&sm[INT_X][INT_X],&sm[INT_X][INT_Y],&sm[INT_X][INT_Z],
		 &sm[INT_Y][INT_Y],&sm[INT_Y][INT_Z],&sm[INT_Z][INT_Z])==9){
      sm[INT_Y][INT_X]=sm[INT_X][INT_Y];
      sm[INT_Z][INT_X]=sm[INT_X][INT_Z];
      sm[INT_Z][INT_Y]=sm[INT_Y][INT_Z];
      resolve_force(fault->normal,sm,trac);
      st = project_vector(trac,rake_vec);
      sn = project_vector(trac,fault->normal);
      fprintf(stdout,"%15.8e %15.8e %15.8e %15.8e %15.8e\n",x,y,z,st,sn);
    }
  }
  exit(0);
}
