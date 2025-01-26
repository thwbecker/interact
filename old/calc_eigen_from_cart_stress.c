/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: calc_eigen_from_cart_stress.c,v 1.4 2003/02/13 22:45:12 becker Exp $

*/
#include "interact.h"
/*

  converts stress matrix into principal stress components

*/
int main(int argc, char **argv)
{
  COMP_PRECISION sm[3][3],x,y,z,eval[3],evec[9],strike,dip;
  my_boolean z_unity = FALSE,p_angle=FALSE;
  int i;
  for(i=1;i<argc;i++){
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){// help
      fehelp(argv);
    }
    else if(strcmp(argv[i],"-u")==0){// normalize to evec z compoenent = 1 
      z_unity = TRUE;
    }
    else if(strcmp(argv[i],"-a")==0){// output vectors as strike and dip 
      p_angle = TRUE;
    }
    else{
      fprintf(stderr,"%s: option \"%s\" not understood\n\n",argv[0],argv[i]);
      fehelp(argv);
    }
  }
  while(fscanf(stdin,NINE_CP_FORMAT,&x,&y,&z,
	       &sm[X][X],&sm[X][Y],&sm[X][Z],
	       &sm[Y][Y],&sm[Y][Z],&sm[Z][Z])==9){
    sm[Y][X] = sm[X][Y];sm[Z][X] = sm[X][Z];sm[Z][Y] = sm[Y][Z];
    // calculate eigenvalues
    eigensystem3d(sm,eval,evec);
    fprintf(stdout,"%g %g %g\t\t",x,y,z);
    for(i=2;i>=0;i--){
      if(p_angle){
	vec_to_angles((evec+i*3),&dip,&strike);
	fprintf(stdout,"%g %g %g       \t\t",eval[i],strike,dip);
      }else{
	if((z_unity)&&(fabs(evec[i*3+Z])>EPS_COMP_PREC))
	  fprintf(stdout,"%g %9.6f %9.6f %9.6f\t\t",
		  eval[i],evec[i*3+X]/evec[i*3+Z],evec[i*3+Y]/evec[i*3+Z],1.0);
	else
	  fprintf(stdout,"%g %9.6f %9.6f %9.6f\t\t",
		  eval[i],evec[i*3+X],evec[i*3+Y],evec[i*3+Z]);
      }
    }
    fprintf(stdout,"\n");
  }
  exit(0);
}
void fehelp(char **argv)
{
  fprintf(stderr,"%s\n",argv[0]);
  fprintf(stderr,"\treads stresses in format\n\t x y z s_xx s_xy s_xz s_yy s_yz s_zz\n\t from stdin\n");
  fprintf(stderr,"\tand converts stresses into principal stress axis\n");
  fprintf(stderr,"\toutput will be to stdout in format\n");
  fprintf(stderr,"\tx y z s1 s1x s1y s1z s2 s2x s2y s2z s3 s3x s3y s3z\n");
  fprintf(stderr,"\twhere s1>s2>s3 are the eigenvalues and six siy siz (i=1,2,3) the respective eigenvector components\n");
  fprintf(stderr,"\n\toptions:\n");
  fprintf(stderr,"\t-u\tnormalize eigenvectors such that z component is unity (else unity length normalization)\n");
  fprintf(stderr,"\t-a\toutput of eigenvectors in strike and dip angles (this ordering, see interact for definition) instead of cartesian\n");
  fprintf(stderr,"\t\tthis overrides -u\n");
  fprintf(stderr,"\n");
  exit(-1);
}
