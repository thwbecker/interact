#include "interact.h"

/*

  reads in patch file and checks for positive 
  feedback between individual patches
  $Id: check_feedback.c,v 1.13 2003/03/02 01:37:41 becker Exp $

*/

int main(int argc,char **argv)
{
  struct med *medium;
  struct flt *fault;
  int evil_pair[2];
  my_boolean bailout=TRUE;
  medium=(struct med *)calloc(1,sizeof(struct med));
  switch(argc){
  case 2:{
    break;
  }
  case 3:{
    sscanf(argv[1],"%i",evil_pair);
    if(evil_pair[0]==0)
      bailout=FALSE;
    break;
  }
  default:{
    fprintf(stderr,"%s file.patch [0/1]\n\treads file.patch file in patch format (like %s)\n\tand checks for positive feedback patches\n",
	    argv[0],GEOMETRY_FILE);
    fprintf(stderr,"\tie., pairs of patches whose self-interaction is smaller than the cross triggering\n");
    fprintf(stderr,"\tif the second argument is set to zero, then program checks all patches\n");
    fprintf(stderr,"\tif it is 1, program will stop at first fatal pair (this is the default\n");
    exit(-1);
  }}
  // read in fault geometry and friction properties
  read_geometry(argv[1],&medium,&fault,TRUE,FALSE,FALSE,FALSE);
  // check all interactions
  if(!check_coulomb_stress_feedback(medium->nrflt,0,
				    fault,medium, TRUE,
				    CALC_I_COEFF_NOW,bailout,
				    evil_pair,-1.0))
    fprintf(stderr,"%s: could not find fatal interactions\n",argv[0]);
  else
    fprintf(stderr,"%s: found at least one fatal interaction\n",argv[0]);
  exit(0);
}
