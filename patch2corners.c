#include "interact.h"
#include "properties.h"
//
// reads in patch format and writes x_1 y_1 z_1 .... x_4 y_4 z_4  to stdout
// where 1...4 are the four corners of the patch
// $Id: patch2corners.c,v 1.11 2003/03/02 01:37:41 becker Exp $ 


int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  COMP_PRECISION *dummy;
  int opmode=CORNEROUT_MODE,i;
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(argc!=1){
    fprintf(stderr,"%s: reads in patch format from stdin and writes x,y,z tripels\n\tfor each corners to stdout\n",
	    argv[0]);
    exit(-1);
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing geom to stdout\n",
	  argv[0]); 
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,FALSE);
   for(i=0;i<medium->nrflt;i++)
     print_patch_geometry_and_bc(0,(fault+i),opmode,0.0,TRUE,stdout,FALSE,dummy);

   return 0;
}

