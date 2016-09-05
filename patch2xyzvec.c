#include "interact.h"
#include "properties.h"

//
// reads in patch format and writes xyz coordinates of corners 
// and base vectors 
//


int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  COMP_PRECISION *dummy;
  int i;
  int opmode=XYZ_AND_VEC_MODE;
  if(argc != 1){
    fprintf(stderr,"%s \n\t reads in patch format from stdin and writes corners and vectors to stdout\n",
	    argv[0]);
    fprintf(stderr,"output format is x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 sx sy sz dx dy dz nx ny nz\n");
    fprintf(stderr,"where 1,2,3,4 are the four corners of the patch (FE ordering) and s,d, and n\n");
    fprintf(stderr,"are the strike, dip, and normal vectors\n");
    exit(-1);
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing coordinates and vectors to stdout\n",
	  argv[0]);
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,FALSE);
  for(i=0;i<medium->nrflt;i++)
    print_patch_geometry_and_bc(0,(fault+i),
				opmode,0.0,TRUE,stdout,FALSE,dummy);

  return 0;
}

