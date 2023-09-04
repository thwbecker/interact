#include "interact.h"
#include "properties.h"

// reads in patch format and guesses activation modes
// for faults in bc.in format



int main(int argc, char **argv)
{
  struct flt *fault;
  struct med *medium;
  int i;
  int opmode=BC_OUT_MODE;
  COMP_PRECISION *dummy;
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(argc!=1){
    fprintf(stderr,"%s: reads in patch format from stdin and writes bc.in modes to stdout\n",
	    argv[0]);
    exit(-1);
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing bc activation modes to stdout\n",
	  argv[0]);
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,FALSE);
  for(i=0;i<medium->nrflt;i++)
    print_patch_geometry_and_bc(0,(fault+i),opmode,0.0,TRUE,stdout,FALSE,dummy);

  return 0;
}

