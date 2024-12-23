#include "interact.h"
#include "properties.h"
//
// reads in patch format and writes all vertices to stdout
// irrespective of connectivity


int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  COMP_PRECISION *dummy=NULL;
  int opmode=VERTEXOUT_MODE,i,nvert;
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(argc!=1){
    fprintf(stderr,"%s: reads in patch format from stdin and writes x,y,z tripels\n\tfor each corners to stdout\n",
	    argv[0]);
    exit(-1);
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing geom to stdout\n",
	  argv[0]); 
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,FALSE);
  nvert = 0;
  for(i=0;i<medium->nrflt;i++)
    nvert += print_patch_geometry_and_bc(0,(fault+i),opmode,0.0,TRUE,stdout,FALSE,dummy);
  fprintf(stderr,"%s: wrote %i path sets of vertices, %i total\n",argv[0],medium->nrflt,nvert);
  return 0;
}

