#include "interact.h"
#include "properties.h"
//
// reads in patch format and writes all vertices to stdout
// irrespective of connectivity


int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  int i;
  COMP_PRECISION tarea;
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(argc!=1){
    fprintf(stderr,"%s: reads in patch format from stdin and writes area for each patch to stdout\n",
	    argv[0]);
    exit(-1);
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing geom to stdout\n",
	  argv[0]); 
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,FALSE);
  tarea=0;
  for(i=0;i<medium->nrflt;i++){
    tarea += fault[i].area;
    if(fabs(patch_area((fault+i))-fault[i].area) > EPS_COMP_PREC)
      fprintf(stderr,"WARNING: fault %12.5e and recomputed %12.5e differ by %e\n",fault[i].area,patch_area((fault+i)),fabs(fault[i].area-patch_area((fault+i))));
    printf("%g\n", fault[i].area);
  }
  fprintf(stderr,"%s: total area: %e\n",argv[0],tarea);
  return 0;
}

