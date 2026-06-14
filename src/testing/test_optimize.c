#include "interact.h"
/*
  reads in geometry and tries to optimize
  the fault patch ordering so that the 
  resulting i matrix bandwidth is reduced
  $Id: test_optimize.c,v 1.10 2003/02/13 22:45:12 becker Exp $
*/

int main(int argc,char **argv)
{
  struct med *medium;
  struct flt *fault;
  int i;
  COMP_PRECISION maxdist,*dummy;
  medium=(struct med *)calloc(1,sizeof(struct med));
  read_geometry(GEOMETRY_FILE,&medium,&fault,FALSE,FALSE,FALSE,FALSE);
  maxdist=max_dist(fault,medium);

  fprintf(stderr,"%s: total distance before: %g\n",
	  argv[0],penalty_dist(fault,medium,maxdist));
  optimize(fault,medium);
  fprintf(stderr,"%s: total distance after: %g\n",
	  argv[0],penalty_dist(fault,medium,maxdist));
  for(i=0;i<medium->nrflt;i++)
    print_patch_geometry_and_bc(i,fault,PATCH_OUT_MODE,0.0,
				FALSE,stdout,FALSE,dummy);
}
