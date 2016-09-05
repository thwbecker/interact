#include "interact.h"

/*
  reads in patch file and determines average quantities
  of the groups of faults that the patches form
  $Id: patch2group.c,v 1.9 2003/03/02 01:37:41 becker Exp $
*/

int main(int argc,char **argv)
{
  struct med *medium;
  struct flt *fault;
  struct geog *grp;
  int i;
  switch(argc){
  case 2:{
    break;
  }
  default:{
    fprintf(stderr,"%s file.patch\n",argv[0]);
    fprintf(stderr,"\treads file.patch file in patch format (like %s)\n\tand determines geometrical quantities of fault groups\n",
	    GEOMETRY_FILE);
    exit(-1);
  }}
  // read in patches
  read_geometry(argv[1],&medium,&fault,FALSE,FALSE,FALSE,FALSE);
  /* init group structs */
  if(!(grp=(struct geog *)calloc(medium->nrgrp,
				 sizeof(struct geog))))
    MEMERROR("main");
  // determine average geometrical quantities of the groups
  calc_group_geometry(medium,fault,grp);
  fprintf(stderr,"%s: writing geometry of groups of patches to stdout\n",argv[0]);
  // output in patch format 
  for(i=0;i<medium->nrgrp;i++)
    fprintf(stdout,"%g %g %g %g %g %g %g %i\n",
	    grp[i].center[X],grp[i].center[Y],grp[i].center[Z],
	    vec_to_strike(grp[i].strike_vec),vec_to_dip(grp[i].dip_vec),
	    grp[i].prange[STRIKE],grp[i].prange[DIP],i);

  exit(0);
}
