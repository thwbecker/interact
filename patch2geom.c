#include "interact.h"
#include "properties.h"
//
// reads in patch format and writes geom format for geomview
//
/* 

   this file can also produce colored slip COFF files using geom.in
   and a headerless flt.nohdr.dat file 


 */

int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  int i;
  FILE *in;COMP_PRECISION *dummy;
  my_boolean shrink_patches=FALSE,read_slip=FALSE;
  COMP_PRECISION fixed_range=0.025;
  char filename[2000];
  if(argc >= 2){
    sscanf(argv[1],"%i",&i);
    shrink_patches = (my_boolean)i;
  }
  if(argc == 3){
    sscanf(argv[2],"%i",&i);
    read_slip = (my_boolean)i;
  }
  if(argc > 3){
    fprintf(stderr,"%s [shrink_patches, 0] [read_slip, 0]\n\t reads in patch format from stdin and writes Geomview OFF format to stdout\n",
	    argv[0]);
    fprintf(stderr,"if shrink_patches is set, will make patches smaller for plotting\n");
    exit(-1);
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing OFF to stdout. shrink: %i read_slip: %i\n",
	  argv[0],shrink_patches,read_slip);
  
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,FALSE);
  if(read_slip){
    fprintf(stderr,"%s: attemping to read slip values from flt.nohdr.dat\n",argv[0]);
    in=myopen("flt.nohdr.dat","r");
    for(i=0;i<medium->nrflt;i++){ /* read in slip values from headerless file */
      if(fscanf(in,"%*f %*f %*f %*f %*f %lf %lf %lf %*f %*f %*f %*i %*i",
		&fault[i].u[STRIKE],&fault[i].u[DIP],&fault[i].u[NORMAL])!=3){
	fprintf(stderr,"%s: read error file  flt.nohdr.dat, patch %i\n",
		argv[0],i+1);
	exit(-1);
      }
    }
    for(i=0;i<medium->nrgrp;i++){
      sprintf(filename,"flt.%i.abs.off",i);
      print_group_data_geom(filename,medium,fault,i,0,fixed_range);
      sprintf(filename,"flt.%i.strike.off",i);
      print_group_data_geom(filename,medium,fault,i,1,fixed_range);
      sprintf(filename,"flt.%i.dip.off",i);
      print_group_data_geom(filename,medium,fault,i,2,fixed_range);
    }
  }else{
    for(i=0;i<medium->nrflt;i++)
      print_patch_geometry_and_bc(0,(fault+i),GEOMVIEW_MODE,
				  0.0,TRUE,stdout,shrink_patches,dummy);
  }
  exit(0);
}

