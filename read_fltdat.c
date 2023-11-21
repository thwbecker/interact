/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: read_fltdat.c,v 1.1 2011/01/07 07:19:58 becker Exp becker $
*/
#include "properties.h"
#include "interact.h"
#include <string.h>

/* 

   read in the fault slip file after the fact, need fault to be
   initialized

   this reads in the version with the two header lines, unless told
   otherwise

   return 1 if file found, zero else
*/

int read_fltdat(char *filename,struct flt *fault,struct med *medium, 
		my_boolean verbose)
{
  int i,grp,igrp,ipatch,iread;
  size_t size_dummy;
  FILE *in;
  char stdin_string[6]="stdin";
  COMP_PRECISION coulomb,mud_normal;
  char *line=NULL;
  if(!medium->geometry_init){
    fprintf(stderr,"read_fltdat: error, initialize fault geometry first\n");
    exit(-1);
  }
  if(!strings_match(filename,stdin_string)){
    if(verbose)
      fprintf(stderr,"read_fltdat: reading flt slip data file \"%s\"\n",
	      filename);
    in = fopen(filename,"r");
    if(!in){
      fprintf(stderr,"read_fltdat: could not open %s, will not have slip information\n",filename);
      return 0;
    }
  }else{
    fprintf(stderr,"read_fltdat: stdin\n");
    in=stdin;
  }

  /* 
     read two header 
  */
  getline(&line,&size_dummy,in);free(line);line=NULL;
  getline(&line,&size_dummy,in);free(line);line=NULL;

  for(grp=0;grp<medium->nrgrp;grp++)
    for(i=0;i<medium->nrflt;i++){
      if(fault[i].group == grp){
	if((iread=fscanf(in,FLTDAT_FORMAT,&fault[i].pos[0],&fault[i].pos[1],&fault[i].area,&coulomb,
		  &mud_normal,&fault[i].u[STRIKE],&fault[i].u[NORMAL],&fault[i].u[DIP],
			 &fault[i].s[STRIKE],&fault[i].s[NORMAL],&fault[i].s[DIP],&ipatch,&igrp))!=13){
	  fprintf(stderr,"read_fltdat: read error, %i items\n",iread);
	  exit(-1);
	}
	if(ipatch != i){
	  fprintf(stderr,"read_fltdat: patch assign error: %i %i\n",ipatch, i);
	  exit(-1);
	}
	if(igrp != grp){
	  fprintf(stderr,"read_fltdat: patch assign error: %i %i\n",igrp, grp);
	  exit(-1);
	}

      }
    }
  fclose(in);
  if(verbose)			/* should really check, but oh well */
    fprintf(stderr,"read_fltdat: read %i patches OK\n",medium->nrflt);
  return 1;

}
