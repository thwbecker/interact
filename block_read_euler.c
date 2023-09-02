/*

read in euler poles for constrained blocks, if flag 
constrain_euler is true. in this case, will open file and try
to read euler poles for the constrained blocks. else, will
set nrbc to zero

if constrain_euler == 1: read from file
   constrain_euler == 2: lock last block in list

$Id: block_read_euler.c,v 1.5 2004/10/05 01:09:46 becker Exp $

*/

#include "interact.h"
#include "blockinvert.h"


void read_constrained_euler_poles(struct bmd *mod,char **argv,
				  int constrain_euler)
{
  int i,n,lock_block;
  FILE *in;
  COMP_PRECISION o[3];
  static COMP_PRECISION gfac = BLOCK_GFAC;
#ifndef BLOCK_SPHERICAL
  fprintf(stderr,"read_constrained_euler_poles only works for spherical\n");
  exit(-1);
#endif
  if(BLOCK_NBASE != 3){
    fprintf(stderr,"read_constrained_euler_poles: NBASE != 3\n");
    exit(-1);
  }
  if((!mod->block_init)||(!mod->nrb)){
    fprintf(stderr,"read_constrained_euler_poles: error: block structure not init or zero blocks\n");
    exit(-1);
  }
  mod->nrbc = 0;
  n = mod->nrb * 3;
  /* 
     make sure that solution is initialized 
  */
  my_vecrealloc(&mod->xsol,n,"xsol");
  for(i=0;i < n;i++)
    mod->xsol[i] = 0.0;
  switch(constrain_euler){
  case 1:
    /* 

    read in constraints from file

    */
    /* array for constrained motions */
    in=myopen(EULER_POLE_FILE,"r");
    /* 
       reads in euler motions in format
       
       block_code[1...nrb] wx wy wz (deg/Myr)
       
    */
    fprintf(stderr,"%s: reading in fixed Euler poles (block[1...nrb] wx wy wz (deg/Myr)) from %s\n",
	    argv[0],EULER_POLE_FILE);
    while(fscanf(in,"%i %lf %lf %lf",&lock_block,o,(o+1),(o+2)) == 4){
      mod->nrbc++;
      if((lock_block < 1) || (lock_block > mod->nrb)){ /* check range for block code */
	fprintf(stderr,"%s: error: code %i out of bounds\n",
		argv[0],lock_block);
	exit(-1);
      }
      lock_block--;			/* move to 0 .. nrb-1 system */
      if(mod->block[lock_block].rot_c){
	fprintf(stderr,"%s: error: block %i was constrained already\n",
		argv[0],lock_block+1);
	exit(-1);
      }
      mod->block[lock_block].rot_c = TRUE;	/* set constrained flag */
      for(i=0;i < 3;i++)	/* assign pole to global solution array,
				   scale from deg/Myr to internal format */
	mod->xsol[lock_block*3+i] = o[i] * gfac;
    }
    fclose(in);
    if(!mod->nrbc){
      fprintf(stderr,"%s: error: no Euler poles were read in\n",
	      argv[0]);
      exit(-1);
    }else{
      fprintf(stderr,"%s: WARNING: %i block motions will be contrained\n",
	      argv[0],mod->nrbc);
    }
    break;
  case 2:
    /* 
       constrain last block 
    */
    lock_block = mod->nrb-1;
    fprintf(stderr,"%s: WARNING: locking last block %c\n",argv[0],bname(lock_block));
    mod->nrbc++;
    mod->block[lock_block].rot_c = TRUE;	/* set constrained flag */
    for(i=0;i < 3;i++)	/* lock last block */
      mod->xsol[lock_block*3 + i] = 0.0;
    break;
  case 0:
    break;
  default:
    fprintf(stderr,"%s: error: constrain_euler %i not defined\n",argv[0],constrain_euler);
    exit(-1);
    break;
  }
  if(mod->nrbc){
    /* 
       check if we have only assigned blocks from the end of the list
       and output 
    */
    mod->first_c = -1;
    for(i=0;i < mod->nrb;i++){
      if(mod->block[i].rot_c){
	if(mod->first_c == -1)
	  mod->first_c = i;
	fprintf(stderr,"%s: block %c motion constrained to: wx: %11g wy: %11g wz: %11g (deg/Myr)\n",
		argv[0],bname(i),mod->xsol[i*3+INT_X]/gfac,mod->xsol[i*3+INT_Y]/gfac,
		mod->xsol[i*3+INT_Z]/gfac);
      }else{
	/* test for continuity */
	if(mod->first_c > 0){
	  fprintf(stderr,"%s: error: block %i was assigned, implying that all block with\n",
		  argv[0],mod->first_c+1);
	  fprintf(stderr,"%s: error: higher code number have to be constrained, too. %i is not.\n",
		  argv[0],i+1);
	  exit(-1);
	}
      }
    }
  }
  /* 
     number of contrained parameters 
  */
  mod->nc = mod->nrbc*BLOCK_NBASE;
}
