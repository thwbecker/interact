#include "interact.h"
#include "properties.h"
/*
 
reads in sets of four points in 3-D and converts to patch format

uses the fit_plane and points2patch subroutines, both are in
fit_plane.c

input format:

p_1^x p_1^y p_1^z 
p_2^x p_2^y p_2^z 
p_3^x p_3^y p_3^z 
p_4^x p_4^y p_4^z 

...
or

p_1^x p_1^y p_1^z code
p_2^x p_2^y p_2^z code
p_3^x p_3^y p_3^z code
p_4^x p_4^y p_4^z code

if code is set



...


repeated for every rectangle. points have to be in finite element
ordering (counterclockwise)

the opmode flag determines the output format, see interact.h

$Id: points2patch.c,v 1.8 2003/07/02 07:57:34 becker Exp $

*/

my_boolean read_points_local(COMP_PRECISION *,int *, my_boolean , FILE *);

int main(int argc, char **argv)
{
  
  struct flt fault[1],*patch;
  struct med *medium;
  COMP_PRECISION x[12],dx;
  COMP_PRECISION *dummy;
  int i,opmode=PATCH_OUT_MODE,nrpatches,seg[2],code;
  long seed = -1;
  my_boolean adjust_area=TRUE,use_code=FALSE;
  medium=(struct med *)calloc(1,sizeof(struct med));

  dx = 1.5;

  if(argc > 4){
    fprintf(stderr,"%s: usage:\n\t%s [adjust_area, %i] [dx, %g] [use_code, %i]\n\tread in four 3-D points per line from stdin\n",
	    argv[0],argv[0],adjust_area,dx,use_code);
    fprintf(stderr,"\tif adjust_area is set, will attempt to make input and output same area\n");
    fprintf(stderr,"\tformat:\n\tx1_x x1_y x1_z x2_x x2_y x2_z ...\n\n");
    fprintf(stderr,"\tto form a regular quad, points have to be in FE (CCW) ordering, starting lower left\n");
    fprintf(stderr,"\tassumes that patch can be described by dip and strike only, no rake!\n");
    fprintf(stderr,"\twrites patch format to stdout\n");
    fprintf(stderr,"\tdx (%g) will subdivide the faults into patches with dx width/length\n",dx);
    fprintf(stderr,"\tif use_code is set, will read in fault number codes\n");
    exit(-1);
  }
  if(argc >= 2){
    sscanf(argv[1],"%i",&i);
    adjust_area=(my_boolean)i;
  }
  if(argc > 2)
    sscanf(argv[2],ONE_CP_FORMAT,&dx);
  if(argc > 3)
    sscanf(argv[3],"%i",&i);
  use_code = (my_boolean)i;
  
  fprintf(stderr,"%s: reading points from %s, writing patch to stdout,",
	  argv[0],"stdin");
  if(adjust_area)
    fprintf(stderr," adjusting area");
  else
    fprintf(stderr," area unadjusted");
  fprintf(stderr,", spacing %g, use code: %i\n",dx,use_code);

  medium->nrflt=0;
  nrpatches=0;
  while(read_points_local(x,&code, use_code, stdin)){
    for(i=0;i<4;i++)
      if(x[i*3+INT_Z] > 0.0){
	fprintf(stderr,"%s: rectangle %i: point %i: z should be <= 0 (%g, %g, %g)\n",
		argv[0],medium->nrflt+1,i+1,x[i*3+INT_X],x[i*3+INT_Y],x[i*3+INT_Z]);
	exit(-1);
      }

    points2patch(fault,x,adjust_area);
    if(use_code)
      fault[0].group = code;
    else
      fault[0].group = medium->nrflt;
    //fprintf(stderr,"%s: set of points: %i code: %i\n",argv[0],medium->nrflt,fault[0].group);
    /* hack for now */
    seg[0] = 1 + fault[0].l/dx;
    seg[1] = 1 + fault[0].w/dx;
    /* subdivide */
    divide_fault_in_patches(0,fault,&patch,&nrpatches,
			    seg,FALSE,TRUE,0,0,&seed);
    medium->nrflt++;
  }
  for(i=0;i<nrpatches;i++){
    print_patch_geometry_and_bc(i,patch,opmode,0.0,
				FALSE,stdout,FALSE,dummy);
  }
  fprintf(stderr,"%s: read %i sets of four points\n",
	  argv[0],medium->nrflt);
  fprintf(stderr,"%s: produced total %i faults with %i patches\n",
	  argv[0],medium->nrflt,nrpatches);
  exit(0);
}
my_boolean read_points_local(COMP_PRECISION *x,int *code, my_boolean read_code, FILE *in)
{
  my_boolean ok=FALSE;
  int i,c=0;
  if(read_code){
    for(i=0;i<4;i++)
      c+= fscanf(stdin,THREE_CPI_FORMAT,&x[INT_X+i*3],&x[INT_Y+i*3],&x[INT_Z+i*3],code);
    if(c == 16)
      ok = TRUE;
  }else{
    for(i=0;i<4;i++)
      c+= fscanf(stdin,THREE_CP_FORMAT,&x[INT_X+i*3],&x[INT_Y+i*3],&x[INT_Z+i*3]);
    if(c == 12)
      ok = TRUE;
  }
  return ok;
}






