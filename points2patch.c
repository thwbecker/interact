#include "interact.h"
#include "properties.h"
/*
 
reads in sets of four points in 3-D and converts to patch format
either by

- fitting the best regular, Okada rectangle, or
- using iquads irregular quads

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



*/

int read_points_local(COMP_PRECISION *,int *, my_boolean , FILE *, int );

int main(int argc, char **argv)
{
  
  struct flt fault[1],*patch;
  struct med *medium;
  COMP_PRECISION x[12],dx,xc[3];
  COMP_PRECISION *dummy=NULL;
  int i,j,opmode=PATCH_OUT_MODE,nrpatches,seg[2],code;
  long seed = -1;
  my_boolean adjust_area=TRUE,use_code=FALSE;
  int iquad=0,np;
  medium=(struct med *)calloc(1,sizeof(struct med));
#ifndef ALLOW_NON_3DQUAD_GEOM
  // simply for ease of storage reasons
  fprintf(stderr,"%s needs ALLOW_NON_3DQUAD_GEOM to be set\n",argv[0]);exit(-1);
#endif

  //dx = 1.5;
  dx=0;

  if((argc > 5)||((argc>1)&&(strcmp(argv[1],"-h")==0)) ){
    fprintf(stderr,"%s: usage:\n\t%s [adjust_area, %i] [dx, %g] [use_code, %i] [iquad, %i]\n\tread in four or three 3-D points per line from stdin\n",
	    argv[0],argv[0],adjust_area,dx,(int)use_code,iquad);
    fprintf(stderr,"\tif adjust_area is set, will attempt to make input and output same area\n");
    fprintf(stderr,"\tformat:\n\tx1_x x1_y x1_z x2_x x2_y x2_z ...\n\n");
    fprintf(stderr,"\tto form a regular quad, points have to be in FE (CCW) ordering, starting lower left\n");
    fprintf(stderr,"\tassumes that patch can be described by dip and strike only, no rake!\n");
    fprintf(stderr,"\twrites patch format to stdout\n");
    fprintf(stderr,"\tdx (%g) will subdivide the faults into patches with dx width/length\n",dx);
    fprintf(stderr,"\t\t a <=0 dx will use no subdivision\n",dx);
    fprintf(stderr,"\tif use_code is set, will read in fault number codes\n");
    fprintf(stderr,"\tif iquad is zero, will try to fit Okada, 1: uses irregular quads and subdivides into triangles\n");
    fprintf(stderr,"\t\t2: expect three points, prints triangular elements\n");
    exit(-1);
  }
  if(argc >= 2){
    sscanf(argv[1],"%i",&i);
    adjust_area=(my_boolean)i;
  }
  if(argc > 2)
    sscanf(argv[2],ONE_CP_FORMAT,&dx);
  if(argc > 3){
    sscanf(argv[3],"%i",&i);
    use_code = (my_boolean)i;
  }
  if(argc > 4){
    sscanf(argv[4],"%i",&iquad);
  }
  if(iquad == 2){
    np=3;
  }else{
    np=4;
  }
  
  fprintf(stderr,"%s: reading sets of %i points from %s, writing patch to stdout,",
	  argv[0],np,"stdin");
  if(adjust_area){
    fprintf(stderr," adjusting area");
  }else{
    fprintf(stderr," area unadjusted");
  }
  fprintf(stderr,", spacing %g, use code %i, iquad: %i\n",dx,use_code,iquad);
  
  medium->nrflt=0;
  nrpatches=0;
  while(read_points_local(x,&code, use_code, stdin,np)==np){
    for(i=0;i < np;i++)
      if(x[i*3+INT_Z] > 0.0){
	fprintf(stderr,"%s: rectangle %i: point %i: z should be <= 0 (%g, %g, %g)\n",
		argv[0],medium->nrflt+1,i+1,x[i*3+INT_X],x[i*3+INT_Y],x[i*3+INT_Z]);
	exit(-1);
      }
    if(use_code)
      fault[0].group = code;
    else
      fault[0].group = medium->nrflt;
    if(iquad==2){		/* triangle output */
      for(i=0;i<3;i++){
	xc[i]=0;
	for(j=0;j < np;j++)
	  xc[i] += x[j*3+i];
	xc[i] /= (double)np;
      }
      /* print triangle patch format */
      fprintf(stdout,"%19.12e %19.12e %19.12e %10.6f %10.6f %19.12e %19.12e %6i ",
	      xc[INT_X],xc[INT_Y],xc[INT_Z],0.,0.,-1.,-1.,fault[0].group);
      fprintf(stdout,"%19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e\n",
	      x[0*3+INT_X],	      x[0*3+INT_Y],	      x[0*3+INT_Z],
	      x[1*3+INT_X],	      x[1*3+INT_Y],	      x[1*3+INT_Z],
	      x[2*3+INT_X],	      x[2*3+INT_Y],	      x[2*3+INT_Z]);
      nrpatches++;
      medium->nrflt++;
    }else if(iquad == 1){	/* irregular quad subdivided by three trianggles 

  2--------------1
  |\            /|
  | \          / |
  |  \   X    /  |
  |   \      /   |
  |    \    /    |
  | N1  \  / N2  |
  |      \/      |
  3 ---- 0 ----- 4
				*/
      for(i=0;i<3;i++){
	xc[i]=0;
	for(j=0;j < np;j++)
	  xc[i] += x[j*3+i];
	xc[i] /= (double)np;
      }
      fprintf(stdout,"%19.12e %19.12e %19.12e %10.6f %10.6f %19.12e %19.12e %6i ",
	      xc[INT_X],xc[INT_Y],xc[INT_Z],0.,0.,1.,-1.,fault[0].group);
      fprintf(stdout,"%19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e\n",
	      x[0*3+INT_X],	      x[0*3+INT_Y],	      x[0*3+INT_Z],
	      x[1*3+INT_X],	      x[1*3+INT_Y],	      x[1*3+INT_Z],
	      x[2*3+INT_X],	      x[2*3+INT_Y],	      x[2*3+INT_Z],
	      x[3*3+INT_X],	      x[3*3+INT_Y],	      x[3*3+INT_Z]);
      nrpatches++;
      medium->nrflt++;
    }else{
      points2patch(fault,x,adjust_area); /* fit best quad (if possible...) */
      medium->nrflt++;      
      if(dx<=0){
	/* no subdivision */
	seg[0]=seg[1]=1;
      }else{
	/* hack for now */
	seg[0] = 1 + fault[0].l/dx;
	seg[1] = 1 + fault[0].w/dx;
      }
      /* subdivide */
      divide_fault_in_patches(0,fault,&patch,&nrpatches,
			      seg,FALSE,TRUE,0,0,&seed,FALSE);

      for(i=0;i<nrpatches;i++){
	print_patch_geometry_and_bc(i,patch,opmode,0.0,
				    FALSE,stdout,FALSE,dummy);
      }
    }
  }
  fprintf(stderr,"%s: read %i sets of %i points\n",argv[0],medium->nrflt,np);
  fprintf(stderr,"%s: produced total %i faults with %i patches\n",
	  argv[0],medium->nrflt,nrpatches);
  exit(0);
}

int read_points_local(COMP_PRECISION *x,int *code, my_boolean read_code, FILE *in, int np)
{
  int i,j;
  for(i=0;i < np;i++){
    for(j=0;j < 3;j++){
      if(fscanf(in,ONE_CP_FORMAT,(x+i*3+j))!=1){
	return FALSE;
      }
    }
    if(read_code){
      if(fscanf(in,"%i",code)!=1){
	return FALSE;
      }
    }
  }
  return i;
}






