/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: makefault.c,v 2.23 2003/03/02 07:34:11 becker Exp $
*/
#include "interact.h"
#include "properties.h"

int main(int argc, char **argv)
{
  int i,j,seg[2],grp,opmode=PATCH_OUT_MODE,nrpatches=0,*list;
  COMP_PRECISION l,w,d,x,y,z,s,a,stop_time,td[2]={-1,-1},srand,drand,*dummy;
  my_boolean circular,randomize_order=FALSE;
  long seed = -1;
#ifdef ALLOW_NON_3DQUAD_GEOM
  my_boolean point_source,twod=FALSE;
#endif 
  struct flt fault[1],*patch;
  a=1.0;/* aspect */
  w=1;
  l=a*w;/* half length */
  seg[0]=1;/* seg in strike dir */
  seg[1]=(int)(((COMP_PRECISION)seg[0])/a);/* seg in dip dir */
  x=0;/* location of centre point of patch */
  y=0;
  srand = drand = 0.0;
  z=-w;
  s=0;/* strike angle in degrees counterclockwise 
	 from west*/
  d=90;/* dip angle, 90 menaing vertical */
  grp=0;/* group indicator */
  circular=FALSE;/* fault patches will form a circular 
		    shape */
#ifdef ALLOW_NON_3DQUAD_GEOM
  point_source=FALSE;// point source instead of rectangular
#endif
  i=1;
  stop_time=10.0;
  while(i<argc){
    if(strcmp(argv[i],"-l")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&l);
      i+=2;
    }else if(strcmp(argv[i],"-w")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&w);
      i+=2;
    }else if(strcmp(argv[i],"-x")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&x);
      i+=2;
    }else if(strcmp(argv[i],"-y")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&y);
      i+=2;
    }else if(strcmp(argv[i],"-z")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&z);
      i+=2;
    }else if(strcmp(argv[i],"-dip")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&d);
      i+=2;
    }else if(strcmp(argv[i],"-strike")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&s);
      i+=2;
    }else if(strcmp(argv[i],"-tw")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&td[1]);
      i+=2;
    }else if(strcmp(argv[i],"-tl")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&td[0]);
      i+=2;
    }else if(strcmp(argv[i],"-srand")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&srand);
      i+=2;
    }else if(strcmp(argv[i],"-drand")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&drand);
      i+=2;
    }else if(strcmp(argv[i],"-seed")==0){
      sscanf(argv[i+1],"%i",&j);
      seed = (long)((j<0)?(j):(-j));
      if(!seed)
	seed=-1;
      i+=2;
    }else if(strcmp(argv[i],"-grp")==0){
      sscanf(argv[i+1],"%i",&grp);
      i+=2;
    }else if(strcmp(argv[i],"-n")==0){
      sscanf(argv[i+1],"%i",&seg[0]);
      i+=2;
    }else if(strcmp(argv[i],"-m")==0){
      sscanf(argv[i+1],"%i",&seg[1]);
      i+=2;
    }else if(strcmp(argv[i],"-stime")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&stop_time);
      i+=2;
    }else if(strcmp(argv[i],"-geom")==0){
      opmode= GEOMVIEW_MODE;
      i++; 
    }else if(strcmp(argv[i],"-xyz")==0){
      opmode= PSXYZ_MODE;
      i++; 
    }else if(strcmp(argv[i],"-rl")==0){
      randomize_order = TRUE;
      i++; 
    }else if(strcmp(argv[i],"-bc")==0){
      opmode= BC_OUT_MODE;
      i++; 
    }else if(strcmp(argv[i],"-circ")==0){
      circular=TRUE;
      i++; 
#ifdef ALLOW_NON_3DQUAD_GEOM
    }else if(strcmp(argv[i],"-point")==0){
      point_source=TRUE;
      i++; 
    }else if(strcmp(argv[i],"-twod")==0){
      twod=TRUE;
      i++; 
#endif
    }else{
      fprintf(stderr,"%s: %s produce patches filling a rectangular fault\n",
	      argv[0],argv[0]);
      fprintf(stderr,"%s: all length scales should be given with fault width set to unity\n",
	      argv[0]);
      fprintf(stderr,"%s: on output, they will be multiplied with the scaling factor %11e\n\n",
	      argv[0],CHAR_FAULT_DIM);
      fprintf(stderr,"%s: the following flags change default values which are given in brackets\n",
	      argv[0]);
      fprintf(stderr,"\tuse -x      value for x (%g)\n",x);
      fprintf(stderr,"\tuse -y      value for y (%g)\n",y);
      fprintf(stderr,"\tuse -z      value for z (%g)\n",z);
      fprintf(stderr,"\tuse -l      value for half-length (%g)\n",l);
      fprintf(stderr,"\tuse -w      value for half-width (%g)\n",w);
      fprintf(stderr,"\tuse -n      value for segments in strike dir (%i)\n",seg[0]);
      fprintf(stderr,"\tuse -m      value for segments in dip dir (%i)\n",seg[1]);
      fprintf(stderr,"\tuse -dip    value for dip angle (degrees down from horiz.)(%g)\n",d);
      fprintf(stderr,"\tuse -strike value for strike angle (deg clockw. from N)(%g)\n",s);
      fprintf(stderr,"\tuse -grp    value for group flag (%i)\n\n",grp);
      fprintf(stderr,"\tuse -tl     value for target length for patches (subdivides faults) (%g)\n",td[0]);
      fprintf(stderr,"\tuse -tw     value for target width  for patches (subdivides faults) (%g)\n",td[1]);
      fprintf(stderr,"\t            -tw and -tl override -n and -m settings if they are set to\n");
      fprintf(stderr,"\t            non-negative values\n\n");
      fprintf(stderr,"\tuse -srand  value for std of strike randomization on patch level  (%g)\n",srand);
      fprintf(stderr,"\tuse -drand  value for std of dip    randomization on patch level  (%g)\n",drand);
      fprintf(stderr,"\tuse -seed   value for random seed in case srand or drand != 0 (%i)\n\n",(int)seed);
      fprintf(stderr,"\tuse -stime  value for stop time in bc file (%g)\n",stop_time);
      fprintf(stderr,"\tuse -circ         for circular patch\n");
      fprintf(stderr,"\n\tuse -geom    for geomview 3D output\n");
      fprintf(stderr,"\tuse -xyz     for psxyz output\n");
      fprintf(stderr,"\tuse -bc      for bc.in output\n");
#ifdef ALLOW_NON_3DQUAD_GEOM
      fprintf(stderr,"\tuse -point   to produce a point source (i.e. write negative length)\n");
      fprintf(stderr,"\t             in this case, w will be the fault area, i.e. 4lw\n");
      fprintf(stderr,"\tuse -twod    to produce 2-D fault segments with zero width, m=1 in this case\n\n");
#endif
      fprintf(stderr,"\tuse -rl      to randomize the order of the patch output \n\n");
      exit(-1);
    }
  }
  if(td[0]>0){
    seg[0]=(int)(l/td[0]+0.5);
    if(seg[0]<1)seg[0]=1;
  }
  if(td[1]>0){
    seg[1]=(int)(w/td[1]+0.5);
    if(seg[1]<1)seg[1]=1;
  }
  if(circular){/* circular fault */
    seg[1]=MAX(seg[0],seg[1]);
    seg[0]=seg[1];
    l=w;
  }
#ifdef ALLOW_NON_3DQUAD_GEOM
  if(twod){			/* switch to 2D geometry */
    if(d != 90){
      fprintf(stderr,"%s: error, dip has to be 90 for 2D\n",
	      argv[0]);
      exit(-1);
    }
    w = 0.0;
    d = 90.0;
    z = 0.0;
    seg[1]=1;
  }
#endif
  fault[0].x[X]=x;
  fault[0].x[Y]=y;
  fault[0].x[Z]=z;
  fault[0].area= l * w * 4.0;
  fault[0].w = w;
  fault[0].l = l;
  fault[0].dip=d;
  fault[0].strike=s;
  fault[0].group=grp;
#ifdef ALLOW_NON_3DQUAD_GEOM
  if(point_source){
    fault[0].type = POINT_SOURCE;
    fault[0].l = -l/w;// - aspect_ratio
    fault[0].w = fault[0].area;
  }else{
    if(!twod)
      fault[0].type = RECTANGULAR_PATCH;
    else
      fault[0].type = TWO_DIM_SEGMENT_PLANE_STRAIN;
  }
#endif
  //
  // check for illegal values
  //
  check_fault_angles(fault);
  divide_fault_in_patches(0,fault,&patch,&nrpatches,seg,circular,TRUE,
			  srand,drand,&seed);
#ifdef ALLOW_NON_3DQUAD_GEOM
  fprintf(stderr,"makefault: x: (%g, %g, %g) l: %g w: %g s: %g d: %g grp: %i n: %5i m: %5i patches: %5i scaling: %11e %s\n",
	  fault[0].x[X],fault[0].x[Y],fault[0].x[Z],
	  fault[0].l,fault[0].w,fault[0].strike,
	  fault[0].dip,fault[0].group,seg[0],seg[1],nrpatches,CHAR_FAULT_DIM,
	  (twod)?("2-D"):((point_source)?("point source"):("")));
  if(circular)
    fprintf(stderr,"circular patches\n");
  else
    fprintf(stderr,"\n");
#else
  fprintf(stderr,"makefault: x: (%g, %g, %g) l: %g w: %g s: %g d: %g grp: %i n: %5i m: %5i patches: %5i scaling: %11e\n",
	  fault[0].x[X],fault[0].x[Y],fault[0].x[Z],
	  fault[0].l,fault[0].w,fault[0].strike,
	  fault[0].dip,fault[0].group,seg[0],seg[1],nrpatches,CHAR_FAULT_DIM);
#endif
  randomize_list(&list,nrpatches,randomize_order);
    
  for(i=0;i<nrpatches;i++){
    j = list[i];
    print_patch_geometry_and_bc(j,patch,opmode,
				stop_time,FALSE,stdout,FALSE,dummy);
  }
  free(list);
  exit(0);
}
