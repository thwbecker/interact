/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: randomflt.c,v 2.29 2003/02/13 22:45:12 becker Exp $


  generates a set of randomly oriented, non-intersecting 
  faults


  /usr/bin/time randomflt -f 100 -m 5 -n 5 > geom.random
  Athlon 1400Mhz    367.790u
  Origin 200        962.209u 
  O2                1007.903u
*/
#include "interact.h"
#include "properties.h"
#include "randomflt.h"


int main(int argc, char **argv)
{
  int i,j,hit,n=1,seg[2]={1,1},random_mode,icnt,iclim,
    opmode,iter,evil_pair[2],nrpatches,old_nrpatches;
  COMP_PRECISION *dummy;
  long int seed=-1;
  COMP_PRECISION 
    w0=1.,/* length bounds, base length */
    dw=0.,// standard deviation in case of Gauss dist
    // otherwise range over which value varies
    dip0=90.,/* dip bounds, base dip */
    ddip=0.,// variation, and so on
    aspect0=2,/* aspect ratio variations */
    daspect=0.,
    x0=0.,y0=0.,/* location variations */
    dx=10.,dy=10.,
    z0=-3,/* depth variations */
    dz=0,
    strike0=0,/* strike variations */
    dstrike=10.0,
    srand = 0,			/* inter-patch variability */
    drand = 0,
    sin_dip,cos_dip,td[2]={-1,-1},asp,nval;
  COMP_PRECISION w1,dip1,stop_time=10.0,tmpdbl,
    aspect1,z1,x1,y1,strike1,bvalue,leeway=1.0;
  double alpha;
  my_boolean circular=FALSE,
    check_for_interactions=TRUE, /* interaction test */
    gr_random= FALSE,gra_random = FALSE,
    check_for_intersections = TRUE; /* intersection test */
  struct med *medium;
  struct flt *fault,*patch;
  // init medium with zeros
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("main: 1");

  /* default values */
  opmode=PATCH_OUT_MODE;

  random_mode=GAUSS_DISTRIBUTED;

  check_input_parameters(argc, argv,&n,&seed, seg,&bvalue, &gr_random,&gra_random,
			 td, &z0, &dz, &strike0, &dstrike,&dip0,&ddip,&w0, 
			 &dw, &aspect0, 
			 &daspect,&x0,  &dx,&y0,  &dy, &opmode, 
			 &check_for_interactions,&random_mode,&leeway,
			 &check_for_intersections,&srand,&drand);
  /* print patches with several segments */
  fprintf(stderr,"%s: n: %i seed: %i segs: %i %i mode: %i\n",
	  argv[0],n,((int)seed),seg[0],seg[1],opmode);

  /* construct base patches */
  if((fault=malloc(sizeof(struct flt)*n))==NULL)
    MEMERROR("main: 2");
  /* 
     assign constants for uniform deviates
  */
  // dip
  dip1 = dip0-ddip/2.0;
  /* width  */
  w1=w0 - dw/2.0;
  /* apsect ratio */
  aspect1=aspect0-daspect/2.0;
  /* location parameters */
  x1=x0-dx/2.0;
  y1=y0-dy/2.;
  /* strike */
  strike1=strike0-dstrike/2.0;
  // depth
  z1 = z0 - dz/2.0;
  //
  nval = -1.0 - bvalue;
  // intialize random number generator
  RAND_GEN(&seed);
  //
  // assign fault dimensions first 
  //
#define _TMP_RNLIM 100
  for(i=0;i<n;i++){
    icnt=0;
    do{
      if(gr_random || gra_random){		/* powerlaw distribution */
	/* assign either intended fault width or area from powerlaw */
	fault[i].w  =  mypower_randnr(w0,dw,nval,0,&seed);
      }else{
	assign_random(&fault[i].w,w0,w1,dw,&seed,random_mode);
      }
      icnt++;
    }while((fault[i].w <= 0.0)&&(icnt < _TMP_RNLIM));
    if(icnt == _TMP_RNLIM){
      fprintf(stderr,"%s: reached iteration limit in fault width assignment loop\n",
	      argv[0]);
      exit(-1);
    }
  }
  qsort(fault,n,sizeof(struct flt),compare_fault_width);
  if(!gra_random)
    if((fault[0].w > -z0)&&(dz==0)){
      fprintf(stderr,"%s: largest fault half width (%g) is larger than z0 (%g), adjusting\n",
	      argv[0],fault[0].w,z0);
      z0 = -fault[0].w;
    }
  for(nrpatches=i=0;i<n;i++){
    // need that for coulomb check 
    fault[i].mu_s=STATIC_MU;
    fault[i].group=i;
#ifdef ALLOW_NON_3DQUAD_GEOM
    fault[i].type=RECTANGULAR_PATCH;
#endif
    /* search for a random location that is far enough 
       from all other faults */
    iter=0;
    do{
      hit=FALSE;
      iter++;
      //
      // start assigning random numbers
      //
      assign_random(&tmpdbl,strike0,strike1,dstrike,&seed,random_mode);
      fault[i].strike = tmpdbl;
      
      if(!gra_random){
	assign_random(&tmpdbl,dip0,dip1,ddip,&seed,random_mode);
	fault[i].dip = tmpdbl;
	my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[i].dip);
	
	// assign depth
	iclim = 300;
	icnt=0;
	tmpdbl = fault[i].w * sin_dip;
	do{
	  assign_random(&fault[i].x[INT_Z],z0,z1,dz,&seed,random_mode);
	  icnt++;
	}while((fault[i].x[INT_Z]  + tmpdbl > 0) && (icnt < iclim));
	if(icnt == iclim){
	  fprintf(stderr,"%s: reached iteration limit of %i in fault depth assignment loop\n",
		  argv[0],iclim);
	  fprintf(stderr,"%s: w: %g dip: %g z0: %g dz: %g\n",argv[0],fault[i].w,fault[i].dip,
		  z0,dz);
	  fprintf(stderr,"%s: ztop: %g\n",argv[0],fault[i].x[INT_Z]  + tmpdbl);
	  exit(-1);
	}
	asp = -1;
	while(asp <= 0)
	  assign_random(&asp,aspect0,aspect1,daspect,&seed,random_mode);
      }else{
	/* don't allow zero dip */
	sin_dip = 0;
	while(fabs(sin_dip) < 1e-5){
	  assign_random(&tmpdbl,dip0,dip1,ddip,&seed,random_mode);
	  fault[i].dip = tmpdbl;
	  my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[i].dip);
	}
	/* 
	   area dominated  
	*/
	fault[i].area = fault[i].w;
	// assign depth
	iclim = 10000;
	icnt=0;
	do{
	  assign_random(&fault[i].x[INT_Z],z0,z1,dz,&seed,random_mode);
	  asp = -1;
	  while(asp <= 0)
	    assign_random(&asp,aspect0,aspect1,daspect,&seed,random_mode);
	  fault[i].w = sqrt(fault[i].area/(4.*asp));
	  tmpdbl = fault[i].w * sin_dip;
	  //fprintf(stderr,"%g a: %g %g %g\n",fault[i].w,asp,fault[i].x[INT_Z],sin_dip);
	  icnt++;
	}while(((fault[i].x[INT_Z] + tmpdbl > 0 )||
		(fault[i].w > (MAX_SEIS_DEPTH+fault[i].x[INT_Z])/sin_dip))
	       &&(icnt < iclim));
	if(icnt == iclim){
	  fprintf(stderr,"%s: -gra reached iteration limit of %i in fault depth assignment loop\n",
		  argv[0],iclim);
	  fprintf(stderr,"%s: area: %g dip: %g z0: %g\n",argv[0],fault[i].area,fault[i].dip,
		  z0);
	  exit(-1);
	}
      }
      // obtain half length and area from half width and aspect ratio
      fault[i].l = fault[i].w * asp;
      fault[i].area = fault[i].l * fault[i].w * 4.0;


      assign_random(&fault[i].x[INT_X],x0,x1,dx,&seed,random_mode);
      assign_random(&fault[i].x[INT_Y],y0,y1,dy,&seed,random_mode);

     
      if(!hit){
	// check for illegal values
	check_fault_angles((fault+i));
	/* derived quantities, need them for base vectors */
	alpha=90.0-(double)fault[i].strike;
	// get angles
	my_sincos_degd(&fault[i].sin_alpha,&fault[i].cos_alpha,alpha);
	my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[i].dip);
	calc_base_vecs(fault[i].t_strike,fault[i].normal,fault[i].t_dip,
		       fault[i].sin_alpha,fault[i].cos_alpha,sin_dip,cos_dip);
	if(check_for_intersections){
	  // compare with all old faults and see if intersects
	  for(j=0;j<i;j++){
	    if(!far_enough((fault+i),(fault+j),leeway)){
	      hit=TRUE;
	      break;
	    }
	  }
	}
	if(!hit){// add to list of patches
	  old_nrpatches=nrpatches;
	  create_patches(i,fault,&patch,&nrpatches,seg,circular,FALSE,td,
			 srand,drand,&seed);
	  if(check_for_interactions){
	    hit=incrementally_check_interaction_coefficients(
		   patch,nrpatches,medium,FALSE,evil_pair,-1.0);
	    if(hit){
	      nrpatches=old_nrpatches;
	    }
	  }
	}
      }
    }while((hit)&&(iter<MAXITER));
    if(iter==MAXITER){
      fprintf(stderr,"%s: too many iterations (%i) at fault %i\n",
	      argv[0],iter,i);
      exit(-1);
    }
  } 
  // patch output
  for(i=0;i<nrpatches;i++)
    print_patch_geometry_and_bc(i,patch,opmode,stop_time,FALSE,
				stdout,FALSE,dummy);
  return 0;
}




void check_input_parameters(int argc, char **argv, int *n,long *seed, int *seg,
			    COMP_PRECISION *bvalue, my_boolean *gr_random,
			    my_boolean *gra_random,
			    COMP_PRECISION *td, 
			    COMP_PRECISION *z0, COMP_PRECISION *dz,
			    COMP_PRECISION *strike0, COMP_PRECISION *dstrike,
			    COMP_PRECISION *dip0,COMP_PRECISION *ddip,
			    COMP_PRECISION *w0, COMP_PRECISION *dw,
			    COMP_PRECISION *aspect0,COMP_PRECISION *daspect,
			    COMP_PRECISION *x0, COMP_PRECISION *dx,
			    COMP_PRECISION *y0, COMP_PRECISION *dy,
			    int *opmode,my_boolean *check_for_interactions,
			    int *random_mode,
			    COMP_PRECISION *leeway,
			    my_boolean *check_for_intersections,
			    COMP_PRECISION *srand, COMP_PRECISION *drand)
{
  int i,j;
  i=1;
  while(i < argc){
    if(strcmp(argv[i],"-f")==0){
      sscanf(argv[i+1],"%i",n);
      i+=2;
    }else if(strcmp(argv[i],"-seed")==0){
      sscanf(argv[i+1],"%i",&j);
      *seed = (long)((j<0)?(j):(-j));
      if(! *seed)
	*seed=-1;
      i+=2;
    }else if(strcmp(argv[i],"-n")==0){
      sscanf(argv[i+1],"%i",seg);
      i+=2;
    }else if(strcmp(argv[i],"-gr")==0){// b-value for GR type width distribution
      sscanf(argv[i+1],ONE_CP_FORMAT,bvalue);
      *gr_random = TRUE;
      i+=2;
    }else if(strcmp(argv[i],"-gra")==0){// b-value for GR type area distribution
      sscanf(argv[i+1],ONE_CP_FORMAT,bvalue);
      *gra_random = TRUE;
      i+=2;
    }else if(strcmp(argv[i],"-m")==0){
      sscanf(argv[i+1],"%i",(seg+1));
      i+=2;
    }else if(strcmp(argv[i],"-tw")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,(td+1));
      i+=2;
    }else if(strcmp(argv[i],"-tl")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,td);
      i+=2;
    }else if(strcmp(argv[i],"-zm")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,z0);
      i+=2;
     }else if(strcmp(argv[i],"-sm")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,strike0);
      i+=2;
    }else if(strcmp(argv[i],"-dm")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,dip0);
      i+=2;
    }else if(strcmp(argv[i],"-srand")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,srand);
      i+=2;
    }else if(strcmp(argv[i],"-drand")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,drand);
      i+=2;
    }else if(strcmp(argv[i],"-wm")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,w0);
      i+=2;
    }else if(strcmp(argv[i],"-am")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,aspect0);
      i+=2;
    }else if(strcmp(argv[i],"-xm")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,x0);
      i+=2;
    }else if(strcmp(argv[i],"-ym")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,y0);
      i+=2;
    }else if(strcmp(argv[i],"-zstd")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,dz);
      i+=2;
    }else if(strcmp(argv[i],"-sstd")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,dstrike);
      i+=2;
    }else if(strcmp(argv[i],"-wstd")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,dw);
      i+=2;
    }else if(strcmp(argv[i],"-dstd")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,ddip);
      i+=2;
    }else if(strcmp(argv[i],"-astd")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,daspect);
      i+=2;
    }else if(strcmp(argv[i],"-xstd")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,dx);
      i+=2;
    }else if(strcmp(argv[i],"-ystd")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,dy);
      i+=2;
    }else if(strcmp(argv[i],"-xyz")==0){
      *opmode=PSXYZ_MODE;
      i++;
    }else if(strcmp(argv[i],"-geom")==0){
      *opmode=GEOMVIEW_MODE;
      i++;
    }else if(strcmp(argv[i],"-nfic")==0){
      *check_for_interactions=FALSE;
      i++;
    }else if(strcmp(argv[i],"-nis")==0){
      *check_for_intersections=FALSE;
      i++;
    }else if(strcmp(argv[i],"-l")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,leeway);
      i+=2;
    }else if(strcmp(argv[i],"-bc")==0){
      *opmode=BC_OUT_MODE;
      i++;
    }else if(strcmp(argv[i],"-gauss")==0){
      *random_mode=GAUSS_DISTRIBUTED;
      i++;
    }else if(strcmp(argv[i],"-uniform")==0){
      *random_mode=UNIFORM_DISTRIBUTED;
      i++;
    }else{
      fprintf(stderr,"%s: %s produces randomly spaced and oriented rectangular faults\n",
	      argv[0],argv[0]);
      fprintf(stderr,"%s: the output is in \"patch\" format as used for the geom.in files\n",argv[0]);
      fprintf(stderr,"%s: for interact\n",argv[0]);
      fprintf(stderr,"options (default values are in brackets)\n\n");
      fprintf(stderr,"\tuse -f    value  for number of faults                             (%i)\n",*n);
      
      fprintf(stderr,"\tuse -sm   value  to change the mean of the strike                 (%g)\n",*strike0);
      fprintf(stderr,"\tuse -dm   value  to change the mean of the dip                    (%g)\n",*dip0);
      fprintf(stderr,"\tuse -am   value  to change the mean of the aspect                 (%g)\n",*aspect0);
      fprintf(stderr,"\tuse -wm   value  to change the mean of the half width             (%g)\n",*w0);
      fprintf(stderr,"\tuse -xm   value  to change the mean of the x coordinate           (%g)\n",*x0);
      fprintf(stderr,"\tuse -xm   value  to change the mean of the y coordinate           (%g)\n",*y0);

      fprintf(stderr,"\tuse -zm   value  to change the mean of the depth                  (%g)\n",*z0);
      
      fprintf(stderr,"\tuse -sstd value  to change the standard deviation of the strike   (%g)\n",*dstrike);
      fprintf(stderr,"\tuse -dstd value  to change the standard deviation of the dip      (%g)\n",*ddip);
      fprintf(stderr,"\tuse -astd value  to change the standard deviation of the aspect   (%g)\n",*daspect);
      fprintf(stderr,"\tuse -wstd value  to change the standard deviation of the width    (%g)\n",*dw);
      fprintf(stderr,"\tuse -xstd value  to change the standard deviation of the x coord. (%g)\n",*dx);
      fprintf(stderr,"\tuse -ystd value  to change the standard deviation of the y coord. (%g)\n",*dy);

      fprintf(stderr,"\tuse -zstd value  to change the standard deviation of the depth    (%g)\n\n",*dz);
       
      fprintf(stderr,"\tuse -l    value  length fraction to use for intersection checks   (%g)\n",*leeway);
      fprintf(stderr,"                   (e.g. l=1.1 will add 10%% to patch length)\n\n");

      fprintf(stderr,"\tuse -n    value  for segments in strike dir                       (%i)\n",seg[0]);
      fprintf(stderr,"\tuse -m    value  for segments in dip dir                          (%i)\n",seg[1]);
      fprintf(stderr,"\tuse -tl   value  for target half-length (divides faults into patches)\n");
      fprintf(stderr,"\tuse -tw   value  for target half-width  (divides faults into patches)\n");
      fprintf(stderr,"\t                 -tw and -tl override -n and -m settings\n");
 
      fprintf(stderr,"\tuse -srand  value for std of strike randomization on patch level  (%g)\n",*srand);
      fprintf(stderr,"\tuse -drand  value for std of strike randomization on patch level  (%g)\n",*drand);


      fprintf(stderr,"\tuse -seed value  for random number seed (<0)                      (%i)\n",
	      (int)*seed);
      fprintf(stderr,"\tuse -geom        for geomview 3D output\n");
      fprintf(stderr,"\tuse -xyz         for psxyz output\n");
      fprintf(stderr,"\tuse -bc          for bc.in output, guessing the activation mode of faults\n\n");

      fprintf(stderr,"\tuse -gauss       to have sll random parameters vary according\n");
      fprintf(stderr,"\t                 to a Gauss distribution (this is the default,\n");
      fprintf(stderr,"\t                 stddev actually means standard deviation)\n");
      fprintf(stderr,"\tuse -uniform     to have random variable vary uniformly. Then, stddev\n");
      fprintf(stderr,"\t                 indicates the range dx over which variables vary from\n");
      fprintf(stderr,"\t                 x_0-dx/2 through x_0+dx/2\n");

      fprintf(stderr,"\tuse -gr value    to have the width vary according to a powerlaw P(x) = C x^(-1-value)\n");
      fprintf(stderr,"\t                 between w1 and w2 extrema (C is determined from normalization),\n");
      fprintf(stderr,"\t                 and b = -value(?) is some sort of b-value (note that area = LW)\n");
      fprintf(stderr,"\t                 In this case, -wm and -wstd will set the w1 and w2 range\n\n");

      fprintf(stderr,"\tuse -gra value   to have the area vary according to a powerlaw P(x) = C x^(-1-value)\n");
      fprintf(stderr,"\t                 between a1 and a2 extrema (C is determined from normalization),\n");
      fprintf(stderr,"\t                 and b = -value(?) is some sort of b-value\n");
      fprintf(stderr,"\t                 In this case, -wm and -wstd will set the a1 and a2 range, dip != 0,\n");
      fprintf(stderr,"\t                 and the aspect ratio will be adjusted for a maximum depth of %g\n\n",
	      MAX_SEIS_DEPTH);

      fprintf(stderr,"\tuse -nfic        to suppress the check for fatal interactions\n\n");
      fprintf(stderr,"\tuse -nis         to suppress the check for intersecting faults\n\n");
      
      exit(-1);
    }
  }
}

