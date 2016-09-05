/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu
  

  $Id: randomize_strike.c,v 1.47 2003/02/13 22:45:12 becker Exp $
*/
#include "interact.h"
#include "properties.h"

//#define SUPER_DEBUG
// what sort of distributions? 
//#define USE_UNIFORM_RANDOM
#define USE_GAUSS_DISTRIBUTED
// max nr of iterations before program gives up
#define MAXITER 100000
// depth limits for z variations
#define UPPER_DEPTH 0.0
#define LOWER_DEPTH -15.0
// max distance for interaction check as a fraction of the 
// max half length of a fault patch
#define MAX_INTER_DIST_FRAC 1.0



int main(int argc, char **argv)
{
  int i,j,hit,nr_inp_flts,seg[2],opmode,iter,iter2,segn[2],
    nrpatches,old_nrpatches,interact_prob,mindi[2]={0,0},
    evil_pair[2];
    COMP_PRECISION *dummy;
  long int seed,titer;
  my_boolean depth_adjust=FALSE,fix_first=FALSE,
    fix_second=FALSE,sort=TRUE,
    circular=FALSE,check_for_interactions=TRUE;
  COMP_PRECISION s01,mindist,tmpdbl,
    s02,ds,alpha,sin_dip,cos_dip,dip,mu_s,
    leeway,oldz,tmpw,depthstd,lower_depth,max_dist,maxl,
    upper_depth,dastdw,dipstddev,first_strike,
    second_strike,td[2],srand=0,drand=0;
  struct flt *fault,*patch;
  struct med *medium;
#ifdef DEBUG
  COMP_PRECISION max_distance;
#endif
#ifdef NO_COMMAND_LINE_PARAMETERS
  FILE *in;
#endif

  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("main: 1");
#ifdef USE_UNIFORM_RANDOM
  COMP_PRECISION s11,s12,s2;
#endif
  /* initialize */
  s01=0;
  s02=0;
  ds=10.0;
  leeway=1.0;
  td[0]=-1;
  td[1]=-1;
  seed=-1;
  seg[0]=1;
  seg[1]=1;
  opmode=PATCH_OUT_MODE;
  lower_depth=LOWER_DEPTH;
  upper_depth=UPPER_DEPTH;
  dastdw=2.0;
  dipstddev=0.0;
  max_dist = MAX_INTER_DIST_FRAC;
  first_strike=second_strike=0.0;
#ifndef   NO_COMMAND_LINE_PARAMETERS
  i=1;
  while(i<argc){
    if(strcmp(argv[i],"-seed")==0){
      sscanf(argv[i+1],"%i",&j);
      seed=(long)(j<0)?(j):(-j);
      if(!seed)
	seed=-1;
      i+=2;
    }else if(strcmp(argv[i],"-n")==0){
      sscanf(argv[i+1],"%i",&seg[0]);
      i+=2;
    }else if(strcmp(argv[i],"-s01")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&s01);
      i+=2;
    }else if(strcmp(argv[i],"-s02")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&s02);
      i+=2;
    }else if(strcmp(argv[i],"-ds")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&ds);
      i+=2;
    }else if(strcmp(argv[i],"-l")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&leeway);
      i+=2;
    }else if(strcmp(argv[i],"-m")==0){
      sscanf(argv[i+1],"%i",&seg[1]);
      i+=2;
    }else if(strcmp(argv[i],"-tw")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&td[1]);
      i+=2;
    }else if(strcmp(argv[i],"-tl")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&td[0]);
      i+=2;
    }else if(strcmp(argv[i],"-xyz")==0){
      opmode=PSXYZ_MODE;
      i++;
    }else if(strcmp(argv[i],"-da")==0){
      depth_adjust=TRUE;
      i++;
    }else if(strcmp(argv[i],"-nsort")==0){
      sort=FALSE;
      i++;
    }else if(strcmp(argv[i],"-nfis")==0){
      check_for_interactions=FALSE;
      i++;
    }else if(strcmp(argv[i],"-ff")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&first_strike);
      fix_first=TRUE;
      i+=2;
    }else if(strcmp(argv[i],"-fs")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&second_strike);
      fix_second=TRUE;
      i+=2;
    }else if(strcmp(argv[i],"-mdf")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&max_dist);
      i+=2;
    }else if(strcmp(argv[i],"-dal")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&lower_depth);
      i+=2;
    }else if(strcmp(argv[i],"-dau")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&upper_depth);
      i+=2;
    }else if(strcmp(argv[i],"-sdip")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&dipstddev);
      i+=2;
    }else if(strcmp(argv[i],"-dasz")==0){
      sscanf(argv[i+1],ONE_CP_FORMAT,&dastdw);
      i+=2;
    }else if(strcmp(argv[i],"-geom")==0){
      opmode=GEOMVIEW_MODE;
      i++;
    }else if(strcmp(argv[i],"-bc")==0){
      opmode=BC_OUT_MODE;
      i++;
    }else{
      fprintf(stderr,"%s: randomize the strike component of input fault patches\n",
	      argv[0]);
      fprintf(stderr,"%s: input is a list of patches in the format\n\n",argv[0]);
      fprintf(stderr,"x y z dip half-length half-width group\n\n\n");
      fprintf(stderr,"%s: program will attempt to arrange strike directions randomly\n",argv[0]);
      fprintf(stderr,"%s: the output is in \"patch\" format as used for the geom.in files\n",
	      argv[0]);
      fprintf(stderr,"%s: for interact. The program randomly select from the mean strike values\n",argv[0]);
      fprintf(stderr,"%s: s01 and s02, makes these values the same if you want only one mean orientation.\n",argv[0]);
      fprintf(stderr,"%s: Will sort patches such that longest faults are placed first and group numbers go\n",argv[0]);
      fprintf(stderr,"%s: from 0 (largest) to N-1 (smallest input patch)\n",argv[0]);
      fprintf(stderr,"%s: options (numbers in brackets are default values):\n\n",argv[0]);
      fprintf(stderr,"\t    -s01   value for first  mean strike direction (%g)\n",s01);
      fprintf(stderr,"\t    -s02   value for second mean strike direction (%g)\n",s02);
      fprintf(stderr,"\t    -ds    value for standard deviation in strike direction (%g)\n",ds);
      fprintf(stderr,"\t    -ff    value to fix the strike (and dip) of patches that belong to group 0\n");
      fprintf(stderr,"\t                strike is then set to value\n");
      fprintf(stderr,"\t    -fs    value to fix the strike (and dip) of patches that belong to group 1\n");
      fprintf(stderr,"\t                strike is then set to value\n");
      fprintf(stderr,"\t    -l     value to prevent intersections around faults that\n");
      fprintf(stderr,"\t                appear to be longer and wider by a factor of value (%g)\n",
	      leeway);
      fprintf(stderr,"\t    -sdip  value to allow for random variations in the dip\n");
      fprintf(stderr,"\t                in this case, the std dev of dip (Gaussian) is given by value (%g)\n",
	      dipstddev);
      fprintf(stderr,"\t    -da         to allow for depth adjustment such that more patches can be fit\n");
      fprintf(stderr,"\t    -dal   value to change the lower bound for depth adjustment to value (%g)\n",
	      lower_depth);
      fprintf(stderr,"\t    -dau   value to change the upper bound for depth adjustment to value (%g)\n",
	      upper_depth);
      fprintf(stderr,"\t    -dasz  value for std dev of depth variations (Gaussian) as multiples of patch width (%g)\n",
	      dastdw);
      fprintf(stderr,"\t    -tl    value for target half-length (divides faults into patches)\n");
      fprintf(stderr,"\t    -tw    value for target half-width  (divides faults into patches)\n");
      fprintf(stderr,"\t                -tw and -tl override -n and -m settings\n");
      fprintf(stderr,"\t    -n     value for segments in strike dir\n");
      fprintf(stderr,"\t    -m     value for segments in dip dir\n");
      fprintf(stderr,"\t    -seed  value for random number seed (has to be smaller than zero) (%i)\n",
	      (int)seed);
      fprintf(stderr,"\t    -geom       for geomview 3D output\n");
      fprintf(stderr,"\t    -xyz        for psxyz output\n");
      fprintf(stderr,"\t    -bc         for bc.in output\n\n");
      fprintf(stderr,"\t    -nsort      to suppress sorting of patches and renaming of group numbers\n\n");
      fprintf(stderr,"\t    -nfis       to suppress searching for fatal interaction situations\n\n");
      fprintf(stderr,"\t    -mdf   value to set the maximum distance between fault patch center points\n");
      fprintf(stderr,"\t                 that is considered for an interaction loop check as a fraction\n");
      fprintf(stderr,"\t                 of the maximum fault group half length.\n");
      fprintf(stderr,"\t                 Default: %g, use -1 to check all patches.\n",
	      MAX_INTER_DIST_FRAC);
      fprintf(stderr,"\tuse -h          for this help page\n\n");
      exit(-1);
    }
  }

#endif
  /* print patches with several segments */
  fprintf(stderr,"%s: mean strike: %g/%g stddev(strike): %g l': %g*l rand_seed: %i segs(s/r): %i/%i mode: %i\n",
	  argv[0],s01,s02,ds,leeway,((int)seed),seg[0],seg[1],opmode);
  if(depth_adjust){
    fprintf(stderr,"%s: adjusting z values to fit more patches, lower/upper bounds: %g/%g\n",
	    argv[0],lower_depth,upper_depth);
    fprintf(stderr,"%s: standard deviation of depth distribution is %g times fault width\n",
	    argv[0],dastdw);
  }
  if(fix_first)
    fprintf(stderr,"%s: fixing strike of first fault group patches to %g\n",
	    argv[0],first_strike);
  if(fix_second)
    fprintf(stderr,"%s: fixing strike of second fault group patches to %g\n",
	    argv[0],second_strike);
  fprintf(stderr,"%s: allowing variation in dip for each patch with std_dev %g\n",
	  argv[0],dipstddev);
  if(check_for_interactions){
    fprintf(stderr,"%s: checking for positive interaction loops while adding patches.\n",
	    argv[0]);
    if(max_dist == -1)
      fprintf(stderr,"%s: checking all patches, no matter what distance\n",argv[0]);
    else
      fprintf(stderr,"%s: max centerpoint distance considered is %g times max fault group half length\n",
	      argv[0],max_dist);
  }else
    fprintf(stderr,"%s: checking for positive interaction loops suppressed\n",
	    argv[0]);

  nr_inp_flts=0;
  /* construct base patches */
  if((fault=malloc(sizeof(struct flt)*(nr_inp_flts+1)))==NULL)
    MEMERROR("main: 2");
#ifndef NO_COMMAND_LINE_PARAMETERS
  fprintf(stderr,"%s: expecting to read fault group locations from stdin\n%s: (use \"%s\" -h for help page)\n",
	  argv[0],argv[0],argv[0]);
#define INPUT_CHANNELL stdin
#else
  depth_adjust=TRUE;
  in=myopen("tmp.cat","r");
  fprintf(stderr,"%s: WARNING: NO PARAMETERS, reading from tmp.cat, allowing for depth adjustment\n",argv[0]);
#define INPUT_CHANNELL in
#endif
  /* 
     input of locations and lenght/widths
  */
  maxl=0;
  while(fscanf(INPUT_CHANNELL,RS_CP_FORMAT,
	       &fault[nr_inp_flts].x[X], &fault[nr_inp_flts].x[Y], 
	       &fault[nr_inp_flts].x[Z], &fault[nr_inp_flts].dip,
	       &fault[nr_inp_flts].l, &fault[nr_inp_flts].w,
	       &fault[nr_inp_flts].group) == 7){
    if((fault[nr_inp_flts].x[Z]<lower_depth)||(fault[nr_inp_flts].x[Z]>upper_depth)){
      fprintf(stderr,"%s: z value fault %i (%g) does not work with preset bounds %g (lower) and %g (upper)\n",
	      argv[0],nr_inp_flts,fault[nr_inp_flts].x[Z],lower_depth,upper_depth);
      exit(-1);
    }
    // for now don't check for variable friction
    fault[nr_inp_flts].mu_s = STATIC_MU;
#ifdef ALLOW_NON_3DQUAD_GEOM
    fault[nr_inp_flts].type=RECTANGULAR_PATCH;
#endif
    if(fault[nr_inp_flts].l > maxl)
      maxl=fault[nr_inp_flts].l;
    fault[nr_inp_flts].area = fault[nr_inp_flts].l * fault[nr_inp_flts].w * 4.0;
    nr_inp_flts++;
    if((fault=realloc(fault,sizeof(struct flt)*(nr_inp_flts+1)))==NULL)
      MEMERROR("main: 3");
  } 
  fprintf(stderr,"%s: read in %i fault locations\n",argv[0],nr_inp_flts);
  if(!nr_inp_flts){
    fprintf(stderr,"%s: that's not enough\n",argv[0]);
    exit(-1);
  }
  // scla max distance
  if(max_dist != -1.0){
    max_dist *= maxl;
    if(check_for_interactions)
      fprintf(stderr,"%s: determined max group length: %g, max dist therefore: %g\n",
	      argv[0],maxl,max_dist);
  }
  /* 
     sorting the patches according to length, longest patch
     first
  */
  if(sort){ 
    fprintf(stderr,"%s: sorting patches according to length and renamed group numbers\n",
	    argv[0]);
    qsort(fault,nr_inp_flts,sizeof(struct flt),compare_fault_length);
    for(i=0;i<nr_inp_flts;i++)
      fault[i].group = i;
    fprintf(stderr,"%s: sorting done\n",argv[0]);
  }else
    fprintf(stderr,"%s: sorting suppressed\n",argv[0]);
  /*
    determine average friction coefficient and estimate 
    total number of patches
  */
  for(mindist=FLT_MAX,mu_s=0.0,nrpatches=i=0;i<nr_inp_flts;i++){
    mu_s += fault[i].mu_s;
    determine_segments(seg,segn,(fault+i),FALSE,td);
    // determine total number of patches
    nrpatches += segn[0]*segn[1];
    for(j=i+1;j<nr_inp_flts;j++)
      if((tmpdbl=distance_3d(fault[i].x,fault[j].x))<mindist){
	mindist=tmpdbl;
	mindi[0]=i;mindi[1]=j;
      }
  } 
  mu_s /= (COMP_PRECISION)nr_inp_flts;
  fprintf(stderr,"%s: average static friction for coulomb check is %g\n",
	  argv[0],mu_s);
  fprintf(stderr,"%s: minimum distance is %12.5e between fault %i and %i\n",
	  argv[0],mindist,mindi[0],mindi[1]);
  fprintf(stderr,"%s: estimated total number of patches: %i\n",
	  argv[0],nrpatches);
  /* strike parameters for uniform directions */
#ifdef USE_UNIFORM_RANDOM
  s11=s01-ds;s2=2.0*ds;
  s12=s02-ds;s2=2.0*ds;
#endif  
  // intialize random number generator
  RAND_GEN(&seed);
  /*
    loop through all input faults to produce subdivided 
    patches
  */
  titer=0;
  nrpatches=0;
  fprintf(stderr,"%s: starting loop through all input faults\n",argv[0]);
  for(i=0;i<nr_inp_flts;i++){
    // save original dip angle
    dip=fault[i].dip;
#ifdef ALLOW_NON_3DQUAD_GEOM
    fault[i].type=RECTANGULAR_PATCH;
#endif
    /* 
       start randomizing loop here
    */
    interact_prob=iter=0;
    do{
      if((dipstddev != 0.0) && ((!fix_first) || (fault[i].group != 0)) && 
	 ((!fix_second) || (fault[i].group != 1)))
	fault[i].dip= dip + mygauss_randnr(dipstddev,&seed);
      /* 
	 assign either of two strike angles by flipping
	 a coin if not suppressed for first patch
      */
      if(((!fix_first)  || (fault[i].group != 0)) && 
	 ((!fix_second) || (fault[i].group != 1))) {
#ifdef  USE_UNIFORM_RANDOM
	if(myrandnr(1.0,&seed) < 0.5)
	  fault[i].strike = s11 + myrandnr(s2,&seed);
	else
	  fault[i].strike = s12 + myrandnr(s2,&seed);
#endif
#ifdef USE_GAUSS_DISTRIBUTED
	if(myrandnr(1.0,&seed) < 0.5)
	  fault[i].strike = s01 + mygauss_randnr(ds,&seed);
	else
	  fault[i].strike = s02 + mygauss_randnr(ds,&seed);
#endif
      }
      if((fix_first) && (fault[i].group==0)){// fix the strike of the first fault
	fault[i].strike=first_strike;
      }
      if((fix_second) && (fault[i].group==1)) {// fix the strike of the second fault
       fault[i].strike = second_strike;
      }
      // check for illegal values
      check_fault_angles((fault+i));
      /* 
	 derived quantities, need them for base vectors 
      */
      my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[i].dip);
      alpha=90.0-(COMP_PRECISION)fault[i].strike;
      // get angles
      my_sincos_deg(&fault[i].sin_alpha,&fault[i].cos_alpha,(COMP_PRECISION)alpha);
      calc_base_vecs(fault[i].t_strike,fault[i].normal,fault[i].t_dip,
		     fault[i].sin_alpha,fault[i].cos_alpha,sin_dip,cos_dip);
      // move depth around to allow for more patches to be inserted
      if(depth_adjust){
	tmpw=fault[i].w*sin_dip;
	oldz=fault[i].x[Z];
	// this is the range that z will be spread over
	depthstd=dastdw*fault[i].w;
	fault[i].x[Z] = oldz ;//+ mygauss_randnr(depthstd,&seed);
	iter2=0;
	while(((fault[i].x[Z]+tmpw > upper_depth)||
	       (fault[i].x[Z]-tmpw < lower_depth))&&
	      (iter2 < MAXITER)){
	  fault[i].x[Z] = oldz + mygauss_randnr(depthstd,&seed);
	  iter2++;
	}
	if(iter2==MAXITER){
	  fprintf(stderr,"%s: max number of iterations (%i) reached in depth loop\n",
		  argv[0],MAXITER);
	  exit(-1);
	}
      }
      // compare with all other faults if the new fault fits
      for(hit=FALSE,j=0;(j<i) && (!hit);j++){
	if(!far_enough((fault+i),(fault+j),leeway)){
	  hit=TRUE;
	  break;
	}
      }
      /*
	now try the subdivision into patches and check for
	evil interactions
      */
      if(!hit){
	// we will use this set of patches if they 
	// are not fatally interacting
	old_nrpatches=nrpatches;
	// create subdivision
	create_patches(i, fault, &patch,&nrpatches,seg,circular,FALSE,td,srand,drand,&seed);
	if(check_for_interactions){
	  hit=incrementally_check_interaction_coefficients(patch,nrpatches,
							   medium,FALSE,evil_pair,
							   max_dist);
	  if(hit){
	    nrpatches=old_nrpatches;
	    interact_prob++;
	  }
	}
      }
      iter ++;
#ifdef DEBUG
      fprintf(stderr,"%s: working on fault %5i, patch: %6i, iter: %5i, ici: %5i, ep: %5i and %5i, dist: %10.4e\r",
	      argv[0],i,nrpatches,iter,interact_prob,evil_pair[0],evil_pair[1],
	      (evil_pair[0]!= -1 && evil_pair[1]!=-1)?
	      (distance_3d(patch[evil_pair[0]].x,patch[evil_pair[1]].x)):(0));
#endif      
    }while((hit) && (iter<MAXITER));
    if(iter==MAXITER){
      fprintf(stderr,"\n%s: fault %i: max number of iterations (%i) exceeded, %i interact problems\n",
	      argv[0],i,MAXITER,interact_prob);
      if(evil_pair[0]!= -1 ){
	printf("%s: last problematic patches: %i and %i\n",
	       argv[0],evil_pair[0],evil_pair[1]);
	print_patch_geometry_and_bc(evil_pair[0],patch,
				    PATCH_OUT_MODE,0.0,FALSE,stderr,FALSE,dummy);
	print_patch_geometry_and_bc(evil_pair[1],patch,
				    PATCH_OUT_MODE,0.0,FALSE,stderr,FALSE,dummy);
      }
      exit(-1);
    }
    titer += iter;
    
  } 
  fprintf(stderr,"\n");
  fprintf(stderr,"%s: total number of randomizing iterations: %i\n",
	  argv[0],(int)titer);
  fprintf(stderr,"%s: total number of actually created patches: %i\n",
	  argv[0],nrpatches);
#ifdef DEBUG
  if(check_for_interactions){
    // check if we really have no positive interaction
    fprintf(stderr,"%s: checking if final set of patches has no feedback loops...\n",
	    argv[0]);
    check_coulomb_stress_feedback(nrpatches,0,patch,medium, 
				  TRUE,CALC_I_COEFF_NOW,TRUE,
				  evil_pair,max_distance);
    fprintf(stderr,"%s: check completed\n",argv[0]);
  }
#endif
  for(i=0;i<nrpatches;i++)
    print_patch_geometry_and_bc(i,patch,opmode,10.0,FALSE,stdout,FALSE,dummy);


  exit(0);
}

