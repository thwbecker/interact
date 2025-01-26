/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: input.c,v 2.38 2001/04/15 00:17:25 becker Exp becker $
*/
#include "properties.h"
#include "interact.h"
#include <string.h>

/* read in the fault geometry */
void read_geometry(char *patch_filename,
		   struct med **medium, struct flt **fault,
		   BOOLEAN read_fault_properties)
{
  int i,j,k,tmpint;
  FILE *in,*in2;
  COMP_PRECISION sin_dip,cos_dip,alpha,mus_avg,mud_avg,d_mu,
    corner[4][3],tmpdbl;
  float wmin,wmax,lmin,lmax;
#ifdef ALLOW_POINT_SOURCE
  int nr_pt_src=0;
#endif
  /* 
     read in fault geometry, documented in help_and_comments.c
  */
  fprintf(stderr,"read_geometry: reading geometry from \"%s\"\n",
	  patch_filename);
  if(read_fault_properties){
    in2=fopen(FAULT_PROP_FILE,"r");
    if(in2)
      fprintf(stderr,"read_geometry: WARNING: reading fault properties from from \"%s\"\n",
	      FAULT_PROP_FILE);
    else{
      fprintf(stderr,"read_geometry: could not open \"%s\" for fault properties, using default values\n",
	      FAULT_PROP_FILE);
      read_fault_properties=FALSE;
    }
  }
  if((*medium=(struct med *)calloc(1,sizeof(struct med)))
     ==NULL)MEMERROR("read_geometry: 1:");
  if(strcmp(patch_filename,"stdin")!=0)
    in=myopen(patch_filename,"r");
  if((*fault=(struct flt *)calloc(1,sizeof(struct flt)))==
     NULL)MEMERROR("read_geometry: 2:");
  for(i=0;i<3;i++){// intialize medium boundaries for plotting
    (*medium)->xmax[i]=FLT_MIN;
    (*medium)->xmin[i]=FLT_MAX;
  } 
  (*medium)->wmean= (*medium)->lmean=0.0;
  lmin=wmin=FLT_MAX;
  lmax=wmax=FLT_MIN;
  i=0;mus_avg=mud_avg=0.0;
 
  while(fscanf(in,PATCH_CP_FORMAT,
	       &(*fault+i)->x[X],
	       &(*fault+i)->x[Y],
	       &(*fault+i)->x[Z],
	       &(*fault+i)->strike,
	       &(*fault+i)->dip,
	       &(*fault+i)->l,&(*fault+i)->w,
	       &tmpint) == 8){
#ifdef ALLOW_POINT_SOURCE
    if(((*fault+i)->w < 0)&&((*fault+i)->l< 0)){
      /*
	
	triangular element
	
      */
      (*fault+i)->type=TRIANGULAR;
      (*fault+i)->xt=(COMP_PRECISION *)
	malloc(sizeof(COMP_PRECISION)*9);
      if(!(*fault+i)->xt)MEMERROR("read_geometry");
      // next nine fields are nodal coordinates
      if(fscanf(in,NINE_CP_FORMAT,
		&(*fault+i)->xt[  X],
		&(*fault+i)->xt[  Y],&(*fault+i)->xt[  Z],
		&(*fault+i)->xt[3+X],
		&(*fault+i)->xt[3+Y],&(*fault+i)->xt[3+Z],
		&(*fault+i)->xt[6+X],&(*fault+i)->xt[6+Y],
		&(*fault+i)->xt[6+Z])!=9)
	READ_ERROR("read_geometry");
      // x will be the centroid, to be calculated
      centroid((*fault+i)->xt,(*fault+i)->x);
      // determine local reference frame by means of the angles
      // w will be the triangle area
      get_alpha_dip_gh((*fault+i)->xt,
		       &(*fault+i)->sin_alpha,
		       &(*fault+i)->cos_alpha,
		       &tmpdbl,&(*fault+i)->w);
      (*fault+i)->dip=(float)tmpdbl;
      (*fault+i)->l=(*fault+i)->w;
      alpha=RAD2DEGF(asin((*fault+i)->sin_alpha));
      (*fault+i)->strike= 90.0 - alpha;
      fprintf(stderr,"read_geometry: fault %i is triangle, x1: (%g, %g, %g) x2: (%g, %g, %g) x3: (%g, %g, %g)\n",
	      i,
	      (*fault+i)->xt[  X],
	      (*fault+i)->xt[  Y],(*fault+i)->xt[  Z],
	      (*fault+i)->xt[3+X],
	      (*fault+i)->xt[3+Y],(*fault+i)->xt[3+Z],
	      (*fault+i)->xt[6+X],
	      (*fault+i)->xt[6+Y],(*fault+i)->xt[6+Z]);
      nr_pt_src++;
    }else if((*fault+i)->l < 0){
      /*

	point source
	
      */
      fprintf(stderr,"read_geometry: fault %i is assumed to be a point source with area %g\n",
	      i,(*fault+i)->w);
      // L and W will hold the `fault' area
      (*fault+i)->l=(*fault+i)->w;
      (*fault+i)->type=POINT_SOURCE;
      nr_pt_src++;
    }else{
      (*fault+i)->type=RECTANGULAR_PATCH; 
    }
#else
    if(((*fault+i)->l <= 0)||((*fault+i)->w <= 0)){
      fprintf(stderr,"read_geometry: fault %i: half length l and width have to be >= 0!\n",i);
      fprintf(stderr,"read_geometry: if point source or triangular elements were\n");
      fprintf(stderr,"what you were looking for,\n");
      fprintf(stderr,"read_geometry: recompile with  ALLOW_POINT_SOURCE flag set\n");
      exit(-1);
    }
#endif
    if((*fault+i)->x[Z] > 0){
      fprintf(stderr,"read_geometry: dfault %i: depth z has to be < 0!\n",i);
      exit(-1);
    }
    if(tmpint < 0){
      fprintf(stderr,"read_geometry: smalles allowed group number is 0, not %i\n",
	      tmpint);
      exit(-1);
    }else{
      (*fault+i)->group=(unsigned int)tmpint;
    }
    // check for illegal angles
    check_fault_angles((*fault+i));
    for(j=0;j<3;j++){
      if((*fault+i)->x[j]<(*medium)->xmin[j])
	(*medium)->xmin[j]=(float)(*fault+i)->x[j];
      if((*fault+i)->x[j]>(*medium)->xmax[j])
	(*medium)->xmax[j]=(float)(*fault+i)->x[j];
    }
    if((*fault+i)->l>lmax)lmax=(*fault+i)->l;
    if((*fault+i)->l<lmin)lmin=(*fault+i)->l;
    if((*fault+i)->w > wmax)wmax=(*fault+i)->w;
    if((*fault+i)->w < wmin)wmin=(*fault+i)->w;
    (*medium)->wmean += (*fault+i)->w;
    (*medium)->lmean += (*fault+i)->l;
    /*
      strike is defined as degrees clockwise from north (azimuth)
      we need the angle alpha, which is 
      counterclockwise from east and used for all rotations 
    */
#ifdef ALLOW_POINT_SOURCE
    if((*fault+i)->type != TRIANGULAR){
#endif
      // if we have different kinds of faults, don't do this
      // calculation if it's a triangular patch since we 
      // have already calculates sin and cos (alpha) above
      alpha= 90.0 - (*fault+i)->strike;
      sincos_deg(&(*fault+i)->sin_alpha,
		 &(*fault+i)->cos_alpha,
		 alpha);
#ifdef ALLOW_POINT_SOURCE
    }
#endif
    sincos_deg(&sin_dip,&cos_dip,
	       (COMP_PRECISION)(*fault+i)->dip);
    // calculate the unity vectors in strike, dip, and normal
    // direction
    calc_base_vecs((*fault+i)->t_strike,(*fault+i)->normal,
		   (*fault+i)->t_dip,
		   (*fault+i)->sin_alpha,
		   (*fault+i)->cos_alpha,
		   sin_dip,cos_dip);
#ifdef SUPER_DEBUG
    fprintf(stderr,"fault %5i: alpha: %10.4e dip: %10.4e sc_alpha:  %10.4e/%10.4e sc_dip: %10.4e/%10.4e\n",
	    i,alpha,(*fault+i)->dip,(*fault+i)->sin_alpha,(*fault+i)->cos_alpha,
	    sin_dip,cos_dip);
    fprintf(stderr," vec: s: (%10.4e,%10.4e,%10.4e) d: (%10.4e,%10.4e,%10.4e) n: (%10.4e,%10.4e,%10.4e)\n",
	    (*fault+i)->t_strike[X],(*fault+i)->t_strike[Y],(*fault+i)->t_strike[Z],
	    (*fault+i)->t_dip[X],(*fault+i)->t_dip[Y],(*fault+i)->t_dip[Z],
	    (*fault+i)->normal[X],(*fault+i)->normal[Y],(*fault+i)->normal[Z]);
#endif
#ifdef ALLOW_POINT_SOURCE
    if((*fault+i)->type != TRIANGULAR){
#endif
      // check depth alignment
      if((*fault+i)->w*sin_dip+(*fault+i)->x[Z]>EPS_COMP_PREC){
	fprintf(stderr,"read_geometry: WARNING: fault %i above surface, z: %20.10e\n",
		i,(*fault+i)->w*sin_dip+(*fault+i)->x[Z]);
      }
      // determine geometrical boundaries for plotting
      calculate_corners(corner,(*fault+i));
      for(j=0;j<4;j++)
	for(k=0;k<3;k++){
	  if((*medium)->xmax[k]<corner[j][k])
	    (*medium)->xmax[k]=corner[j][k];
	  if((*medium)->xmin[k]>corner[j][k])
	    (*medium)->xmin[k]=corner[j][k];
	}
#ifdef ALLOW_POINT_SOURCE
    }
#endif
    if(read_fault_properties){
      /* 
	 frictional properties, static and dynamic 
	 either read from file (-f switch)
      */
      if(fscanf(in2,"%f %f",&(*fault+i)->mu_s,
		&(*fault+i)->mu_d)!=2){
	fprintf(stderr,"read_geometry: read error: properties file: %s\n",
		FAULT_PROP_FILE);
	fprintf(stderr,"read_geometry: could not read two parameters for fault %i\n",
		i);
	exit(-1);
      }
      d_mu=(*fault+i)->mu_s - (*fault+i)->mu_d;
      if(d_mu>1 || d_mu<0){
	fprintf(stderr,"read_geometry: WARNING: Delta mu is %g for fault %i\n",
		d_mu,i);
      }
    }else{
      // or use constant values
      (*fault+i)->mu_s=(float)STATIC_MU;
      (*fault+i)->mu_d=(float)(STATIC_MU-DELTA_MU);
    }
    mus_avg += (COMP_PRECISION)(*fault+i)->mu_s;
    mud_avg += (COMP_PRECISION)(*fault+i)->mu_d;
    /* initialize some of the arrays with zeros */
    (*fault+i)->u[X]=(*fault+i)->u[Y]=(*fault+i)->u[Z]=0.0;
    if(((*fault+i)->mode=(MODE_TYPE *)malloc(sizeof(MODE_TYPE)))
       ==NULL)MEMERROR("read_geometry: 3:");
    (*fault+i)->mode[0]=INACTIVE;
    /* add one more fault to the list */
    i++;
    if((*fault=(struct flt *)
	realloc(*fault,sizeof(struct flt)*(i+1)))==NULL)
      MEMERROR("read_geometry: 4:");
  }
  /* 
     end fault input loop
  */
  if(strcmp(patch_filename,"stdin")!=0)
    fclose(in);
  if(read_fault_properties)
    fclose(in2);
  // now we have the number of faults
  (*medium)->nrflt = i;
  if(!(*medium)->nrflt){
    fprintf(stderr,"read_geometry: did not read in any faults from geometry file\n");
    exit(-1);
  }
  if((*medium)->nrflt>99998)
    fprintf(stderr,"read_geometry: WARNING: parts of the program I/O might not work for more patches than 99999!\n");
  *fault=(struct flt *)
    realloc(*fault,sizeof(struct flt)*(*medium)->nrflt);
  //
  // determine the number of groups that were assigned
  //
  for((*medium)->nrgrp=i=0;i<(*medium)->nrflt;i++){
    if((*fault+i)->group > (*medium)->nrflt-1){
      fprintf(stderr,"read_geometry: fault %i has group %i, should be between 0 and %i (nr of patches - 1)\n",
	      i,(*fault+i)->group,(*medium)->nrflt-1);
      exit(-1);
    }
    if((*fault+i)->group+1 > (*medium)->nrgrp)
      (*medium)->nrgrp=(*fault+i)->group+1;
  }
  fprintf(stderr,"read_geometry: working with %i group(s) of patches\n",(*medium)->nrgrp);
  // allocate space for the fault group arrays
  if(((*medium)->fault_group=
      (struct fgrp *)calloc((*medium)->nrgrp,
			   sizeof(struct fgrp)))==NULL)
    MEMERROR("read_geometry: 5:");
  // calculate average friction values and geometries
  mus_avg /= (COMP_PRECISION)(*medium)->nrflt;
  mud_avg /= (COMP_PRECISION)(*medium)->nrflt;
  (*medium)->wmean/=(COMP_PRECISION)(*medium)->nrflt;
  (*medium)->lmean/=(COMP_PRECISION)(*medium)->nrflt;
  
  if((*medium)->nrflt<=(*medium)->max_nr_flt_files)
    for(i=0;i<(*medium)->nrflt;i++)
      fprintf(stderr,"read_geometry: fault %4i: x: (%9.6g, %9.6g, %9.6g) strike: %6.2f dip: %6.2f half_length: %6.3g aspect: %9.6g\n",
	      i,(*fault+i)->x[X],(*fault+i)->x[Y],
	      (*fault+i)->x[Z],(*fault+i)->strike,
	      (*fault+i)->dip,(*fault+i)->l,
	      (*fault+i)->l/(*fault+i)->w);
  fprintf(stderr,"read_geometry: read in %i fault patch(es)\nread_geometry: half length: min/mean/max: %g/%g/%g\nread_geometry: half width:  min/mean/max: %g/%g/%g\n",
	  (*medium)->nrflt,
	  lmin,(*medium)->lmean,lmax,
	  wmin,(*medium)->wmean,wmax);
  fprintf(stderr,"read_geometry: average values for mu_s/mu_d: %g/%g\n",
	  mus_avg,mud_avg);
  //
  // adjust min/max of medium (for plotting only)
  //
  for(i=0;i<3;i++){
    (*medium)->xmin[i]=
      ((*medium)->xmin[i] < -MIN_GEOM_RANGE) ? 
      (*medium)->xmin[i] : -MIN_GEOM_RANGE;
    (*medium)->xmax[i]=
      (*medium)->xmax[i] >  MIN_GEOM_RANGE ?  
      (*medium)->xmax[i] : MIN_GEOM_RANGE;
  }
#ifdef ALLOW_POINT_SOURCE
  if(nr_pt_src == 0){
    fprintf(stderr,"read_geometry: no point source \"faults\" were read in, recompile without\n");
    fprintf(stderr,"read_geometry: ALLOW_POINT_SOURCE flag set for improved speed\n");
  }
#endif
  (*medium)->geometry_init=TRUE;
}


/* 
   read in boundary conditions for faults 
*/
void read_boundary_conditions(struct med *medium,
			      struct flt *fault)
{
  FILE *in;
  int flt_nr,bc_code,n=0,inc,start_flt,stop_flt,i,j;
  BOOLEAN printevery;
  in=fopen(BC_FILE,"r");
  if(in != NULL){
    fprintf(stderr,"read_boundary_conditions: reading boundary conditions\n");
    fscanf(in,"%i",&medium->op_mode);
    switch(medium->op_mode){
    case ONE_STEP_CALCULATION:{
      /*
	resize the fault mode array to allow 
	for three different activations
      */
      medium->nr_flt_mode=3;
      for(i=0;i<medium->nrflt;i++){
	if((fault[i].mode=(MODE_TYPE *)
	    realloc(fault[i].mode,medium->nr_flt_mode*sizeof(MODE_TYPE)))
	   ==NULL)MEMERROR("read_boundary_conditions: 1:");
	for(j=0;j<medium->nr_flt_mode;j++)
	  fault[i].mode[j]=INACTIVE;
      }
      read_one_step_bc(in,medium,fault);
      break;
    }
    /*

      SIMULATE LOADING EXPERIMENT

    */
    case SIMULATE_LOADING:
    case SIMULATE_LOADING_AND_PLOT:{
      medium->nr_flt_mode=1;
      fprintf(stderr,"read_boundary_conditions: simulating loading\n");
      /* timing issues */
      fscanf(in,THREE_CP_FORMAT,
	     &medium->dt0,&medium->print_interval,
	     &medium->stop_time);
#ifndef VARIABLE_TIME_STEP
      fprintf(stderr,"read_boundary_conditions: until time %g in fixed steps of %g\n",
	      medium->stop_time,medium->dt0);
#else
      fprintf(stderr,"read_boundary_conditions: until time %g in variable time steps\n",
	      medium->stop_time);
      fprintf(stderr,"read_boundary_conditions: the timestep is calculated from MIN(time_to_failure, print_interval)\n");
      /*
	importantly, we will set the base time step to 
	unity to facilitate computations of the adequate 
	timestep based on the standard increment
      */
      medium->dt0=1.0;
#endif
      if(medium->dt0 <= 0){
	fprintf(stderr,"read_boundary_conditions: dt0 has to be larger than zero, dt0: %g\n",
		medium->dt0);
	exit(-1);
      }
      medium->dt=medium->dt0;

      if(medium->op_mode == SIMULATE_LOADING_AND_PLOT){
#ifdef USE_PGPLOT
	/* plotting interval */
	fscanf(in,TWO_CP_FORMAT,
	       &medium->x_plot_interval,
	       &medium->x_scroll_interval);
	fprintf(stderr,"read_boundary_conditions: x plotting every %g on window 1, scroll every %g\n",
		medium->x_plot_interval,
		medium->x_scroll_interval);
	medium->x_scroll_inc=
	  medium->x_scroll_interval;
#else
	fprintf(stderr,"read_boundary_conditions: you requested X output but program was compiled without PGPLOT support\n");
	fscanf(in,"%*f %*f");
#endif
      }
      /* 
	 read in fault slip mode codes 
	 for loading simulation
      */
      n=0;
      while(fscanf(in,"%i %i",&flt_nr,&bc_code) == 2){
	n++;
	if(flt_nr < 0){// we will set several faults to same boundary code
	  start_flt=bc_code;
	  fscanf(in,"%i %i",&stop_flt,&bc_code);
	  if(stop_flt < 0){
	    start_flt= 0;
	    stop_flt=medium->nrflt-1;
	  }
	  inc=-flt_nr;
	  fprintf(stderr,"read_boundary_conditions: const. %5i: flts %5i to %5i: activation mode %3i, %s\n",
		  n,start_flt,stop_flt,bc_code,comment_on_code(bc_code));
	  printevery=FALSE;
	}else{// single fault mode
	  start_flt=flt_nr;
	  stop_flt=flt_nr;
	  inc=1;
	  printevery=TRUE;
	}
	if(abs(bc_code)>=SHRT_MAX){
	  fprintf(stderr,"read_boundary_conditions: bc code %i won't fit in a short int variable\n",
		  bc_code);
	  exit(-1);
	}
	for(flt_nr=start_flt;flt_nr<=stop_flt;flt_nr += inc){
	  if(flt_nr > (medium->nrflt-1) || flt_nr<0){
	    fprintf(stderr,"read_boundary_conditions: file %s, line %i, fault number %i out of range (max nr flt: %i)\n",
		    BC_FILE,n,flt_nr,medium->nrflt);
	    exit(-1);
	  }
	  switch(bc_code){
	  case INACTIVE:
	  case STRIKE_SLIP:
	  case STRIKE_SLIP_LEFTLATERAL:
	  case STRIKE_SLIP_RIGHTLATERAL:
	  case DIP_SLIP:
	  case DIP_SLIP_UPWARD:
	  case DIP_SLIP_DOWNWARD:
	  case COULOMB_STRIKE_SLIP:
	  case COULOMB_STRIKE_SLIP_LEFTLATERAL:
	  case COULOMB_STRIKE_SLIP_RIGHTLATERAL:
	  case COULOMB_DIP_SLIP:
	  case COULOMB_DIP_SLIP_UPWARD:
	  case COULOMB_DIP_SLIP_DOWNWARD:
#ifndef NO_OPENING_MODES
	  case NORMAL_SLIP:
	  case NORMAL_SLIP_OUTWARD:
	  case NORMAL_SLIP_INWARD:
#endif
	  case COULOMB_MAXSDIR_SLIP:
	  case MAXSDIR_SLIP:{
	    if(printevery)
	      fprintf(stderr,"read_boundary_conditions: const. %5i fault %5i activation mode %3i, %s\n",
		      n,flt_nr,bc_code,comment_on_code(bc_code));
	    if(fault[flt_nr].mode[0] != INACTIVE){
	      fprintf(stderr,"read_boundary_conditions: fault %i has already been assigned mode %i\n",
		      flt_nr,fault[flt_nr].mode[0]);
	      fprintf(stderr,"read_boundary_conditions: loading simulation run can only cope with one mode\n");
	      exit(-1);
	    }
	    fault[flt_nr].mode[0]=(MODE_TYPE)bc_code;
#ifdef SUPER_DEBUG
	    fprintf(stderr,"read_boundary_conditions: fault %5i: code: %i\n",
		    flt_nr,fault[flt_nr].mode[0]);
#endif
	    break;
	  }
#ifdef NO_OPENING_MODES
	  case NORMAL_SLIP:
	  case NORMAL_SLIP_OUTWARD:
	  case NORMAL_SLIP_INWARD:{
	    fprintf(stderr,"read_boundary_conditions: normal (opening) slip modes have been switched off, recompile\n");
	    fprintf(stderr,"read_boundary_conditions: without NO_OPENING_MODES switch\n");
	    exit(-1);
	    break;
	  }
#endif
	  default:{
	    fprintf(stderr,"read_boundary_conditions: fault activation mode %i undefined\n",
		    bc_code);
	    exit(-1);
	    break;
	  }}
	}
      }
      break;
    }
    default:{
      fprintf(stderr,"read_boundary_conditions: operational mode %i unddefined\n",medium->op_mode);
      exit(-1);
      break;
    }}
    fclose(in);
  }else{
    fprintf(stderr,"read_boundary_conditions: no bc's specified\n");
  }
  medium->bc_init=TRUE;
}     



void read_one_step_bc(FILE *in,struct med *medium,
		      struct flt *fault)
{
  /* 
     read in boundary conditions for one step calculation 
     mode
  */
  BOOLEAN *sma,printevery,added_to_con,added_to_uncon;
  int flt_nr,bc_code,n,inc,start_flt,stop_flt,i,j,nrflt3;
  COMP_PRECISION bc_value,*rhs;
  

  if(medium->nr_flt_mode != 3){
    fprintf(stderr,"read_one_step_bc: have to allow for %i modes, only %i set\n",
	    3,medium->nr_flt_mode);
    exit(-1);
  }
  fscanf(in,"%i",&i);
  medium->print_bulk_fields=(BOOLEAN)i;
  if(medium->print_bulk_fields){
    fprintf(stderr,"read_boundary_conditions: will print bulk fields\n");
    /* 
       geometry for field outputs 
    */
    fscanf(in,"%f %f %i %f %f %i %f %f %i",
	   &medium->pxmin[X],&medium->pxmax[X],
	   &medium->n[X],
	   &medium->pxmin[Y],&medium->pxmax[Y],
	   &medium->n[Y],
	   &medium->pxmin[Z],&medium->pxmax[Z],
	   &medium->n[Z]);
    if(medium->n[Z]>0){
      fprintf(stderr,"read_boundary_conditions: xmin: %g xmax: %g n: %i dx: %g\n",
	      medium->pxmin[X],medium->pxmax[X],
	      medium->n[X],medium->n[X]!=1?(medium->pxmax[X]-medium->pxmin[X])/
	      ((COMP_PRECISION)(medium->n[X]-1)):0);
      fprintf(stderr,"read_boundary_conditions: ymin: %g ymax: %g m: %i dy: %g\n",
	      medium->pxmin[Y],medium->pxmax[Y],
	      medium->n[Y],medium->n[Y]!=1?(medium->pxmax[Y]-medium->pxmin[Y])/
	      ((COMP_PRECISION)(medium->n[Y]-1)):0.0);
      fprintf(stderr,"read_boundary_conditions: zmin: %g zmax: %g o: %i dz: %g\n",
	      medium->pxmin[Z],medium->pxmax[Z],
	      medium->n[Z],medium->n[Z]!=1?(medium->pxmax[Z]-medium->pxmin[Z])/
	      ((COMP_PRECISION)(medium->n[Z]-1)):0.0);
    }else{
      // create a field that is in the plane of the average strike and dip directions 
      // of fault group 0
      fprintf(stderr,"read_boundary_conditions: uses field in line with the plane of fault group 0\n");
      fprintf(stderr,"read_boundary_conditions: strike extent min: %g max: %g n: %i d_strike: %g\n",
	      medium->pxmin[X],medium->pxmax[X],
	      medium->n[X],medium->n[X]!=1?(medium->pxmax[X]-medium->pxmin[X])/
	      ((COMP_PRECISION)(medium->n[X]-1)):0);
      fprintf(stderr,"read_boundary_conditions: dip extent:   min: %g max: %g m: %i d_dip:    %g\n",
	      medium->pxmin[Y],medium->pxmax[Y],
	      medium->n[Y],medium->n[Y]!=1?(medium->pxmax[Y]-medium->pxmin[Y])/
	      ((COMP_PRECISION)(medium->n[Y]-1)):0.0);
    }
  }else{
    fprintf(stderr,"read_boundary_conditions: will not print bulk fields\n");
  }
  fprintf(stderr,"read_boundary_conditions: running one step calculation\n");
  /* 
     allocate temporary arrays for all faults, will be deleted further down 
     if these get too big, have to do more bookkeeping,
     now, for simplicity keep like that (otherwise fault constraints
     would have to be sorted)
  */
  nrflt3=medium->nrflt*3;
  if((sma=(BOOLEAN *)
      malloc(nrflt3*sizeof(BOOLEAN))) == NULL)
    MEMERROR("read_boundary_conditions: 2:");
  for(i=0;i<nrflt3;i++)sma[i]=INACTIVE;
  if((rhs=(COMP_PRECISION *)
      calloc(nrflt3,sizeof(COMP_PRECISION)))
     ==NULL)MEMERROR("read_boundary_conditions: 3:");
  n=0;
  /* read in contraints */
  while(fscanf(in,IIF_CP_FORMAT,&flt_nr,&bc_code,
	       &bc_value) == 3){
    n++;
    if(flt_nr < 0){/* the fault indicator was negative,
		      meaning we will assign the next 
		      values to a sequence of faults */
      if((bc_code < 0) && (bc_value <0)){
	start_flt=0;
	stop_flt=medium->nrflt-1;
      }else{
	start_flt=bc_code;
	stop_flt=(int)bc_value;
      }
      inc=-flt_nr;
      fscanf(in,IF_CP_FORMAT,&bc_code,&bc_value);
      printevery=FALSE;
      fprintf(stderr,"read_boundary_conditions: const: %5i, flts %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
	      n,start_flt,stop_flt,bc_code,bc_value,comment_on_code_bc(bc_code,bc_value));
    } else {
      printevery=TRUE;
      start_flt=flt_nr;
      stop_flt=flt_nr;
      inc=1;
    }
    /*
      now assign to temporary arrays since we have sorted 
      by faults and added multiple
      activations
    */ 
    for(flt_nr=start_flt;flt_nr <= stop_flt;flt_nr+=inc){
      if(flt_nr > (medium->nrflt-1)){
	fprintf(stderr,"read_boundary_conditions: file %s, line %i, fault number %i out of range (max nr flt: %i)\n",
		BC_FILE,n,flt_nr,medium->nrflt);
	exit(-1);
      }
      switch(bc_code){
	/* displacements */
      case STRIKE:
      case DIP:
      case NORMAL:{/* can be added to fault values */
	fault[flt_nr].u[bc_code] += bc_value;
	if(printevery)
	  fprintf(stderr,"read_boundary_conditions: const: %5i, flts %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
		  n,start_flt,stop_flt,bc_code,bc_value,comment_on_code_bc(bc_code,bc_value));
	break;
      }
      /* 
	 stress boundary conditions will be assembled in 
	 a huge array first, then assigned to rhs. later
	 so we can keep up with a general ordering scheme 
	 as in rupture.c 
	 
	 we will also have to assign fault slip modes
      */
      case COULOMB_STRIKE_SLIP_LEFTLATERAL:
      case COULOMB_STRIKE_SLIP_RIGHTLATERAL:
      case COULOMB_STRIKE_SLIP:
      case STRIKE_SLIP_LEFTLATERAL:
      case STRIKE_SLIP_RIGHTLATERAL:
      case STRIKE_SLIP:{
	if(sma[flt_nr*3+STRIKE]){
	  fprintf(stderr,"read_boundary_conditions: stress component %i on fault %i was already constrained\n",
		  bc_code,flt_nr);
	  exit(-1);
	} 
	rhs[flt_nr*3+STRIKE]=bc_value;
	if(bc_code>OS_C_OFFSET){// need to correct for normal stress changes
	  if(fabs(rhs[flt_nr*3+STRIKE])<EPS_COMP_PREC)
	    sma[flt_nr*3+STRIKE]=ACTIVATED;
	  else if(rhs[flt_nr*3+STRIKE]>0)
	    sma[flt_nr*3+STRIKE]=
	      ACTIVATED_AND_POSITIVE_NORMAL_CORRECTION;
	  else
	    sma[flt_nr*3+STRIKE]=
	      ACTIVATED_AND_NEGATIVE_NORMAL_CORRECTION;
	}else{
	  sma[flt_nr*3+STRIKE]=ACTIVATED;
	}
	fault[flt_nr].mode[STRIKE]=(MODE_TYPE)bc_code;
	if(printevery)
	  fprintf(stderr,"read_boundary_conditions: const: %5i, flts %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
		 n,start_flt,stop_flt,bc_code,bc_value,comment_on_code_bc(bc_code,bc_value));
	break;
      }
      case COULOMB_DIP_SLIP_UPWARD:
      case COULOMB_DIP_SLIP_DOWNWARD:
      case COULOMB_DIP_SLIP:
      case DIP_SLIP_UPWARD:
      case DIP_SLIP_DOWNWARD:
      case DIP_SLIP:{
	if(sma[flt_nr*3+DIP]){
	  fprintf(stderr,"read_boundary_conditions: stress component %i on fault %i was already constrained\n",
		  bc_code-3,flt_nr);
	  exit(-1);
	}
	rhs[flt_nr*3+DIP]=bc_value;
	if(bc_code > OS_C_OFFSET){// need to correct for normal stress changes
	  if(fabs(rhs[flt_nr*3+DIP])<EPS_COMP_PREC)
	    sma[flt_nr*3+DIP]=ACTIVATED;
	  else if(rhs[flt_nr*3+DIP]>0)
	    sma[flt_nr*3+DIP]=
	      ACTIVATED_AND_POSITIVE_NORMAL_CORRECTION;
	  else
	    sma[flt_nr*3+DIP]=
	      ACTIVATED_AND_NEGATIVE_NORMAL_CORRECTION;
	}else{
	  sma[flt_nr*3+DIP]=ACTIVATED;
	}
	fault[flt_nr].mode[DIP]=(MODE_TYPE)bc_code;
	if(printevery)
	  fprintf(stderr,"read_boundary_conditions: const: %5i, flts %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
		  n,start_flt,stop_flt,bc_code,bc_value,comment_on_code_bc(bc_code,bc_value));
	break;
      }
      case NORMAL_SLIP_OUTWARD:
      case NORMAL_SLIP_INWARD:
      case NORMAL_SLIP:{
#ifdef NO_OPENING_MODES
	fprintf(stderr,"read_boundary_conditions: normal mode slipping off, how did we get this far?\n");
	exit(-1);
#else
	if(sma[flt_nr*3+NORMAL]){
	  fprintf(stderr,"read_boundary_conditions: stress component %i on fault %i was already constrained\n",
		  bc_code-3,flt_nr);
	  exit(-1);
	}
	sma[flt_nr*3+NORMAL]=ACTIVATED;
	rhs[flt_nr*3+NORMAL]=bc_value;
	fault[flt_nr].mode[NORMAL]=(MODE_TYPE)bc_code;
	if(printevery)
	  fprintf(stderr,"read_boundary_conditions: const: %5i, flts %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
		 n,start_flt,stop_flt,bc_code,bc_value,comment_on_code_bc(bc_code,bc_value));
#endif
	break;
      }
      default:{
	fprintf(stderr,"read_boundary_conditions: file %s, line %i, bc code %i undefined \n",
		BC_FILE,n,bc_code);
	exit(-1);
      }}
    }
  }
  /* 
     ASSIGN STRESS CONDITIONS TO MATRIX in the same 
     fashion as in rupture.c, first initialize the 
     real arrays
  */
  init_equation_system(medium,fault);
  /*
    
    ASSIGN BOUNDARY CONDITIONS and decide if we want a 
    constrained or unconstrained (positive only or positive and negative) 
    solution for slip
    
  */
  for(i=0;i<medium->nrflt;i++){
    added_to_con=  FALSE;
    added_to_uncon=FALSE;
    for(j=0;j<3;j++){// loop through possible modes
      switch(fault[i].mode[j]){
	/* 
	   first the unconstrained modes, 
	   slip can go either way
	*/
      case COULOMB_STRIKE_SLIP:
      case COULOMB_DIP_SLIP:
      case STRIKE_SLIP:
      case DIP_SLIP:
      case NORMAL_SLIP:{
	/* 
	   add a stress equation for each active 
	   slipping mode
	*/
	if(sma[i*3+j]){
	  medium->sma[medium->naflt*3+j]=sma[i*3+j];
	  add_to_right_hand_side(rhs[i*3+j],
				 &medium->b,
				 &medium->nreq);
	  added_to_uncon=TRUE;/*
				have to increment fault 
				contraint counter
			      */
	}
	break;
      }
      /* 
	 now CONSTRAINED equations, slip of specific direction 
	 was specified 
      */
      default:{
	if(sma[i*3+j]){
	  medium->sma_con[medium->naflt_con*3+j]=sma[i*3+j];
	  /* 

	     check for strange settings that would lead to funny 
	     stuff  without interaction-modified stresses 

	  */
	  if((((fault[i].mode[j]==STRIKE_SLIP_LEFTLATERAL) ||  
	       (fault[i].mode[j]==COULOMB_STRIKE_SLIP_LEFTLATERAL))
	      && (rhs[i*3+j]>0))||
	     (((fault[i].mode[j]==STRIKE_SLIP_RIGHTLATERAL) || 
	       (fault[i].mode[j]==COULOMB_STRIKE_SLIP_RIGHTLATERAL))
	      && (rhs[i*3+j]<0))||
	     (((fault[i].mode[j]==DIP_SLIP_UPWARD ) ||         
	       (fault[i].mode[j]==COULOMB_DIP_SLIP_UPWARD))
	      && (rhs[i*3+j]>0))||
	     (((fault[i].mode[j]==DIP_SLIP_DOWNWARD) ||
	       (fault[i].mode[j]==COULOMB_DIP_SLIP_DOWNWARD))
	      && (rhs[i*3+j]<0))||
	     ((fault[i].mode[j] == NORMAL_SLIP_OUTWARD)      && (rhs[i*3+j] < 0))||
	     ((fault[i].mode[j] == NORMAL_SLIP_INWARD)       && (rhs[i*3+j] > 0))){
	    fprintf(stderr,
		    "read_boundary_conditions: WARNING: fault: %5i bc_code: %3i dir: %1i might be incompatible with stress: %g\n",
		    i,fault[i].mode[j],j,rhs[i*3+j]);
	  }
	  if(fault[i].mode[j]==STRIKE_SLIP_RIGHTLATERAL || 
	     fault[i].mode[j]==COULOMB_STRIKE_SLIP_RIGHTLATERAL ||
	     fault[i].mode[j]==DIP_SLIP_DOWNWARD || 
	     fault[i].mode[j]==COULOMB_DIP_SLIP_DOWNWARD  ||
	     fault[i].mode[j]==NORMAL_SLIP_INWARD)
	    /* 
	       assign negative b for constrained, 'anormal' modes 
	    */
	    add_to_right_hand_side(-rhs[i*3+j],
				   &medium->b_con,
				   &medium->nreq_con);
	  else
	    add_to_right_hand_side(rhs[i*3+j],
				   &medium->b_con,
				   &medium->nreq_con);
	  // have to increment counter
	  added_to_con=TRUE;
	}
	break;
      }}
    }
#ifdef SUPER_DEBUG
    fprintf(stderr,"read_boundary_condition: fault %i: codes: %i %i %i\n",
	    i,sma[i*3],sma[i*3+1],sma[i*3+2]);
    fprintf(stderr,"read_boundary_condition: uc flt %5i medium codes: %i %i %i\n",
	    medium->naflt,
	    medium->sma[medium->naflt*3],
	    medium->sma[medium->naflt*3+1],
	    medium->sma[medium->naflt*3+2]);
    fprintf(stderr,"read_boundary_condition:  c flt %5i medium codes: %i %i %i\n",
	    medium->naflt_con,
	    medium->sma_con[medium->naflt_con*3],
	    medium->sma_con[medium->naflt_con*3+1],
	    medium->sma_con[medium->naflt_con*3+2]);
    
#endif    
    // step up the active fault counters and sma arrays
    if(added_to_uncon)
      add_to_active_fault_list(i,
			       &medium->nameaf,
			       &medium->naflt,
			       &medium->sma);
    if(added_to_con)
      add_to_active_fault_list(i,
			       &medium->nameaf_con,
			       &medium->naflt_con,
			       &medium->sma_con);
  }
  if(medium->nreq_con + medium->nreq > 0){
    fprintf(stderr,
	    "read_boundary_condition: total of %i constrained and %i unconstrained equations\n",
	    medium->nreq_con,medium->nreq);
    fprintf(stderr,
	    "read_boundary_condition: with     %i constrained and %i unconstrained faults\n",
	    medium->naflt_con,medium->naflt);
  }
  free(rhs);free(sma);
}
/*
  
  assign the stress matrix loading factors 
  s[i]=a[i]+b[i]*time

*/
void read_stress_fac(BOOLEAN read_stress_relation_factors,
		     COMP_PRECISION *a,COMP_PRECISION *b)
{
  int i;
  FILE *in;
  for(i=0;i<6;i++){
    a[i]=0.0;
    b[i]=0.0;
  }
  // default is simple shear, \sigma_xy
  b[1]=STRESSING_RATE;
  if(!read_stress_relation_factors)
     fprintf(stderr,"read_stress_fac: using default simple shear stress relation for loading\n");
  else if(!(in=fopen(STRESS_RELATION_FILE,"r"))){
    fprintf(stderr,"read_stress_fac: can not open \"%s\"\n",STRESS_RELATION_FILE);
    fprintf(stderr,"read_stress_fac: using default simple shear stress relation for loading\n");
  }else{
    
    for(i=0;i<6;i++)
      if(fscanf(in,TWO_CP_FORMAT,(a+i),(b+i))!=2){
	fprintf(stderr,"read_stress_fac: error reading, need a_i b_i, 0<=i<=6\n");
	exit(-1);
      }
    fclose(in);
    fprintf(stderr,"read_stress_fac: read stress relation from \"%s\"\n",
	    STRESS_RELATION_FILE);
  }
}

void init_parameters(char **argv, int argc, 
		     BOOLEAN *read_fault_properties,
		     BOOLEAN *read_initial_fault_stress,
		     BOOLEAN *suppress_interactions,
		     BOOLEAN *whole_fault_mode,
		     COMP_PRECISION *med_cohesion,
		     int *max_nr_flt_files,
		     BOOLEAN *read_stress_relation_factors,
		     BOOLEAN *use_slip_files,
		     BOOLEAN *whole_fault_deactivations,
		     COMP_PRECISION *min_stress_drop,
		     BOOLEAN *use_sparse_storage,
		     I_MATRIX_PREC *i_mat_cutoff,
		     BOOLEAN *use_old_imat,
		     BOOLEAN *save_imat,
		     BOOLEAN *check_for_interaction_feedback)
{
  int i;
  /* 
     assign default values 
  */
  *read_fault_properties=READ_FAULT_PROPERTIES_DEF;
  *read_initial_fault_stress=READ_INITIAL_FAULT_STRESS_DEF;
  *suppress_interactions=SUPPRESS_INTERACTIONS_DEF;
  *whole_fault_mode=WHOLE_FAULT_MODE_DEF;
  *med_cohesion=COHESION_DEF;
  *max_nr_flt_files=MAX_NR_FLT_FILES_DEF;
  *read_stress_relation_factors=READ_STRESS_RELATION_FACTORS_DEF;
  *use_slip_files=PRINT_SLIPLINE_DEF;
  *whole_fault_deactivations=WHOLE_FAULT_DEACTIVATIONS_DEF;
  *use_sparse_storage=USE_SPARSE_STORAGE_DEF;
  *min_stress_drop=0.0;
  *i_mat_cutoff=I_MAT_CUTOFF_DEF;
  *use_old_imat=USE_OLD_IMAT_DEF;
  *save_imat=SAVE_IMAT_DEF;
  *check_for_interaction_feedback=CHECK_FOR_INTERACTION_FEEDBACK_DEF;
  /* 
     check for input options 
  */
  for(i=1;i<argc;i++){
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){// help
      phelp();
      exit(-1);
    }else if(strcmp(argv[i],"-f")==0){// fault prop file
      toggle(read_fault_properties);
    }else if(strcmp(argv[i],"-s")==0){// fault stress init file
      toggle(read_initial_fault_stress);
    }else if(strcmp(argv[i],"-b")==0){// no individual fault group files
      *max_nr_flt_files=0;
    }else if(strcmp(argv[i],"-i")==0){// no fault interactions
      toggle(suppress_interactions);
    }else if(strcmp(argv[i],"-w")==0){// whole fault activation mode
      toggle(whole_fault_mode);
    }else if(strcmp(argv[i],"-nd")==0){// sparse storage mode
      toggle(use_sparse_storage);
    }else if(strcmp(argv[i],"-si")==0){// save I matrix on file?
      toggle(save_imat);
    }else if(strcmp(argv[i],"-ci")==0){// check for feedback interaction?
      toggle(check_for_interaction_feedback);
    }else if(strcmp(argv[i],"-v")==0){// whole fault de-activation mode
      toggle(whole_fault_deactivations);
    }else if(strcmp(argv[i],"-l")==0){// read in a/b factors for stress relation
      toggle(read_stress_relation_factors);
    }else if(strcmp(argv[i],"-p")==0){// print to the slipline files
      toggle(use_slip_files);
    }else if(strcmp(argv[i],"-c")==0){// cohesion
      i++;
      sscanf(argv[i],ONE_CP_FORMAT,med_cohesion);
#ifdef NO_COHESION
      fprintf(stderr,"init_parameters program was compiled with NO_COHESION option\n");
      fprintf(stderr,"init_parameters therefore, cohesion is always zero and -c does not make sense\n");
      exit(-1);
#endif
    }else if(strcmp(argv[i],"-ms")==0){// minimum stress drop
      i++;
      sscanf(argv[i],ONE_CP_FORMAT,min_stress_drop);
    }else if(strcmp(argv[i],"-oi")==0){// use old i matrix
      toggle(use_old_imat);
    }else if(strcmp(argv[i],"-ei")==0){// i matrix cutoff value
      i++;
      sscanf(argv[i],ONE_IP_FORMAT,i_mat_cutoff);
    }else{
      fprintf(stderr,"init_parameters can not use parameter %s, use -h for help page\n",
	      argv[i]);
      exit(-1);
    }
  }
}

// deal with boolean values/switches

char *name_boolean(BOOLEAN value)
{
  if(value)
    return("ON");
  else
    return("OFF");
}

BOOLEAN toggle(BOOLEAN *variable)
{
  if(*variable){
    *variable=FALSE;
    return(FALSE);
  }else{
    *variable=TRUE;
    return(TRUE);
  }
}


