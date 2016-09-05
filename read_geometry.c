/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: read_geometry.c,v 1.28 2003/03/02 07:33:49 becker Exp $
*/
#include "properties.h"
#include "interact.h"
#include <string.h>

/* 

   read in the fault patch geometries

   initialize the fault and medium structures

*/

void read_geometry(char *patch_filename,struct med **medium, 
		   struct flt **fault,
		   my_boolean read_fault_properties,
		   my_boolean twod_approx_is_plane_stress,
		   my_boolean half_plane,
		   my_boolean verbose)
{
  int i,j,k,tmpint,jlim,ic;
  FILE *in,*in2=NULL;
  COMP_PRECISION sin_dip,cos_dip,alpha,mus_avg,mud_avg,d_mu,corner[4][3],
    lloc,wloc,eps_for_z = EPS_COMP_PREC * 100.0;
  float wmin,wmax,lmin,lmax;
  static my_boolean init=FALSE;
#ifdef ALLOW_NON_3DQUAD_GEOM
  int nr_pt_src=0,nr_triangle=0,nr_2d=0;
  COMP_PRECISION tmpdbl;
#endif
  /* 

     READ IN FAULT GEOMETRY in the patch format as 
     documented in help_and_comments.c

  */ 
  if(strcmp(patch_filename,"stdin")!=0){
    if(verbose)
      fprintf(stderr,"read_geometry: reading geometry from file \"%s\"\n",
	      patch_filename);
    in=myopen(patch_filename,"r");
  }else{
    if(verbose)
      fprintf(stderr,"read_geometry: reading geometry from stdin\n");
    in=stdin;
  }
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
  /*
    
    initialize the medium structure


  */
  if(init){
    fprintf(stderr,"read_geometry: error: routine should only be called once\n");
    exit(-1);
  }
  if((*medium=(struct med *)calloc(1,sizeof(struct med)))
       ==NULL)MEMERROR("read_geometry: 1:");
#ifdef ALLOW_NON_3DQUAD_GEOM
  // kind of elastic approximation for 2D segments
  (*medium)->twod_approx_is_plane_stress = twod_approx_is_plane_stress;
  if(((*medium)->twod_approx_is_plane_stress)&&(half_plane)){
    fprintf(stderr,"read_geometry: error: half-plane only implemented for plane strain\n");
    exit(-1);
  }
#endif
  if((*fault=(struct flt *)calloc(1,sizeof(struct flt)))==
     NULL)MEMERROR("read_geometry: 2:");
  for(i=0;i<3;i++){// intialize medium boundaries for plotting
    (*medium)->xmax[i]=FLT_MIN;
    (*medium)->xmin[i]=FLT_MAX;
  } 
  (*medium)->nan = sqrt(-1.0);
  (*medium)->wmean= (*medium)->lmean=0.0;
  lmin = wmin =  FLT_MAX;
  lmax = wmax = -FLT_MAX;// L and W are positive quantities
  i=0;mus_avg=mud_avg=0.0;
  /*
    
    start geometry line input loop

  */
  while(fscanf(in,PATCH_CP_FORMAT,&(*fault+i)->x[X],&(*fault+i)->x[Y],
	       &(*fault+i)->x[Z],&(*fault+i)->strike,&(*fault+i)->dip,
	       &(*fault+i)->l,&(*fault+i)->w,&tmpint) == 8){

#ifdef ALLOW_NON_3DQUAD_GEOM
    if((*fault+i)->w == 0.0){
      /*

	2-D element
	
      */
      if(((*fault+i)->dip != 90.0)||((*fault+i)->x[Z] != 0.0)){
	fprintf(stderr,"read_geometry: width 0 selects: 2D mode: dip, z should be 90, and 0, respectively\n");
	exit(-1);
      }
      if((*medium)->twod_approx_is_plane_stress)
	(*fault+i)->type = TWO_DIM_SEGMENT_PLANE_STRESS;
      else{
	if(!half_plane)
	  (*fault+i)->type = TWO_DIM_SEGMENT_PLANE_STRAIN;
	else
	  (*fault+i)->type = TWO_DIM_HALFPLANE_PLANE_STRAIN;
      }
      nr_2d++;
    }else if(((*fault+i)->w < 0)&&((*fault+i)->l< 0)){
      /*
	
	TRIANGULAR ELEMENT
	
      */
      (*fault+i)->type = TRIANGULAR;
      (*fault+i)->xt=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*9);
      if(!(*fault+i)->xt)MEMERROR("read_geometry");
      // next nine fields are nodal coordinates
      if(fscanf(in,NINE_CP_FORMAT,
		&(*fault+i)->xt[  X],&(*fault+i)->xt[  Y],
		&(*fault+i)->xt[  Z],
		&(*fault+i)->xt[3+X],&(*fault+i)->xt[3+Y],
		&(*fault+i)->xt[3+Z],
		&(*fault+i)->xt[6+X],&(*fault+i)->xt[6+Y],
		&(*fault+i)->xt[6+Z])!=9)
	READ_ERROR("read_geometry");
      // x will be the centroid, to be calculated
      calc_centroid_tri((*fault+i)->xt,(*fault+i)->x);
      // determine local reference frame by means of the angles
      // w, l, and area will be the triangle area
      get_alpha_dip_tri_gh((*fault+i)->xt,&(*fault+i)->sin_alpha,
			   &(*fault+i)->cos_alpha,&tmpdbl,&(*fault+i)->area);
      (*fault+i)->dip=(float)RAD2DEGF(tmpdbl);
      (*fault+i)->l = (*fault+i)->w = (*fault+i)->area;
      alpha=RAD2DEGF(asin((*fault+i)->sin_alpha));
      (*fault+i)->strike= 90.0 - alpha;
#ifdef DEBUG
      fprintf(stderr,"read_geometry: fault %i is triangular, x1: (%g, %g, %g) x2: (%g, %g, %g) x3: (%g, %g, %g), area: %g\n",
	      i,(*fault+i)->xt[  X],(*fault+i)->xt[  Y],(*fault+i)->xt[  Z],
	      (*fault+i)->xt[3+X],(*fault+i)->xt[3+Y],(*fault+i)->xt[3+Z],
	      (*fault+i)->xt[6+X],(*fault+i)->xt[6+Y],(*fault+i)->xt[6+Z],
	      (*fault+i)->w);
#endif
      nr_triangle++;
    }else if((*fault+i)->l < 0){
      /*

	point source
	
      */
      fprintf(stderr,"read_geometry: fault %i: point source: area %g and \"aspect ratio\": %g\n",
	      i,(*fault+i)->w,-(*fault+i)->l);
      // W will hold the `fault' area, read in as w'
      (*fault+i)->area = (*fault+i)->w;
      (*fault+i)->type=POINT_SOURCE;
      nr_pt_src++;
    }else{
      /*

	regular, rectangular patch

      */
      (*fault+i)->type=RECTANGULAR_PATCH; 
    }
#else
    if(((*fault+i)->l <= 0)||((*fault+i)->w <= 0)){
      fprintf(stderr,"read_geometry: fault %i: half length l and width have to be >= 0 (%g/%g)!\n",
	      i,(*fault+i)->l,(*fault+i)->w);
      fprintf(stderr,"read_geometry: if 2-D, point source, or triangular elements were\n");
      fprintf(stderr,"read_geometry: what you were looking for,\n");
      fprintf(stderr,"read_geometry: recompile with  ALLOW_NON_3DQUAD_GEOM flag set\n");
      exit(-1);
    }
#endif
    if((*fault+i)->x[Z] > 0){
      fprintf(stderr,"read_geometry: patch %i: z has to be < 0, mid point z: %g\n",
	      i,(*fault+i)->x[Z]);
      exit(-1);
    }
    if(tmpint < 0){
      fprintf(stderr,"read_geometry: smallest allowed group number is 0, not %i\n",
	      tmpint);
      exit(-1);
    }else{
      (*fault+i)->group=(unsigned int)tmpint;
    }
    // check for illegal angles
    check_fault_angles((*fault+i));
    // do some stats for length and position
    for(j=0;j<3;j++){
      if((*fault+i)->x[j]<(*medium)->xmin[j])
	(*medium)->xmin[j] = (*fault+i)->x[j];
      if((*fault+i)->x[j]>(*medium)->xmax[j])
	(*medium)->xmax[j] = (*fault+i)->x[j];
    }
    if((*fault+i)->l > lmax)lmax=(*fault+i)->l;
    if((*fault+i)->l < lmin)lmin=(*fault+i)->l;
    if((*fault+i)->w > wmax)wmax=(*fault+i)->w;
    if((*fault+i)->w < wmin)wmin=(*fault+i)->w;
    (*medium)->wmean += (*fault+i)->w;
    (*medium)->lmean += (*fault+i)->l;
    (*fault+i)->area = (*fault+i)->l * (*fault+i)->w * 4.0;
    /*
      strike is defined as degrees clockwise from north (azimuth)
      we need the angle alpha, which is 
      counterclockwise from east and used for all rotations 
    */
#ifdef ALLOW_NON_3DQUAD_GEOM
    if((*fault+i)->type != TRIANGULAR){
#endif
      // if we have different kinds of faults, don't do this
      // calculation if it's a triangular patch since we 
      // have already calculates sin and cos (alpha) above
      alpha= 90.0 - (*fault+i)->strike;
      my_sincos_deg(&(*fault+i)->sin_alpha,&(*fault+i)->cos_alpha,
		    (COMP_PRECISION)alpha);
      if(fabs((*fault+i)->sin_alpha)< EPS_COMP_PREC)
	(*fault+i)->sin_alpha=0.0;
      if(fabs((*fault+i)->cos_alpha)< EPS_COMP_PREC)
	(*fault+i)->cos_alpha=0.0;
#ifdef ALLOW_NON_3DQUAD_GEOM
    }
#endif
    my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)(*fault+i)->dip);
    //
    // calculate the unity vectors in strike, dip, and normal
    // direction
    //
    calc_base_vecs((*fault+i)->t_strike,
		   (*fault+i)->normal,(*fault+i)->t_dip,
		   (*fault+i)->sin_alpha,
		   (*fault+i)->cos_alpha,sin_dip,cos_dip);
#ifdef SUPER_DEBUG
    fprintf(stderr,"fault %5i: strike: %g dip: %g sc_alpha:  %10.3e/%10.3e sc_dip: %10.3e/%10.3e\n",
	    i,90.0-alpha,(*fault+i)->dip,(*fault+i)->sin_alpha,
	    (*fault+i)->cos_alpha,sin_dip,cos_dip);
    fprintf(stderr," vec: s: (%10.3e,%10.3e,%10.3e) d: (%10.3e,%10.3e,%10.3e) n: (%10.3e,%10.3e,%10.3e)\n",
	    (*fault+i)->t_strike[X],(*fault+i)->t_strike[Y],(*fault+i)->t_strike[Z],
	    (*fault+i)->t_dip[X],(*fault+i)->t_dip[Y],(*fault+i)->t_dip[Z],
	    (*fault+i)->normal[X],(*fault+i)->normal[Y],(*fault+i)->normal[Z]);
#endif
#ifdef ALLOW_NON_3DQUAD_GEOM
    if((*fault+i)->type != TRIANGULAR){
#endif
      // determine geometrical boundaries for plotting
      calculate_corners(corner,(*fault+i),&lloc,&wloc);
#ifdef ALLOW_NON_3DQUAD_GEOM// select the number of corners
      jlim = (patch_is_2d((*fault+i)->type))?(2):(4);
#else
      jlim = 4;
#endif
      for(j=0;j<jlim;j++){
	// check depth alignment
	if(corner[j][Z] > eps_for_z){
	  fprintf(stderr,"read_geometry: patch %i, corner %i above surface, z: %20.10e (eps: %g)\n",
		  i,j,corner[j][Z],eps_for_z);
	  fprintf(stderr,"z: %g w: %g dip: %g\n",(*fault+i)->x[Z],(*fault+i)->w,(*fault+i)->dip);
	  fprintf(stderr,"read_geometry: exiting\n");
	  exit(-1);
	}
#ifdef ALLOW_NON_3DQUAD_GEOM
	if((*fault+i)->type == TWO_DIM_HALFPLANE_PLANE_STRAIN){
	  /* check if segment is sticking out into the air */
	  if((corner[0][Y] > 0 )||(corner[1][Y] > 0)){
	    fprintf(stderr,"read_geometry: error, half-plane segment %i endpoints: %g,%g and %g,%g\n",
		    i,corner[0][X] ,corner[0][Y],corner[1][X] ,corner[1][Y]);
	    exit(-1);
	  }
	}
#endif
	for(k=0;k<3;k++){
	  if(((*medium)->xmax[k]) < corner[j][k])
	    (*medium)->xmax[k] = corner[j][k];
	  if(((*medium)->xmin[k]) > corner[j][k])
	    (*medium)->xmin[k] = corner[j][k];
	}
      }
#ifdef ALLOW_NON_3DQUAD_GEOM
    }
#endif
    if(read_fault_properties){
      /* 
	 frictional properties, static and dynamic 
	 either read from file (-f switch)
      */
      if(fscanf(in2,"%f %f",&(*fault+i)->mu_s,&(*fault+i)->mu_d)!=2){
	fprintf(stderr,"read_geometry: read error: properties file: %s\n",
		FAULT_PROP_FILE);
	fprintf(stderr,"read_geometry: could not read two parameters for fault %i\n",
		i);
	exit(-1);
      }
      d_mu=(*fault+i)->mu_s - (*fault+i)->mu_d;
      if(d_mu > 1 || d_mu < 0){
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
    /* 
       initialize some of the fault arrays with zeros 
    */
    // coulomb stress correction factors
    (*fault+i)->cf[0] = (*fault+i)->cf[1] = 0.0;
    for(ic=0;ic<3;ic++){
      (*fault+i)->u[ic] = (*fault+i)->s[ic] = 0.0;
      (*fault+i)->sinc[ic] = 0.0;
    }
    (*fault+i)->mode[0] = INACTIVE;
    /* 
       add one more fault to the fault structure array 
    */
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
    fprintf(stderr,"read_geometry: %s %i\n",patch_filename,(*medium)->nrflt);
    exit(-1);
  }
  if((*medium)->nrflt > 99998)
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
  if(verbose)
  fprintf(stderr,"read_geometry: working with %i group(s) of patches\n",
	  (*medium)->nrgrp);
  // allocate space for the fault group arrays and initialize with zero
  if(((*medium)->fault_group=(struct fgrp *)calloc((*medium)->nrgrp,
						   sizeof(struct fgrp)))==NULL)
    MEMERROR("read_geometry: 5:");
  // calculate average friction values and geometries
  mus_avg /= (COMP_PRECISION)(*medium)->nrflt;
  mud_avg /= (COMP_PRECISION)(*medium)->nrflt;
  (*medium)->wmean /= (COMP_PRECISION)(*medium)->nrflt;
  (*medium)->lmean /= (COMP_PRECISION)(*medium)->nrflt;
  
  if((*medium)->nrflt<=(*medium)->max_nr_flt_files)
    for(i=0;i<(*medium)->nrflt;i++)
      fprintf(stderr,"read_geometry: fault %4i: x: (%9.6g, %9.6g, %9.6g) strike: %6.2f dip: %6.2f half_length: %6.3g aspect: %9.6g\n",
	      i,(*fault+i)->x[X],(*fault+i)->x[Y],
	      (*fault+i)->x[Z],(*fault+i)->strike,(*fault+i)->dip,(*fault+i)->l,
	      (*fault+i)->l/(*fault+i)->w);
  if(verbose){
    fprintf(stderr,"read_geometry: read in %i fault patch(es)\nread_geometry: half length: min/mean/max: %g/%g/%g\nread_geometry: half width:  min/mean/max: %g/%g/%g\n",
	    (*medium)->nrflt,lmin,(*medium)->lmean,lmax,wmin,(*medium)->wmean,wmax);
    fprintf(stderr,"read_geometry: average values for friction coefficients: mu_s/mu_d: %g/%g\n",
	    mus_avg,mud_avg);
  }
  //
  // adjust min/max of medium (for plotting only)
  //
  for(i=0;i<3;i++){
    (*medium)->xmin[i]=
      ((*medium)->xmin[i] < -MIN_GEOM_RANGE) ? 
      ((*medium)->xmin[i]) : (-MIN_GEOM_RANGE);
    (*medium)->xmax[i]=
      ((*medium)->xmax[i] >  MIN_GEOM_RANGE ) ?  
      ((*medium)->xmax[i] ): (MIN_GEOM_RANGE);
  }
#ifdef ALLOW_NON_3DQUAD_GEOM
  if(nr_pt_src + nr_triangle + nr_2d == 0){
    if(verbose){
      fprintf(stderr,"read_geometry: no special geometry patches were read in, recompiling without ALLOW_NON_3DQUAD_GEOM flag\n");
      fprintf(stderr,"read_geometry: would improve speed and size requirements of interact\n");
    }
  }else{
    fprintf(stderr,"read_geometry: WARNING: read in %i 2D elements, %i points, %i triangles, and %i rectangles\n",
	    nr_2d,nr_pt_src,nr_triangle,(*medium)->nrflt - nr_triangle - nr_pt_src - nr_2d);
    if(nr_2d)
      fprintf(stderr,"read_geometry: two dimensional approximation: plane %s %s\n",
	      ((*medium)->twod_approx_is_plane_stress)?("stress"):("strain"),
	      (half_plane)?("(half plane)"):("(full plane)"));
  }
#endif
  //
  // determine the position of each patch in a group local coordinate system
  // and initialize the group structure
  calculate_position_of_patch(*medium,*fault);
  (*medium)->geometry_init = TRUE;
  init = TRUE;
}

