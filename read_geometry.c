/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (c) Thorsten Becker, thwbecker@post.harvard.edu

*/
#include "properties.h"
#include "interact.h"


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
  int i,j,k,tmpint,ic,l,off;
  FILE *in,*in2=NULL;
  COMP_PRECISION sin_dip,cos_dip,mus_avg,mud_avg,d_mu,vertex[MAX_NR_EL_VERTICES*3],four_points[12],
    lloc,wloc,eps_for_z = EPS_COMP_PREC * 100.0,t_strike[3],t_dip[3],normal[3],area,earea;
#ifdef DEBUG
  COMP_PRECISION p[3];
#endif
  float wmin,wmax,lmin,lmax;
  double alpha;
  static my_boolean init=FALSE;
#ifdef ALLOW_NON_3DQUAD_GEOM
  int nr_pt_src=0,nr_triangle=0,nr_2d=0,nr_iquad=0;
#endif
  if(init){
    if((*medium)->comm_rank == 0)
      fprintf(stderr,"read_geometry: error: routine should only be called once\n");
    exit(-1);
  }
  if((*medium)->comm_rank != 0)
    verbose = 0;
  /* 

     READ IN FAULT GEOMETRY in the patch format as 
     documented in help_and_comments.c

  */ 
  if(strcmp(patch_filename,"stdin")!=0){
    if(verbose){
#ifdef ALLOW_NON_3DQUAD_GEOM
      fprintf(stderr,"read_geometry: reading geometry from \"%s\", non quads allowed\n",patch_filename);
#else
      fprintf(stderr,"read_geometry: reading geometry from \"%s\", non quads not allowed\n",patch_filename);
#endif
    }
    in = myopen(patch_filename,"r");
  }else{
    if(verbose){
#ifdef ALLOW_NON_3DQUAD_GEOM
      fprintf(stderr,"read_geometry: reading geometry from stdin, non quads allowed\n");
#else
      fprintf(stderr,"read_geometry: reading geometry from stdin, non quads not allowed\n");
#endif
    }
    in = stdin;
  }
  if(read_fault_properties){
    in2=fopen(FAULT_PROP_FILE,"r");
    if(in2){
      if((*medium)->comm_rank == 0)
	fprintf(stderr,"read_geometry: WARNING: reading fault properties from \"%s\"\n",
		FAULT_PROP_FILE);
    }else{
#ifdef DEBUG
      if((*medium)->comm_rank == 0)
	fprintf(stderr,"read_geometry: could not open \"%s\" for fault properties, using defaults\n",
		FAULT_PROP_FILE);
#endif
      read_fault_properties=FALSE;
    }
  }

 
#ifdef ALLOW_NON_3DQUAD_GEOM
  //
  // kind of elastic approximation for 2D segments
  //
  (*medium)->twod_approx_is_plane_stress = twod_approx_is_plane_stress;
  if(((*medium)->twod_approx_is_plane_stress)&&(half_plane)){
    if((*medium)->comm_rank == 0)
      fprintf(stderr,"read_geometry: error: half-plane only implemented for plane strain\n");
    exit(-1);
  }
#endif
  if((*fault=(struct flt *)calloc(1,sizeof(struct flt)))==NULL)MEMERROR("read_geometry: 2:");
  for(i=0;i<3;i++){// intialize medium boundaries for plotting
    (*medium)->xmax[i] = FLT_MIN;
    (*medium)->xmin[i] = FLT_MAX;
  } 
  (*medium)->nan = NAN;
  (*medium)->wmean= (*medium)->lmean=0.0;
  lmin = wmin =  FLT_MAX;
  lmax = wmax = -FLT_MAX;// L and W are positive quantities
  i=0;mus_avg=mud_avg=0.0;
  /*
    
    start geometry line input loop

  */
  while(fscanf(in,PATCH_CP_FORMAT,&(*fault+i)->x[INT_X],&(*fault+i)->x[INT_Y],
	       &(*fault+i)->x[INT_Z],&(*fault+i)->strike,&(*fault+i)->dip,
	       &(*fault+i)->l,&(*fault+i)->w,&tmpint) == 8){

#ifdef ALLOW_NON_3DQUAD_GEOM
    if((*fault+i)->w == 0.0){
      /*
	
	W == 0 --> 2-D

	2-D element
	
      */
      if(((*fault+i)->dip != 90.0)||((*fault+i)->x[INT_Z] != 0.0)){
	 if((*medium)->comm_rank == 0)
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
    }else if(((*fault+i)->w < 0)&&((*fault+i)->l< 0)){ /* length and
							  width
							  negative
							  will switch
							  to
							  triangle */
      /*
	L,W < 0 --> triangle
	
	TRIANGULAR ELEMENT
	
      */
      (*fault+i)->type = TRIANGULAR;
      (*fault+i)->xn=(COMP_PRECISION *)
	malloc(sizeof(COMP_PRECISION)*9);
      if(!(*fault+i)->xn)MEMERROR("read_geometry");
      // next nine fields are nodal coordinates
      if(fscanf(in,NINE_CP_FORMAT,
		&(*fault+i)->xn[  INT_X],&(*fault+i)->xn[  INT_Y],
		&(*fault+i)->xn[  INT_Z],
		&(*fault+i)->xn[3+INT_X],&(*fault+i)->xn[3+INT_Y],
		&(*fault+i)->xn[3+INT_Z],
		&(*fault+i)->xn[6+INT_X],&(*fault+i)->xn[6+INT_Y],
		&(*fault+i)->xn[6+INT_Z])!=9)
	READ_ERROR("read_geometry");
      /* check */
      for(j=0;j<3;j++){
	if((*fault+i)->xn[j*3+INT_Z] > 0){
	  fprintf(stderr,"read_geometry: triangular patch %i node %i above ground, error\n",
		  i,j);
	  exit(-1);
	}
      }
      /* 
	 compute all the triangular properties , basis vectors and
	 centroid
      */
      get_tri_prop_based_on_gh((*fault+i));
#ifdef DEBUG
      check_fault_normal_vectors((*fault+i));
      {			/* check for consistency */
	/*  */

	calc_global_strike_dip_from_local((*fault+i),t_strike,normal,t_dip);
	p[0] = dotp_3d((*fault+i)->t_strike,t_strike);
	p[1] = dotp_3d((*fault+i)->t_dip,t_dip);
	p[2] = dotp_3d((*fault+i)->normal,normal);
	if((p[0] < .95)||(p[1] < .95)||(p[2] < .95)){
	  fprintf(stderr,"tri: strike %5.2f dip %5.2f\n",(*fault+i)->strike,(*fault+i)->dip);
	  fprintf(stderr,"tri*g: s %5.2f d %5.2f n %5.2f\n",p[0],p[1],p[2]);
	  fprintf(stderr,"tri s: %5.2f %5.2f %5.2f (%.4e) d: %5.2f %5.2f %5.2f (%.4e) n: %5.2f %5.2f %5.2f (%.4e)\n qbsd: %5.2f %5.2f %5.2f (%.4e) d: %5.2f %5.2f %5.2f (%.4e) n: %5.2f %5.2f %5.2f (%.4e)\n",
		  (*fault+i)->t_strike[INT_X],(*fault+i)->t_strike[INT_Y],(*fault+i)->t_strike[INT_Z],norm_3d((*fault+i)->t_strike),
		  (*fault+i)->t_dip[INT_X],(*fault+i)->t_dip[INT_Y],(*fault+i)->t_dip[INT_Z],norm_3d((*fault+i)->t_dip),
		  (*fault+i)->normal[INT_X],(*fault+i)->normal[INT_Y],(*fault+i)->normal[INT_Z],norm_3d((*fault+i)->normal),
		  t_strike[INT_X],t_strike[INT_Y],t_strike[INT_Z],norm_3d(t_strike),
		  t_dip[INT_X],t_dip[INT_Y],t_dip[INT_Z],norm_3d(t_dip),
		  normal[INT_X],normal[INT_Y],normal[INT_Z],norm_3d(normal));
	}

      }
      if((*medium)->comm_rank == 0)
	fprintf(stderr,"read_geometry: fault %5i is triangular, x1: (%10.3e, %10.3e, %10.3e) x2: (%10.3e, %10.3e, %10.3e) x3: (%10.3e, %10.3e, %10.3e), area: %10.3e\n",
		i,(*fault+i)->xn[  INT_X],(*fault+i)->xn[  INT_Y],(*fault+i)->xn[  INT_Z],
		(*fault+i)->xn[3+INT_X],(*fault+i)->xn[3+INT_Y],(*fault+i)->xn[3+INT_Z],
		(*fault+i)->xn[6+INT_X],(*fault+i)->xn[6+INT_Y],(*fault+i)->xn[6+INT_Z],
		(*fault+i)->area);
#endif
      nr_triangle++;
   }else if(((*fault+i)->w < 0)){
      /*
	W < 0 --> iquad
	
      */
      (*fault+i)->type = IQUAD;
      /* 

	 five nodes in element, we read in only four (A = 3, B = 4, C = 1, D = 2) for
	 now, and assign auxiliary base node from average from A and B, call that 0
	 
	 there is a master triangle (012) and two auxiliary ones (302
	 and 041), and the stresses are evaluated within the master
	 triangle, weighting 

	 vertices 0, 1, and 2 by 0.5, 0.25, and 0.25, respectively (version I), or
	 vertices 3, 4, 1, 2 by 0.25 each (center of iquad (version II))

	 
	 D              C
	 2--------------1
         |\            /|
	 | \          / |
         |  \        /  |
	 |   \  X   /   |
	 |    \    /    |
	 | N1  \  / N2  |
	 |      \/      |
	 3 ---- 0 ----- 4
	 A              B
      */
      /* fault holds five nodes, plus two auxiliary fault plane angle
	 projection sets 3*(5+2*3)

      */
      (*fault+i)->xn=(COMP_PRECISION *)
	malloc(sizeof(COMP_PRECISION)*(MAX_NR_EL_VERTICES+3*2)*3);
      if(!(*fault+i)->xn)MEMERROR("read_geometry");
      // next 12 fields are nodal coordinates, read in four points
      for(j=0;j < 12;j++)
	if(fscanf(in,ONE_CP_FORMAT,(four_points+j))!=1)
	  READ_ERROR("read_geometry: iquad points cannot be read");
      /* check */
      for(j=0;j < 4;j++){
	if(four_points[j*3+INT_Z] > 0){
	  fprintf(stderr,"read_geometry: iquad patch %i node %i above ground, %.20e, error\n",
		  i,j,four_points[j*3+INT_Z] );
	  exit(-1);
	}
      }
      
      for(j=0;j < 3;j++){
	/* bottom middle point */
	(*fault+i)->xn[0*3+j] = (four_points[0*3+j] + four_points[1*3+j])/2.0;
	(*fault+i)->xn[1*3+j] = four_points[2*3+j]; /* top right */
	(*fault+i)->xn[2*3+j] = four_points[3*3+j]; /* top left */
	/*  */
	(*fault+i)->xn[3*3+j] = four_points[0*3+j]; /* bottom left */
	(*fault+i)->xn[4*3+j] = four_points[1*3+j]; /* bottom right */
      }
      /* compute all the triangular properties for the main triangle */
      get_tri_prop_based_on_gh((*fault+i));
#ifdef INTERACT_IQUAD_XC_VERSION_I
      /* compute the stress evaluation point (~ centroid) */
      for(j=0;j<3;j++){
	(*fault+i)->x[j]  = 0.5 * (*fault+i)->xn[0*3+j];
	(*fault+i)->x[j] += 0.25* (*fault+i)->xn[1*3+j];
	(*fault+i)->x[j] += 0.25* (*fault+i)->xn[2*3+j];
      }
#else
      for(j=0;j<3;j++){
	(*fault+i)->x[j]  = 0.25 * (*fault+i)->xn[3*3+j];
	(*fault+i)->x[j] += 0.25 * (*fault+i)->xn[4*3+j];
	(*fault+i)->x[j] += 0.25 * (*fault+i)->xn[1*3+j];
	(*fault+i)->x[j] += 0.25 * (*fault+i)->xn[2*3+j];
      }

#endif

      
      area = (*fault+i)->l * (*fault+i)->l;

      for(off=5,l=1;l<3;l++){
	/* compute auxiliary triagnle one and two */
	get_sub_normal_vectors((*fault+i),l,t_strike,t_dip, normal,&earea);
	/* aseemble projection matrix used to go from main triangle slip
	   to local slip systems */
	(*fault+i)->xn[off*3+0] = dotp_3d((*fault+i)->t_strike,t_strike); /* 5 and 8 for N1 and N2 */
	(*fault+i)->xn[off*3+1] = dotp_3d((*fault+i)->t_dip,   t_strike);
	(*fault+i)->xn[off*3+2] = dotp_3d((*fault+i)->normal,  t_strike);
	off++;
	(*fault+i)->xn[off*3+0] = dotp_3d((*fault+i)->t_strike,t_dip); /* 6 and 9 for N1 and N2 */
	(*fault+i)->xn[off*3+1] = dotp_3d((*fault+i)->t_dip,   t_dip);
	(*fault+i)->xn[off*3+2] = dotp_3d((*fault+i)->normal,  t_dip);
	off++;
	(*fault+i)->xn[off*3+0] = dotp_3d((*fault+i)->t_strike,normal); /* 7 and 10 for N1 and N2 */
	(*fault+i)->xn[off*3+1] = dotp_3d((*fault+i)->t_dip,   normal);
	(*fault+i)->xn[off*3+2] = dotp_3d((*fault+i)->normal,  normal);
	off++;
	area += earea;
      }
    
      (*fault+i)->l = (*fault+i)->w = sqrt(area);
      
#ifdef DEBUG
      if((*medium)->comm_rank == 0){
	fprintf(stderr,"read_geometry: fault %5i is iquad, x1: (%10.2e, %10.2e, %10.2e) x2: (%10.2e, %10.2e, %10.2e) x3: (%10.2e, %10.2e, %10.2e) x4: (%10.2e, %10.2e, %10.2e), area: %10.2e\n",
		i,
		(*fault+i)->xn[3*3+INT_X],(*fault+i)->xn[3*3+INT_Y],(*fault+i)->xn[3*3+INT_Z],
		(*fault+i)->xn[4*3+INT_X],(*fault+i)->xn[4*3+INT_Y],(*fault+i)->xn[4*3+INT_Z],
		(*fault+i)->xn[1*3+INT_X],(*fault+i)->xn[1*3+INT_Y],(*fault+i)->xn[1*3+INT_Z],
		(*fault+i)->xn[2*3+INT_X],(*fault+i)->xn[2*3+INT_Y],(*fault+i)->xn[2*3+INT_Z],
		(*fault+i)->area);
	fprintf(stderr,"read_geometry: main triangle: angles strike: %10.2e dip: %10.2e, l: %10.2e w: %10.2e\n",
		(*fault+i)->strike,		(*fault+i)->dip,(*fault+i)->l,(*fault+i)->w);
	fprintf(stderr,"read_geometry: main triangle:   strike: (%10.2e, %10.2e, %10.2e) dip:  (%10.2e, %10.2e, %10.2e) normal: (%10.2e, %10.2e, %10.2e)\n",
		(*fault+i)->t_strike[INT_X],		(*fault+i)->t_strike[INT_Y],		(*fault+i)->t_strike[INT_Z],
		(*fault+i)->t_dip[INT_X],		(*fault+i)->t_dip[INT_Y],		(*fault+i)->t_dip[INT_Z],
		(*fault+i)->normal[INT_X],		(*fault+i)->normal[INT_Y],		(*fault+i)->normal[INT_Z]);
	fprintf(stderr,"read_geometry: aux1 projection: strike: (%10.2e, %10.2e, %10.2e) dip:  (%10.2e, %10.2e, %10.2e) normal: (%10.2e, %10.2e, %10.2e)\n",
		(*fault+i)->xn[5*3+0],(*fault+i)->xn[5*3+1],(*fault+i)->xn[5*3+2],
		(*fault+i)->xn[6*3+0],(*fault+i)->xn[6*3+1],(*fault+i)->xn[6*3+2],
		(*fault+i)->xn[7*3+0],(*fault+i)->xn[7*3+1],(*fault+i)->xn[7*3+2]);
	fprintf(stderr,"read_geometry: aux2 projection: strike: (%10.2e, %10.2e, %10.2e) dip:  (%10.2e, %10.2e, %10.2e) normal: (%10.2e, %10.2e, %10.2e)\n",
		(*fault+i)->xn[8*3+0],(*fault+i)->xn[8*3+1],(*fault+i)->xn[8*3+2],
		(*fault+i)->xn[9*3+0],(*fault+i)->xn[9*3+1],(*fault+i)->xn[9*3+2],
		(*fault+i)->xn[10*3+0],(*fault+i)->xn[10*3+1],(*fault+i)->xn[10*3+2]);
      }
#endif
      nr_iquad++;
    }else if((*fault+i)->l < 0){
      /*

	L < 0 --> point source 
	
	POINT SOURCE
	
      */
      if((*medium)->comm_rank == 0)
	fprintf(stderr,"read_geometry: fault %i: point source: area %g and \"aspect ratio\": %g\n",
		i,(*fault+i)->w,-(*fault+i)->l);
      // W will hold the `fault' area, read in as w'
      (*fault+i)->area = (*fault+i)->w;
      (*fault+i)->type=POINT_SOURCE;
      nr_pt_src++;
    }else{
      /*

	regular, rectangular Okada type patch

      */
      (*fault+i)->type = OKADA_PATCH; 
    }
#else
    if(((*fault+i)->l <= 0)||((*fault+i)->w <= 0)){
      if((*medium)->comm_rank == 0){
	fprintf(stderr,"read_geometry: fault %i: half length l and width have to be >= 0 (%g/%g)!\n",
		i,(*fault+i)->l,(*fault+i)->w);
	fprintf(stderr,"read_geometry: if 2-D, point source, iquads, or triangular elements were\n");
	fprintf(stderr,"read_geometry: what you were looking for,\n");
	fprintf(stderr,"read_geometry: recompile with  ALLOW_NON_3DQUAD_GEOM flag set\n");
      }
      exit(-1);
    }
#endif
    if((*fault+i)->x[INT_Z] > 0){
      if((*medium)->comm_rank == 0)
	fprintf(stderr,"read_geometry: patch %03i: z has to be < 0, mid point z: %g\n",
		i,(*fault+i)->x[INT_Z]);
      exit(-1);
    }
    if(tmpint < 0){
      if((*medium)->comm_rank == 0)
	fprintf(stderr,"read_geometry: smallest allowed group number is 0, not %i\n",
		tmpint);
      exit(-1);
    }else{
      (*fault+i)->group=(unsigned int)tmpint;
    }
    //
    // do some stats for length and position
    //
    for(j=0;j<3;j++){
      if((*fault+i)->x[j]<(*medium)->xmin[j])
	(*medium)->xmin[j] = (*fault+i)->x[j];
      if((*fault+i)->x[j]>(*medium)->xmax[j])
	(*medium)->xmax[j] = (*fault+i)->x[j];
    }
#ifdef ALLOW_NON_3DQUAD_GEOM
    if( ( (*fault+i)->type != TRIANGULAR) || ( (*fault+i)->type != IQUAD) ){
#endif
      // check for illegal angles
      check_fault_angles((*fault+i));
      
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
      // if we have different kinds of faults, don't do this
      // calculation if it's a triangular patch since we 
      // have already calculates sin and cos (alpha) above
      alpha= 90.0 - (*fault+i)->strike;
      my_sincos_degd(&(*fault+i)->sin_alpha,&(*fault+i)->cos_alpha,
		     alpha);
      /* if(fabs((*fault+i)->sin_alpha)< EPS_COMP_PREC) */
      /* 	(*fault+i)->sin_alpha=0.0; */
      /* if(fabs((*fault+i)->cos_alpha)< EPS_COMP_PREC) */
      /* 	(*fault+i)->cos_alpha=0.0; */
      my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)(*fault+i)->dip);
      //
      // calculate the unity vectors in strike, dip, and normal
      // direction
      //
      /* 
	 executed for quads here 
      */
      calc_quad_base_vecs((*fault+i)->t_strike,
			  (*fault+i)->normal,(*fault+i)->t_dip,
			  (*fault+i)->sin_alpha,
			  (*fault+i)->cos_alpha,sin_dip,cos_dip);
      //fprintf(stderr,"tdq %g %g %g\n",(*fault+i)->t_dip[0],(*fault+i)->t_dip[1],(*fault+i)->t_dip[2]);
#ifdef SUPER_DUPER_DEBUG
      check_fault_normal_vectors((*fault+i));
      if((*medium)->comm_rank == 0){
	fprintf(stderr,"fault %5i: strike: %g dip: %g sc_alpha:  %10.2e/%10.2e sc_dip: %10.2e/%10.2e\n",
		i,90.0-alpha,(*fault+i)->dip,(*fault+i)->sin_alpha,
		(*fault+i)->cos_alpha,sin_dip,cos_dip);
	fprintf(stderr," vec: s: (%10.2e,%10.2e,%10.2e) d: (%10.2e,%10.2e,%10.2e) n: (%10.2e,%10.2e,%10.2e)\n",
		(*fault+i)->t_strike[INT_X],(*fault+i)->t_strike[INT_Y],(*fault+i)->t_strike[INT_Z],
		(*fault+i)->t_dip[INT_X],(*fault+i)->t_dip[INT_Y],(*fault+i)->t_dip[INT_Z],
		(*fault+i)->normal[INT_X],(*fault+i)->normal[INT_Y],(*fault+i)->normal[INT_Z]);
      }
#endif
      //
      // determine geometrical boundaries for plotting
      //
      calculate_vertices(vertex,(*fault+i),&lloc,&wloc); /* lloc and
							   wloc will
							   be full l
							   and w, not
							   half */
      for(j=0; j < nvert_of_patch((*fault+i));j++){ /* loop through all vertices */
	//
	// check depth alignment
	//
	if(vertex[j*3+INT_Z] > eps_for_z){
	  if((*medium)->comm_rank == 0){
	    fprintf(stderr,"read_geometry: patch %i, vertex %i above surface, z: %20.10e (eps: %g)\n",
		    i,j,vertex[j*3+INT_Z],eps_for_z);
	    fprintf(stderr,"z: %g w: %g dip: %g\n",(*fault+i)->x[INT_Z],(*fault+i)->w,(*fault+i)->dip);
	    fprintf(stderr,"read_geometry: exiting\n");
	  }
	  exit(-1);
	}
#ifdef ALLOW_NON_3DQUAD_GEOM
	if((*fault+i)->type == TWO_DIM_HALFPLANE_PLANE_STRAIN){
	  /* check if segment is sticking out into the air */
	  if((vertex[0*3+INT_Y] > 0 )||(vertex[1*3+INT_Y] > 0)){
	    if((*medium)->comm_rank == 0){
	      fprintf(stderr,"read_geometry: error, half-plane segment %i endpoints: %g,%g and %g,%g\n",
		      i,vertex[0*3+INT_X] ,vertex[0*3+INT_Y],vertex[1*3+INT_X] ,vertex[1*3+INT_Y]);
	    }
	    exit(-1);
	  }
	}
#endif
	for(k=0;k<3;k++){
	  if(((*medium)->xmax[k]) < vertex[j*3+k])
	    (*medium)->xmax[k] = vertex[j*3+k];
	  if(((*medium)->xmin[k]) > vertex[j*3+k])
	    (*medium)->xmin[k] = vertex[j*3+k];
	}
      }
#ifdef ALLOW_NON_3DQUAD_GEOM
    }else{			/* triangular */
      lloc = sqrt((*fault+i)->area)/2;
      (*medium)->wmean += lloc;
      (*medium)->lmean += lloc;
         
      if(lloc > lmax)lmax=lloc;
      if(lloc < lmin)lmin=lloc;
      if(lloc > wmax)wmax=lloc;
      if(lloc < wmin)wmin=lloc;
  
    }
#endif
    if(read_fault_properties){
      if(!in2){
	if((*medium)->comm_rank == 0)
	  fprintf(stderr,"read_geometry: in2 stream not open?!\n");
	exit(-1);
      }
      /* 
	 frictional properties, static and dynamic 
	 either read from file (-f switch)
      */
      if(fscanf(in2,"%f %f",&(*fault+i)->mu_s,&(*fault+i)->mu_d)!=2){
	if((*medium)->comm_rank == 0){
	  fprintf(stderr,"read_geometry: read error: properties file: %s\n",
		  FAULT_PROP_FILE);
	  fprintf(stderr,"read_geometry: could not read two parameters for fault %i\n",
		  i);
	}
	exit(-1);
      }
      d_mu=(*fault+i)->mu_s - (*fault+i)->mu_d;
      if(d_mu > 1 || d_mu < 0){
	if((*medium)->comm_rank == 0)
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
    if((*medium)->comm_rank == 0){
      fprintf(stderr,"read_geometry: did not read in any faults from geometry file\n");
      fprintf(stderr,"read_geometry: %s %i\n",patch_filename,(*medium)->nrflt);
    }
    exit(-1);
  }
  if((*medium)->nrflt > 99998)
    if((*medium)->comm_rank == 0)
      fprintf(stderr,"read_geometry: WARNING: parts of the program I/O might not work for more patches than 99999!\n");
  *fault=(struct flt *)
    realloc(*fault,sizeof(struct flt)*(*medium)->nrflt);
  //
  // determine the number of groups that were assigned
  //
  for((*medium)->nrgrp=i=0;i<(*medium)->nrflt;i++){
    if((*fault+i)->group > (*medium)->nrflt-1){
      if((*medium)->comm_rank == 0)
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
      if((*medium)->comm_rank == 0)
	fprintf(stderr,"read_geometry: fault %4i: x: (%9.6g, %9.6g, %9.6g) strike: %6.2f dip: %6.2f half_length: %6.3g aspect: %9.6g\n",
		i,(*fault+i)->x[INT_X],(*fault+i)->x[INT_Y],
		(*fault+i)->x[INT_Z],(*fault+i)->strike,(*fault+i)->dip,(*fault+i)->l,
		(*fault+i)->l/(*fault+i)->w);
  if(verbose){
    fprintf(stderr,"read_geometry: read in %i fault patch(es)\nread_geometry: half length: min/mean/max: %8.3e/%8.3e/%8.3e\nread_geometry: half width:  min/mean/max: %8.3e/%8.3e/%8.3e\n",
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
  if(nr_pt_src + nr_triangle + nr_2d + nr_iquad == 0){
    if(verbose){
      fprintf(stderr,"read_geometry: no non-quad patches were read in, recompiling without ALLOW_NON_3DQUAD_GEOM flag\n");
      fprintf(stderr,"read_geometry: might possibly improve speed and size requirements of interact\n");
    }
  }else{
    if((*medium)->comm_rank == 0){
      fprintf(stderr,"read_geometry: %i 2D elements, %i points, %i triangles, %i iquads, and %i regular quads\n",
	      nr_2d,nr_pt_src,nr_triangle,nr_iquad,(*medium)->nrflt - nr_triangle - nr_pt_src - nr_2d - nr_iquad);
      if(nr_2d)
	fprintf(stderr,"read_geometry: two dimensional approximation: plane %s %s\n",
		((*medium)->twod_approx_is_plane_stress)?("stress"):("strain"),
		(half_plane)?("(half plane)"):("(full plane)"));
    }
  }
#endif
  //
  // determine the position of each patch in a group local coordinate system
  // and initialize the group structure
  calculate_position_of_patch(*medium,*fault);
  
  (*medium)->geometry_init = TRUE;
  init = TRUE;
}

