/*

read in faults for the blockinvert program and assign displacements
resolved at the observational points gx[nrp*DIM]

if calc_intcoeff is false, will just read in faults

$Id: block_read_bflt.c,v 1.7 2011/01/07 07:19:58 becker Exp $


*/
#include "interact.h"
#include "blockinvert.h"
#include "properties.h"

void read_bflt(struct bmd *mod,COMP_PRECISION global_locking_depth,
	       struct prj projection,my_boolean echo_input,
	       COMP_PRECISION dl,
	       my_boolean rigid,my_boolean calc_intcoeff,
	       my_boolean check_doubles,
	       my_boolean override_locking_depth,
	       my_boolean override_dip,my_boolean fgeo_out,
	       my_boolean damp_nslip,my_boolean mod_codes_out,
	       my_boolean constrain_slip_direction)
{
  int i,j,n,nnf,norig,lblock[2],ncs,norig4;
  int *cs=NULL,cs_loc[2];
  COMP_PRECISION fac,df,newloc[2],oldloc[2],eps2,lfac,dip,locking_depth,
    rlon[2],rlat[2],dummy,*exloc,xctest1[6],xctest2[6],master_azi,
    sdamp[2];
  FILE *in,*gout;
  my_boolean dipwarned = FALSE,use_damp[2];
  char string[100];
#ifdef SUPERDUPER_DEBUG
  FILE *out;
  out=myopen("tmp.flt.dat","w");
  fprintf(stderr,"block_read_bflt: writing projected fault displacements to tmp.flt.dat\n");
#endif
  if(override_dip)		/* warn for dips */
    fprintf(stderr,"block_read_bflt: WARNING: all dips overwritten with 90 degree (vertical)\n");
#ifdef DEBUG
  fprintf(stderr,"block_read_bflt: attempting to read fault related data\n");
#endif
  /* number of constrained slip faults */
  ncs = 0;
  if(constrain_slip_direction){
    /* 
       
    read in constraints for slip directions
    
    */
    in = myopen(CSD_FILE,"r");
    my_ivecrealloc(&cs,(ncs+1)*4,"read_bflt: cs");
    while(fscanf(in,"%i %i %i",(cs+ncs*4),(cs+ncs*4+1),(cs+ncs*4+2))==3){
      if((cs[ncs*4] < 1)){ /* test for code in range */
	fprintf(stderr,"block_read_bflt: error: line %i: code %i: should be between 1 and total orig nflt\n",
		ncs+1,cs[ncs*4]);
	exit(-1);
      }
      cs[ncs*4]--;		/* this is the real code, 0... norig - 1 */
      for(j=1;j<3;j++){	/* test for range of strike and normal codes 
			   and assign
			*/
	if((cs[ncs*4+j] < -1)||(cs[ncs*4+j] > 1)){
	  fprintf(stderr,"block_read_bflt: error: line %i: slip code %i should be -1, 0, or 1 but is %i\n",
		  ncs+1,j,cs[ncs*4+j]);
	  exit(-1);
	}
      }
      cs[ncs*4+3] = 0;	/* this is the 'assigned' flag for later */
      /* incremeant */
      ncs++;
      my_ivecrealloc(&cs,(ncs+1)*4,"read_bflt: cs");
    }
    fclose(in);
    fprintf(stderr,"block_read_bflt: WARNING: read %i slip direction constraints from %s\n",
	    ncs,CSD_FILE);
  }
  /* 
     
     read in fault geometry from flt.block

  */
  in = fopen(FLTBLOCK_FILE,"r");
  if(in){
    eps2 = SQUARE(EPS_COMP_PREC);
				/* get fault structure */


  

    my_vecalloc(&exloc,4,"read_bflt: ex");
    /* 
       read in original faults in format:
       
       lon1 lat1 lon2 lat2 dip block_1 block_2 locking_factor locking_depth
       
       dip is from horizontal, ie. 90 degree is vertical.
       if dip is positive, this means that the fault is the hanging wall, ie.
       moving counterclockwise around a block, the dip is to your right.
       
       lon1,lat1 and lon2,lat2 define the surface trace

       format (lon, lat in degrees)
       
       n: real number of produced faults
       norig: number of faults read in (n >= norig)
       locking_factor: 1: completely locked, 0: freely slipping

       will be subdivided accoding to llim flag
       
    */
    n = norig = norig4 = 0;
    new_fault(&mod->fault,n);    
    while(fscanf(in,"%lf %lf %lf %lf %lf %i %i %lf %lf ",
		 (exloc+norig4+  INT_X),(exloc+norig4  +INT_Y),
		 (exloc+norig4+2+INT_X),(exloc+norig4+2+INT_Y),
		 &dip,(lblock),(lblock+1),&lfac,&locking_depth) == 9){
      for(i=0;i < 2;i++){
	if((lblock[i] <= 0) && (lblock[i] > mod->nrb)){
	  fprintf(stderr,"block_read_bflt: error: fault %i block codes not betweeen 1 and nrb: %i\n",
		  norig+1,mod->nrb);
	  exit(-1);
	}
	lblock[i]--;	/* internal codes from 0 .. nb-1 */
      }
      if(lfac < 0){	/* check for lfac values */
	fprintf(stderr,"block_read_bflt: lfac: %g doesn't make sense, should be >= 0\n",
		lfac);
	exit(-1);
      }
      if(lfac > 1){	/* check some more */
	fprintf(stderr,"block_read_bflt: lfac: WARNING: %g probably doesn't make sense, should likely be 0 < lfac < 1\n",
		lfac);
      }
      if(override_dip){
	/* overriding dip with 90 degree vertical faults */
	dip = 90.0;
      }
      if(fabs(dip - 90.0) > EPS_COMP_PREC){
	if(!dipwarned){
	  fprintf(stderr,"block_read_bflt: WARNING: at least one fault has dip != 90\n");
	  dipwarned=TRUE;
	}
	if(mod->nslip == 1){
	  fprintf(stderr,"block_read_bflt: WARNING: nslip was set to unity, resetting dip for fault %i to 90!\n",
		  norig+1);
	  dip = 90;
	}
      }
      if(locking_depth < 0){
	fprintf(stderr,"block_read_bflt: error, negative locking depth: %g flt: %i\n",
		locking_depth,norig+1);
	exit(-1);
      }
      if(check_doubles){
	//
	// check if this geometry path has been read in already
	//
	lonlat2xyz(exloc[norig4  +INT_X],exloc[norig4  +INT_Y],(xctest1));
	lonlat2xyz(exloc[norig4+2+INT_X],exloc[norig4+2+INT_Y],(xctest1+3));
	j = norig - 1;
	for(i=0;i < j;i++){
	  lonlat2xyz(exloc[i*4  +INT_X],exloc[i*4  +INT_Y],(xctest2));
	  lonlat2xyz(exloc[i*4+2+INT_X],exloc[i*4+2+INT_Y],(xctest2+3));
	  if(((distance_squared_3d(xctest1,xctest2) < eps2) &&
	      (distance_squared_3d((xctest1+3),(xctest2+3)) < eps2))||
	     ((distance_squared_3d((xctest1+3),xctest2) < eps2) &&
	      (distance_squared_3d(xctest1,(xctest2+3)) < eps2))){
	    fprintf(stderr,"block_read_bflt: error, fault geometry %i: (%g, %g) - (%g, %g) already read in as orig code %i\n",
		    norig+1,exloc[norig4+INT_X],exloc[norig4+INT_Y],
		    exloc[norig4+2+INT_X],exloc[norig4+2+INT_Y],i+1);
	    exit(-1);
	  }
	}
      }
      //
      // assign coordinates to test fault
      //
      for(i=0;i < 2;i++)
	for(j=0;j < 2;j++)
	  (mod->fault+n)->ex[i][j] = exloc[norig4 + i*2 + j];
      /* 
	 assign damping values for normal slip motion
      */
      if(damp_nslip){
	use_damp[0] = FALSE;use_damp[1] = TRUE; /* use only normal */
	sdamp[0] = 0.0; /* strike */
	sdamp[1] = mod->nslip_damp; /* normal slip damping value */
      }else{
	use_damp[0] = use_damp[1] = FALSE;
	sdamp[0] = sdamp[1] = 0.0;
      }
      /* init with no slip constraints */
      cs_loc[0] = cs_loc[1] = 0; 
      /* 
	 loop through the sloop direction 
	 constrained faults if we have any
	 
      */
      for(i=0;i < ncs;i++){	/* loop through all constraints */
	if(!cs[i*4+3]){	/* this constraint has not been 
			   assigned yet */
	  if(cs[i*4] == norig){ /* this fault is constrained */
	    cs[i*4+3]++;	/* set assigned flag */
	    if(cs[i*4+3] > 1){
	      fprintf(stderr,"block_read_bflt: error, orig fault %i was listed more than once as constrained in %s\n",
		      norig+1,CSD_FILE);
	      exit(-1);
	    }
	    for(j=1;j < 3;j++){
	      /* 
		 assign slip constraints to fault local array
	      */
	      cs_loc[j-1] = cs[i*4+j];
	    }
	  }
	}
      }	/* end loop through constraints */
      /*
	
      generate a new fault and increment the fault 
      number counter n
      
      this is a master segment, and we should determine the 
      master azimuth
      
      */
      generate_new_fault(&mod->fault,lblock,&n,mod->gx,mod->nrgp,
			 mod->sx,mod->nrsp,(override_locking_depth)?
			 (global_locking_depth):(locking_depth),
			 projection,mod->stress_depths,rigid,lfac,
			 calc_intcoeff,TRUE,&master_azi,dip,
			 norig,sdamp,use_damp,cs_loc);
      /*
	
      create subdividing faults if segment too long 
      
      */
      if((mod->fault+(n-1))->l > dl){	/* fault length > desired length */
	nnf = (int)((mod->fault+(n-1))->l/dl) + 1; /* nr of subdivisisons */
	if(nnf > 1){
	  n--;			
	  df = 1.0/((COMP_PRECISION)(nnf));	
	  for(i=0;i < 2;i++){	/* initialize great circle */
	    rlon[i] = DEG2RADF((mod->fault+n)->ex[i][0]);
	    rlat[i] = DEG2RADF((mod->fault+n)->ex[i][1]);
	    oldloc[i] = (mod->fault+n)->ex[0][i];
	  }
	  for(fac=df,i=0;i<nnf;i++,fac += df){
	    get_point_on_gc(rlon[0],rlat[0],rlon[1],rlat[1],
			    fac,(newloc+INT_X),(newloc+INT_Y),&dummy);
	    for(j=0;j<2;j++)
	      newloc[j] = RAD2DEGF(newloc[j]);
	    for(j=0;j<2;j++){	/* assign to new fault */
	      (mod->fault+n)->ex[0][j] = oldloc[j];
	      (mod->fault+n)->ex[1][j] = newloc[j];
	    }
	    generate_new_fault(&mod->fault,lblock,&n,mod->gx,mod->nrgp,
			       mod->sx,mod->nrsp,(override_locking_depth)?
			       (global_locking_depth):(locking_depth),
			       projection,mod->stress_depths,rigid,lfac,
			       calc_intcoeff,FALSE,&master_azi,dip,
			       norig,sdamp,use_damp,cs_loc);
	    a_equals_b_vector(oldloc,newloc,2);
	  }
	}
      }
      // make more room for original fault coordinates
      norig++;
      norig4 += 4;
      sprintf(string,"read_bflt: exloc: re: size; %i",4*(norig+1));
      my_vecrealloc(&exloc,4*(norig+1),string);
    }
    fprintf(stderr,"block_read_bflt: %s read in %i faults, produced %i\n",
	    (norig==0)?("WARNING:"):(""),norig,n);
    mod->nflt = n;			/* real number of faults */
    free(exloc);
    /* 
       
      check if assignment worked out OK for constrained slip directions
      
    */
    for(i=0;i<ncs;i++){
      if(cs[i*4] > norig){
	fprintf(stderr,"block_read_bflt: error, line %i of %s attempted to constrain fault %i but only %i original faults\n",
		i+1,CSD_FILE,cs[i*4]+1,norig);
	exit(-1);
      }
      if(cs[i*4+3] != 1){
	fprintf(stderr,"block_read_bflt: error, line %i of %s attempted to constrain fault %i but not found\n",
		i+1,CSD_FILE,cs[i*4]+1);
	  exit(-1);
      }
    }
    /* 
       
       check to see if we are damping slip modes 

    */
    mod->nfdamp=0;
    for(i=0;i < mod->nflt;i++)
      for(j=0;j < 2;j++)
	if((mod->fault+i)->use_damp[j])
	  mod->nfdamp++;
    if(mod->nfdamp)
      fprintf(stderr,"block_read_bflt: WARNING: damping normal slip motion with factor %g\n",
	      mod->nslip_damp);
    /* 
       make sure longitudes are in 0 ... 360 range 
    */
    for(i=0;i < mod->nflt;i++)
      for(j=0;j<2;j++)
	if((mod->fault+i)->ex[j][INT_X] < 0)
	  (mod->fault+i)->ex[j][INT_X] += 360.0;
    /* output? */
    if(echo_input)
      for(i=0;i < mod->nflt;i++){
	fprintf(stderr,"flt %4i: (%12g, %12g) - (%12g, %12g) l: %11g w: %11g bcode: %3i %3i azi: %11g dip: %g lfac: %g ld: %g\n",
		i+1,(mod->fault+i)->ex[0][INT_X],(mod->fault+i)->ex[0][INT_Y],
		(mod->fault+i)->ex[1][INT_X],(mod->fault+i)->ex[1][INT_Y],
		(mod->fault+i)->l,(mod->fault+i)->w,
		(mod->fault+i)->block[0]+1,
		(mod->fault+i)->block[1]+1,
		(mod->fault+i)->azi,(mod->fault+i)->dip,
		(mod->fault+i)->lfac,(mod->fault+i)->ld);
      }
    if(fgeo_out){
      fprintf(stderr,"block_read_bflt: writing fault geometry for testing to %s\n",
	      FGEOOUT_FILE);
      gout=myopen( FGEOOUT_FILE,"w");
      print_fault_geometry((mod->fault),mod->nflt,gout);
      fclose(gout);
    }
    if(mod_codes_out){		/* for output of all new files to code.llim.dat file */
      sprintf(string,"code.%g.dat",dl);
      gout=myopen(string,"w");
      fprintf(stderr,"block_read_bflt: writing modified flt.block type file to %s\n",
	      string);
      for(i=0;i < mod->nflt;i++){
	fprintf(gout,"%.9f %.9f %.9f %.9f %11g %5i %5i %11g %11g\n",
		(mod->fault+i)->ex[0][INT_X],(mod->fault+i)->ex[0][INT_Y],
		(mod->fault+i)->ex[1][INT_X],(mod->fault+i)->ex[1][INT_Y],
		(mod->fault+i)->dip,
		(mod->fault+i)->block[0]+1,
		(mod->fault+i)->block[1]+1,(mod->fault+i)->lfac,
		(mod->fault+i)->ld);
      }
      fclose(gout);
    }
#ifdef SUPERDUPER_DEBUG
    fclose(out);
#endif
  }else{			/* couldn't open file */
    mod->nflt = 0;
  }
  free(cs);
}

/*


  create a new fault structure with interaction coefficients etc.


*/
void generate_new_fault(struct bflt **fault,int *lblock, int *n,
			COMP_PRECISION *gx,int nrgp,
			COMP_PRECISION *sx,int nrsp,
			COMP_PRECISION locking_depth,
			struct prj projection,
			COMP_PRECISION *stress_depths,
			my_boolean rigid, COMP_PRECISION lfac,
			my_boolean calc_intcoeff,
			my_boolean is_master_fault,
			COMP_PRECISION *master_azi,
			COMP_PRECISION master_dip,
			int code, COMP_PRECISION *sdamp, 
			my_boolean *use_damp,
			int *cs_loc)

{
  COMP_PRECISION depth;
  int i;
  /* assign code to fault */
  (*fault+(*n))->orig_code = code;
   /* if dip = 90, set the vertical flag to true, else false */
  assign_bflt_dip_mode((*fault+(*n)),master_dip);
  /*

    get projection and fault parameters, depth gets overwritten
    
    input: ex, depth, rest output

    l, w are in the projected frame!

    gx is the center of the fault in geographical coordinates
    
  */
  get_projected_fault_parameters((*fault+(*n))->ex,// endpoints
				 locking_depth,
				 (*fault+(*n))->sx, /* lon, lat 
						       surface trace center */
				 &(*fault+(*n))->azi,// azimuth
				 &(*fault+(*n))->dip, /* dip */
				 &(*fault+(*n))->l, /* half length */
				 &(*fault+(*n))->w, /* half width  */
				 &depth);	/* center depth, > 0 */
  /* 
     assign preliminary center of patch coordinates
     those are OK, if dip == 90, else have to readjust
  */

  a_equals_b_vector((*fault+(*n))->x,(*fault+(*n))->sx,2); /* lon and lat */
  (*fault+(*n))->x[INT_Z] = -depth;	/* depth z is negative below ground */
  if((*fault+(*n))->l < 1.0e-6){
    fprintf(stderr,"error: block_read_bflt: error: fault half length: %g\n",
	    (*fault+(*n))->l);
    exit(-1);
  }
  //
  // check the azimuth with respect to the master azimuth
  //
  if(is_master_fault)
    *master_azi = (*fault+(*n))->azi;
  else{
    if(fabs(*master_azi - (*fault+(*n))->azi) > 10.0){// 180 flip
      (*fault+(*n))->azi -= 180.;
      (*fault+(*n))->dip = - (*fault+(*n))->dip;
      if(fabs(*master_azi - (*fault+(*n))->azi) > 10.0){
	fprintf(stderr,"block_read_bflt: error: master azi: %g slave: %g\n",
		*master_azi,(*fault+(*n))->azi);
	exit(-1);
      }
    }
  }
  //
  // get angles
  // alpha = 90 - azimuth
  //
  (*fault+(*n))->alpha = 90.0 - (*fault+(*n))->azi;
  my_sincos_deg(&(*fault+(*n))->sa,&(*fault+(*n))->ca,
		(*fault+(*n))->alpha);
  my_sincos_deg(&(*fault+(*n))->sd,&(*fault+(*n))->cd,
		(*fault+(*n))->dip);
  /*
    
    calculate the unity basis vectors (normalized) in strike,
    normal, and dip direction
    
  */

  calc_base_vecs(&(*fault+(*n))->evec[STRIKE*3], 
		 &(*fault+(*n))->evec[NORMAL*3], 
		 &(*fault+(*n))->evec[DIP*3], 
		 (*fault+(*n))->sa,(*fault+(*n))->ca,
		 (*fault+(*n))->sd,(*fault+(*n))->cd);
#ifdef BLOCK_SPHERICAL
  /* 
     
  convert the coordinates into spherical system, also prepare polar basis
  vector

  */
  // coordinates 
  lonlat2xyz_deg((*fault+(*n))->x[INT_X],(*fault+(*n))->x[INT_Y],
		 (*fault+(*n))->xc);
  // polar base
  calculate_polar_base((*fault+(*n))->x[INT_X],(*fault+(*n))->x[INT_Y],
		       (*fault+(*n))->pbase);
#endif

  /*  
      assign all parameters that depend on the locking depth,
      this includes the projected coordinates of the patch 
      center

      there is some redundancy in here for first time inits,
      leave in for now to make more modulare
  */
  assign_fault_locking_depth_parameters((*fault+(*n)),
					locking_depth,
					projection,TRUE,*n);
#ifdef SUPERDUPER_DEBUG  
  fprintf(stderr,"flt %i: x: %g, %g px: %g, %g (proj: %i: %g/%g/%g)\n", 
   	  *n,(*fault+(*n))->x[INT_X],(*fault+(*n))->x[INT_Y], 
	  (*fault+(*n))->px[INT_X],(*fault+(*n))->px[INT_Y], 
	  projection.type,projection.clon,projection.clat, 
	  projection.azi); 
#endif
  //
  // assign lblock codes to this fault structure
  for(i=0;i < 2;i++)
    (*fault+(*n))->block[i] = lblock[i];
  //
  // assign locking factor
  (*fault+(*n))->lfac = lfac;
  if(calc_intcoeff)
    // get the v and s interaction coefficients
    get_bflt_intcoeff(fault,*n,gx,nrgp,sx,nrsp,stress_depths,
		      rigid);
  /* 
     assign all damping values, but only activate if flag set 
  */
  for(i=0;i < 2;i++){
    (*fault+(*n))->sdamp[i] =    sdamp[i];
    (*fault+(*n))->use_damp[i] = use_damp[i];
  }
  /* 
     constrained slip
  */
  for(i=0;i < 2;i++){
    (*fault+(*n))->sc[i] = cs_loc[i];
    if(i == 0){
      /* strike */
      if(cs_loc[i])
	fprintf(stderr,"block_read_bflt: fault %5i (orig: %4i) strike slip constraint: %s\n",
		*n+1,code+1,(cs_loc[i] > 0)?("left-lateral"):("right-lateral"));
    }else{
      /* normal/thrust */
      if(cs_loc[i]){
	if(fabs((*fault+(*n))->dip - 90) > EPS_COMP_PREC){ /* thrust/normal */
	  fprintf(stderr,"block_read_bflt: fault %5i  (orig: %4i) (dipping normal) constraint: %s\n",
		  *n+1,code+1,(cs_loc[i] > 0)?("thrust, up-dip"):("normal, down-dip"));
	}else{	/* opening/closing */
	  fprintf(stderr,"block_read_bflt: fault %5i  (orig: %4i) vertical dip, normal constraint: %s\n",
		  *n+1,code+1,(cs_loc[i] > 0)?("opening"):("closing"));
	}
      }
    }
  }
  /* 
     increment fault number counter 
  */
  (*n) += 1;
  /* and make more room for new faults */
  new_fault(fault,*n);
}
/*

  initialize or reassign the v structure that holds the interaction
  coefficients for displacements and the s structure which holds the
  stress matrix. v is evaluated at the nrgp gx points and s at the
  nrsp sx points
  
*/
void get_bflt_intcoeff(struct bflt **fault, int iflt, 
		       COMP_PRECISION *gx, int nrgp,
		       COMP_PRECISION *sx,int nrsp,
		       COMP_PRECISION *stress_depths,
		       my_boolean rigid)
{
  COMP_PRECISION disp[3],pu[3],u[3],ps[3][3],s[3][3],x[3],px[3],
    dummy=0,sfac;
  int i,j,k,iret;
  /* 
     reallocate and clear the displacement matrix [3][3]
  */
  (*fault+iflt)->v = (struct vcm *)
    realloc((*fault+iflt)->v,nrgp*sizeof(struct vcm));
  if(!(*fault+iflt)->v)
    MEMERROR("block_read_bflt,2b1");
  for(i=0;i<nrgp;i++)		/* clear */
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	(*fault+iflt)->v[i].vc[j][k] = 0.0;
  if(nrsp){			
    /* 
       there are stress observations 
    */
    /* 
       reallocate stress interaction coefficients 
    */
    (*fault+iflt)->s = (struct scm *)
      realloc((*fault+iflt)->s,nrsp * sizeof(struct scm));
    if(!(*fault+iflt)->s)
      MEMERROR("block_read_bflt,2b2");
    /* clear */
    for(i=0;i<nrsp;i++)
      for(j=0;j<3;j++)
	for(k=0;k<6;k++)
	  (*fault+iflt)->s[i].sc[j][k] =  0.0;
  }
  if(!rigid){			
    /* 

    only if not rigid, will the fault displacement make any
    difference to strain or stress
    
    first assign strains to all points that have velocity
    observations
    

    */
    for(i=0;i < nrgp;i++){/* loop through all observational
			     points for velocities */
      /* 
	 
	 project observational point  
	 in the fault local oblique Mercator projection where
	 the fault itself is in the center
	 
      */
      a_equals_b_vector(x,(gx+i*BLOCK_DIM),2); 
      /* evaluate velocities at surface */
      x[INT_Z] = 0.0;			
      geoproject(x, px,FLT_ROT_PROJECTION ,(*fault+iflt)->x[INT_X],
		 (*fault+iflt)->x[INT_Y],(*fault+iflt)->azi,
		 dummy, dummy, dummy, dummy, (int)FALSE);
      for(j=0;j < 3;j++){	/* assign all slip modes */
	/* set up unity slip vector */
	get_right_slip(disp,j,(*fault+iflt)->lfac);
	// evaluate displacements in rotated frame
	eval_rectangle_basic(px,(*fault+iflt)->l,
			     (*fault+iflt)->w,
			     (*fault+iflt)->dip,
			     -(*fault+iflt)->x[INT_Z],
			     disp,pu,ps,&iret);
	if((iret)&&(j==0)){
	  fprintf(stderr,"block_read_bflt: error: singular at fault %i GPS vel obs point %i: %g, %g\n",
		  iflt+1,i+1,*(gx+i*BLOCK_DIM+INT_X),*(gx+i*BLOCK_DIM+INT_Y));
	  exit(-1);
	}
	// rotate pu back into the observational frame
	rotate_vec(pu, u, (*fault+iflt)->ca, -(*fault+iflt)->sa);
	for(k=0;k < 3;k++)	/* assign to interaction coeff 
				   and scale with locking factor
				*/
	  (*fault+iflt)->v[i].vc[j][k] = u[k];
#ifdef SUPERDUPER_DEBUG
	fprintf(out,"flt: %i op: %i  gx: %g %g ogx: %g %g u(%1i): %g %g\n",
		iflt+1,i+1,(*fault+iflt)->x[INT_X],(*fault+iflt)->x[INT_Y],
		gx[i*BLOCK_DIM+INT_X],gx[i*BLOCK_DIM+INT_Y],j+1,
		(*fault+iflt)->v[i].vc[j][INT_X],
		(*fault+iflt)->v[i].vc[j][INT_Y]);
#endif
      }
    }
    /*

    now assign stress interaction coefficients for stress
    observation points

    */
    //
    // multiply with slip factor
    sfac = (*fault+iflt)->lfac ;
    for(i=0;i < nrsp;i++){/* loop through all observational
			     points for stresses */
      a_equals_b_vector(x,(sx+i*BLOCK_DIM),2); 
      /* the observational frame counts  z depths negative */
      x[INT_Z] = -stress_depths[i];
      // stress observational point in the fault local
      // projection
      geoproject(x,px,FLT_ROT_PROJECTION,(*fault+iflt)->x[INT_X],
		 (*fault+iflt)->x[INT_Y],(*fault+iflt)->azi,
		 dummy, dummy, dummy, dummy, (int)FALSE);
      for(j=0;j < 3;j++){	/* use all slip modes */
	// get a displacement vector with fac slip
	// in the appropriate component (typically unity)
	get_right_slip(disp,j,sfac);	
	//
	eval_rectangle_basic(px,(*fault+iflt)->l,
			     (*fault+iflt)->w,
			     (*fault+iflt)->dip,
			     -(*fault+iflt)->x[INT_Z],
			     disp,pu,ps,&iret);
	if((iret)&&(j==0)){
	  fprintf(stderr,"block_read_bflt: error: singular at fault %i stress obs point %i: %g, %g\n",
		  iflt+1,i+1,*(sx+i*BLOCK_DIM+INT_X),*(sx+i*BLOCK_DIM+INT_Y));
	  exit(-1);
	}
	// rotate ps back into the observational frame
	rotate_mat_z(ps,s,(*fault+iflt)->ca, -(*fault+iflt)->sa);
	// assign stress tensor
	/* 
	   stress tensor s(k,l) for slip in the j direction at stress
	   observation i for fault iflt 
	*/
	(*fault+iflt)->s[i].sc[j][0] = s[INT_X][INT_X] / SHEAR_MODULUS;
	(*fault+iflt)->s[i].sc[j][1] = s[INT_X][INT_Y] / SHEAR_MODULUS;
	(*fault+iflt)->s[i].sc[j][2] = s[INT_X][INT_Z] / SHEAR_MODULUS;
	(*fault+iflt)->s[i].sc[j][3] = s[INT_Y][INT_Y] / SHEAR_MODULUS;
	(*fault+iflt)->s[i].sc[j][4] = s[INT_Y][INT_Z] / SHEAR_MODULUS;
	(*fault+iflt)->s[i].sc[j][5] = s[INT_Z][INT_Z] / SHEAR_MODULUS;
      }
    }
  }
}


/*

trash fault structure

*/
void free_bflt(struct bflt **fault, int nflt,
	       int nrsp)
{
  int i;
  for(i=0;i < nflt;i++){
    free((*fault+i)->v);
    if(nrsp)
      free((*fault+i)->s);
#ifdef BLOCK_SPHERICAL
    free((*fault+i)->xc);
    free((*fault+i)->pbase);
#endif    
  }
  free(*fault);
}      
/* determine fault dip mode, ie. is it vertical? */
void assign_bflt_dip_mode(struct bflt *fault,COMP_PRECISION dip)
{
  fault->dip = dip;
  if(fabs(fabs(dip) - 90.0) < EPS_COMP_PREC)
    fault->vertical = TRUE;
  else
    fault->vertical = FALSE;
}
/* used by block_checkflt */
void flip_block_code(struct bflt *fault, COMP_PRECISION *dir)
{
  int a;
  COMP_PRECISION x;
  a = fault->block[0];
  fault->block[0] = fault->block[1];
  fault->block[1] = a;
  x = dir[0];
  dir[0] = dir[1];
  dir[1] = x;
}
/* 
   print fault surface end points and center for testing 
   purposes. format:

   lon1 lat1 lon2 lat2   clon clat cz(>0 km)   azi dip l w
   
*/
void print_fault_geometry(struct bflt *fault,int nflt, FILE *gout)
{
  int i;
  
  for(i=0;i<nflt;i++)
    fprintf(gout,"%g %g %g %g\t%g %g %g\t%g %g\t%g %g\n",
	    fault[i].ex[0][INT_X], fault[i].ex[0][INT_Y],
	    fault[i].ex[1][INT_X], fault[i].ex[1][INT_Y],
	    fault[i].x[INT_X],fault[i].x[INT_Y],fault[i].x[INT_Z],
	    fault[i].azi, fault[i].dip,
	    fault[i].l,fault[i].w);
}
/* 
   given an initialized fault structure, (re)assign all
   variables that depend on the locking depth ld

*/
void assign_fault_locking_depth_parameters(struct bflt *fault,
					   COMP_PRECISION ld,
					   struct prj projection,
					   my_boolean verbose, /* those
							       two only 
							       for output 
							    */
					   int nflt)
{
  COMP_PRECISION x[3],px[3],dummy=0;
  int i;
  fault->ld = ld;
#ifdef DEBUG
  if(fabs(fault->dip) < EPS_COMP_PREC){
    fprintf(stderr,"reassign_fault_locking_depth: error: dip: %g\n",
	    fault->dip);
    exit(-1);
  }
#endif
  /* depth of center of patch */
  fault->x[INT_Z] = -ld/2.0;	/* x[INT_Z] is < 0 */
  fault->w = -fault->x[INT_Z]/sin((fault->dip)*DEG2RAD);
  /* 

     move the mid point coordinates according to the fault dip

  */
  if(!fault->vertical){
    //
    // mid point of fault trace at surface
    a_equals_b_vector(x,fault->sx,2);x[INT_Z]=0.0; 
    // go to projected system
    geoproject(x,px,projection.type,projection.clon,
	       projection.clat,projection.azi, dummy, dummy,
	       projection.lat1,projection.lat2,(int)FALSE);
    //
    // in the projected system, we can add in km
    // x_c = x_h - w * d
    //
    for(i=0;i < 3;i++)
      px[i] -= fault->w * fault->evec[DIP*3+i];
    //
    // do inverse projection, now for the patch center coordinates
    //
    geoproject(px,fault->x,
	       projection.type,projection.clon,
	       projection.clat,projection.azi, dummy, dummy,
	       projection.lat1,projection.lat2,(int)TRUE);
    if(verbose)
      fprintf(stderr,"flt: %i dip: %g smp: (%g, %g) ctr: (%g, %g, %g)\n",
	      nflt,fault->dip,x[INT_X],x[INT_Y],fault->x[INT_X],fault->x[INT_Y], 
	      fault->x[INT_Z]);
  }
  //
  // convert the fault patch center to general projected 
  // coordinates
  //
  geoproject(fault->x,fault->px,
	     projection.type,projection.clon,
	     projection.clat,projection.azi, dummy, dummy,
	     projection.lat1,projection.lat2,(int)FALSE);
}

/* 

add another fault to list and return increment number, don't change 
current number of faults as input

*/
int new_fault(struct bflt **flt, int n)
{
  *flt = (struct bflt *)realloc(*flt,sizeof(struct bflt) * (n + 1));
  if(! *flt){
    fprintf(stderr,"new_fault: memory error: current n: %i\n",n);
    exit(-1);
  }
  init_bflt((*flt+n)); 
  return n+1;
}

/*

  initialize bflt structure

*/
void init_bflt(struct bflt *fault)
{
  fault->block[0] = fault->block[1] = 0;
  fault->x[INT_X]=fault->x[INT_Y]=fault->x[INT_Z] = 0.0;
  fault->px[INT_X]=fault->px[INT_Y]=fault->px[INT_Z] = 0.0;
  // so that we can use realloc
  fault->s = NULL;
  fault->v = NULL;
#ifdef BLOCK_SPHERICAL
  /* for the computation of  location in cartesian systems */
  my_vecalloc(&fault->xc,3,"init_bflt: xc");
  my_vecalloc(&fault->pbase,9,"init_bflt: pbase");
#else
  fault->xc = NULL;
  fault->pbase = NULL;
#endif
}
