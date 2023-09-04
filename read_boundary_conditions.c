/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu

*/
#include "properties.h"
#include "interact.h"



/* 
   read in boundary conditions for faults,
   a and b are the background stressing values

   if init_system is FALSE, will only read BC conditions,
   for regular computation, has to be TRUE
*/
void read_boundary_conditions(struct med *medium, 
			      struct flt *fault,
			      COMP_PRECISION *a, 
			      COMP_PRECISION *b,
			      my_boolean init_system)
{
  FILE *in;
  int patch_nr,bc_code,n=0,inc,start_patch,stop_patch,i,j;
  my_boolean printevery;
  in=fopen(BC_FILE,"r");
  if(in != NULL){
    HEADNODE
      fprintf(stderr,"read_boundary_conditions: reading boundary conditions from \"%s\"\n",
	      BC_FILE);
    fscanf(in,"%i",&medium->op_mode);
    switch(medium->op_mode){
    case ONE_STEP_CALCULATION:{
      /*
	resize the fault mode array to allow 
	for three different activations
      */
      medium->nr_flt_mode=3;
      for(i=0;i<medium->nrflt;i++){
	for(j=0;j<medium->nr_flt_mode;j++)
	  fault[i].mode[j]=INACTIVE;
      }
      read_one_step_bc(in,medium,fault,a,b,init_system);
      break;
    }
    /*

      SIMULATE LOADING EXPERIMENT part

    */
    case SIMULATE_LOADING:
    case SIMULATE_LOADING_AND_PLOT:{
      medium->nr_flt_mode=1;
      HEADNODE
	fprintf(stderr,"read_boundary_conditions: simulating loading\n");
      /* timing issues */
      fscanf(in,THREE_CP_FORMAT,&medium->dt0,&medium->print_interval,
	     &medium->stop_time);
      if(!medium->variable_time_step){
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: until time %g in fixed steps of %g\n",
		  medium->stop_time,medium->dt0);
      }else{
	HEADNODE{
	  fprintf(stderr,"read_boundary_conditions: until time %g in variable time steps\n",
		  medium->stop_time);
	  fprintf(stderr,"read_boundary_conditions: the timestep is calculated from MIN(time_to_failure, print_interval)\n");
	}
	/*
	  importantly, we will set the base time step to 
	  unity to facilitate computations of the adequate 
	  timestep based on the standard increment
	*/
	medium->dt0=1.0;
      }
      if(medium->dt0 <= 0){
	HEADNODE
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
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: x plotting every %g on window 1, scroll every %g\n",
		  medium->x_plot_interval,
		  medium->x_scroll_interval);
	medium->x_scroll_inc=
	  medium->x_scroll_interval;
#else
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: you requested X output but program was compiled without PGPLOT support\n");
	fscanf(in,"%*f %*f");
#endif
      }
      /* 
	 read in fault slip mode codes for loading simulation
      */
      n=0;
      while(fscanf(in,"%i %i",&patch_nr,&bc_code) == 2){
	n++;
	if(patch_nr < 0){// we will set several faults to same boundary code
	  start_patch=bc_code;
	  if(fscanf(in,"%i %i",&stop_patch,&bc_code)!=2){
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: read error for multiple assignment\n");
	    exit(-1);
	  }
	  if(stop_patch < 0){
	    start_patch= 0;
	    stop_patch = medium->nrflt - 1;
	  }
	  inc = -patch_nr;
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: const. %5i: patches %5i to %5i: activation mode %3i, %s\n",
		    n,start_patch,stop_patch,bc_code,comment_on_code(bc_code));
	  printevery=FALSE;
	}else{// single fault mode
	  start_patch = patch_nr; 
	  stop_patch  = patch_nr;
	  inc=1;
	  printevery=TRUE;
	}
	if(abs(bc_code)>=SHRT_MAX){
	  HEADNODE{
	    fprintf(stderr,"read_boundary_conditions: bc code %i won't fit in a short int variable\n",
		    bc_code);
	    fprintf(stderr,"read_boundary_conditions: this means either an input or an implementation error\n");
	  }
	  exit(-1);
	}

	for(patch_nr=start_patch;patch_nr <= stop_patch;patch_nr += inc){
	  // check for out of bounds flt nr
	  if((patch_nr > (medium->nrflt-1)) || (patch_nr < 0)){
	    HEADNODE{
	      fprintf(stderr,"read_boundary_conditions: file %s, line %i, patch number %i out of range (max nr patch: %i)\n",
		      BC_FILE,n,patch_nr,medium->nrflt-1);
	      fprintf(stderr,"read_boundary_conditions: numbering for patches goes from 0 ... N-1, where N is the number of patches\n");
	    }
	    exit(-1);
	  }
#ifdef ALLOW_NON_3DQUAD_GEOM
	  // check for 2-D disallowed codes
	  if(patch_is_2d(fault[patch_nr].type)){
	    if((bc_code == DIP_SLIP)||(bc_code == DIP_SLIP_UPWARD)||
	       (bc_code == DIP_SLIP_DOWNWARD) || 
	       (bc_code == COULOMB_DIP_SLIP_UPWARD) ||
	       (bc_code == COULOMB_DIP_SLIP_DOWNWARD) ||
	       (bc_code == COULOMB_MAXSDIR_SLIP) || (bc_code == MAXSDIR_SLIP))
	      HEADNODE
		fprintf(stderr,"read_boundary_conditions: patch  %i is 2-D, dip motion is not allowed\n",
			patch_nr);
	    exit(-1);
	  }
#endif
	  // select the modes
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
	      HEADNODE
		fprintf(stderr,"read_boundary_conditions: const. %5i fault %5i activation mode %3i, %s\n",
			n,patch_nr,bc_code,comment_on_code(bc_code));
	    if(fault[patch_nr].mode[0] != INACTIVE){
	      HEADNODE{
		fprintf(stderr,"read_boundary_conditions: fault %i has already been assigned mode %i\n",
			patch_nr,fault[patch_nr].mode[0]);
		fprintf(stderr,"read_boundary_conditions: loading simulation run can only cope with one mode\n");
	      }
	      exit(-1);
	    }
	    fault[patch_nr].mode[0]=(MODE_TYPE)bc_code;
#ifdef SUPER_DEBUG
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: fault %5i: code: %i\n",
		      patch_nr,fault[patch_nr].mode[0]);
#endif
	    break;
	  }
#ifdef NO_OPENING_MODES
	  case NORMAL_SLIP:
	  case NORMAL_SLIP_OUTWARD:
	  case NORMAL_SLIP_INWARD:{
	    HEADNODE{
	      fprintf(stderr,"read_boundary_conditions: normal (opening) slip modes have been switched off, recompile\n");
	      fprintf(stderr,"read_boundary_conditions: without NO_OPENING_MODES switch\n");
	    }
	    exit(-1);
	    break;
	  }
#endif
	  default:{
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: fault activation mode %i undefined\n",
		      bc_code);
	    exit(-1);
	    break;
	  }}
	}
      }
      break;
    } /* end loading simulation part */
    default:{
      HEADNODE
	fprintf(stderr,"read_boundary_conditions: operational mode %i undefined\n",medium->op_mode);
      exit(-1);
      break;
    }}
    fclose(in);
  }else{
    HEADNODE
      fprintf(stderr,"read_boundary_conditions: can not open BC file \"%s\", exiting.\n",
	      BC_FILE);
    exit(-1);
  }
  medium->bc_init=TRUE;
}   

/* 

   read in boundary conditions for one step calculation mode

   
   init_equation_system should be TRUE for normal operation
   only set to FALSE if you are just interested in reading
   in the boundary conditions, but not performing any calculations
   


*/
void read_one_step_bc(FILE *in,struct med *medium,struct flt *fault,
		      COMP_PRECISION *a,COMP_PRECISION *b,
		      my_boolean init_system)
{
  my_boolean *sma,printevery,added_to_con,added_to_uncon,
    print_res_stress_screen=FALSE,
    rotate_to_local=FALSE,check_for_res_stress_output=FALSE,use_dip,
    sma_local[3],slip_bc_assigned=FALSE,stress_bc_assigned=FALSE;
  int patch_nr,bc_code,n,inc,start_patch,stop_patch,i,j,i3,nrflt3,
    patch_nr3;
  COMP_PRECISION bc_value,*rhs,global_strike,global_dip,slip[3],
    global_dip_rad,global_alpha_rad,gstrike[3],gnormal[3],gdip[3],
    sm[3][3],bstress[3],stress_drop,sin_global_alpha_rad,
    cos_global_alpha_rad, sin_global_dip_rad, cos_global_dip_rad;
  FILE *tmp_in,*rsout=NULL;
  if(medium->nr_flt_mode != 3){
    HEADNODE
      fprintf(stderr,"read_one_step_bc: have to allow for %i modes, only %i set\n",
	      3,medium->nr_flt_mode);
    exit(-1);
  }
  HEADNODE
    fprintf(stderr,"read_boundary_conditions: running one step calculation\n");
  //
  // read in the output field option
  //
  if(fscanf(in,"%i",&i)!=1){
    HEADNODE
      fprintf(stderr,"read_boundary_conditions: read error, output field option\n");
    exit(-1);
  }
  if(i == 1)// use a grid for output locations
    medium->print_bulk_fields = TRUE;
  else if(i == 2)// read output locations from file
    medium->read_oloc_from_file = TRUE;
  else if(i != 0){
    HEADNODE
      fprintf(stderr,"read_boundary_conditions: error, output field mode %i undefined\n",
	      i);
    exit(-1);
  }
  if(medium->print_bulk_fields){
    HEADNODE
      fprintf(stderr,"read_boundary_conditions: will print bulk fields on grid\n");
    /* 
       geometry for field outputs 
    */
    if(fscanf(in,FIELD_CP_FORMAT ,
	      &medium->pxmin[INT_X],&medium->pxmax[INT_X],&medium->n[INT_X],
	      &medium->pxmin[INT_Y],&medium->pxmax[INT_Y],&medium->n[INT_Y],
	      &medium->pxmin[INT_Z],&medium->pxmax[INT_Z],&medium->n[INT_Z])!=9){
      fprintf(stderr,"read_boundary_conditions: read error grid bounds\n");
      exit(-1);
    }
    if(medium->n[INT_Z] >= 0){// normal grid
      HEADNODE{
	fprintf(stderr,"read_boundary_conditions: xmin: %11g xmax: %11g n: %5i dx: %11g\n",
		medium->pxmin[INT_X],medium->pxmax[INT_X],
		medium->n[INT_X],medium->n[INT_X]!=1?(medium->pxmax[INT_X]-medium->pxmin[INT_X])/
		((COMP_PRECISION)(medium->n[INT_X]-1)):0);
	fprintf(stderr,"read_boundary_conditions: ymin: %11g ymax: %11g m: %5i dy: %11g\n",
		medium->pxmin[INT_Y],medium->pxmax[INT_Y],
		medium->n[INT_Y],medium->n[INT_Y]!=1?(medium->pxmax[INT_Y]-medium->pxmin[INT_Y])/
		((COMP_PRECISION)(medium->n[INT_Y]-1)):0.0);
	fprintf(stderr,"read_boundary_conditions: zmin: %11g zmax: %11g o: %5i dz: %11g\n",
		medium->pxmin[INT_Z],medium->pxmax[INT_Z],
		medium->n[INT_Z],medium->n[INT_Z]!=1?(medium->pxmax[INT_Z]-medium->pxmin[INT_Z])/
		((COMP_PRECISION)(medium->n[INT_Z]-1)):0.0);
      }
    }else{
      if(medium->n[INT_Z] < -2){
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: error: nz: %i\n",
		  medium->n[INT_Z]);
	exit(-1);
      }
      //
      // create a field that is in the plane of the average strike and
      // dip or strike and normal directions of fault group 0
      //
      HEADNODE{
	fprintf(stderr,"read_boundary_conditions: uses field output within the average fault plane of group 0\n");
	fprintf(stderr,"read_boundary_conditions: strike extent: min: %g max: %g n: %i d_strike: %g\n",
		medium->pxmin[INT_X],medium->pxmax[INT_X],
		medium->n[INT_X],
		(medium->n[INT_X]!=1)?(medium->pxmax[INT_X]-medium->pxmin[INT_X])/
		((COMP_PRECISION)(medium->n[INT_X]-1)):0);
	fprintf(stderr,"read_boundary_conditions: %s extent:   min: %g max: %g m: %i d_%s:    %g\n",
		(medium->n[INT_Z] == -1)?("dip"):("normal"),medium->pxmin[INT_Y],medium->pxmax[INT_Y],
		medium->n[INT_Y],(medium->n[INT_Z] == -1)?("dip"):("normal"),
		(medium->n[INT_Y]!=1)?(medium->pxmax[INT_Y]-medium->pxmin[INT_Y])/
		((COMP_PRECISION)(medium->n[INT_Y]-1)):0.0);
      }
    }
  }else{
    HEADNODE
      fprintf(stderr,"read_boundary_conditions: will not print bulk fields\n");
  }
  if(medium->read_oloc_from_file){
    /*

      read in the irregularly spaced output locations from file which has to be in 

      x y z 

      format

    */
    HEADNODE
      fprintf(stderr,"read_boundary_conditions: reading irregular output locations from \"%s\"\n",
	      OLOC_FILE);
    i = medium->olocnr = 0;
    medium->xoloc=(float *)malloc((i+3)*sizeof(float));
    if(!medium->xoloc)MEMERROR("");
    tmp_in=myopen(OLOC_FILE,"r");
    while(fscanf(tmp_in,"%f %f %f",(medium->xoloc+i+INT_X),
		 (medium->xoloc+i+INT_Y),(medium->xoloc+i+INT_Z))==3){
      medium->olocnr++;i+=3;
      medium->xoloc=(float *)realloc(medium->xoloc,(i+3)*sizeof(float));
      if(!medium->xoloc)MEMERROR("");
    }
    fclose(tmp_in);
    if(medium->olocnr){
      HEADNODE
	fprintf(stderr,"read_boundary_conditions: read %i locations for output\n",
		medium->olocnr);
    }else{
      HEADNODE
	fprintf(stderr,"read_boundary_conditions: error, no locations read from \"%s\"\n",
		OLOC_FILE);
      exit(-1);
    }
  }
  /* 

     allocate temporary arrays for all faults, will be deleted further
     down if these get too big, have to do more bookkeeping, now, for
     simplicity keep like that (otherwise fault constraints would have
     to be sorted)

  */
  nrflt3 = medium->nrflt * 3;
  if((sma=(my_boolean *)malloc(nrflt3*sizeof(my_boolean))) == NULL)
    PMEMERROR("read_boundary_conditions: 2:");
  for(i=0;i < nrflt3;i++)
    sma[i]=INACTIVE;
  if((rhs=(COMP_PRECISION *)calloc(nrflt3,sizeof(COMP_PRECISION)))==NULL)
    PMEMERROR("read_boundary_conditions: 3:");
  n=0;
  /* 

     ----------------------------------------------------------------------

     read in constraints in format 

     fault_patch_number boundary_condition_code boundary_condition_value


     ----------------------------------------------------------------------
     
  */
  while(fscanf(in,IIF_CP_FORMAT,&patch_nr,&bc_code,&bc_value) == 3){
    n++;
    if(patch_nr < 0){/* the fault patch indicator was negative, meaning
		      we will assign the next values to a sequence of
		      fault patches */
      if((bc_code < 0) && (bc_value < 0)){
	/* both bc_code and bc_values are negative,
	   assign this boundary code to all patches */
	start_patch = 0;
	stop_patch = medium->nrflt-1;
      }else{
	/* bc_code indicates the first patch */
	start_patch = bc_code;
	/* bc_values the last patch for assignment */
	stop_patch = (int)bc_value;
      }
      /* increment */
      inc = -patch_nr;
      /* read in the boundary code and bc_value for real */
      if(fscanf(in,IF_CP_FORMAT,&bc_code,&bc_value) != 2){
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: read error for assigning codes from %i to %i in %i increments\n",
		  start_patch,stop_patch,inc);
	exit(-1);
      }
      printevery = FALSE;
      HEADNODE
	fprintf(stderr,"read_boundary_conditions: const: %5i, patches %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
		n,start_patch,stop_patch,bc_code,bc_value,
		comment_on_code_bc(bc_code,bc_value));
#ifdef ALLOW_NON_3DQUAD_GEOM
      // check if we have triangular elements and prescribed slip
      if(slip_type_bc(bc_code)){
	for(rotate_to_local=FALSE,
	      i=start_patch;i <= stop_patch;i++)
	  if((i > medium->nrflt - 1)||(i<0)){
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: error, trying to address patch %i out of %i total\n",
		      i+1,medium->nrflt);
	    exit(-1);
	  }
	  if(fault[i].type == TRIANGULAR){
	    rotate_to_local = TRUE;
	    break;
	  }
      }else{
	rotate_to_local=FALSE;
      }
#endif      
    } else {			/* single patch assignment */
      printevery  = TRUE;
      start_patch = patch_nr;
      stop_patch  = patch_nr;
      inc = 1;
      if((patch_nr > medium->nrflt - 1)||(patch_nr<0)){
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: error, trying to address patch %i out of %i total\n",
		  i+1,medium->nrflt);
	exit(-1);
      }
#ifdef ALLOW_NON_3DQUAD_GEOM
      /* only for triangular elements do we have to rotate into
	 a local cooridnate system */
      if(slip_type_bc(bc_code) && 
	 (fault[patch_nr].type==TRIANGULAR))
	rotate_to_local=TRUE;
      else
	rotate_to_local=FALSE;	
#endif      
    }
    // message for mixed BCs
    if((!slip_type_bc(bc_code)) && (slip_bc_assigned)){
      HEADNODE
	fprintf(stderr,"read_boundary_conditions: mixed bcs, final stress value of %g will take pre-stressing by slip into account\n",
		bc_value);
    }
    if(rotate_to_local){
      if(fscanf(in,TWO_CP_FORMAT,&global_strike,&global_dip) != 2){
	HEADNODE{
	  fprintf(stderr,"read_boundary_conditions: expected additional strike and dip specification for slip of triangular patch(es)\n");
	  fprintf(stderr,"read_boundary_conditions: and could not read it\n");
	}
	exit(-1);
      }
      HEADNODE{
	fprintf(stderr,"read_boundary_conditions: since at least one patch is triangular, expected strike and dip\n");
	fprintf(stderr,"read_boundary_conditions: for global orientation for slip, those were read as s: %g and d: %g\n",
		global_strike,global_dip);
      }
      // this is now alpha
      global_alpha_rad = DEG2RADF((90.0-global_strike));
      global_dip_rad   = DEG2RADF(global_dip);
      // 
      my_sincos(&sin_global_alpha_rad,&cos_global_alpha_rad,global_alpha_rad);
      my_sincos(&sin_global_dip_rad,&cos_global_dip_rad,global_dip_rad);
      //
      calc_base_vecs(gstrike, gnormal, gdip,
		     sin_global_alpha_rad, cos_global_alpha_rad,
		     sin_global_dip_rad,   cos_global_dip_rad);
#ifdef SUPER_DEBUG
      HEADNODE
	fprintf(stderr," vec: s: (%10.3e,%10.3e,%10.3e) d: (%10.3e,%10.3e,%10.3e) n: (%10.3e,%10.3e,%10.3e)\n",
		gstrike[INT_X],gstrike[INT_Y],gstrike[INT_Z],gdip[INT_X],gdip[INT_Y],gdip[INT_Z],
		gnormal[INT_X],gnormal[INT_Y],gnormal[INT_Z]);
#endif
    } /* end rotate_to_local part */
    /* 

       now assign to temporary arrays since we have sorted by faults
       and added multiple activations

    */
    for(patch_nr3 = start_patch*3,patch_nr=start_patch;
	patch_nr <= stop_patch;
	patch_nr += inc, patch_nr3 += 3){
#ifdef ALLOW_NON_3DQUAD_GEOM
      if(patch_is_2d(fault[patch_nr].type)){
	//
	// check for disallowed codes in 2-D
	//
	if((bc_code == DIP)||(bc_code == COULOMB_DIP_SLIP_UPWARD)||
	   (bc_code == COULOMB_DIP_SLIP_DOWNWARD) || 
	   (bc_code == COULOMB_DIP_SLIP) || (bc_code == DIP_SLIP_UPWARD) ||
	   (bc_code == DIP_SLIP_DOWNWARD) || (bc_code == DIP_SLIP) ||
	   (bc_code == SHEAR_STRESS_FREE) || (bc_code == SHEAR_STRESS_FREE_STRIKE_ONLY) 
	   || (bc_code == SHEAR_STRESS_FRICTION)){
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: segment %i is 2-D, dip motion not allowed\n",
		    patch_nr);
	  exit(-1);
	}
      }
#endif
      if(patch_nr > (medium->nrflt - 1)){
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: file %s, line %i, fault number %i out of range (max nr flt: %i)\n",
		  BC_FILE,n,patch_nr,medium->nrflt);
	exit(-1);
      }
      switch(bc_code){
	/* 
	   specified slip boundary conditions
	*/
      case DIP:
#ifdef ALLOW_NON_3DQUAD_GEOM
	if(patch_is_2d(fault[patch_nr].type)){
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: error, dip slip assigned to 2-D patch %i (use normal)\n",
		    patch_nr);
	  exit(-1);
	}
#endif
      case STRIKE:
      case NORMAL:{
	if(!slip_bc_assigned)
	  slip_bc_assigned = TRUE;
	if(stress_bc_assigned){
	  HEADNODE{
	    fprintf(stderr,"read_boundary_conditions: error patch %i: slip BC assigned after stress BCs\n",
		    patch_nr);
	    fprintf(stderr,"read_boundary_conditions: slip has to be fully specified first for mixed boundary conditions\n");
	  }
	  exit(-1);
	}
	if(!rotate_to_local){
	  //
	  // simply assign the boundary condition value and mark it as a
	  // quake, ie. add to fault's own slip and calculate the stress
	  // effect on the slipping fault and all other faults
	  //
	  // initialize 3-D slip vectors and activation flags
	  for(i=0;i < 3;i++)
	    if(i == bc_code){
	      slip[i] = bc_value;
	      sma_local[i] = TRUE;
	    }else{
	      slip[i] = 0.0;
	      sma_local[i] = FALSE;
	    }
	  if(init_system)	/* this should be TRUE */
	    quake(sma_local,slip,patch_nr,fault,medium,TRUE,FALSE);
	  else			/* else, just add to slip array */
	    add_b_to_a_vector(fault[patch_nr].u,slip,3);
	  if(printevery)
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: const: %5i, flts %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
		      n,start_patch,stop_patch,bc_code,bc_value,
		      comment_on_code_bc(bc_code,bc_value));
	}else{
	  // have to project to local coordinate system, this is the
	  // case for triangular elements
	  if(bc_code == STRIKE){// global strike component
	    slip[STRIKE] = 
	      project_vector(fault[patch_nr].t_strike,gstrike)*bc_value;
	    slip[DIP]    = 
	      project_vector(fault[patch_nr].t_dip,   gstrike)*bc_value;
	    slip[NORMAL] = 
	      project_vector(fault[patch_nr].normal,  gstrike)*bc_value;
	  }else if(bc_code == DIP){// global dip component
	    slip[STRIKE] = 
	      project_vector(fault[patch_nr].t_strike,gdip)*bc_value;
	    slip[DIP]    = 
	      project_vector(fault[patch_nr].t_dip,   gdip)*bc_value;
	    slip[NORMAL] = 
	      project_vector(fault[patch_nr].normal,  gdip)*bc_value;
	  }else{// global normal component
	    slip[STRIKE] = 
	      project_vector(fault[patch_nr].t_strike,gnormal)*bc_value;
	    slip[DIP]    = 
	      project_vector(fault[patch_nr].t_dip,   gnormal)*bc_value;
	    slip[NORMAL] = 
	      project_vector(fault[patch_nr].normal,  gnormal)*bc_value;
	  }
	  for(i=0;i<3;i++)
	    if(fabs(slip[i]) > EPS_COMP_PREC)
	      sma_local[i] = TRUE;
	    else
	      sma_local[i] = FALSE;
	  if(init_system)	/* should be TRUE */
	    quake(sma_local,slip,patch_nr,fault,medium,TRUE,FALSE);
	  else			/* just add to slip array */
	    add_b_to_a_vector(fault[patch_nr].u,slip,3);
#ifdef DEBUG
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions:  patch %i (s: %g, d: %g) cum. rot. u: s:%g d:%g n:%g\n",
		    patch_nr,90.0-RAD2DEGF(asin(fault[patch_nr].sin_alpha)),
		    fault[patch_nr].dip,
		    fault[patch_nr].u[STRIKE],fault[patch_nr].u[DIP],
		    fault[patch_nr].u[NORMAL]);
#endif
	}
	break;
      }
      /* 

	 STRESS BOUNDARY CONDITIONS 


	 will be assembled in a huge array
	 first, then assigned to rhs. later so we can keep up with a
	 general ordering scheme as in rupture.c
	 
	 we will also have to assign fault slip modes 

      */
      case COULOMB_STRIKE_SLIP_LEFTLATERAL:
      case COULOMB_STRIKE_SLIP_RIGHTLATERAL:
      case COULOMB_STRIKE_SLIP:
      case STRIKE_SLIP_LEFTLATERAL:
      case STRIKE_SLIP_RIGHTLATERAL:
      case STRIKE_SLIP:{
	if(!stress_bc_assigned)
	  stress_bc_assigned = TRUE;
	if(sma[patch_nr3+STRIKE]){
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: stress component %i on fault %i was already constrained\n",
		    bc_code,patch_nr);
	  exit(-1);
	} 
	// assign reight hand side and activate slipping mode
	rhs[patch_nr3+STRIKE] = bc_value - fault[patch_nr].s[STRIKE];
	sma[patch_nr3+STRIKE] = ACTIVATED;
#ifdef SUPER_DEBUG
	fprintf(stderr,"read_boundary_conditions: patch: %5i s_strike bc: %12.4e background: %12.4e target: %12.4e\n",
		patch_nr,bc_value,fault[patch_nr].s[STRIKE],rhs[patch_nr3+STRIKE]);
#endif
	//
	if(bc_code > OS_C_OFFSET){// need to correct for normal stress changes
	  if(fault[patch_nr].mu_s == 0.0){
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: for friction adjustment, need mu_s first\n");
	    exit(-1);
	  }
	  if(rhs[patch_nr3+STRIKE]>0)
	    fault[patch_nr].cf[STRIKE]=  (COMP_PRECISION)fault[patch_nr].mu_d;
	  else
	    fault[patch_nr].cf[STRIKE]= -(COMP_PRECISION)fault[patch_nr].mu_d;
	}
	fault[patch_nr].mode[STRIKE]=(MODE_TYPE)bc_code;
	if(printevery)
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: const: %5i, flts %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
		    n,start_patch,stop_patch,bc_code,bc_value,comment_on_code_bc(bc_code,bc_value));
	break;
      }
      case COULOMB_DIP_SLIP_UPWARD:
      case COULOMB_DIP_SLIP_DOWNWARD:
      case COULOMB_DIP_SLIP:
      case DIP_SLIP_UPWARD:
      case DIP_SLIP_DOWNWARD:
      case DIP_SLIP:{
	if(!stress_bc_assigned)
	  stress_bc_assigned = TRUE;
	if(sma[patch_nr3+DIP]){
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: stress component %i on fault %i was already constrained\n",
		    bc_code-3,patch_nr);
	  exit(-1);
	}
#ifdef ALLOW_NON_3DQUAD_GEOM
	if(patch_is_2d(fault[patch_nr].type)){
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: error, dip stress assigned to 2-D patch %i\n",patch_nr);
	  exit(-1);
	}
#endif
	// assign right hand side and activate slipping mode
	rhs[patch_nr3+DIP] = bc_value - fault[patch_nr].s[DIP];
	sma[patch_nr3+DIP] = ACTIVATED;
#ifdef SUPER_DEBUG
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: patch: %5i s_dip    bc: %12.4e background: %12.4e target: %12.4e\n",
		  patch_nr,bc_value,fault[patch_nr].s[DIP],rhs[patch_nr3+DIP]);
#endif
	//
	if(bc_code > OS_C_OFFSET){// need to correct for normal stress changes
	  if(fault[patch_nr].mu_s == 0.0){
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: for friction adjustment, need mu_s first\n");
	    exit(-1);
	  }
	  if(rhs[patch_nr3+DIP]>0)
	    fault[patch_nr].cf[DIP]=  (COMP_PRECISION)fault[patch_nr].mu_d;
	  else
	    fault[patch_nr].cf[DIP]= -(COMP_PRECISION)fault[patch_nr].mu_d;
	}
	fault[patch_nr].mode[DIP]=(MODE_TYPE)bc_code;
	if(printevery)
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: const: %5i, flts %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
		    n,start_patch,stop_patch,bc_code,bc_value,comment_on_code_bc(bc_code,bc_value));
	break;
      }
      case NORMAL_SLIP_OUTWARD:
      case NORMAL_SLIP_INWARD:
      case NORMAL_SLIP:{
	if(!stress_bc_assigned)
	  stress_bc_assigned = TRUE;
#ifdef NO_OPENING_MODES
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: normal mode slipping off, how did we get this far?\n");
	exit(-1);
#else
	if(sma[patch_nr3+NORMAL]){
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: stress component %i on fault %i was already constrained\n",
		    bc_code-3,patch_nr);
	  exit(-1);
	}
	// assign right hand side and activate slipping mode
	sma[patch_nr3+NORMAL] = ACTIVATED;
	rhs[patch_nr3+NORMAL] = bc_value - fault[patch_nr].s[NORMAL];
#ifdef SUPER_DEBUG
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: patch: %5i s_normal bc: %12.4e background: %12.4e target: %12.4e\n",
		  patch_nr,bc_value,fault[patch_nr].s[NORMAL],rhs[patch_nr3+NORMAL]);
#endif
	//
	fault[patch_nr].mode[NORMAL]=(MODE_TYPE)bc_code;
	if(printevery)
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: const: %5i, flts %5i to %5i, bcode: %4i, value: %9.6f, %s\n",
		    n,start_patch,stop_patch,bc_code,bc_value,
		    comment_on_code_bc(bc_code,bc_value));
#endif
	break;
      }
      case SHEAR_STRESS_FREE:
      case SHEAR_STRESS_FREE_STRIKE_ONLY:
      case SHEAR_STRESS_FRICTION:{
	/*

	  stress boundary conditions that are based on the resolved background stress

	  (we adjust for rpeexisting stresses (from slip) at the end!)

	*/
	if(!stress_bc_assigned)
	  stress_bc_assigned = TRUE;
	
	if(sma[patch_nr3+DIP]){
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: patch %i: dip component was already constrained, bc_code: %i does not work\n",
		    patch_nr,bc_code);
	  exit(-1);
	}
	if(sma[patch_nr3+STRIKE]){
	  HEADNODE
	    fprintf(stderr,"read_boundary_conditions: patch %i: strike component was already constrained, bc_code: %i does not work\n",
		    patch_nr,bc_code);
	  exit(-1);
	} 
	//
	// decide on the way that resolved stresses are output for debugging
	//
	if(!check_for_res_stress_output){
	  if(stop_patch - start_patch > 20){
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: too many patches for stderr, resolved stresses (strike, dip, normal) written to \"%s\" only\n",
		      RES_STRESS_FILE);
	    print_res_stress_screen = FALSE;
	  }else{
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: output of resolved stresses (strike, dip, normal) to stderr and \"%s\"\n",
		      RES_STRESS_FILE);
	    print_res_stress_screen = TRUE;
	  }
	  rsout = myopen(RES_STRESS_FILE,"w");
	  fprintf(rsout,"# resolved stresses from BC code %i for each patch: s_strike s_dip s_normal\n",bc_code);
	  check_for_res_stress_output = TRUE;
	}
	/*
	  
	  achieve shear stress free or friction stress by slipping 
	  in dip and strike direction so that the resolved background stress is 
	  reduced to zero
	  
	*/
	// time is 0
	background_stress(sm,fault[patch_nr].x,0.0,a,b,medium->pressure);
	// determine stress on patch
	calc_three_stress_components(sm, fault[patch_nr].normal,
				     fault[patch_nr].t_strike,
				     fault[patch_nr].t_dip,fault[patch_nr].normal,
				     &bstress[STRIKE],&bstress[DIP],
				     &bstress[NORMAL]);
	//
	// activate the strike and dip slip modes for 3-D, only strike
	// for 2-D
	//
	if(bc_code == SHEAR_STRESS_FREE_STRIKE_ONLY){
	  use_dip = FALSE;	/* only strike slip motion */
	}else{
#ifdef ALLOW_NON_3DQUAD_GEOM
	  if(patch_is_2d(fault[patch_nr].type))
	    use_dip = FALSE;
	  else
	    use_dip = TRUE;
#else
	  use_dip = TRUE;
#endif
	}
	sma[patch_nr3 + STRIKE] = ACTIVATED;
	fault[patch_nr].mode[STRIKE] = (MODE_TYPE) STRIKE_SLIP;
	if(use_dip){
	  sma[patch_nr3 + DIP]   = ACTIVATED;
	  fault[patch_nr].mode[DIP]   = (MODE_TYPE) DIP_SLIP;
	}
	if((bc_code ==  SHEAR_STRESS_FREE)||(bc_code == SHEAR_STRESS_FREE_STRIKE_ONLY)){
	  //
	  // activate both dip and strike direction (for 3-D, if allowed) and 
	  // slip to shear stress free boundary condition
	  //
	  rhs[patch_nr3+STRIKE]= -bstress[STRIKE];
	  if(use_dip)
	    rhs[patch_nr3+DIP]=  -bstress[DIP];
	}else{
	  //
	  // slip according to friction criterion
	  // 
	  if(bstress[NORMAL] > 0){
	    //
	    // the fault is under extension
	    //
	    HEADNODE
	      fprintf(stderr,"read_boundary_conditions: WARNING: patch %i: under tension, shear stress free conditions used instead of %i\n",
		      patch_nr,bc_code);
	    rhs[patch_nr3+STRIKE]= -bstress[STRIKE];
	    if(use_dip)
	      rhs[patch_nr3+DIP]=  -bstress[DIP];
	  }else{
	    //
	    // under compression, obtain strike and dip stress values
	    //

	    stress_drop=((COMP_PRECISION)fault[patch_nr].mu_s)*
			 (-bstress[NORMAL]);

	    if(use_dip){/* decide how to distribute stress in striek and dip 
			   directions */
	      get_maxsdir_stress_drops2(bstress,stress_drop,(rhs+patch_nr3));
	    }else{// no dip motion, that is equivalent to fi = -1
	      rhs[patch_nr3+STRIKE] = (bstress[STRIKE] > 0)?
		(-stress_drop):(stress_drop);
	    }
	  }
	}
	//
	// adjust for preexisting stress
	//
	if(sma[patch_nr3+STRIKE])
	  rhs[patch_nr3+STRIKE] -= fault[patch_nr].s[STRIKE];
	if(sma[patch_nr3+DIP])
	  rhs[patch_nr3+DIP] -= fault[patch_nr].s[DIP];
	//
#ifdef SUPER_DEBUG
	HEADNODE{
	  fprintf(stderr,"read_boundary_conditions: patch: %5i res stress strike:  background: %12.4e target: %12.4e\n",
		  patch_nr,fault[patch_nr].s[STRIKE],rhs[patch_nr3+STRIKE]);
	  fprintf(stderr,"read_boundary_conditions: patch: %5i res stress dip:     background: %12.4e target: %12.4e\n",
		  patch_nr,fault[patch_nr].s[DIP],rhs[patch_nr3+DIP]);
	}
#endif
	HEADNODE{
	  if(print_res_stress_screen){
	    if(use_dip){
	      fprintf(stderr,"read_boundary_conditions: bc: %3i resolved s: s/d/n: %9.2e/%9.2e/%9.2e target s: s/d: %9.2e/%9.2e on patch %3i %s\n",
		      bc_code,bstress[STRIKE],bstress[DIP],bstress[NORMAL],
		      rhs[patch_nr3+STRIKE],rhs[patch_nr3+DIP],patch_nr,
		      ((hypot(rhs[patch_nr3+STRIKE],rhs[patch_nr3+DIP])<EPS_COMP_PREC)?
		       ("(amp = zero!)"):("")));
	    }else{
	      fprintf(stderr,"read_boundary_conditions: bc: %3i resolved s: s/d/n: %9.2e/%9.2e/%9.2e target s: strike: %9.2e on patch %3i %s\n",
		      bc_code,bstress[STRIKE],bstress[DIP],bstress[NORMAL],
		      rhs[patch_nr3+STRIKE],patch_nr,
		      ((hypot(rhs[patch_nr3+STRIKE],rhs[patch_nr3+DIP])<EPS_COMP_PREC)?
		       ("(amp = zero!)"):("")));
	    }
	  }
	  if(check_for_res_stress_output)
	    fprintf(rsout,"%12.5e %12.5e %12.5e\n",bstress[STRIKE],
		    bstress[DIP],bstress[NORMAL]);
	}
	break;
      }	/* end of shear stress friction */
      default:{
	HEADNODE
	  fprintf(stderr,"read_boundary_conditions: file %s, line %i, bc code %i undefined \n",
		  BC_FILE,n,bc_code);
	exit(-1);
      }}			/* end of stress BC cases */
    }
  }
  if(!n){
    HEADNODE
      fprintf(stderr,"read_boundary_conditions: error, couldn't find any boundary conditions specified\n");
    exit(-1);
  }
  if(check_for_res_stress_output)// close the resolved stress file
    fclose(rsout);


  /* 

     ASSIGN STRESS CONDITIONS TO MATRIX in the same 
     fashion as in rupture.c, first initialize the real arrays


  */
  if(init_system){
    init_equation_system(medium,fault);
    /*
      
    ASSIGN BOUNDARY CONDITIONS and decide if we want a constrained or
    unconstrained (positive only or positive and negative) solution
    for slip
    
    */
    for(i=i3=0;i<medium->nrflt;i++,i3+=3){
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
	  if(sma[i3+j]){
	    medium->sma[medium->naflt*3+j] = sma[i3+j];
	    // the x vector is in here just so that it gets resized at the same time
	    // as b
	    add_to_right_hand_side(rhs[i3+j],&medium->b,&medium->xsol,
				   &medium->nreq);
	    added_to_uncon=TRUE;/* have to increment fault contraint counter */
	  }
	  break;
	}
	  /* 
	     now CONSTRAINED equations, slip of specific direction 
	     was specified 
	  */
	default:{
	  if(sma[i3+j]){
	    medium->sma_con[medium->naflt_con*3+j]=sma[i3+j];
	    /* 
	       
	    check for strange settings that would lead to funny 
	    stuff  without interaction-modified stresses 
	    
	    */
	    if((((fault[i].mode[j]==STRIKE_SLIP_LEFTLATERAL) ||  
		 (fault[i].mode[j]==COULOMB_STRIKE_SLIP_LEFTLATERAL))
		&& (rhs[i3+j]>0))||
	       (((fault[i].mode[j]==STRIKE_SLIP_RIGHTLATERAL) || 
		 (fault[i].mode[j]==COULOMB_STRIKE_SLIP_RIGHTLATERAL))
		&& (rhs[i3+j]<0))||
	       (((fault[i].mode[j]==DIP_SLIP_UPWARD ) ||         
		 (fault[i].mode[j]==COULOMB_DIP_SLIP_UPWARD))
		&& (rhs[i3+j]>0))||
	       (((fault[i].mode[j]==DIP_SLIP_DOWNWARD) ||
		 (fault[i].mode[j]==COULOMB_DIP_SLIP_DOWNWARD))
		&& (rhs[i3+j]<0))||
	       ((fault[i].mode[j] == NORMAL_SLIP_OUTWARD)      
		&& (rhs[i3+j] < 0))||
	       ((fault[i].mode[j] == NORMAL_SLIP_INWARD)       
		&& (rhs[i3+j] > 0))){
	      HEADNODE
		fprintf(stderr,
			"read_boundary_conditions: WARNING: fault: %5i bc_code: %3i dir: %1i might be incompatible with stress: %g\n",
			i,fault[i].mode[j],j,rhs[i3+j]);
	    }
	    if(fault[i].mode[j]==STRIKE_SLIP_RIGHTLATERAL || 
	       fault[i].mode[j]==COULOMB_STRIKE_SLIP_RIGHTLATERAL ||
	       fault[i].mode[j]==DIP_SLIP_DOWNWARD || 
	       fault[i].mode[j]==COULOMB_DIP_SLIP_DOWNWARD  ||
	       fault[i].mode[j]==NORMAL_SLIP_INWARD)
	      /* 
		 assign negative b for constrained, 'anormal' modes 
	      */
	      // add to RHS and resize b and x
	      add_to_right_hand_side(-rhs[i3+j],&medium->b_con,
				     &medium->xsol_con,&medium->nreq_con);
	    else
	      add_to_right_hand_side( rhs[i3+j],&medium->b_con,
				      &medium->xsol_con,&medium->nreq_con);
	    // have to increment counter
	    added_to_con=TRUE;
	  }
	  break;
	}}
      }
#ifdef SUPER_DEBUG
      HEADNODE{
	if(sma[i3] || sma[i3+1] || sma[i3+2])
	  fprintf(stderr,"read_boundary_conditions:    flt %5i: activ. codes: %i %i %i rhs: %11.4e %11.4e %11.4e\n",
		  i,sma[i3+0],sma[i3+1],sma[i3+2],rhs[i3+0],rhs[i3+1],rhs[i3+2]);
	if(medium->sma[medium->naflt*3] || 
	   medium->sma[medium->naflt*3+1] || 
	   medium->sma[medium->naflt*3+2])
	  fprintf(stderr,"read_boundary_conditions: uc flt %5i: medium a.cd.: %i %i %i\n",
		  medium->naflt,
		  medium->sma[medium->naflt*3],
		  medium->sma[medium->naflt*3+1],
		  medium->sma[medium->naflt*3+2]);
	if(medium->sma_con[medium->naflt_con*3] ||
	   medium->sma_con[medium->naflt_con*3+1] ||
	   medium->sma_con[medium->naflt_con*3+2])
	  fprintf(stderr,"read_boundary_conditions:  c flt %5i: medium a.cd.: %i %i %i\n",
		  medium->naflt_con,
		  medium->sma_con[medium->naflt_con*3],
		  medium->sma_con[medium->naflt_con*3+1],
		  medium->sma_con[medium->naflt_con*3+2]);
      }
#endif    
      // step up the active fault counters and sma arrays
      if(added_to_uncon)
	add_to_active_fault_list(i,&medium->nameaf,&medium->naflt,
				 &medium->sma);
      if(added_to_con)
	add_to_active_fault_list(i,&medium->nameaf_con,&medium->naflt_con,
				 &medium->sma_con);
    }
    HEADNODE{
      if(medium->nreq_con + medium->nreq > 0){
	fprintf(stderr,
		"read_boundary_conditions: total of %i constrained and %i unconstrained equations\n",
		medium->nreq_con,medium->nreq);
	fprintf(stderr,
		"read_boundary_conditions: with     %i constrained and %i unconstrained faults\n",
		medium->naflt_con,medium->naflt);
      }
    }

  }else{
    HEADNODE
      fprintf(stderr,"initialize: read_one_step_bc: WARNING: skipping EQ system init, assigning RHS to fault[].s[]\n");
    /* instead, assign the right hand sides to the stress components */
    for(i=0;i<medium->nrflt;i++)
      for(j=0;j<3;j++)
	fault[i].s[j] = rhs[i*3+j];
  }
  free(rhs);free(sma);
}
/*

  decide if one-step boundary condition is of
  prescribed slip type (as opposed to stress)

 */
my_boolean slip_type_bc(int bc_code)
{
  if((bc_code == STRIKE) || (bc_code == DIP) || 
     (bc_code == NORMAL))
    return(TRUE);
  else
    return(FALSE);
}
