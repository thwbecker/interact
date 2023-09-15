/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu


  main initialization routine, sets defaults and optional parameters 
  from the command line
  
  calls initialize, which reads in geometry and boundary conditions

*/
#include "interact.h"
#include "properties.h"

void check_parameters_and_init(int argc, char **argv,
			       struct med **medium,
			       struct flt **fault,
			       my_boolean *read_initial_fault_stress,
			       COMP_PRECISION *a,COMP_PRECISION *b)
{
  int max_nr_flt_files;
  my_boolean read_fault_properties,suppress_interactions,whole_fault_mode,
    read_stress_relation_factors,use_slip_files,whole_fault_deactivations,
    use_sparse_storage,use_old_imat,use_old_amat,save_imat,save_amat,
    check_for_interaction_feedback,keep_slipping,attempt_restart,suppress_nan_output,
    geomview_output,twod_approx_is_plane_stress,print_plane_coord,half_plane,
    variable_time_step,debug,no_interactions,force_petsc;
  short int solver_mode;
  COMP_PRECISION pressure,med_cohesion,min_stress_drop,wcutoff;
  I_MATRIX_PREC i_mat_cutoff;
  if((*medium)->comm_rank==0){
    if((*medium)->comm_size==1)
      fprintf(stderr,"init: compiled on %s %s, running in serial\n",__DATE__,__TIME__);
    else
      fprintf(stderr,"init: compiled on %s %s, running on %i cores\n",__DATE__,__TIME__,(*medium)->comm_size);
  }
  // initialization phase, get parameters from command
  // line options
  init_parameters(argv,argc, &read_fault_properties,
		  read_initial_fault_stress,
		  &suppress_interactions,&whole_fault_mode,
		  &med_cohesion,&max_nr_flt_files,
		  &read_stress_relation_factors,&use_slip_files,
		  &whole_fault_deactivations,
		  &min_stress_drop,&use_sparse_storage,
		  &i_mat_cutoff,&use_old_imat,&use_old_amat,
		  &save_imat,&save_amat,&check_for_interaction_feedback,
		  &keep_slipping,&attempt_restart,&solver_mode,
		  &suppress_nan_output,&geomview_output,&pressure,
		  &twod_approx_is_plane_stress,&print_plane_coord,
		  &half_plane,&variable_time_step,&debug,&wcutoff,
		  &no_interactions,&force_petsc,(*medium)->comm_rank);
  // load files, etc
  initialize(medium,fault,read_fault_properties,max_nr_flt_files,
	     suppress_interactions,whole_fault_mode,med_cohesion,a,b,
	     read_stress_relation_factors,use_slip_files,whole_fault_deactivations,
	     min_stress_drop,use_sparse_storage,i_mat_cutoff,use_old_imat,
	     use_old_amat,save_imat,save_amat,check_for_interaction_feedback,
	     keep_slipping,attempt_restart,solver_mode,suppress_nan_output,
	     geomview_output,pressure,twod_approx_is_plane_stress,
	     print_plane_coord,half_plane,variable_time_step,debug,TRUE,wcutoff,
	     no_interactions,force_petsc);
}
/*

  read geometry and boundary conditions

  write default parameters to the model structure

  open files

  call the interaction matrix part


  init_system should normally be TRUE

*/
void initialize(struct med **medium, struct flt **fault,
		my_boolean read_fault_properties,int max_nr_flt_files,
		my_boolean suppress_interactions,my_boolean whole_fault_mode,
		COMP_PRECISION med_cohesion,COMP_PRECISION *a,COMP_PRECISION *b,
		my_boolean read_stress_relation_factors,my_boolean use_slip_files,
		my_boolean whole_fault_deactivations,COMP_PRECISION min_stress_drop,
		my_boolean use_sparse_storage,I_MATRIX_PREC i_mat_cutoff,
		my_boolean use_old_imat,my_boolean use_old_amat,
		my_boolean save_imat,my_boolean save_amat,
		my_boolean check_for_interaction_feedback,my_boolean keep_slipping,
		my_boolean attempt_restart,short int solver_mode,
		my_boolean suppress_nan_output,my_boolean geomview_output,
		COMP_PRECISION pressure,my_boolean twod_approx_is_plane_stress,
		my_boolean print_plane_coord,my_boolean half_plane,
		my_boolean variable_time_step,my_boolean debug,
		my_boolean init_system,COMP_PRECISION wcutoff,
		my_boolean no_interactions,my_boolean force_petsc)
{
  int serr;
  char tmpstring[STRLEN];
  //
  // read in the main model description files
  //
  // first, need the FAULT GEOMETRY, and calculate the position of patches
  // also, allocate space and initialize the medium structure
  //
  
  read_geometry(GEOMETRY_FILE,medium,fault,read_fault_properties,
		twod_approx_is_plane_stress,half_plane,TRUE);
  if((*medium)->comm_rank==0)
    fprintf(stderr,"initialize: all stress values are based on a shear modulus of %g\n",
	    SHEAR_MODULUS);
  // assign the background pressure 
  (*medium)->pressure = pressure;
  // read in the background stressing factors, 
  read_stress_fac(read_stress_relation_factors,a,b,pressure,*medium);
  /* check for time stepping before reading in BCs */
  (*medium)->variable_time_step = variable_time_step;
  if(!(*medium)->variable_time_step)
    fprintf(stderr,"initialize: WARNING: using constant time steps\n");

  if(((*medium)->debug = debug)){
    fprintf(stderr,"initialize: running in debugging mode, output in directory %s\n",INT_TMP_DIR);
    snprintf(tmpstring,STRLEN,"mkdir -p %s",INT_TMP_DIR);
    serr = system(tmpstring);
    if(serr){
      fprintf(stderr,"initialize: could not make temporary output dir %s\n",INT_TMP_DIR);
      exit(-1);
    }
  }
  //
  // read in boundary conditions, e.g. what kind of simulations, if
  // one step, prescribed slip or stress etc.
  //
  read_boundary_conditions(*medium,*fault,a,b,init_system);  
  //
  //
  (*medium)->no_interactions = no_interactions;
  if((*medium)->no_interactions)
    fprintf(stderr,"initialize: WARNING: suppressing all fault-to-fault interactions!\n");
  (*medium)->solver_mode = solver_mode;
  (*medium)->max_nr_flt_files=max_nr_flt_files;
  (*medium)->suppress_interactions = suppress_interactions;
  (*medium)->print_plane_coord = print_plane_coord;
  if((*medium)->suppress_interactions){
    if((*medium)->op_mode == ONE_STEP_CALCULATION){
      fprintf(stderr,"initialize: suppressed interactions are only meant for loading simulation!\n");
      exit(-1);
    }
    else
      fprintf(stderr,"initialize: WARNING: suppressing interactions\n");
  }
  (*medium)->whole_fault_activations=whole_fault_mode;
  (*medium)->whole_fault_deactivations=whole_fault_deactivations;
  if((*medium)->whole_fault_activations){
    /* 

       when we are in whole fault mode, this switch 
       will turn the whole fault deactivations OFF, if set 

    */
    if((*medium)->whole_fault_deactivations)
      (*medium)->whole_fault_deactivations = FALSE;
    else
      (*medium)->whole_fault_deactivations = TRUE;
  }
  (*medium)->keep_slipping = keep_slipping;
  if((*medium)->op_mode != ONE_STEP_CALCULATION){
    // some messages only if we are really in a loading simulation
    if((*medium)->whole_fault_activations)
      fprintf(stderr,"initialize: whole groups are activated when patches are critical\n");
    if((*medium)->whole_fault_deactivations)
      fprintf(stderr,"initialize: whole fault gets deactivated when patch deactivated\n");
    if(EXHAUSTIVE_CRITICAL_STRESS_EPS != CRITICAL_STRESS_EPS){
      fprintf(stderr,"initialize: note that the critical stress threshold is different for exhaustive and iterative approaches\n");
    }
    if((*medium)->keep_slipping)
      fprintf(stderr,"initialize: initially slipping patches remain active\n");
    else
      fprintf(stderr,"initialize: initially slipping patches are exhausted by triggering\n");
  }
  (*medium)->use_slip_files=use_slip_files;
  (*medium)->cohesion=med_cohesion;
#ifdef NO_COHESION
  if((*medium)->cohesion != 0.0){
    fprintf(stderr,"initialize: program was compiled assuming that cohesion is zero\n");
    fprintf(stderr,"initialize: therefore cohesion = %g will not work. Exiting.\n",
	    (*medium->cohesion));
    exit(-1);
  }
#else
  if((*medium)->cohesion == 0.0)
    if((*medium)->op_mode != ONE_STEP_CALCULATION)
      fprintf(stderr,"initialize: you could use the -DNO_COHESION option to recompile for C=0 at all times\n");
#endif

  (*medium)->min_stress_drop = min_stress_drop;
  if((*medium)->min_stress_drop != 0.0){
    if((*medium)->op_mode != ONE_STEP_CALCULATION)
      fprintf(stderr,"initialize: minimum stress drop for loading simulation set to %g\n",
	      (*medium)->min_stress_drop);
  }
  (*medium)->use_sparse_storage = use_sparse_storage;
  if(((*medium)->use_sparse_storage) && ((*medium)->op_mode == ONE_STEP_CALCULATION)){
    fprintf(stderr,"initialize: sparse storage only works for loading simulation!\n");
    exit(-1);
  }
  (*medium)->i_mat_cutoff = i_mat_cutoff;
  (*medium)->use_old_imat = use_old_imat;
  (*medium)->save_imat = save_imat;
  (*medium)->use_old_amat = use_old_amat;
  (*medium)->save_amat = save_amat;
  (*medium)->wcutoff = wcutoff;
  if(((*medium)->op_mode != ONE_STEP_CALCULATION)&&
     ((*medium)->save_amat || (*medium)->use_old_amat)){
    fprintf(stderr,"initialize: saving/restoring the current equations (A) matrix only makes sense for one step calculations.\n");
    fprintf(stderr,"initialize: if you want to store the interaction (I) matrix, use -si and -oi, respectively.");
    exit(-1);
  }
  (*medium)->check_for_interaction_feedback =
    check_for_interaction_feedback;
  /* intervals for slip line printing */
  (*medium)->slip_line_dt = (*medium)->stop_time/10.0;
  // for spatial correlations
  (*medium)->spcorr_interval = SPCORR_INTERVAL_DEF;
  // restart?
  (*medium)->attempt_restart = attempt_restart;
  // output of NaNs in stress or displacement?
  (*medium)->suppress_nan_output = suppress_nan_output;
  // Geomview COFF file output?
  (*medium)->geomview_output = geomview_output;
  if(force_petsc){
#ifdef USE_PETSC
    (*medium)->force_petsc = force_petsc;
#else
    HEADNODE
      fprintf(stderr,"initialize: force_petsc is true, but no Petsc support compiled in\n");
    
#endif
  }
  /*
     initialize files, possibly plotting windows,
     and interaction matrix 
  */
  switch((*medium)->op_mode){
  case SIMULATE_LOADING_AND_PLOT:{// loading simulation with plotting
    init_files(medium,fault);
#ifdef USE_PGPLOT
    (*medium)->nr_xwindows=3; 
    (*medium)->time_tic_interval=1;
    (*medium)->l_time_tic_interval=5*
      (*medium)->time_tic_interval;
    init_plot_window(*medium,*fault); 
#endif
    calc_interaction_matrix(*medium,*fault,FALSE); 
    break;
  }
  case SIMULATE_LOADING:{// loading simulation
    init_files(medium,fault);
    calc_interaction_matrix(*medium,*fault,FALSE); 
    break;
  }
  case ONE_STEP_CALCULATION:{// one step doesn't need I matrix
    break;
  }
  default:{
    fprintf(stderr,"init: opmode %i undefined\n",(*medium)->op_mode);
    exit(-1);
  }}
}


/* 
   initialize files for output 
*/
void init_files(struct med **medium,struct flt **fault)
{
  int i;
  char tmpstring[STRLEN],tmpstring2[STRLEN];
  /* 
     initialize files 
  */
  // first flt.*.dat files
  if((*medium)->nrgrp <= (*medium)->max_nr_flt_files){
    (*medium)->flt_stress_out=
      malloc(sizeof(FILE *)*(*medium)->nrgrp);
    for(i=0;i<(*medium)->nrgrp;i++){
      snprintf(tmpstring,STRLEN,"%s.%i.dat",FAULT_DATA_PREFIX,i);
      (*medium)->flt_stress_out[i]=myopen(tmpstring,"w");
      fprintf((*medium)->flt_stress_out[i],
	      "# fault patch group %i\n",i);
      fprintf((*medium)->flt_stress_out[i],
	      "# time <s_c> <d_m*s_n> <u_s> <u_n> <u_d> <s_s> <s_n> <s_d> min(s_c) max(s_c)\n");
    }
    sprintf(tmpstring,"%s.%i.dat",FAULT_DATA_PREFIX,0);
    sprintf(tmpstring2,"%s.%i.dat",FAULT_DATA_PREFIX,(*medium)->nrgrp-1);
    fprintf(stderr,"init_files: writing to fault group data files \"%s\" ... \"%s\" every %g timesteps\n",
	    tmpstring,tmpstring2,(*medium)->print_interval);
    if((*medium)->variable_time_step)
      fprintf(stderr,"init_files: will use %g as upper limit for variable dt\n",
	      (*medium)->print_interval);
    (*medium)->flt_stress_init=TRUE;
  }
  // slip event files
#ifdef BINARY_PATCH_EVENT_FILE
  fprintf(stderr,"init_files: writing to \"%s\", format: time nr_iter patch_nr slip[s,d,n] moment",EVENT_FILE_BINARY);
  fprintf(stderr,",\ninit_files: in binary format (flt/int/int/flt/flt/flt/flt)\n");
  (*medium)->events_out= myopen(EVENT_FILE_BINARY,"w");
#else
  fprintf(stderr,"init_files: writing to \"%s\", format: time nr_iter patch_nr slip[s,d,n] moment",EVENT_FILE_ASCII);
  fprintf(stderr,", ascii format\n");
  (*medium)->events_out= myopen(EVENT_FILE_ASCII,"w");
#endif
  // cevents.dat file for events per quake
  fprintf(stderr,"init_files: writing to \"%s\", format: time group_nr moment total_moment m_slip_s m_slip_d m_slip_n m_slip_x m_slip_y m_slip_z m_x m_y m_z\n",
	  CEVENT_FILE);
  (*medium)->cevents_out=myopen(CEVENT_FILE,"w");
  (*medium)->events_init=TRUE;

  fprintf(stderr,"init_files: writing to \"%s\", format: time group_nr mead/std_s,d,n\n",
	  STRESS_STAT_FILE);
  (*medium)->stress_stat_out=myopen(STRESS_STAT_FILE,"w");
  
  // slipline files
  if((*medium)->nrgrp <= (*medium)->max_nr_flt_files){
    if((*medium)->use_slip_files){
      if((*medium)->slip_line_dt != 0.0){
	(*medium)->slip_line_out=malloc(sizeof(FILE *)*
					(*medium)->nrgrp);
	for(i=0;i<(*medium)->nrgrp;i++){
	  sprintf(tmpstring,"%s.%i.dat",SLIP_LINE_PREFIX,i);
	  (*medium)->slip_line_out[i]=myopen(tmpstring,"w");
	}
	sprintf(tmpstring,"%s.%i.dat",SLIP_LINE_PREFIX,0);
	sprintf(tmpstring2,"%s.%i.dat",SLIP_LINE_PREFIX,(*medium)->nrgrp-1);
	fprintf(stderr,"init_files: writing to slipline files \"%s\" ... \"%s\" every %g timesteps\n",
		tmpstring,tmpstring2,(*medium)->slip_line_dt);
	(*medium)->slip_line_init=TRUE;
      }
    }else{
      fprintf(stderr,"init_files: slipline.*.dat files suppressed\n");
    }
  }
}


void init_parameters(char **argv, int argc, my_boolean *read_fault_properties,
		     my_boolean *read_initial_fault_stress,
		     my_boolean *suppress_interactions,my_boolean *whole_fault_mode,
		     COMP_PRECISION *med_cohesion,int *max_nr_flt_files,
		     my_boolean *read_stress_relation_factors,
		     my_boolean *use_slip_files,
		     my_boolean *whole_fault_deactivations,
		     COMP_PRECISION *min_stress_drop,
		     my_boolean *use_sparse_storage,
		     I_MATRIX_PREC *i_mat_cutoff,my_boolean *use_old_imat,
		     my_boolean *use_old_amat,my_boolean *save_imat,
		     my_boolean *save_amat,
		     my_boolean *check_for_interaction_feedback,
		     my_boolean *keep_slipping,
		     my_boolean *attempt_restart,short int *solver_mode,
		     my_boolean *suppress_nan_output,my_boolean *geomview_output,
		     COMP_PRECISION *pressure,
		     my_boolean *twod_approx_is_plane_stress,
		     my_boolean *print_plane_coord,
		     my_boolean *half_plane,
		     my_boolean *variable_time_step,
		     my_boolean *debug,
		     COMP_PRECISION *wcutoff,
		     my_boolean *no_interactions,
		     my_boolean *force_petsc,
		     int rank)
{
  int i;
  my_boolean warned = FALSE;
  /* 
     assign default values 
  */
  *read_fault_properties = READ_FAULT_PROPERTIES_DEF;
  *read_initial_fault_stress = READ_INITIAL_FAULT_STRESS_DEF;
  *suppress_interactions = SUPPRESS_INTERACTIONS_DEF;
  *whole_fault_mode = WHOLE_FAULT_MODE_DEF;
  *med_cohesion = COHESION_DEF;
  *max_nr_flt_files = MAX_NR_FLT_FILES_DEF;
  *read_stress_relation_factors = READ_STRESS_RELATION_FACTORS_DEF;
  *use_slip_files = PRINT_SLIPLINE_DEF;
  *whole_fault_deactivations = WHOLE_FAULT_DEACTIVATIONS_DEF;
  *use_sparse_storage = USE_SPARSE_STORAGE_DEF;
  *min_stress_drop = 0.0;
  *i_mat_cutoff = I_MAT_CUTOFF_DEF;
  *use_old_imat = USE_OLD_IMAT_DEF;
  *save_imat = SAVE_IMAT_DEF;
  *use_old_amat = USE_OLD_AMAT_DEF;
  *save_amat = SAVE_AMAT_DEF;
  *check_for_interaction_feedback = CHECK_FOR_INTERACTION_FEEDBACK_DEF;
  *keep_slipping = KEEP_SLIPPING_DEF;
  *attempt_restart = FALSE;
  *solver_mode = LU_SOLVER;// default solver
  *suppress_nan_output = SUPPRESS_NAN_OUTPUT_DEF;
  *geomview_output = GEOMVIEW_OUTPUT_DEF;
  *pressure = PRESSURE_DEF;
  *twod_approx_is_plane_stress = TWO_DIM_APPROX_IS_PLANE_STRESS_DEF;
  *print_plane_coord = PRINT_PLANE_COORD_DEF;
  *half_plane = HALF_PLANE_DEF;
  *variable_time_step = !(CONSTANT_TIME_STEP_DEF);
  *debug = DEBUG_DEF;
  *wcutoff = SVD_THRESHOLD;
  *no_interactions = FALSE;
  *force_petsc = FALSE;
  /* 
     check for input options 
  */
  for(i=1;i<argc;i++){
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){// help
      if(rank==0)
	phelp();
      exit(-1);
    }else if(strcmp(argv[i],"-f")==0){// fault prop file
      toggle(read_fault_properties);
    }else if(strcmp(argv[i],"-r")==0){// attempt a restart
      *attempt_restart = TRUE;
    }else if(strcmp(argv[i],"-ni")==0){// no interactions
      *no_interactions = TRUE;
    }else if(strcmp(argv[i],"-s")==0){// fault stress init file
      toggle(read_initial_fault_stress);
    }else if(strcmp(argv[i],"-ps")==0){// plane strain/stress approximation
      toggle(twod_approx_is_plane_stress);
    }else if(strcmp(argv[i],"-hp")==0){	/* half-plane for 2-D */
      toggle(half_plane);
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
    }else if(strcmp(argv[i],"-sa")==0){// save A matrix on file?
      toggle(save_amat);
    }else if(strcmp(argv[i],"-sn")==0){// output of NaNs?
      toggle(suppress_nan_output);
    }else if(strcmp(argv[i],"-ci")==0){// check for feedback interaction?
      toggle(check_for_interaction_feedback);
    }else if(strcmp(argv[i],"-v")==0){// whole fault de-activation mode
      toggle(whole_fault_deactivations);
    }else if(strcmp(argv[i],"-d")==0){// debugging mode
      toggle(debug);
    }else if(strcmp(argv[i],"-l")==0){// read in a/b factors for stress relation
      toggle(read_stress_relation_factors);
    }else if(strcmp(argv[i],"-p")==0){// print to the slipline files
      toggle(use_slip_files);
    }else if(strcmp(argv[i],"-pc")==0){// print plane coordinates
      toggle(print_plane_coord);
    }else if(strcmp(argv[i],"-gv")==0){// Geomview output
      toggle(geomview_output);
    }else if(strcmp(argv[i],"-ct")==0){// time step
      toggle(variable_time_step);
    }else if(strcmp(argv[i],"-sv")==0){// use SVD solver for Ax=b
      *solver_mode = SVD_SOLVER;
    }else if(strcmp(argv[i],"-sl")==0){// use LU solver for Ax=b
      *solver_mode = LU_SOLVER;
    }else if(strcmp(argv[i],"-ss")==0){// use sparse solver for Ax=b
      *solver_mode = SPARSE_SOLVER;
    }else if(strcmp(argv[i],"-k")==0){// keep initial patches slipping during rupture
      toggle(keep_slipping);
    }else if(strcmp(argv[i],"-c")==0){// cohesion
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,med_cohesion);
#ifdef NO_COHESION
      if(rank == 0){
	fprintf(stderr,"init_parameters program was compiled with NO_COHESION option\n");
	fprintf(stderr,"init_parameters therefore, cohesion is always zero and -c does not make sense\n");
      }
      exit(-1);
#endif
    }else if(strcmp(argv[i],"-ms")==0){// minimum stress drop
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,min_stress_drop);
    }else if(strcmp(argv[i],"-pr")==0){// pressure
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,pressure);
    }else if(strcmp(argv[i],"-oi")==0){// use old i matrix
      toggle(use_old_imat);
    }else if(strcmp(argv[i],"-oa")==0){// use old A matrix
      toggle(use_old_amat);
    }else if(strcmp(argv[i],"-ei")==0){// i matrix cutoff value
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_IP_FORMAT,i_mat_cutoff);
    }else if(strcmp(argv[i],"-wc")==0){// SVD wmax
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,wcutoff);
    }else if(strcmp(argv[i],"-fpetsc")==0){// 
      toggle(force_petsc);
    }else{
      if((rank == 0)&&(!warned)){
	fprintf(stderr,"init_parameters encountered at least one parameter which cannot be interpreted by interact\n");
	warned = TRUE;
      }
    }
  }
}

// check, if we can read a value for the option flag
void advance_argument(int *i,int argc, char **argv)
{
  if(argc <= *i + 1){// no arguments left
    fprintf(stderr,"%s: input parameters: error: option \"%s\" needs a value\n",
	    argv[0],argv[*i]);
    exit(-1);
  }
  *i += 1;
}
//
// deal with boolean values/switches
char *name_boolean(my_boolean value)
{
  if(value)
    return("ON");
  else
    return("OFF");
}

my_boolean toggle(my_boolean *variable)
{
  if(*variable){
    *variable=FALSE;
    return(FALSE);
  }else{
    *variable=TRUE;
    return(TRUE);
  }
}

/*
  
  assign the stress matrix loading factors 

  s[i]=a[i]+b[i]*time

  or read them from file, if read_stress_relation_factors is set

  also add a background pressure bpressure

*/
void read_stress_fac(my_boolean read_stress_relation_factors,COMP_PRECISION *a,
		     COMP_PRECISION *b, COMP_PRECISION bpressure,
		     struct med *medium)
{
  int i;
  FILE *in;
  COMP_PRECISION x[3],sm[3][3],dm[3][3],loc_pressure,s2;
  for(i=0;i<6;i++){// initialize stressing constants ands rates with zero
    a[i]=0.0;
    b[i]=0.0;
  }
  //
  // default is simple shear, \sigma_xy
  b[1]=STRESSING_RATE;
  //
  if(!read_stress_relation_factors){
    HEADNODE
      fprintf(stderr,"read_stress_fac: using default simple shear stress relation for loading\n");
  }else if(!(in=fopen(STRESS_RELATION_FILE,"r"))){
    HEADNODE
      fprintf(stderr,"read_stress_fac: can not open \"%s\", will use default background stressing\n",
	      STRESS_RELATION_FILE);
  }else{
    for(i=0;i < 6;i++)// read in in sxx sxy sxz syy syz szz format for A and B
      if(fscanf(in,TWO_CP_FORMAT,(a+i),(b+i))!=2){
	HEADNODE
	  fprintf(stderr,"read_stress_fac: error reading, need six a_i b_i pairs\n");
	exit(-1);
      }
    fclose(in);
    HEADNODE
      fprintf(stderr,"read_stress_fac: read stress relation in vec[6] format from \"%s\"\n",
	      STRESS_RELATION_FILE);
  }
  // default location to evaluate stress
  // (doesn't matter if stress is homogeneous)
  x[INT_X] = x[INT_Y] = x[INT_Z] = 0.0;
  if(bpressure != 0.0){
#ifdef HYDROSTATIC_PRESSURE
    HEADNODE
      fprintf(stderr,"read_stress_fac: to this a depth dependent background pressure of %g will be added\n",
	      bpressure);
#else
    HEADNODE
      fprintf(stderr,"read_stress_fac: to this a constant background pressure of %g will be added\n",
	      bpressure);
#endif
  }
  // stress at time t=0
  HEADNODE
    fprintf(stderr,"read_stress_fac: at x:(%10.3e, %10.3e, %10.3e) and time t: 0\n",
	    x[INT_X],x[INT_Y],x[INT_Z]);
  background_stress(sm,x,0.0,a,b,bpressure);
  HEADNODE{
    fprintf(stderr,"read_stress_fac: stress matrix: ((%10.3e,%10.3e,%10.3e),\n",
	    sm[INT_X][INT_X],sm[INT_X][INT_Y],sm[INT_X][INT_Z]);
    fprintf(stderr,"read_stress_fac: stress matrix:  (%10.3e,%10.3e,%10.3e),\n",
	    sm[INT_Y][INT_X],sm[INT_Y][INT_Y],sm[INT_Y][INT_Z]);	 
    fprintf(stderr,"read_stress_fac: stress matrix:  (%10.3e,%10.3e,%10.3e))\n",
	    sm[INT_Z][INT_X],sm[INT_Z][INT_Y],sm[INT_Z][INT_Z]);
  }
  // deviatoric
  calc_deviatoric_stress(sm,dm,&loc_pressure,&s2);
  HEADNODE
    fprintf(stderr,"read_stress_fac: dev. stress vec: (%10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e) pressure: %10.3e sII: %10.3e\n",
	    dm[INT_X][INT_X],dm[INT_X][INT_Y],dm[INT_X][INT_Z],
	    dm[INT_Y][INT_Y],dm[INT_Y][INT_Z],dm[INT_Z][INT_Z],loc_pressure,s2);
  
  // stress at time t=1
  background_stress(sm,x,1.0,a,b,bpressure);
  calc_deviatoric_stress(sm,dm,&loc_pressure,&s2);
  HEADNODE{
    fprintf(stderr,"read_stress_fac: at x:(%10.3e, %10.3e, %10.3e) and time t: 1\n",x[INT_X],x[INT_Y],x[INT_Z]);
    fprintf(stderr,"read_stress_fac: dev. stress vec: (%10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e) pressure: %10.3e sII: %10.3e\n",
	    dm[INT_X][INT_X],dm[INT_X][INT_Y],dm[INT_X][INT_Z],dm[INT_Y][INT_Y],dm[INT_Y][INT_Z],dm[INT_Z][INT_Z],
	    loc_pressure,s2);
  }
  if(fabs(loc_pressure) < EPS_COMP_PREC)
    HEADNODE
      fprintf(stderr,"read_stress_fac: WARNING: total background pressure is zero\n");
}


