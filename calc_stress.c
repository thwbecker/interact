/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: calc_stress.c,v 2.36 2003/01/13 06:43:35 becker Exp becker $
*/
#include "interact.h"
#include "properties.h"

void initialize_stress_state(struct flt *fault,struct med *medium,
			     my_boolean read_initial_fault_stress,
			     COMP_PRECISION *a,COMP_PRECISION *b)
{
  COMP_PRECISION sm[3][3],trac[3],as[3];
#ifdef SUPER_DEBUG
  COMP_PRECISION abs_tau,abs_tau1;
#endif  
  int i;
  FILE *in=NULL;
  /* 
     initialize incremental loading algorithm 
  */
  if(!medium->stress_state_init){
    /* 
       resolve the loading function as incremental increase on faults 
    */
    HEADNODE{
      fprintf(stderr,"initialize_stress_state: assuming loading is linear with time\n");
#ifdef HYDROSTATIC_PRESSURE
      fprintf(stderr,"initialize_stress_state: assuming hydrostatic pressure scales as %g/%g\n",
	      medium->pressure,CHAR_FAULT_DIM);
#else
      fprintf(stderr,"initialize_stress_state: pressure is held constant at %g at all depths\n",
	      medium->pressure);
#endif
    }
    if((medium->cohesion != 0.0) || (medium->cohesion!=COHESION_DEF))
      HEADNODE
	fprintf(stderr,"initialize_stress_state: WARNING: the cohesion was set to %g\n",
		medium->cohesion);
    if(read_initial_fault_stress){
      in=fopen(FAULT_STRESS_INIT_FILE,"r");
      if(in){
	HEADNODE
	  fprintf(stderr,"initialize_stress_state: WARNING: reading initial non-hydrostatic stresses for faults from \"%s\"\n",
		  FAULT_STRESS_INIT_FILE);
      }else{
	HEADNODE
	  fprintf(stderr,"initialize_stress_state: could not find initial stress file \"%s\", using homogeneous values\n",
		  FAULT_STRESS_INIT_FILE);
	read_initial_fault_stress=FALSE;
      }
    }
#ifdef DEBUG
    HEADNODE
      fprintf(stderr,"initialize_stress_state: initial stresses (including hydrostatic part):\n");
#endif
    if(medium->variable_time_step)
      if(medium->dt0 != 1.0){
	HEADNODE
	  fprintf(stderr,"initialize_stress_state: error, for adjustable time step to work, need dt0=1 (%g)\n",
		  medium->dt0);
	exit(-1);
      }
    for(i=0;i<medium->nrflt;i++){
      /* 
	 use background stress at position x and 
	 initial time to 
	 initialize fault stresses
      */
      background_stress(sm,fault[i].x,medium->time,a,b,medium->pressure);
      calc_three_stress_components(sm,fault[i].normal,fault[i].t_strike,
				   fault[i].t_dip,fault[i].normal,
				   &fault[i].s[STRIKE],&fault[i].s[DIP],
				   &fault[i].s[NORMAL]);
      /* 
	 use background stress at time dt0 to 
	 determine the incremental stress change due
	 to future loading but don't add increment yet
      */
      background_stress(sm,fault[i].x,(medium->time + medium->dt0),a,b,
			medium->pressure);
      resolve_force(fault[i].normal,sm,trac);
      //
      // these are the increments of stress for each dt0 step
      //
      fault[i].sinc[STRIKE] = project_vector(trac,fault[i].t_strike) -
	fault[i].s[STRIKE];
      fault[i].sinc[DIP] =    project_vector(trac,fault[i].t_dip) -
	fault[i].s[DIP];
      fault[i].sinc[NORMAL] = project_vector(trac,fault[i].normal) -
	fault[i].s[NORMAL];
      // read in additional initial stress from file
      if(read_initial_fault_stress){
	if(fscanf(in,THREE_CP_FORMAT,(as+STRIKE),(as+DIP),(as+NORMAL)) != 3){
	  HEADNODE
	    fprintf(stderr,"initialize_stress_state: could not read 3 components of init stress for fault %i\n",
		    i+1);
	  exit(-1);
	}else{
	  fault[i].s[STRIKE] += as[STRIKE];
	  fault[i].s[DIP]    += as[DIP];
	  fault[i].s[NORMAL] += as[NORMAL];
  	}
      }
#ifdef SUPER_DEBUG
      calc_absolute_shear_stress_and_inc(&abs_tau,&abs_tau1,i,fault);
      HEADNODE{
	fprintf(stderr,"s(t0): flt %5i: s: %25.16e d: %25.16e n: %25.16e C: %25.16e\n",
		i,fault[i].s[STRIKE],fault[i].s[DIP],fault[i].s[NORMAL],
		coulomb_stress(abs_tau,(COMP_PRECISION)fault[i].mu_s,
			       fault[i].s[NORMAL],medium->cohesion));
	fprintf(stderr,"s(t0+dt):         s: %25.16e d: %25.16e n: %25.16e C: %25.16e\n",
		fault[i].s[STRIKE]+fault[i].sinc[STRIKE],
		fault[i].s[DIP]   +fault[i].sinc[DIP],
		fault[i].s[NORMAL]+fault[i].sinc[NORMAL],
		coulomb_stress(abs_tau1,(COMP_PRECISION)fault[i].mu_s,
			       fault[i].s[NORMAL]+fault[i].sinc[NORMAL],
			       medium->cohesion));
      }
      
      if(!FINITE_TEST(fault[i].s[STRIKE])|| !FINITE_TEST(fault[i].s[DIP]) || !FINITE_TEST(fault[i].s[NORMAL])){
	HEADNODE
	  fprintf(stderr,"initialize_stress_state: stresses not finite, error!\n");
	exit(-1);
      }
#endif
    }
    medium->stress_state_init=TRUE;
    if(read_initial_fault_stress)
      fclose(in);
    // read in possibly previous activations
    if(medium->attempt_restart)
      adjust_medium_for_restart(medium,fault);
  }
}
/* 
   
   adds incremental load to fault stresses accoding to medium->dt
   
*/
void update_stress_state(struct flt *fault,struct med *medium)
{
  int i;
  if(!medium->variable_time_step){
    /* 
       using base time step dt0 at all times 
    */
#ifdef DEBUG
    if(medium->dt != medium->dt0){
      HEADNODE
	fprintf(stderr,"update_stress_state: expecting dt = dt0\n");
      exit(-1);
    }
#endif
    for(i=0;i<medium->nrflt;i++){
      fault[i].s[STRIKE] += fault[i].sinc[STRIKE];
      fault[i].s[DIP]    += fault[i].sinc[DIP];
      fault[i].s[NORMAL] += fault[i].sinc[NORMAL];
    }
  }else{
    /*
      
    for variable time step calculations, we have set
    dt0 to unity so that we can simply multiply dt*sinc
    
    */
#ifdef DEBUG
    if(medium->dt0 != 1.0){
      HEADNODE
	fprintf(stderr,"update_stress_state: dt0 should be 1\n");
      exit(-1);
    }
#endif 
    for(i=0;i<medium->nrflt;i++){
      fault[i].s[STRIKE] += medium->dt * fault[i].sinc[STRIKE];
      fault[i].s[DIP]    += medium->dt * fault[i].sinc[DIP];
      fault[i].s[NORMAL] += medium->dt * fault[i].sinc[NORMAL];
    }
  }
} 


/* 

   calculate stresses and displacements due to displacements on all
   faults and background stress, either for a bulk field (gridded
   output locations) or for spotted observational locations as read in
   from file

*/
void calc_fields(struct med *medium,struct flt *fault,
		 my_boolean include_background_stress,
		 my_boolean include_background_displacement,
		 COMP_PRECISION *a,COMP_PRECISION *b)
{
  int i,j,k,o,o1,iret,p1,p2,singular_count,not_ok;
  long int nxy,nxyz;
  int nz;
  my_boolean use_fault_plane;
  COMP_PRECISION dx[3],x[3],xl[3],u[3],sm[3][3],vec_1[3],vec_2[3],
    flt_mean_x[3],s1,s2,d1,d2;
  SUM_ARR_PREC *local_u,*local_s;
  FILE *out=NULL;
  singular_count = 0;
  if(medium->print_bulk_fields){
    //
    // grid output locations
    //
    // check if we want a projection plane or 3d field
    fiddle_with_limits_for_plot(medium,&nz,&use_fault_plane,dx,FALSE);
    HEADNODE {
      fprintf(stderr,"calc_fields: calculating bulk stress and displacement fields...\n");
      fprintf(stderr,"calc_fields: field dimensions: %i %i %i %i\n",
	      nz,medium->n[INT_Y],medium->n[INT_X],6);
    }
  }else if(medium->read_oloc_from_file){
    // output on xyz tripels as read from file
    medium->n[INT_X] = medium->n[INT_Y] = 1;
    medium->n[INT_Z] = nz = medium->olocnr;
    HEADNODE
      fprintf(stderr,"calc_fields: calculating fields for %i spotted observations\n",
	      medium->n[INT_Z]);
  }
  nxy = (long int)medium->n[INT_X] * (long int) medium->n[INT_Y];
  nxyz = (long int)nz * nxy;
  if(!nxyz){
    HEADNODE
      fprintf(stderr,"calc_fields: bound of stress/displacement array error: nz %i nxy %li\n",
	      nz,nxy);
    exit(-1);
  }
  HEADNODE{		/* global arrays */
    if((medium->s=(SUM_ARR_PREC *)calloc(nxyz*6,sizeof(SUM_ARR_PREC)))==NULL)
      PMEMERROR("calc_fields: stress array");
    if((medium->u=(SUM_ARR_PREC *)calloc(nxyz*3,sizeof(SUM_ARR_PREC)))==NULL)
      PMEMERROR("calc_fields: displacement array");
  }else{
    medium->u = NULL;
    medium->s = NULL;
  }
  if(medium->comm_size > 1){
    //fprintf(stderr,"calc_fields: allocating local array for core %i\n",medium->comm_rank);
    if((local_s=(SUM_ARR_PREC *)calloc(nxyz*6,sizeof(SUM_ARR_PREC)))==NULL)
      PMEMERROR("calc_fields: local stress array");
    if((local_u=(SUM_ARR_PREC *)calloc(nxyz*3,sizeof(SUM_ARR_PREC)))==NULL)
      PMEMERROR("calc_fields: local displacement array");
  }else{
    local_u = medium->u;
    local_s = medium->s;
  }
  /*
    
    background stress and displacement first
    
  */
  if(medium->print_bulk_fields){
    //
    // ouytput on grid
    //
    //
    // determine if any of the locations are above zero
    // x[INT_X]/x[INT_Y] is used for direction along strike and dip
    if(use_fault_plane){
      //
      // obtain the average strike and dip or normal vectors depending on n[INT_Z]
      //
      get_fault_plane_basevec(flt_mean_x,vec_1,vec_2,fault,medium);
      vec_to_angles(vec_1, &d1, &s1);
      vec_to_angles(vec_2, &d2, &s2);
      HEADNODE
	fprintf(stderr,"calc_fields: avg. fault plane: x:(%g,%g,%g), v_s: strike: %g dip: %g, v_%s: strike: %g dip: %g\n",
		flt_mean_x[INT_X],flt_mean_x[INT_Y],flt_mean_x[INT_Z],
		s1,d1,(medium->n[INT_Z]==-1)?("d"):("n"),s2,d2);
      if(!(medium->ok=(my_boolean *)malloc(nxy*sizeof(my_boolean))))
	PMEMERROR("calc_fields: 3");
      // check if the output is OK
      for(o1=i=not_ok=0,x[INT_X]=medium->pxmin[INT_X];
	  i < medium->n[INT_X];
	  x[INT_X] += dx[INT_X], i++,o1 += medium->n[INT_Y])
	for(j=0,x[INT_Y]=medium->pxmin[INT_Y];j < medium->n[INT_Y];
	    x[INT_Y] += dx[INT_Y],j++){
	  get_local_x_on_plane(xl,x,flt_mean_x,vec_1,vec_2);
	  if(xl[INT_Z] <= 0){
	    medium->ok[o1+j]=TRUE;
	  }else{
	    medium->ok[o1+j]=FALSE;
	    not_ok++;
	  }
	}
      if(not_ok)
	HEADNODE
	  fprintf(stderr,"calc_fields: %i out of %i points were not appropriate for plane mode\n",
		  not_ok,medium->n[INT_Y]*medium->n[INT_X]);
    }
    if(include_background_stress || include_background_displacement){
      HEADNODE
	fprintf(stderr,"calc_fields: WARNING: including background stress for grid\n");
      /* include background stress, else zero */
      for(k=0,x[INT_Z]=medium->pxmin[INT_Z];k<nz;x[INT_Z]+=dx[INT_Z],k++){
	for(i=0,x[INT_X]=medium->pxmin[INT_X];i<medium->n[INT_X];x[INT_X]+=dx[INT_X],i++){
	  for(j=0,x[INT_Y]=medium->pxmin[INT_Y];j<medium->n[INT_Y];x[INT_Y]+=dx[INT_Y],j++){
	    if(!use_fault_plane || medium->ok[i*medium->n[INT_Y]+j]){
	      if(use_fault_plane){
		// determine position along the fault plane
		get_local_x_on_plane(xl,x,flt_mean_x,vec_1,vec_2);
	      }else{
		xl[INT_X]=x[INT_X];xl[INT_Y]=x[INT_Y];xl[INT_Z]=x[INT_Z];
	      }
	      if(include_background_displacement){
		background_disp(u,xl,medium,a,b);
		p1=POSU(i,j,k,INT_X);
		local_u[p1++] = (SUM_ARR_PREC)u[INT_X];
		local_u[p1++] = (SUM_ARR_PREC)u[INT_Y];
		local_u[p1]   = (SUM_ARR_PREC)u[INT_Z];
	      }
	      if(include_background_stress){
		background_stress(sm,xl,medium->time,a,b,medium->pressure);
		p2=POSS(i,j,k,0);
		local_s[p2++] = (SUM_ARR_PREC)sm[INT_X][INT_X];
		local_s[p2++] = (SUM_ARR_PREC)sm[INT_X][INT_Y];
		local_s[p2++] = (SUM_ARR_PREC)sm[INT_X][INT_Z];
		local_s[p2++] = (SUM_ARR_PREC)sm[INT_Y][INT_Y];
		local_s[p2++] = (SUM_ARR_PREC)sm[INT_Y][INT_Z];
		local_s[p2]   = (SUM_ARR_PREC)sm[INT_Z][INT_Z];
	      }
	    }
	  }
	}
      }
    }
  }else if(medium->read_oloc_from_file){
    //
    // output locations are given on xyz tripels
    //
    if(include_background_stress || include_background_displacement){
      HEADNODE
	fprintf(stderr,"calc_fields: WARNING: including background stress for spotted\n");
      for(i=j=0;i<medium->n[INT_Z];i++,j+=3){
	for(k=0;k<3;k++)
	  xl[k] = (COMP_PRECISION)medium->xoloc[j+k];
	if(include_background_displacement){
	  background_disp(u,xl,medium,a,b);
	  p1 = j;
	  local_u[p1++] = (SUM_ARR_PREC)u[INT_X];
	  local_u[p1++] = (SUM_ARR_PREC)u[INT_Y];
	  local_u[p1]   = (SUM_ARR_PREC)u[INT_Z];
	}
	if(include_background_stress){
	  background_stress(sm,xl,medium->time,a,b,medium->pressure);
	  p2= i * 6;
	  local_s[p2++] = (SUM_ARR_PREC)sm[INT_X][INT_X];
	  local_s[p2++] = (SUM_ARR_PREC)sm[INT_X][INT_Y];
	  local_s[p2++] = (SUM_ARR_PREC)sm[INT_X][INT_Z];
	  local_s[p2++] = (SUM_ARR_PREC)sm[INT_Y][INT_Y];
	  local_s[p2++] = (SUM_ARR_PREC)sm[INT_Y][INT_Z];
	  local_s[p2]   = (SUM_ARR_PREC)sm[INT_Z][INT_Z];
	}
      }
    }
  }
  /* 

     add up contributions from all faults 

  */
  if(medium->print_bulk_fields){
    HEADNODE
      if(medium->print_plane_coord){
	out=myopen(PLANE_COORD_FILE,"w");
	fprintf(out,"# format:\n#\txg[X] xg[Y] xg[Z] xl[X] xl[Y]\n#\twith v_s: (%g,%g,%g) v_%s: (%g,%g,%g)\n#\n",
		vec_1[INT_X],vec_1[INT_Y],vec_1[INT_Z],(medium->n[INT_Z]==-1)?("d"):("n"),
		vec_2[INT_X],vec_2[INT_Y],vec_2[INT_Z]);
	fprintf(stderr,"calc_fields: writing fault plane coordinates to \"%s\"\n",
		PLANE_COORD_FILE);
      }

    if(medium->comm_size > 1)
      fprintf(stderr,"calc_stress: core %03i/%03i computing grid for %05i to %05i\n",
	      medium->comm_rank,medium->comm_size,
	      medium->myfault0,medium->myfaultn);
    
    for(k=0,x[INT_Z]=medium->pxmin[INT_Z];k<nz;x[INT_Z]+=dx[INT_Z],k++)
      for(i=0,x[INT_X]=medium->pxmin[INT_X];i<medium->n[INT_X];x[INT_X]+=dx[INT_X],i++)
	for(j=0,x[INT_Y]=medium->pxmin[INT_Y];j<medium->n[INT_Y];x[INT_Y]+=dx[INT_Y],j++){
	  //fprintf(stderr,"calc_fields: working on %04i/%04i/%04i\r",k,i,j);
	  if(!use_fault_plane || medium->ok[i*medium->n[INT_Y]+j]){
	    if(use_fault_plane){
	      // determine position along the fault plane
	      get_local_x_on_plane(xl,x,flt_mean_x,vec_1,vec_2);
	    }else{
	      xl[INT_X]=x[INT_X];xl[INT_Y]=x[INT_Y];xl[INT_Z]=x[INT_Z];
	    }
	    HEADNODE
	      if(medium->print_plane_coord)
		fprintf(out,"%g %g %g %g %g\n",xl[INT_X],xl[INT_Y],xl[INT_Z],x[INT_X],x[INT_Y]);
	    if(xl[INT_Z] > 0.0){
	      if(xl[INT_Z] > 1e-10){
		fprintf(stderr,"calc_fields: positive depth in loop, kij: %i %i %i x: (%g, %g, %g)\n",
			k,i,j,xl[INT_X],xl[INT_Y],xl[INT_Z]);
		exit(-1);
	      }else{
		xl[INT_Z]=0.0;
	      }
	    }
	    p1 = POSU(i,j,k,INT_X);
	    p2 = POSS(i,j,k,0);

	    /* possibly executed only for each core */
	    for(o = medium->myfault0;o < medium->myfaultn;o++){
	      
	      if(norm_3d(fault[o].u) >= EPS_COMP_PREC){
		//
		// actual fault contribution is accounted for HERE
		//
		eval_green(xl,(fault+o),fault[o].u,u,sm,&iret);
		if(!iret){
		  local_u[p1]   += (SUM_ARR_PREC)u[INT_X];
		  local_u[p1+1] += (SUM_ARR_PREC)u[INT_Y];
		  local_u[p1+2] += (SUM_ARR_PREC)u[INT_Z];
		  local_s[p2]   += (SUM_ARR_PREC)sm[INT_X][INT_X];
		  local_s[p2+1] += (SUM_ARR_PREC)sm[INT_X][INT_Y];
		  local_s[p2+2] += (SUM_ARR_PREC)sm[INT_X][INT_Z];
		  local_s[p2+3] += (SUM_ARR_PREC)sm[INT_Y][INT_Y];
		  local_s[p2+4] += (SUM_ARR_PREC)sm[INT_Y][INT_Z];
		  local_s[p2+5] += (SUM_ARR_PREC)sm[INT_Z][INT_Z];
		}else{
		  singular_count++;
		  local_u[p1]   = medium->nan;
		  local_u[p1+1] = medium->nan;
		  local_u[p1+2] = medium->nan;
		  local_s[p2]   = medium->nan;
		  local_s[p2+1] = medium->nan;
		  local_s[p2+2] = medium->nan;
		  local_s[p2+3] = medium->nan;
		  local_s[p2+4] = medium->nan;
		  local_s[p2+5] = medium->nan;
		}
	      }
	    }
	  }
	}
    HEADNODE
      if(medium->print_plane_coord)
	fclose(out);
  }else  if(medium->read_oloc_from_file){
    if(medium->comm_size > 1)
      fprintf(stderr,"calc_stress: core %03i/%03i computing spotted for %05i to %05i\n",
	      medium->comm_rank,medium->comm_size,
	      medium->myfault0,medium->myfaultn);
    //
    // output given on spotted locations
    //
    for(i=j=0;i < medium->olocnr;i++,j+=3){
      for(k=0;k<3;k++)
	xl[k]=(COMP_PRECISION)medium->xoloc[j+k];
      if(xl[INT_Z] > 0.0){
	if(xl[INT_Z] > 1e-10){
	  fprintf(stderr,"calc_fields: positive depth for location %i x: (%g, %g, %g)\n",
		  k,xl[INT_X],xl[INT_Y],xl[INT_Z]);
	  exit(-1);
	}else{
	  xl[INT_Z]=0.0;
	}
      }
      p1 = j;
      p2 = i * 6;
      for(o=medium->myfault0;o < medium->myfaultn;o++){
	
	if(norm_3d(fault[o].u) >= EPS_COMP_PREC){
	  //
	  // actual fault contribution is accounted for HERE
	  //
	  eval_green(xl,(fault+o),fault[o].u,u,sm,&iret);
	  if(!iret){
	    local_u[p1]   += (SUM_ARR_PREC)u[INT_X];
	    local_u[p1+1] += (SUM_ARR_PREC)u[INT_Y];
	    local_u[p1+2] += (SUM_ARR_PREC)u[INT_Z];
	    local_s[p2]   += (SUM_ARR_PREC)sm[INT_X][INT_X];
	    local_s[p2+1] += (SUM_ARR_PREC)sm[INT_X][INT_Y];
	    local_s[p2+2] += (SUM_ARR_PREC)sm[INT_X][INT_Z];
	    local_s[p2+3] += (SUM_ARR_PREC)sm[INT_Y][INT_Y];
	    local_s[p2+4] += (SUM_ARR_PREC)sm[INT_Y][INT_Z];
	    local_s[p2+5] += (SUM_ARR_PREC)sm[INT_Z][INT_Z];
	  }else{
	    singular_count++;
	    local_u[p1]   = medium->nan;
	    local_u[p1+1] = medium->nan;
	    local_u[p1+2] = medium->nan;
	    local_s[p2]   = medium->nan;
	    local_s[p2+1] = medium->nan;
	    local_s[p2+2] = medium->nan;
	    local_s[p2+3] = medium->nan;
	    local_s[p2+4] = medium->nan;
	    local_s[p2+5] = medium->nan;
	  }
	}
      }
    }
  }
  if(singular_count)
    fprintf(stderr,"calc_fields: core %03i WARNING: there were %i singular entries in the field\n",
	    medium->comm_rank,singular_count);
  else
    fprintf(stderr,"calc_fields: core %03i done, no singular entries in the field\n",
	    medium->comm_rank);
  //for(i=0;i<10;i++)fprintf(stderr,"%10g ",local_u[i]);fprintf(stderr,"\n");
#ifdef USE_PETSC
  if(medium->comm_size > 1){

    MPI_Reduce(local_u, medium->u, (int)nxyz*3, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_s, medium->s, (int)nxyz*6, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    free(local_u);free(local_s);
  }
#endif
  medium->bulk_field_init=TRUE;
}
