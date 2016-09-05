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
    fprintf(stderr,"initialize_stress_state: assuming loading is linear with time\n");
#ifdef HYDROSTATIC_PRESSURE
    fprintf(stderr,"initialize_stress_state: assuming hydrostatic pressure scales as %g/%g\n",
	    medium->pressure,CHAR_FAULT_DIM);
#else
    fprintf(stderr,"initialize_stress_state: pressure is held constant at %g at all depths\n",
	    medium->pressure);
#endif
    if((medium->cohesion != 0.0) || (medium->cohesion!=COHESION_DEF))
      fprintf(stderr,"initialize_stress_state: WARNING: the cohesion was set to %g\n",
	      medium->cohesion);
    if(read_initial_fault_stress){
      in=fopen(FAULT_STRESS_INIT_FILE,"r");
      if(in)
	fprintf(stderr,"initialize_stress_state: WARNING: reading initial non-hydrostatic stresses for faults from \"%s\"\n",
		FAULT_STRESS_INIT_FILE);
      else{
	fprintf(stderr,"initialize_stress_state: could not find initial stress file \"%s\", using homogeneous values\n",
		FAULT_STRESS_INIT_FILE);
	read_initial_fault_stress=FALSE;
      }
    }
#ifdef DEBUG
    fprintf(stderr,"initialize_stress_state: initial stresses (including hydrostatic part):\n");
#endif
    if(medium->variable_time_step)
      if(medium->dt0 != 1.0){
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
      
      if(!FINITE_TEST(fault[i].s[STRIKE])|| !FINITE_TEST(fault[i].s[DIP]) || !FINITE_TEST(fault[i].s[NORMAL])){
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
  int i,j,k,o,o1,iret,p1,p2,singular_count,nz,not_ok,nxy;
  my_boolean use_fault_plane;
  COMP_PRECISION dx[3],x[3],xl[3],u[3],sm[3][3],vec_1[3],vec_2[3],
    flt_mean_x[3],s1,s2,d1,d2;
  FILE *out=NULL;
  singular_count = 0;
  if(medium->print_bulk_fields){
    //
    // grid output locations
    //
    // check if we want a projection plane or 3d field
    fiddle_with_limits_for_plot(medium,&nz,&use_fault_plane,dx,FALSE);
    fprintf(stderr,"calc_fields: calculating bulk stress and displacement fields...\n");
    fprintf(stderr,"calc_fields: field dimensions: %i %i %i %i\n",
	    nz,medium->n[Y],medium->n[X],6);
  }
  if(medium->read_oloc_from_file){
    // output on xyz tripels as read from file
    medium->n[X] = medium->n[Y] = 1;
    medium->n[Z] = nz = medium->olocnr;
    fprintf(stderr,"calc_fields: calculating fields for %i spotted observations\n",
	    medium->n[Z]);
  }
  nxy = medium->n[X] * medium->n[Y];
  if(!(nz*nxy*6)){
    fprintf(stderr,"calc_fields: bound of stress/displacement array error: %i\n",
	    nz*nxy*6);
    exit(-1);
  }
  if((medium->s=(float *)calloc(nz*nxy*6,sizeof(float)))==NULL)
    MEMERROR("calc_fields: stress array");
  if((medium->u=(float *)calloc(nz*nxy*3,sizeof(float)))==NULL)
    MEMERROR("calc_fields: displacement array");
  /*
    
    background stress and displacement first
    
  */
  if(medium->print_bulk_fields){
    //
    // ouytput on grid
    //
    //
    // determine if any of the locations are above zero
    // x[X]/x[Y] is used for direction along strike and dip
    if(use_fault_plane){
      //
      // obtain the average strike and dip or normal vectors depending on n[Z]
      //
      get_fault_plane_basevec(flt_mean_x,vec_1,vec_2,fault,medium);
      vec_to_angles(vec_1, &d1, &s1);
      vec_to_angles(vec_2, &d2, &s2);
      fprintf(stderr,"calc_fields: avg. fault plane: x:(%g,%g,%g), v_s: strike: %g dip: %g, v_%s: strike: %g dip: %g\n",
	      flt_mean_x[X],flt_mean_x[Y],flt_mean_x[Z],
	      s1,d1,(medium->n[Z]==-1)?("d"):("n"),s2,d2);
      if(!(medium->ok=(my_boolean *)malloc(nxy*sizeof(my_boolean))))
	MEMERROR("calc_fields: 3");
      // check if the output is OK
      for(o1=i=not_ok=0,x[X]=medium->pxmin[X];
	  i<medium->n[X];
	  x[X]+=dx[X], i++,o1 += medium->n[Y])
	for(j=0,x[Y]=medium->pxmin[Y];j<medium->n[Y];x[Y]+=dx[Y],j++){
	  get_local_x_on_plane(xl,x,flt_mean_x,vec_1,vec_2);
	  if(xl[Z]<=0.0){
	    medium->ok[o1+j]=TRUE;
	  }else{
	    medium->ok[o1+j]=FALSE;
	    not_ok++;
	  }
	}
      if(not_ok)
	fprintf(stderr,"calc_fields: %i out of %i points were not appropriate for plane mode\n",
		not_ok,medium->n[Y]*medium->n[X]);
    }
    if(include_background_stress || include_background_displacement){
      fprintf(stderr,"calc_fields: WARNING: including background stress\n");
      /* include background stress, else zero */
      for(k=0,x[Z]=medium->pxmin[Z];k<nz;x[Z]+=dx[Z],k++)
	for(i=0,x[X]=medium->pxmin[X];i<medium->n[X];x[X]+=dx[X],i++)
	  for(j=0,x[Y]=medium->pxmin[Y];j<medium->n[Y];x[Y]+=dx[Y],j++){
	    if(!use_fault_plane || medium->ok[i*medium->n[Y]+j]){
	      if(use_fault_plane){
		// determine position along the fault plane
		get_local_x_on_plane(xl,x,flt_mean_x,vec_1,vec_2);
	      }else{
		xl[X]=x[X];xl[Y]=x[Y];xl[Z]=x[Z];
	      }
	      if(include_background_displacement){
		background_disp(u,xl,medium,a,b);
		p1=POSU(i,j,k,X);
		medium->u[p1++] = (float)u[X];
		medium->u[p1++] = (float)u[Y];
		medium->u[p1]   = (float)u[Z];
	      }
	      if(include_background_stress){
		background_stress(sm,xl,medium->time,a,b,medium->pressure);
		p2=POSS(i,j,k,0);
		medium->s[p2++] = (float)sm[X][X];
		medium->s[p2++] = (float)sm[X][Y];
		medium->s[p2++] = (float)sm[X][Z];
		medium->s[p2++] = (float)sm[Y][Y];
		medium->s[p2++] = (float)sm[Y][Z];
		medium->s[p2]   = (float)sm[Z][Z];
	      }
	    }
	  }
    }
  }
  if(medium->read_oloc_from_file){
    //
    // output locations are given on xyz tripels
    //
    if(include_background_stress || include_background_displacement){
      fprintf(stderr,"calc_fields: WARNING: including background stress\n");
      for(i=j=0;i<medium->n[Z];i++,j+=3){
	for(k=0;k<3;k++)
	  xl[k] = (COMP_PRECISION)medium->xoloc[j+k];
	if(include_background_displacement){
	  background_disp(u,xl,medium,a,b);
	  p1 = j;
	  medium->u[p1++] = (float)u[X];
	  medium->u[p1++] = (float)u[Y];
	  medium->u[p1]   = (float)u[Z];
	}
	if(include_background_stress){
	  background_stress(sm,xl,medium->time,a,b,medium->pressure);
	  p2= i * 6;
	  medium->s[p2++] = (float)sm[X][X];
	  medium->s[p2++] = (float)sm[X][Y];
	  medium->s[p2++] = (float)sm[X][Z];
	  medium->s[p2++] = (float)sm[Y][Y];
	  medium->s[p2++] = (float)sm[Y][Z];
	  medium->s[p2]   = (float)sm[Z][Z];
	}
      }
    }
  }
  /* 
     add up contributions from all faults 
  */
  if(medium->print_bulk_fields){
    if(medium->print_plane_coord){
      out=myopen(PLANE_COORD_FILE,"w");
      fprintf(out,"# format:\n#\txg[X] xg[Y] xg[Z] xl[X] xl[Y]\n#\twith v_s: (%g,%g,%g) v_%s: (%g,%g,%g)\n#\n",
	      vec_1[X],vec_1[Y],vec_1[Z],(medium->n[Z]==-1)?("d"):("n"),
	      vec_2[X],vec_2[Y],vec_2[Z]);
      fprintf(stderr,"calc_fields: writing fault plane coordinates to \"%s\"\n",
	      PLANE_COORD_FILE);
    }
    for(k=0,x[Z]=medium->pxmin[Z];k<nz;x[Z]+=dx[Z],k++)
      for(i=0,x[X]=medium->pxmin[X];i<medium->n[X];x[X]+=dx[X],i++)
	for(j=0,x[Y]=medium->pxmin[Y];j<medium->n[Y];x[Y]+=dx[Y],j++){
	  //fprintf(stderr,"calc_fields: working on %04i/%04i/%04i\r",k,i,j);
	  if(!use_fault_plane || medium->ok[i*medium->n[Y]+j]){
	    if(use_fault_plane){
	      // determine position along the fault plane
	      get_local_x_on_plane(xl,x,flt_mean_x,vec_1,vec_2);
	    }else{
	      xl[X]=x[X];xl[Y]=x[Y];xl[Z]=x[Z];
	    }
	    if(medium->print_plane_coord)
	      fprintf(out,"%g %g %g %g %g\n",xl[X],xl[Y],xl[Z],x[X],x[Y]);
	    if(xl[Z] > 0.0){
	      if(xl[Z] > 1.0e-10){
		fprintf(stderr,"calc_fields: positive depth in loop, kij: %i %i %i x: (%g, %g, %g)\n",
			k,i,j,xl[X],xl[Y],xl[Z]);
		exit(-1);
	      }else{
		xl[Z]=0.0;
	      }
	    }
	    p1=POSU(i,j,k,X);
	    p2=POSS(i,j,k,0);
	    for(o=0;o<medium->nrflt;o++){
	      if(norm_3d(fault[o].u) >= EPS_COMP_PREC){
		//
		// actual fault contribution is accounted for HERE
		//
		eval_green(xl,(fault+o),fault[o].u,u,sm,&iret);
		if(!iret){
		  medium->u[p1]   += (float)u[X];
		  medium->u[p1+1] += (float)u[Y];
		  medium->u[p1+2] += (float)u[Z];
		  medium->s[p2]   += (float)sm[X][X];
		  medium->s[p2+1] += (float)sm[X][Y];
		  medium->s[p2+2] += (float)sm[X][Z];
		  medium->s[p2+3] += (float)sm[Y][Y];
		  medium->s[p2+4] += (float)sm[Y][Z];
		  medium->s[p2+5] += (float)sm[Z][Z];
		}else{
		  singular_count++;
		  medium->u[p1] = medium->nan;
		  medium->u[p1+1] = medium->nan;
		  medium->u[p1+2] = medium->nan;
		  medium->s[p2] = medium->nan;
		  medium->s[p2+1] = medium->nan;
		  medium->s[p2+2] = medium->nan;
		  medium->s[p2+3] = medium->nan;
		  medium->s[p2+4] = medium->nan;
		  medium->s[p2+5] = medium->nan;
		}
	      }
	    }
	  }
	}
    if(medium->print_plane_coord)
      fclose(out);
  }
  singular_count=0;
  if(medium->read_oloc_from_file){
    // output given on spotted locations
    for(i=j=0;i<medium->olocnr;i++,j+=3){
      for(k=0;k<3;k++)
	xl[k]=(COMP_PRECISION)medium->xoloc[j+k];
      if(xl[Z] > 0.0){
	if(xl[Z] > 1.0e-10){
	  fprintf(stderr,"calc_fields: positive depth for location %i x: (%g, %g, %g)\n",
		  k,xl[X],xl[Y],xl[Z]);
	  exit(-1);
	}else{
	  xl[Z]=0.0;
	}
      }
      p1 = j;
      p2 = i * 6;
      for(o=0;o<medium->nrflt;o++){
	if(norm_3d(fault[o].u) >= EPS_COMP_PREC){
	  //
	  // actual fault contribution is accounted for HERE
	  //
	  eval_green(xl,(fault+o),fault[o].u,u,sm,&iret);
	  if(!iret){
	    medium->u[p1]   += (float)u[X];
	    medium->u[p1+1] += (float)u[Y];
	    medium->u[p1+2] += (float)u[Z];
	    medium->s[p2]   += (float)sm[X][X];
	    medium->s[p2+1] += (float)sm[X][Y];
	    medium->s[p2+2] += (float)sm[X][Z];
	    medium->s[p2+3] += (float)sm[Y][Y];
	    medium->s[p2+4] += (float)sm[Y][Z];
	    medium->s[p2+5] += (float)sm[Z][Z];
	  }else{
	    singular_count++;
	    medium->u[p1]   = medium->nan;
	    medium->u[p1+1] = medium->nan;
	    medium->u[p1+2] = medium->nan;
	    medium->s[p2]   = medium->nan;
	    medium->s[p2+1] = medium->nan;
	    medium->s[p2+2] = medium->nan;
	    medium->s[p2+3] = medium->nan;
	    medium->s[p2+4] = medium->nan;
	    medium->s[p2+5] = medium->nan;
	  }
	}
      }
    }
  }
  if(singular_count)
    fprintf(stderr,"calc_fields: WARNING: there were %i singular entries in the field\n",
	    singular_count);
  else
    fprintf(stderr,"calc_fields: done, no singular entries in the field\n");
  
  medium->bulk_field_init=TRUE;
}

/*

  calculate the background stress (sm[3][3] matrix) at location x[3]
  and time "time" given the constant stress matrix factors a[6] and the 
  loading rates b[6] as well as the pressure "pressure"

  

*/

void background_stress(COMP_PRECISION sm[3][3], COMP_PRECISION *x, 
		       COMP_PRECISION time,COMP_PRECISION *a,
		       COMP_PRECISION *b,COMP_PRECISION pressure)
{
  COMP_PRECISION locp;
#ifdef HYDROSTATIC_PRESSURE
  locp = -(x[Z]/HYDROSTATIC_PRESSURE) * pressure; 
#else
  locp = pressure; 
#endif
  /* isotropic elements, compression negative */
  sm[X][X] = a[0] + time * b[0] - locp;
  sm[Y][Y] = a[3] + time * b[3] - locp;
  sm[Z][Z] = a[5] + time * b[5] - locp;
  /* off diagonal elements  */  
  sm[X][Y]=sm[Y][X] = a[1] + time * b[1];
  sm[X][Z]=sm[Z][X] = a[2] + time * b[2];
  sm[Y][Z]=sm[Z][Y] = a[4] + time * b[4];
}
void background_disp(COMP_PRECISION *u, COMP_PRECISION *x, 
		     struct med *medium,COMP_PRECISION *a,
		     COMP_PRECISION *b)
{
  /* the characteristic strain rate for simple shear 
     is given by 

     characteristic stressing rate
     --------------------------
     2 mu

     integration gives the characteristic strain

     characteristic stressing rate
     ----------------------------- y_location
     mu
     
  */
  int i;
  my_boolean hit=FALSE;
  for(i=0;i<6;i++)
    if(a[i]!=0.0||((i!=1)&&(b[i]!=0.0))){hit=TRUE;break;}
  if(hit){
    fprintf(stderr,"background_disp: EXITING: background displacement is inaccurate since no simple shear stressing\n");
    exit(-1);
  }
  u[X]=medium->time * (b[1]/SHEAR_MODULUS)*u[Y];
  u[Y]=u[Z]=0.0;
}
/*

  obtain the local coordinates given the base vectors vec_1 and vec_2


*/
void get_local_x_on_plane(COMP_PRECISION *xl,COMP_PRECISION *x,
			  COMP_PRECISION *flt_mean_x,COMP_PRECISION *vec_1,
			  COMP_PRECISION *vec_2)
{
  int i;
  for(i=0;i<3;i++){
    xl[i]  = flt_mean_x[i];
    xl[i] += vec_1[i] * x[X];
    xl[i] += vec_2[i] * x[Y];
  }
}
/*


  obtain average faulkt plane vectors and location
  on return flt_mean_x will hold the mean location and vec_1 and vec_2
  the mean strike and dip or the mean strike and normal vectors, depending
  on the n[Z] flag, -1 or -2 


 */
void get_fault_plane_basevec(COMP_PRECISION *flt_mean_x,
			     COMP_PRECISION *vec_1,COMP_PRECISION *vec_2,
			     struct flt *fault,struct med *medium)
{
  int n,i,j;
  // get average fault plane vectors
  // and mean location of patches
  for(i=0;i<3;i++)
    vec_1[i]=vec_2[i]=flt_mean_x[i]=0.0;
  if(medium->n[Z] == -1){
    fprintf(stderr,"get_fault_plane_basevec: base vectors are average strike and dip of fault group 0\n");
    for(n=i=0;i<medium->nrflt;i++)
      if(fault[i].group == 0){
	n++;
	for(j=0;j<3;j++){
	  flt_mean_x[j]   += fault[i].x[j];
	  vec_1[j]        += fault[i].t_strike[j];
	  vec_2[j]        += fault[i].t_dip[j];
	}
      }
  }else if(medium->n[Z] == -2){
    fprintf(stderr,"get_fault_plane_basevec: base vectors are average strike and normal of fault group 0\n");
    for(n=i=0;i<medium->nrflt;i++)
      if(fault[i].group == 0){
	n++;
	for(j=0;j<3;j++){
	  flt_mean_x[j]   += fault[i].x[j];
	  vec_1[j]        += fault[i].t_strike[j];
	  vec_2[j]        += fault[i].normal[j];
	}
      }
  }else{
    fprintf(stderr,"get_fault_plane_basevec: medium->n[Z] has to be -1 or -2 but is %i\n",
	    medium->n[Z]);
    exit(-1);
  }
  if(n)
    for(i=0;i<3;i++){
      flt_mean_x[i]  /=(COMP_PRECISION)n;
      if(fabs(flt_mean_x[i])<EPS_COMP_PREC)
	flt_mean_x[i]=0.0;
      vec_1[i]  /=(COMP_PRECISION)n;
      if(fabs(vec_1[i])<EPS_COMP_PREC)
	vec_1[i]=0.0;
      vec_2[i]     /=(COMP_PRECISION)n;
      if(fabs(vec_2[i])<EPS_COMP_PREC)
	vec_2[i]=0.0;
    }
  normalize_3d(vec_1);normalize_3d(vec_2);
}
/*

  calculate the deviator stress and pressure
  given a stress matrix sm
  
  output is dm and pressure (of original tensor), and second invariant
  of deviatoric tensor

*/

void calc_deviatoric_stress(COMP_PRECISION sm[3][3],COMP_PRECISION dm[3][3],
			    COMP_PRECISION *pressure, COMP_PRECISION *s2)
{
  int i,j;
  *pressure= -(sm[X][X] + sm[Y][Y] + sm[Z][Z])/3.0;
 
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      dm[i][j] = sm[i][j] + ((i==j)?(*pressure):(0.0));
  /* second invariant of deviatoric tensor */
  *s2 = sqrt(0.5* (dm[X][X] * dm[X][X] + 
		   dm[X][Y] * dm[X][Y] * 2.0 + 
		   dm[Y][Y] * dm[Y][Y] + 
		   dm[Y][Z] * dm[Y][Z] * 2.0 + 
		   dm[Z][Z] * dm[Z][Z] + 
		   dm[X][Z] * dm[X][Z] * 2.0));
}
