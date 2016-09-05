/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: rupture.c,v 2.25 2002/10/08 19:24:44 tbecker Exp $

*/
#include "interact.h"
/* 
   check if rupture criterion fulfilled at any faults
   at given time
   if so, solve for the solution so that all faults have
   the target stress 
*/


my_boolean activate_faults(struct flt *fault,struct med *medium)
{
  int i,old_nr_of_active_groups,nrsweeps,old_active[2];
  my_boolean faults_have_slipped,active_modified,mark_quakes;
  // did faults slip at all (not only in last iteration)
  faults_have_slipped=FALSE;
  //
  // start solution loop, iter belongs to medium structure since we want
  // to access it in the quake subroutine, else it could be local
  //
  medium->iter=0;
  do{
    medium->iter++;
    if((medium->keep_slipping && (medium->iter == 1))||
       (!medium->keep_slipping)){
      //
      // INITIALIZATION
      // 
      // first iteration of a keep_slipping run or restart from iteration
      // with exhaustive slip
      //
      // memorize the stress drop that would result from initial 
      // shear stress value and dynamical friction coefficient
      // (this is important for possible later 'dynamic' activation)
      // as well as the initial shear stress configuration for mixed mode 
      // activations
      //
      for(i=0;i<medium->nrflt;i++){
	if(fault[i].mode[0] != INACTIVE){
	  // initialize stress ratio factor
	  if(fabs(fault[i].s[DIP]) > 0)
	    fault[i].f_initial = fabs(fault[i].s[STRIKE]/fault[i].s[DIP]);
	  else// no dip stresses
	    fault[i].f_initial = -1.0;
	  // initialize the dynamic friction stress drop
	  // as |\tau|^i + \mu_d \sigma_n^i
	  calc_absolute_shear_stress(&fault[i].taud,i,fault);
	  fault[i].taud += (COMP_PRECISION)fault[i].mu_d * fault[i].s[NORMAL];
	}
      }
      /* 
	 reset the arrays that deal with equations for the
	 constraint system. also switch all group activation flags off
	 we do this only for the first iterations since we don't want to loose
	 the patches that are already active
      */
      init_equation_system(medium,fault);
      old_active[0] = old_active[1] = 0;
    }else{
      old_active[0] = medium->naflt;
      old_active[1] = medium->naflt_con;
    }
    /* 
       check for all fault patches that will be active at this time step 
       and determine the target stress drop values 
    */
    nrsweeps=0;
    do{
      // this loop is only executed more than once if  whole fault activations
      // are selected. them, we might have to run it twice to get all patches 
      nrsweeps++;
      if(nrsweeps>2)
	fprintf(stderr,"activate_faults: WARNING: %ith loop seeking for activations\n",
		nrsweeps);
      old_nr_of_active_groups = medium->nr_active_groups;
      //
      // check all faults for activation
      //
      for(i=0;i<medium->nrflt;i++)
	fault_criterion(i,fault,medium);
      // if the number of active groups has changed in 
      // this sweep through all patches and the whole
      // fault activation mode is selected, we have to 
      // go back and check again so that we are not missing
      // possibly slipping patches further down
    }while((medium->whole_fault_activations)&&
	   (medium->nr_active_groups != old_nr_of_active_groups));
    //
    // execute the solver routine if we have active faults
    //
    if((medium->naflt)||(medium->naflt_con)){
      faults_have_slipped=TRUE;
      //
      // set flag if number of active patches changed from last iteration
      //
      if((old_active[0]  != medium->naflt) ||(old_active[1] != medium->naflt_con)){
	if((old_active[0] > medium->naflt)||(old_active[1]  > medium->naflt_con)){
	  fprintf(stderr,"rupture: error, number of active patches decreased\n");
	  fprintf(stderr,"rupture: from %i/%i to %i/%i\n",
		  old_active[0],old_active[1],medium->naflt,medium->naflt_con);
	  exit(-1);
	}
	active_modified = TRUE;
      }else
	active_modified = FALSE;
      //
      // decide on actions
      //
      if(medium->keep_slipping && (medium->iter > 1)){
	// we are keepign the faults slip and are in a second or later iteration
	if(active_modified){// number of active patches changed
	  //
	  // remove the old solution from the stress field to form the
	  // new one that includes other possibly activated patch
	  //
	  // therefore: don't mark quake but calc effect
	  //
	  add_solution(old_active[0],medium->sma,medium->xsol,medium->nameaf,
		       medium,fault,FALSE,TRUE,-1.0);
	  add_solution(old_active[1],medium->sma_con,medium->xsol_con,
		       medium->nameaf_con,medium,fault,FALSE,TRUE,-1.0);
	  if(medium->debug)
	    fprintf(stderr,"rupture: time: %11g it: %2i na: %3i nac: %3i removing test solution\n",
		    medium->time,medium->iter,old_active[0],old_active[1]);
	}
      }
      if(active_modified){
	/* 
	   at least one new fault pacth was activated at this timestep, 
	   number of constraints is given by nreq and nreq_con
	   solve the equation system at least once
	*/
	solve(medium,fault);
	// if not in keep_slipping mode, add the solution and mark quakes, else
	// add only for testing purposes but do not print to file yet
	mark_quakes = (medium->keep_slipping ? FALSE : TRUE);
	add_solution(medium->naflt,medium->sma,medium->xsol,medium->nameaf,
		     medium,fault,mark_quakes,TRUE,1.0);
	add_solution(medium->naflt_con,medium->sma_con,medium->xsol_con,
		     medium->nameaf_con,medium,fault,mark_quakes,TRUE,1.0);
	if(medium->debug){
	  fprintf(stderr,"rupture: time: %11g it: %2i na: %3i nac: %3i adding test solution\n",
		  medium->time,medium->iter,medium->naflt,medium->naflt_con);
	  print_solutions(medium->naflt,medium->nameaf,
			  fault,medium,"unconstrained");
	  print_solutions(medium->naflt_con,medium->nameaf_con,
			  fault,medium,"  constrained");
	}
      }
    }else
      active_modified = FALSE;
    /*
      
      re-start and check for activations since now 
      other patches might have been triggered. bailout when the number of active 
      patches has not changed since the previous activation
      
    */
  }while(active_modified && (medium->iter < RUPTURE_ITER_LIM));
  if(medium->iter >= RUPTURE_ITER_LIM){// check if we are above iteration limit
    fprintf(stderr,"solve: t: %g, maximum iter (%i) reached, %i/%i act. flts: ",
	    medium->time,medium->iter,medium->naflt,medium->naflt_con);
    for(i=0;i<medium->nrflt;i++)
      if(fault[i].active)fprintf(stderr,"%i ",i);
    fprintf(stderr,"\n");
#ifdef DEBUG
    // only bail-out if we are debugging, else continue program
    terminate(medium,fault);
#endif
  }
  medium->total_iter += medium->iter;
  // remember highest number of iterations
  if(medium->iter > medium->max_iter_realized)
    medium->max_iter_realized = medium->iter;
  if(medium->keep_slipping && faults_have_slipped){
    //
    // for keep slipping mode, mark the real solution to the system, ie add to 
    // slip and moment lists. else, this has been done at each iteration
    //
    if(medium->debug)
      fprintf(stderr,"rupture: time: %11g it: %2i na: %3i nac: %3i marking final solution\n",
	      medium->time,medium->iter,medium->naflt,medium->naflt_con);
    add_solution(medium->naflt,medium->sma,medium->xsol,medium->nameaf,
		 medium,fault,TRUE,FALSE,1.0);
    add_solution(medium->naflt_con,medium->sma_con,medium->xsol_con,
		 medium->nameaf_con,medium,fault,TRUE,FALSE,1.0);
  }
  return(faults_have_slipped);
}
/* 
   this routine checks the fracture criterion for a
   fault and adds the target stress_drop, if active, to
   the list of equations in the constrained or 
   unconstrained system 
*/

void fault_criterion(int flt,struct flt *fault,struct med *medium)
{
  int i,ip;
  COMP_PRECISION tstress_drop[3];
  //
  // if we activated this fault already, return to 
  // caller to avoid adding the fault twice to the system
  if(fault[flt].active)
    return;
  // decide on what to do depending on the fault slip
  // mode (constrained or unconstrained motion)
  switch(fault[flt].mode[0]){
  case COULOMB_STRIKE_SLIP:
  case STRIKE_SLIP:
  case COULOMB_DIP_SLIP:
  case DIP_SLIP:
  case COULOMB_MAXSDIR_SLIP:
  case MAXSDIR_SLIP:
  case NORMAL_SLIP:{
    /* 
       unconstrained slip modes, patches can move both
       ways
    */
    ip = medium->naflt * 3;
    if(fracture_criterion(flt,fault,tstress_drop,&medium->sma[ip],medium)){
      /* add an equation for each active slipping mode*/
      for(i=0;i<3;i++){
	if(medium->sma[ip+i]){// add to RHS and grow b and xsol
	  add_to_right_hand_side(tstress_drop[i],&medium->b,&medium->xsol,
				 &medium->nreq);
	}
      }
      /* add to list and increment counter */
      add_to_active_fault_list(flt,&medium->nameaf,&medium->naflt,&medium->sma);
    }
    break;
  }
  default:{
    /* 
       assume we have constrained slip modes that have
       to be solved with the non-negative solver
    */
    ip = medium->naflt_con * 3;
    if(fracture_criterion(flt,fault,tstress_drop,&medium->sma_con[ip],medium)){ 
      for(i=0;i<3;i++)
	if(medium->sma_con[ip+i]){
	  // add to RHS and grow xsol and b vectors
	  /* assign negative stress drop here */
	  if(fault[flt].mode[0]==STRIKE_SLIP_RIGHTLATERAL || 
	     fault[flt].mode[0]==COULOMB_STRIKE_SLIP_RIGHTLATERAL ||
	     fault[flt].mode[0]==DIP_SLIP_DOWNWARD || 
	     fault[flt].mode[0]==COULOMB_DIP_SLIP_DOWNWARD ||
	     fault[flt].mode[0]==NORMAL_SLIP_INWARD)
	    add_to_right_hand_side(-tstress_drop[i],&medium->b_con,
				   &medium->xsol_con,&medium->nreq_con);
	  else
	    add_to_right_hand_side(tstress_drop[i],&medium->b_con,
				   &medium->xsol_con,&medium->nreq_con);
	}
      add_to_active_fault_list(flt,&medium->nameaf_con,&medium->naflt_con,
			       &medium->sma_con);
    }
    break;
  }}
  return;
}
