/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: adjust_time_step.c,v 1.17 2003/01/13 06:41:01 becker Exp $
*/
#include "interact.h"
/*

  adjust the time step based on our knowledge
  about the next time of failure

*/


//#define NOTIFY_WHEN_TIMESTEP_CHANGED


void adjust_time_step(struct flt *fault, struct med *medium)
{
  int i,slip_mode;
  COMP_PRECISION normal_dt0,abs_tau,abs_tau1,fac,
    ct,cdt0,loc_time_to_failure,loc_time_to_sign_change;
  COMP_PRECISION maxct=0;
  static COMP_PRECISION old_time_step;
  static my_boolean init=FALSE;
  if(!init){
    old_time_step=medium->dt;
    init=TRUE;
  }
  // determine total time to failure for the medium
  medium->time_to_failure=FLT_MAX;
  // check all active (ie. still participating faults) faults
  for(i=0;i<medium->nrflt;i++){
    if(fault[i].mode[0]==INACTIVE)
      continue;
    // obtain the shear stress values at time t and t+dt0
    slip_mode=calc_absolute_shear_stress_and_inc(&abs_tau,&abs_tau1,i,fault);
    // coulomb stress at time t
    ct=coulomb_stress(abs_tau,(COMP_PRECISION)fault[i].mu_s,
		      fault[i].s[NORMAL],medium->cohesion);
    // positive coulomb stress
    if(ct >= 0.0){
      // check for the stress drop being below limit
      // if, so continue with loop and don't use this
      // predicted time to failure
      if(cstress_drop(ct,(COMP_PRECISION)(fault[i].mu_s-fault[i].mu_d),
		      fault[i].s[NORMAL],medium->cohesion) < 
	 medium->min_stress_drop)
	continue;
      if(medium->time == 0){
	fprintf(stderr,"adjust_time_step: WARNING: patch %i at time 0 should have already ruptured\n",
		i);
	fprintf(stderr,"adjust_time_step: Coulomb stress: %g\n",ct);
	loc_time_to_failure = MIN_TIME_STEP;
      }else{
	if((ct > 1.0e-8)&&(medium->dt>1.0e-8)){
	  fprintf(stderr,"adjust_time_step: WARNING: encountered positive Coulomb stress, t: %g dt: %g fault: %i cs: %g\n",
		  medium->time,medium->dt,i,ct);
	}
	// this to allow rupture really soon
	loc_time_to_failure = MIN_TIME_STEP;
      }
    }else{// coulomb stress negative at this time
      // normal stress at time t + dt0
      normal_dt0=fault[i].s[NORMAL]+fault[i].sinc[NORMAL];
      // coulomb stress at time t + dt0
      cdt0=coulomb_stress(abs_tau1,(COMP_PRECISION)fault[i].mu_s,
			  normal_dt0,medium->cohesion);
      /*
	the coulomb stress on this fault will be decreasing 
	or stay the same. however, we still want to check if one
	of the shear stresses changes the sign
      */
      if(cdt0 <= ct){
	loc_time_to_failure = FLT_MAX;
      }else{			
	/* 
	   calculate time to failure based on \sigma_c == 0
	   FOR THIS TO WORK medium->dt0 HAS TO BE SET TO UNITY!
	   (formula is actually  ttf=(-ct * medium->dt0)/
	                             (cdt0 - ct)
	*/
	loc_time_to_failure  = -ct/(cdt0 - ct)+EPS_COMP_PREC;
      }
      /*
	check if the appropriate shear stress would change 
	sign in the predicted time to failure. if so, use
	the time of the sign change

	HERE, AGAIN, WE HAVE ASSUMED THAT dt0 IS UNITY
      */
      loc_time_to_sign_change=FLT_MAX;
      if((slip_mode == STRIKE_SLIP) || (slip_mode == MAXSDIR_SLIP))
	if((fault[i].s[STRIKE] != 0.0) && (((fac=fault[i].sinc[STRIKE]/
	      fault[i].s[STRIKE])*loc_time_to_failure) < -1.0))
	  loc_time_to_sign_change= -1.0/fac;
      if((slip_mode == DIP_SLIP)||(slip_mode == MAXSDIR_SLIP))
	if((fault[i].s[DIP] != 0.0) && (((fac=fault[i].sinc[DIP]/
	      fault[i].s[DIP])*loc_time_to_failure) < -1.0))
	  loc_time_to_sign_change= -1.0/fac;
      if(loc_time_to_sign_change < loc_time_to_failure)
	loc_time_to_failure=loc_time_to_sign_change;
    }
    if(loc_time_to_failure < medium->time_to_failure){
      medium->time_to_failure = loc_time_to_failure;
      if(medium->debug)
	maxct=ct;
    }
  }
  //
  // adjust time stepping for max and min time steps
  // as given by MIN_TIME_STEP and medium->print_interval
  //
  if(medium->time_to_failure < medium->print_interval)
    medium->dt=(medium->time_to_failure >= MIN_TIME_STEP)?(medium->time_to_failure):(MIN_TIME_STEP);
  else
    medium->dt=medium->print_interval;

  if(medium->debug)
    if(old_time_step != medium->dt){
      fprintf(stderr,"ats: t: %15.6e, iter: %3i ttf: %15.6e cs: %20.14e dtnew %20.10e Ddt %20.10e, next t: %20.15e\n",
	      medium->time,medium->nr_timesteps,
	      medium->time_to_failure,maxct,
	      medium->dt,medium->dt-old_time_step,medium->time+medium->dt);
      old_time_step=medium->dt;
    }
  if(medium->time + medium->dt == medium->time){
    fprintf(stderr,"adjust_time_step: WARNING: time increment of %20.16e (min: %20.16e) at time %g too small\n",
	    medium->dt,MIN_TIME_STEP,medium->time);
    medium->dt += MIN_TIME_STEP;
  }
}


