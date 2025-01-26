/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: fracture_criterion.c,v 2.40 2003/01/13 06:44:24 becker Exp $
*/
#include "interact.h"
#include "properties.h"

#ifdef DEBUG
#define SAFETY_CHECK
#define CHECK_STRESS_DROP_ROUTINE_INPUT
#endif

my_boolean fracture_criterion(int flt, struct flt *fault,COMP_PRECISION *tstress_drop,
			   my_boolean *sma,struct med *medium)
{
  my_boolean activated,slipmode[2];
#ifdef SAFETY_CHECK
  if(fault[flt].active){
    fprintf(stderr,"fracture_criterion: patch %i is already active, avoid calling this routine twice for one fault\n",
	    flt);
    exit(-1);
  }
#endif
  /* 
     initialize with zeros  
  */
  fault[flt].active=FALSE;
  sma[STRIKE]=sma[DIP]=sma[NORMAL]=INACTIVE;
  tstress_drop[STRIKE]=tstress_drop[DIP]=
    tstress_drop[NORMAL]=0.0;
#ifdef LATENCY
  /*
    check for latency, ie. is the fault patch 
    ready to rupture again? if not, do not proceed
    with fracture criterion checks. has to have 
    ruptured at least once, ie. |u| != 0, for the check
    to make sense
  */
  if((norm_3d(fault[flt].u) != 0.0)&&
     ((medium->time-(COMP_PRECISION)fault[flt].last_activation_time) < LATENCY))
    return(FALSE);
#endif
  /* 
     determine if fault should rupture, and obtain target stress drop values 
     depending on the operational mode

     (Coulomb mode modifications will modify the target stress drop by the 
     product of dynamic stress and normal stress changes during rupture)

     WARNING: remember to change adjust_time_step, if variable timesteps are 
     wanted if you change the of Coulomb stress modes
  */
  switch(fault[flt].mode[0]){
  case COULOMB_STRIKE_SLIP_LEFTLATERAL:
  case COULOMB_STRIKE_SLIP_RIGHTLATERAL:
  case COULOMB_STRIKE_SLIP:
  case STRIKE_SLIP_LEFTLATERAL:
  case STRIKE_SLIP_RIGHTLATERAL:
  case STRIKE_SLIP:{
    /*
      
      first the pure strike-slip modes that don't allow any slip in dip 
      direction at the same time, but possible Coulomb correction
      
    */
    if(in_coulomb_compress_regime(fault[flt].s[NORMAL])){
      slipmode[STRIKE]=TRUE;slipmode[DIP]=FALSE;// toggle the slip modes
      two_dir_slip_check(sma,tstress_drop,slipmode,flt,fault,medium);
    }else{// extensional 
      if(medium->iter == 1)// deactivate if this is the first iteration
	// ie. if stresses are still 'static'
	deactivate_patch(flt,fault,medium);
    }
    break;
  }   
  /*
    
    pure dip-slip modes, similar to above

  */
  case COULOMB_DIP_SLIP_UPWARD:
  case COULOMB_DIP_SLIP_DOWNWARD:
  case COULOMB_DIP_SLIP:
  case DIP_SLIP_UPWARD:
  case DIP_SLIP_DOWNWARD:
  case DIP_SLIP:{
    if(in_coulomb_compress_regime(fault[flt].s[NORMAL])){
      slipmode[STRIKE]=FALSE;slipmode[DIP]=TRUE;
      two_dir_slip_check(sma,tstress_drop,slipmode,flt, fault,medium);
    }else{
      if(medium->iter == 1)
	deactivate_patch(flt,fault,medium);
    }
    break;
  }  
  /* 
     fault is active in the max shear stress (strike + dip) 
     mode, this means all slips in the plane of the 
     fault are allowed. otherwise similar to the 
     Coulomb part above. However, we will have to make a choice if zero
     stress drop is enforced or not (ie. is dip stress allowed to build 
     up, if we are slipping in strike direction, eg.)
  */
  case COULOMB_MAXSDIR_SLIP:
  case MAXSDIR_SLIP:{
    if(in_coulomb_compress_regime(fault[flt].s[NORMAL])){
      slipmode[STRIKE]=TRUE;slipmode[DIP]=TRUE;
      two_dir_slip_check(sma,tstress_drop,slipmode,flt,fault,medium);
    }else{ 
      if(medium->iter == 1)
	deactivate_patch(flt,fault,medium);
    }
    break;
  }
  /* 
     insert other operational modes for simulations here
  */
  /*

    EMPTY, for now

  */
  case INACTIVE:{
    /*
      this here just in case, might want to reactivate faults,
      e.g. when they become compressive again
    */
    break;
  }
  default:{
    fprintf(stderr,
	    "patch %3i, dir %i has operational mode %i which is undefined in fracture criterion\n",
	    flt,0,fault[flt].mode[0]);
    exit(-1);
    break;
  }}
  /* 
     set some important activation flags and stuff 
  */
  if(sma[STRIKE] || sma[DIP] || sma[NORMAL]){
    // the fault has been activated
    fault[flt].active=TRUE;
    activated=TRUE;
    // switch on activation flag of group, if not active yet
    if(!medium->fault_group[fault[flt].group].active){
      medium->fault_group[fault[flt].group].active=TRUE;
      // count up number of active groups (need this to decide if we should 
      // recheck the patches)
      medium->nr_active_groups++;
    }
#ifdef LATENCY
    fault[flt].last_activation_time=(float)medium->time;
#endif    
  }else{
    activated=FALSE;
  }
  return(activated);
}

/*
  routine to calculate the coulomb stress based on strike (slipmode[STRIKE]==1) 
  and dip (slipmode[DIP]==1) components of the shear stress. 
  if these slipmode are set, will allow slip in the correspondind directions

  returns the activation codes in array sma, and the stress drops in array 
  tstress_drop, IF THE FAULT IS CRITICAL 

  input is flt, the patch's number

  WARNING: ASSUMES THAT NORMAL STRESS IS NEGATIVE (IE COMPRESSIVE)

*/

void two_dir_slip_check(my_boolean *sma,COMP_PRECISION *tstress_drop,
			my_boolean *slipmode,int flt, struct flt *fault, 
			struct med *medium)
{
  COMP_PRECISION abs_shear_stress,stress_drop,overshoot,f1,f2;
  my_boolean critical;
  static COMP_PRECISION critical_stress;
  static my_boolean init=FALSE;
  // 
  // determine the critical shear stress level for activation (theory: 0)
  // based on the kind of simulation we are running
  //
  if(!init){
    if(medium->keep_slipping)
      critical_stress = CRITICAL_STRESS_EPS;
    else
      critical_stress = EXHAUSTIVE_CRITICAL_STRESS_EPS;
    fprintf(stderr,"fracture_criterion: tdsc: critical stress: %e * stress drop\n",
	    critical_stress/STRESS_DROP);
    if((critical_stress > 0) || (fabs(critical_stress) > 0.5*STRESS_DROP)){
      fprintf(stderr,"fracture_criterion: range inappropriate (stress_drop: %g)\n",
	      STRESS_DROP);
      exit(-1);
    }
    init = TRUE;
  }
  //
  // weighing allows us to reduce the slipping modes 
  // to strike or dip only
  // slipmode[DIP] or slipmode[STRIKE] (or both) 
  // should be non-zero
  // we don't use calc_abs_shear stress since this would involve extra 
  // calcs for stress at time t+dt
  //
  if(!slipmode[DIP]){// consider no dip stress
    abs_shear_stress=fabs(fault[flt].s[STRIKE]);
  }else if(!slipmode[STRIKE]){// no strike stress
    abs_shear_stress=fabs(fault[flt].s[DIP]);
  }else{// both non-zero, use vector length
    abs_shear_stress  = SQUARE(fault[flt].s[STRIKE]);
    abs_shear_stress += SQUARE(fault[flt].s[DIP]);
    abs_shear_stress  = sqrt(abs_shear_stress);
  }
  //
  // use vector to determine critical shear stress for Coulomb criterion
  //
  overshoot = coulomb_stress(abs_shear_stress,(COMP_PRECISION)fault[flt].mu_s,
			     fault[flt].s[NORMAL],medium->cohesion);
#ifdef SUPER_DEBUG
  fprintf(stderr,"fs: t: %11g it: %2i flt: %2i cs: %20.15e\n",
	  medium->time,medium->iter,flt,overshoot);
#endif
  //
  // stress drop as such is by definition positive
  //
  //
  // determine if super critical
  //
  if(overshoot >= critical_stress){
    if(medium->iter == 1){
      // normal activation
      stress_drop = cstress_drop(overshoot,
				 (COMP_PRECISION)(fault[flt].mu_s-fault[flt].mu_d),
				 fault[flt].s[NORMAL],medium->cohesion);

    }else{/* this fault became only critical after other patches
	     ruptured. hence, we use the old shear stress and drop to the 
	     dynamic value 
	  */
      stress_drop = fault[flt].taud;
    }
    critical = TRUE;
  }else if(medium->whole_fault_activations){
    // if whole fault activation mode is on
    // check if other patches of the fault are slipping
    // if so, activate this patch, too
    if(medium->fault_group[fault[flt].group].active){
      stress_drop = fault[flt].taud;
      critical=TRUE;
    }else{
      stress_drop = 0;
      critical = FALSE;
    }
  }else{
    critical = FALSE;
    stress_drop = 0.0;
  }
  /*
    
    activate the fault patch if
    a) the critical stress is > 0 and the resulting stress drop is larger than 
       the minimum stress drop allowed, or
    b) the fault is in whole fault activation mode, other patches have been 
       activated, and the shear stress is higher than the lower friction value
       times the normal stress (as checked by min_stress_drop comparison below)
  */
  if((critical) && (stress_drop > medium->min_stress_drop)){
    //
    // reset possible Coulomb corrections
    fault[flt].cf[STRIKE]=fault[flt].cf[DIP]=0.0;
    // obtain stress drop resolved on strike and
    // dip directions
    if(!slipmode[DIP]){// no dip stress change required 
      //                  (stress[DIP] might still be zero)
      //
      sma[STRIKE]=ACTIVATED;
      if(fault[flt].s[STRIKE] > 0){// need negative stress drop
	tstress_drop[STRIKE] = -stress_drop;
	if(fault[flt].mode[0] >= OS_C_OFFSET)// Coulomb correction
	  fault[flt].cf[STRIKE]=  (COMP_PRECISION)fault[flt].mu_d;
      }else{// need the positive stress drop 
	tstress_drop[STRIKE] = stress_drop;
	if(fault[flt].mode[0] >= OS_C_OFFSET)// negative Coulomb correction
	  fault[flt].cf[STRIKE]= -(COMP_PRECISION)fault[flt].mu_d;
      }
#ifdef DEBUG
      fprintf(stderr,"fc: strike: t: %11g, %4i/%4i: s_s/d/n/c:%10.6e/%10.6e/%10.6e/%10.6e sd_s: %10.6e\n",
	      medium->time,flt,fault[flt].group,
	      fault[flt].s[STRIKE],fault[flt].s[DIP],fault[flt].s[NORMAL],
	      coulomb_stress(abs_shear_stress,(COMP_PRECISION)fault[flt].mu_s,
			     fault[flt].s[NORMAL],medium->cohesion),
	      tstress_drop[STRIKE]);
#endif
    }else if(!slipmode[STRIKE]){// no strike modes needed, repeat from above
      sma[DIP]=ACTIVATED;
      if(fault[flt].s[DIP] > 0){
	tstress_drop[DIP] = -stress_drop;
	if(fault[flt].mode[0] >= OS_C_OFFSET)
	  fault[flt].cf[DIP]=  (COMP_PRECISION)fault[flt].mu_d;
      }else{
	tstress_drop[DIP] = stress_drop;
	if(fault[flt].mode[0] >= OS_C_OFFSET)
	  fault[flt].cf[DIP]= -(COMP_PRECISION)fault[flt].mu_d;
      }
#ifdef DEBUG
      fprintf(stderr,"fc: dip : t: %11g, %4i/%4i s_s/d/c:%10.6e/%10.6e/%10.6e sd_d: %10.6e\n",
	      medium->time,flt,fault[flt].group,fault[flt].s[STRIKE],
	      fault[flt].s[DIP],coulomb_stress(abs_shear_stress,
					       (COMP_PRECISION)fault[flt].mu_s,
					       fault[flt].s[NORMAL],
					       medium->cohesion),tstress_drop[DIP]);
#endif
    }else{
      //
      // get the absolute target stress drop in both ways
      // (even if the other stress drop is set to zero, we will incorporate that
      // here)
      //
      get_maxsdir_stress_drops(fault[flt].s,stress_drop,tstress_drop,
			       fault[flt].f_initial,&f1,&f2);
      sma[STRIKE] = sma[DIP] = ACTIVATED;
      if(fault[flt].mode[0] >= OS_C_OFFSET){
	if(fault[flt].s[STRIKE] > 0)
	  fault[flt].cf[STRIKE] =  (COMP_PRECISION)fault[flt].mu_d * f1;
	else
	  fault[flt].cf[STRIKE] = -(COMP_PRECISION)fault[flt].mu_d * f1;
	if(fault[flt].s[DIP] > 0)
	  fault[flt].cf[DIP] =     (COMP_PRECISION)fault[flt].mu_d * f2;
	else
	  fault[flt].cf[DIP] =    -(COMP_PRECISION)fault[flt].mu_d * f2;
      }
#ifdef DEBUG
    fprintf(stderr,"fc: mixed: t: %11g, %4i/%4i s_s/d/c:%12.5e/%12.5e/%12.5e sd_s/d/f:%12.5e/%12.5e/%5.2f\n",
	    medium->time,flt,fault[flt].group,fault[flt].s[STRIKE],
	    fault[flt].s[DIP],coulomb_stress(abs_shear_stress,
					     (COMP_PRECISION)fault[flt].mu_s,
					     fault[flt].s[NORMAL],medium->cohesion),
	    tstress_drop[STRIKE],tstress_drop[DIP],
	    tstress_drop[STRIKE]/tstress_drop[DIP]);
#endif
    }
  }
#ifdef DEBUG
  //fprintf(stderr,"fc: done\n");
#endif
}

void deactivate_patch(int flt, struct flt *fault, struct med *medium)
{
  /* 
     the fault is in the extensional regime, deactivate it 
  */
  fault[flt].mode[0]=INACTIVE;
  fprintf(stderr,"fracture_criterion: stress on patch %3i is tensile, deactivated at time %10.5g\n",
	  flt,medium->time);
  if(medium->whole_fault_deactivations){
    fprintf(stderr,"fracture_criterion: deactivating group %3i as a consequence\n",fault[flt].group);
    deactivate_group(flt,fault[flt].group, fault,medium);
  }
}

void deactivate_group(int patch, int group, struct flt *fault, struct med *medium)
{
  int i;
  for(i=0;i<medium->nrflt;i++){
    if(fault[i].group==group && i != patch){
      if(fault[i].active)
	fprintf(stderr,"deactivate_group: patch %5i is active at time %11g, will deactivate at next step\n",
		i,medium->time);
      else
	fprintf(stderr,"deactivate_group: patch %5i is deactivated since in group %i\n",
		i,fault[i].group);
      fault[i].mode[0]=INACTIVE;
    }
  }
}

/* 
   return the absolute value of the shear stress
   and the absolute value of the shear stress plus
   the dt0 increment

   shear stress can be strike or dip component or 
   a mixture of both depending on the rupture mode
   of the faults
*/
int calc_absolute_shear_stress_and_inc(COMP_PRECISION *abs_tau,
				       COMP_PRECISION *abs_tau1,
				       int flt,struct flt *fault)
{
  COMP_PRECISION tmpd1,tmpd2;
  switch(fault[flt].mode[0]){
  case COULOMB_STRIKE_SLIP_LEFTLATERAL:
  case COULOMB_STRIKE_SLIP_RIGHTLATERAL:
  case COULOMB_STRIKE_SLIP:
  case STRIKE_SLIP_LEFTLATERAL:
  case STRIKE_SLIP_RIGHTLATERAL:
  case STRIKE_SLIP:{
    *abs_tau =fabs(fault[flt].s[STRIKE]);
    *abs_tau1=fabs(fault[flt].sinc[STRIKE]+fault[flt].s[STRIKE]);
    return(STRIKE_SLIP);
  }
  case COULOMB_DIP_SLIP_UPWARD:
  case COULOMB_DIP_SLIP_DOWNWARD:
  case COULOMB_DIP_SLIP:
  case DIP_SLIP_UPWARD:
  case DIP_SLIP_DOWNWARD:
  case DIP_SLIP:{
    *abs_tau =fabs(fault[flt].s[DIP]);
    *abs_tau1=fabs(fault[flt].sinc[DIP]+fault[flt].s[DIP]);
    return(DIP_SLIP);
  }
  case COULOMB_MAXSDIR_SLIP:
  case MAXSDIR_SLIP:{
    *abs_tau =sqrt(SQUARE(fault[flt].s[DIP])+SQUARE(fault[flt].s[STRIKE]));
    tmpd1=fault[flt].sinc[DIP]   +fault[flt].s[DIP];
    tmpd2=fault[flt].sinc[STRIKE]+fault[flt].s[STRIKE];
    *abs_tau1=sqrt(SQUARE(tmpd1)+SQUARE(tmpd2));
    return(MAXSDIR_SLIP);
  }
  default:{
    *abs_tau=0.0;*abs_tau1=0.0;
    return(INACTIVE);
  }}
}
// same as above without incremental stress 
int calc_absolute_shear_stress(COMP_PRECISION *abs_tau,int flt,struct flt *fault)
{
  switch(fault[flt].mode[0]){
  case COULOMB_STRIKE_SLIP_LEFTLATERAL:
  case COULOMB_STRIKE_SLIP_RIGHTLATERAL:
  case COULOMB_STRIKE_SLIP:
  case STRIKE_SLIP_LEFTLATERAL:
  case STRIKE_SLIP_RIGHTLATERAL:
  case STRIKE_SLIP:{
    *abs_tau =fabs(fault[flt].s[STRIKE]);
    return(STRIKE_SLIP);
  }
  case COULOMB_DIP_SLIP_UPWARD:
  case COULOMB_DIP_SLIP_DOWNWARD:
  case COULOMB_DIP_SLIP:
  case DIP_SLIP_UPWARD:
  case DIP_SLIP_DOWNWARD:
  case DIP_SLIP:{
    *abs_tau =fabs(fault[flt].s[DIP]);
    return(DIP_SLIP);
  }
  case COULOMB_MAXSDIR_SLIP:
  case MAXSDIR_SLIP:{
    *abs_tau =sqrt(SQUARE(fault[flt].s[DIP])+SQUARE(fault[flt].s[STRIKE]));
    return(MAXSDIR_SLIP);
  }
  default:{
    *abs_tau=0.0;
    return(INACTIVE);
  }}
}






