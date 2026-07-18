/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: check_interaction.c,v 1.10 2002/10/08 19:24:44 tbecker Exp $
*/
#include "interact.h"

/*

  calculate interaction coefficients and check if 
  the Coulomb stress change on any fault is larger than
  the Coulomb stress change on the fault itself for
  any given mode of rupture

  these routines here are basically only wrappers,
  the main work is done in solve_mode_dependend.c

*/

my_boolean incrementally_check_interaction_coefficients(
        struct flt *fault,int nrflt,struct med *medium,
	my_boolean call_out_loud, int *evil_pair,
	COMP_PRECISION max_distance)
{
  static int nrflts_to_stay=0;
#ifdef DEBUG
  call_out_loud=TRUE;
#endif
  if(check_coulomb_stress_feedback(nrflt,nrflts_to_stay,
				   fault,medium, call_out_loud,
				   CALC_I_COEFF_NOW,TRUE,
				   evil_pair,max_distance))
    return(TRUE);
  else{
    nrflts_to_stay = nrflt;
   
    return(FALSE);
  }
}

/*
  
  calculate the actual coulomb stress changes due to 
  self-interaction and interaction and return TRUE if
  a positive feedback loop is detected

  the real routines are in solve_mode_dependend.c

*/

my_boolean check_coulomb_stress_feedback(int nrflt,
				      int nrflts_to_stay,
				      struct flt *fault,
				      struct med *medium,
				      my_boolean call_out_loud,
				      int imatmode,
				      my_boolean bailout,
				      int *evil_pair,
				      COMP_PRECISION max_distance)
{
  my_boolean hit;
  switch(imatmode){
  case I_MAT_IN_MEMORY:{
    hit=check_coulomb_stress_feedback_1(nrflt,nrflts_to_stay,fault,
					medium,call_out_loud,
					bailout,evil_pair,
					max_distance);
    break;
  }
  case I_MAT_ON_FILE:{
    hit=check_coulomb_stress_feedback_2(nrflt,nrflts_to_stay,fault,
					medium,call_out_loud,
					bailout,evil_pair,
					max_distance);
    break;
  }
  case CALC_I_COEFF_NOW:{
    hit=check_coulomb_stress_feedback_3(nrflt,nrflts_to_stay,fault,
					medium,call_out_loud,
					bailout,evil_pair,
					max_distance);
    break;
  }
  case SPARSE_I_MAT:{
    hit=check_coulomb_stress_feedback_4(nrflt,nrflts_to_stay,fault,
					medium,call_out_loud,
					bailout,evil_pair,
					max_distance);
    break;
  }
  default:{
    fprintf(stderr,"check_coulomb_stress_feedback: I matrix access code %i not defined\n",
	    imatmode);
    exit(-1);
  }}
  return(hit);
}




