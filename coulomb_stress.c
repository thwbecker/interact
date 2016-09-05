#include "interact.h"
#include "properties.h"

//
// calculate Coulomb stress assuming that the shear stress value is positive
// and that the normal stress value is negative, ie. under compression
//

COMP_PRECISION coulomb_stress(COMP_PRECISION absolute_val_shear_stress,
			      COMP_PRECISION mu_s,
			      COMP_PRECISION compressive_normal_stress,
			      COMP_PRECISION cohesion)
{
  COMP_PRECISION cs;
#ifdef ADD_COULOMB_STRESS_NOISE
  static long seed=-1;
#endif
#ifdef DEBUG
  if(absolute_val_shear_stress < 0){
    fprintf(stderr,"coulomb_stress: we shouldn't call with negative shear stress\n");
    exit(-1);
  }
#endif
  cs  = absolute_val_shear_stress;
  cs += mu_s * compressive_normal_stress;
#ifndef NO_COHESION
  cs -= cohesion;
#endif
#ifdef ADD_COULOMB_STRESS_NOISE
  /*
    add a uniformly distributed random "noise"
    with peak to trough amplitude ADD_COULOMB_STRESS_NOISE
  */
  cs += (-0.5+myrandnr(1.0,&seed))*ADD_COULOMB_STRESS_NOISE;
#endif
  return(cs);
}
/*
  
  calculate the stress drop that would result
  for a given set of parameters
*/
COMP_PRECISION cstress_drop(COMP_PRECISION overshoot,COMP_PRECISION delta_mu,
			    COMP_PRECISION normal_stress,COMP_PRECISION cohesion)
{
  COMP_PRECISION sd;
#ifdef CHECK_STRESS_DROP_ROUTINE_INPUT
  if(overshoot < -1e-7){
    fprintf(stderr,"cstress_drop: WARNING: overshoot < 1.0e-7, %g\n",overshoot);
  }
  if(normal_stress >= EPS_COMP_PREC){
    fprintf(stderr,"cstress_drop: WARNING: normal stress is not negative (%g)\n",
	    normal_stress);
  }
#endif
  sd  = overshoot;
  sd -= delta_mu * normal_stress;// normal stress should be 0
  sd += cohesion;
#ifdef CHECK_STRESS_DROP_ROUTINE_INPUT
  if(sd<=0.0){
    fprintf(stderr,"cstress_drop: error, stress drop negative, %g\n",
	    sd);
    exit(-1);
  }
#endif
  return sd;
}

// routine to check if the fault is in the normal
// stress regime requeired for a coulomb stress
// fracture criterion, ie. at least slighly
// compressive, normal stress smaller than zero
//
my_boolean in_coulomb_compress_regime(COMP_PRECISION 
				   normal_stress)
{
  if(normal_stress <= -TENSILE_RANGE)
    return(TRUE);
  else
    return(FALSE);
}
