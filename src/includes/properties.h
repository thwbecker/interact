/* 
   constant material properties for interact

   individual values might be modified by reading in, eg,
   fault properties when using the -f switch

   REMEMBER TO RECOMPILE ALL FILES AFTER MAKING
   CHANGES


   REMEMBER to check the whole file since changing, say, the 
   shear modulus once might not affect some of the quantities further 
   down

*/
/* STRESS_DROP_NORM flag: if defined, will normalize stresses by
 stress drop if not, shear modulus will be set to unity */

#define STRESS_DROP_NORM 	/* old default */

#ifdef STRESS_DROP_NORM
/* set by stress drop */
#define SHEAR_MODULUS_DEF 1.0e4
/* 
   define stress drop in consistent units
   for Earth, \Delta \sigma \sim 3e6 Pa, \mu= \sim 3e10 Pa,
   therefore, if we set the characteristic stress by the shear modulus 
   to unity, stress_drop should be 1e-4
*/
#define STRESS_DROP_DEF 1.0
#else
/*   characteristic stress set by shear modulus  */
#define SHEAR_MODULUS_DEF  1.0
#define STRESS_DROP_DEF 1.0e-4
#endif
#define POISSON_NU_DEF 0.25

/* 


   static coeff of friction for all faults if
   not read in by separate file


*/
#define STATIC_MU 0.6
/*
  difference between the static and dynamic 
  value of friction. again, can be overwritten by 
  file
 */
#define DELTA_MU 0.1
/* 
   cahracteristic dimension of faults, ie. length scale
   should be on the order of the fault's width if things
   are to work out reasonably
*/
#define CHAR_FAULT_DIM 1.0	/* best to leave that at unity */
/* 
   the pressure scale is given via the characteristic
   stress drop such that in reality
              stress drop               3 x 10^6 Pa
   pressure = -------------------- = ------------------ = 10^-3
              shear modul*delta mu   3 x 10^10 Pa 10^-1

   this will be attained at charcteristic depth, ie. at z=w

   pressure is positive in compression and can be changed
   by a command line option
   
   Note that PRESSURE_DEF is also used by create_random_stress_file
   to determine if patches would start out with positive Coulomb
   stress in a loading simulation

*/
#define PRESSURE_DEF (STRESS_DROP_DEF/DELTA_MU)
/* 
   define HYDROSTATIC_PRESSURE if you want depth dependent pressure. 
   if undefined, pressure will be constant
   if defined, pressure will be given by 

   p = pressure * (- x[Z])/HYDROSTATIC_PRESSURE

   that is, the depth of the fault patch over the value
   of HYDROSTATIC_PRESSURE gives the fraction of pressure
   applying at this fault patch. in this case,
   HYDROSTATIC_PRESSURE gives the length scale of the depth
   dependence

   suitable choice would be CHAR_FAULT_DIM
*/


/* 
   for Coulomb criterion, as defined from Byerlee's law
   
   60 Mpa were found to be appropriate for normal stresses 
   over 200MPa, ie depths larger than 

   correspond to 2 10^2 in non-dimensionalzied units,
   ie. if the shear modulus determines the stress scale

   COHESION is defined POSITIVE, ie. >= 0 !!!

   WARNING: this can be overwritten by a runtime
   switch

*/
#define COHESION_DEF (5.0*STRESS_DROP_DEF)
/* #define COHESION_DEF 0.0 */

/* we consider a compressive stress of amplitude TENSILE_RANGE 
   to be needed for strike-slip (give positive TENSILE_RANGE as compressive!)
   faulting. if you want faulting up to 0 compressive stress, set to zero  */
#define TENSILE_RANGE EPS_COMP_PREC
/* 
   stressing rate for calc_stress.c, routine 
   background_stress()

   sigma_xy =  (PRE_STRESSING_TIME_OFFSET+time)*stressing_rate

   if stressing rate is set to STRESS_DROP_DEF, then
   the characteristic time is unity, ie. one seismic 
   cycle

   WARNING: defaults can be overwritten by smlin.in!

*/
#define STRESSING_RATE STRESS_DROP_DEF




/* 
   visco-elastic hereditary kernel machinery (rsf_ve.c) compile-time
   settings.  These are fixed here rather than run-time options
   because good values are geometry independent and unlikely to need
   changing:

   VE_NIMG_DEF: image-family depth of the Laplace-domain sampling of
   the layered antiplane kernel.  Each sample at s_k truncates the
   reflection series at Gamma(s_k)^n; the smallest sample s_min =
   b/100 has Gamma ~ 0.99, so n = 600 leaves a relative truncation
   of exp(-6) ~ 2.5e-3 in the SLOWEST kernel content (the relaxed /
   late-interseismic part); the fast and coseismic parts converge
   after a few images.  Halving it doubles that bias; assembly cost
   is linear in it.

   VE_FIT_TOL_DEF: abort threshold on the worst per-pair held-out
   residual of the Prony fit.  A pure safety gate: it changes no
   numbers, it only refuses to run when the rate ladder cannot
   represent the sampled kernel (typical healthy values are 3e-5,
   two orders below the gate).
*/
#define VE_NIMG_DEF 600
#define VE_FIT_TOL_DEF 1e-3
