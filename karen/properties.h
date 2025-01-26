/* 
   constant material properties for interact

   individual values might be modified by reading in, eg,
   fault properties when using the -f switch

*/

/* characteristic stress set by shear modul */
#define SHEAR_MODULUS  1.0

/* 
   define stress drop as multiples of the shear modul
*/

#define STRESS_DROP (SHEAR_MODULUS*1.0e-4*1.5)

#define SHEAR_MODULUS2 (2.0*SHEAR_MODULUS)
/* 
   lambda = mu (2 nu)/(1-2*nu) 
   for nu = 1/4:
   lambda = mu
*/
#define LAMBDA SHEAR_MODULUS
/* 
   
   alpha = (lambda+mu)/(lambda+2 mu) 

   where lambda and mu are lame constants 
   (mu is the shear modulus)
*/
#define ALPHA ((LAMBDA+SHEAR_MODULUS)/(LAMBDA+2.0*SHEAR_MODULUS))


/* 
   static coeff of friction for all faults if
   not read in by sparate file
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
#define CHAR_FAULT_DIM 1.0

/* 
   the pressure scale is given via the characteristic
   stress drop such that in reality
              stress drop               3 x 10^6 Pa
   pressure = -------------------- = ------------------ = 10^-3
              shear modul*delta mu   3 x 10^10 Pa 10^-1

   this will be attained at charcteristic depth, ie. at z=w
   
*/
#define PRESSURE (STRESS_DROP/DELTA_MU)

/* 
   define for depth dependent pressure. 
   if undefined, pressure will be constant
   if defined, pressure will be given by 

   p = PRESSURE * (- x[Z])/HYDROSTATIC_PRESSURE

   that is, the depth of the fault patch over the value
   of HYDROSTATIC_PRESSURE gives the fraction of PRESSURE
   applying at this fault patch

   suitable choice would be CHAR_FAULT_DIM
*/
//#define HYDROSTATIC_PRESSURE  CHAR_FAULT_DIM


/* 
   for Coulomb criterion, as defined from Byerlee's law
   
   60 Mpa were found to be appropriate for normal stresses 
   over 200MPa, ie depths larger than 

   correspond to 2 10^2 in non-dimensionalzied units,
   ie. if the shear modulus determines the stress scale

   COHESION is defined POSITIVE, ie. >= 0 !!!
*/
//#define COHESION_DEF (5.0*STRESS_DROP)
#define COHESION_DEF 0.0 

// we consider a compressive stress of 
// amplitude TENSILE_RANGE to be needed for strike-slip
// faulting. if you want faulting up to 0 compressive 
// stress, set to zero 
#define TENSILE_RANGE (STRESS_DROP)



/* 
   stressing rate for calc_stress.c, routine 
   background_stress()

   sigma_xy =  (PRE_STRESSING_TIME_OFFSET+time)*stressing_rate

   if stressing rate is set to STRESS_DROP, then
   the characteristic time is unity, ie. one seismic 
   cycle

   WARNING: defaults can be overwritten by smlin.in!

*/
#define STRESSING_RATE STRESS_DROP
#define PRE_STRESSING_TIME_OFFSET 0.0

