/*

  constants for the blockinvert programs

  $Id: blockinvert.h,v 1.12 2004/10/05 01:10:04 becker Exp $


*/


#define BLOCK_NBASE 3			/* basic sub-col entries for A */
#define BLOCK_DIM 2			/* dimensions for solution,
					   in plane x-y as well as fault-local
					   where it's strike and normal 
					*/
/* 
   projection to use for fault local computations

   project_azi is basically the only one that works properly

*/
//#define FLT_ROT_PROJECTION PROJECT_AZI
#define FLT_ROT_PROJECTION OMERC_AZI

/* maximum number of blocks */
#define MAX_NR_BLOCK 20

/* 

this factor will ensure that omega Euler vectors for

 vel = omega \cross r 

 will be in deg/Myr if the input velocities are in mm/yr 
 and r is normalized to unity length

 gfac = 111.194926644559
*/
#define BLOCK_GFAC (RADIUS_EARTH/180.0*PI)

/* for assign_additional_sol_values, modes for assignment */
#define INIT_ADD_SOL 0
#define RANDOM_ADD_SOL 1
/* 
   filenames 
*/
/* input */
#define GPS_VEL_FILE "velc.data" /* GPS data */
#define FLTBLOCK_FILE "flt.block" /* geometry file */
#define STRESSDATA_FILE "stress.data" /* stress data */
#define STRESSDEPTH_FILE "stress.depth"	/* depth for stress eval for 
					   each point */
#define RIGIDBLOCKSITES_FILE "rsites.dat" /* file for the sites which 
					     are used to define a reference 
					     block */
/* output */
#define VELFITOUT_FILE "vel.fit" /* predicted velocities */
#define SLIPOUT_FILE "slip.fit" /* slip rates */
#define SOLOUT_FILE "solution.bin"
#define OMEGA_OUT "omega.dat"	/* solution vector and uncertainties in ASCII */
#define SFITOUT_FILE "stress.fit" /* predicted stresses */
#define SCSTRESSOUT_FILE "stress.dsca" /* scaled input (data) stresses */
#define FGEOOUT_FILE "fgeo.gmt" 	/* fault geometry for testing */
#define LDOUT_FILE "ld.fit"	/* locking depths */
#define CFACOUT_FILE "cfac.fit" /* coupling factors */
#define COVOUT_FILE "cov.asc"	/* covariance */
#define VELCOR_FILE "vcorr.dat"	/* velocities corrected for net
				   rotation */
#define VEL_GR_FILE "vsub.dat" /* net rotation used for correction */
#define EULER_POLE_FILE "epc.dat" /* file with euler poles
				     to constrain
				  */
#define SV_FILE "sv.dat"	/* singular value output */
#define CSD_FILE "csd.dat"	/* constraints on fault slip motion */
#define BLOCK_CENTROID_FILE "centroids.dat"
