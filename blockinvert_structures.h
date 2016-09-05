/*


  fault structure as used by blockinvert



*/
// these are building block matrices
struct vcm{
  COMP_PRECISION vc[3][3];/* 
				  displacement (rate) vs[i][j]
				  where j indicates the velocity 
				  component in X,Y,Z space 
				  and i the slip component in 
				  STRIKE, DIP, NORMAL space
			  */
};
struct scm{			/* stress interaction matrix */
  COMP_PRECISION sc[3][6];	/* sc[i][j] is the s[j] 
				   stress matrix for slip in i
				   direction. s[j] is in short (6)
				   component format
				*/
};
struct bflt{
  COMP_PRECISION x[3];		/* 
				   geographical coordinates
				   of center of patch

				   X: lon Y: lat 
				   Z: is -depth [km]

				*/
  COMP_PRECISION *xc,*pbase;	/* 
				   cartesian coordinates 
				   of fault and polar base vector 
				   (only used if spherical is set)
				*/
  COMP_PRECISION px[3];		/* 
				   projected coordinates of x[]
				   in mercator frame
				*/

  COMP_PRECISION sx[2];		/* 
				   geographical coordinates 
				   of midpoint at surface (can
				   be different from x[X,Y] if
				   fault has dip != 90
				*/
  COMP_PRECISION ex[2][2];	/* 
				   geographic coordinates of
				   left and right end point 
				   ex[][X,Y] = lon,lat
				   of the surface trace
				   
				*/
  COMP_PRECISION evec[9]; /* strike (0,1,2) normal (3,4,5) and 
			     dip(6,7,8) vectors, normalized */
  COMP_PRECISION azi,alpha,sa,ca,dip,sd,cd; /* 
					       azimuth in degrees,
					       alpha = 90 - azi
					       sin and cos of alpha 
					       dip in degrees,
					       sin and cos of dip
					    */
  COMP_PRECISION l,w;		/* half length (IN THE PROJECTED FRAME) 
				   and width */
  int block[2];			/* blocks on which this fault
				   borders */
  struct vcm *v;		/* 
				   for each observational 
				   point, we will construct a 
				   vcm matrix 
				*/
  struct scm *s;		/* stress interaction matrix */
  COMP_PRECISION lfac;		/* locking factor, 1: completely 
				   locked 0: freely slipping. 
				   gets multiplied to vcm and scm's
				*/
  COMP_PRECISION ld;		/* locking depth */
  my_boolean vertical;		/* if fault->dip == 90, will 
				   use opening mode, else dip mode
				*/
  int orig_code;		/* 
				   original number of the input fault 
				   from which this subsegment might be 
				   generated
				*/
  COMP_PRECISION sdamp[2];	/* slip damping, strike & normal */
  my_boolean use_damp[2];		/* actually use this damping value */
  short int sc[2];		/* constrained slip directions

				for strike:
				positive values of slip mean only left-lateral fault,
				negative only right-lateral; 
				For dip: positive values of slip mean only thrust (up-dip motion)
				fault, negative normal (down-dip motion) fault;
				For normal: positive values of slip mean explosive source,
				negative implosive;

				*/

};
/* 
   general projection structure 
*/
struct prj{
  int type;
  COMP_PRECISION clon,clat;	/* projection center in deg */
  COMP_PRECISION azi;		/* azimuth in degrees */
  COMP_PRECISION lat1,lat2;	/* standard parallels, deg */
};

/* 

block structure (IF YOU CHANGE THIS, MAKE SURE TO ADJUST COPY_BLOCK IN
BLOCK_MATRIX.C)

*/
struct bck{
  COMP_PRECISION center[3];	/* average block coordinates in cartesian 
				   reference frame */
  COMP_PRECISION xrigid[3];	/* 
				   rigid body solution from input velocities,
				   possibly in mofied reference frame 
				*/
  my_boolean fixed;		/* if block is fixed */
  my_boolean rot_c;		/* if block euler motion is constrained */
  int pcnt;			/* number of GPS points in block */
};
/* 

define the model structure for blockinvert type programs

*/

struct bmd{
  COMP_PRECISION *v,*sigv,*vmod,*rho; /* spherical velocities
					 v:(ve, vn) and stresses [6],
					 velocity uncertainties, and
					 model predictions rho:
					 correlation
				      */
  COMP_PRECISION *vc,*vmodc,*sigvc;	/* 
					   cartesian velocities and 
					   velocity uncertainties
					*/
  my_boolean cvel_init,cvelmod_init;/* initialized cartesian observed
				       and model velocities
				*/
  COMP_PRECISION *pbase;	/* polar basis vectors for all GPS
				   velocities */
  my_boolean pbase_init;		/* TRUE, if polar basis
					   initialized
					 */
  COMP_PRECISION *gx,*sx;	/* locations of GPS and stress
				   observations in lon/lat */

  COMP_PRECISION *gpx;		/* projected velocity locations in
				   projection space
				 */
  COMP_PRECISION *gcx;		/* GPS locations in cartesian space
				 */
  /* reference sites for reference block motion */
  COMP_PRECISION *fblock_sites;
  int fblock_nsites;		/* number of ref sites used */
  my_boolean select_fblock_sites; /* switch */
  /* 
     matrices for the solution, covariance and such
  */
  COMP_PRECISION *a,*f,*gf,*e,
    *kmat,*d,*g,*cov,*imat,*k2mat,*sigma;	
  COMP_PRECISION *xsol;	/* solution vector */
  COMP_PRECISION *saved_stress;	/* saved stress vector */
  COMP_PRECISION *stress_depths; /* depth at which to evaluate 
				    stresses */
  COMP_PRECISION    nslip_damp;	/* global damping factor */
  int nfdamp;			/* 
				   nr of damping parameters,
				   between 0, nflt and nflt*nslip
				*/
  COMP_PRECISION xdamp;		/* damping factor for solution */
  int nxdamp;			/* 
				   normal damping of solution vector 0
				   or n
				*/
  int nrgp,nrsp,nrb,nflt,nslip,mgd,nsnf,na,m1,m2,n,m,ncov,ms,
    nrbc,nc;
  /* 
     nrgp:  nr of velocity points
     nrsp:  nr of stress points
     nrb:   nr of blocks
     nrcb:  nr of blocks with constrained Euler poles
     nrf:   nr of faults 
     nslip: nr of slip directions,
     nsnf:  nslip * nrf
     mgd:   nr of GPS observations * BLOCK_DIM
     m1:    nr of rows of A (mgd for cartesian, nrgp*3 for spherical)
     m2:    nr of stress observations * 6 
     m:     m1 + m2
     ms:    m + nsnf (if slip damping activated), else m
     nc:    number of contrained parameters = nrbc * 3
  */
  int *bcode;			/* boundary code of each vel
				   observation */
  struct bflt *fault;		/* fault structure */
  struct bck *block;		/* block structure */
  my_boolean block_init;

  int first_c;	/* first  constrained block in 0..nrb-1 notation */
  /* 
     for changing the reference frame 
  */
  COMP_PRECISION *vcorp;	/* polar coordinate vel corr field  */
  my_boolean changed_reference_frame; /* flag */
  COMP_PRECISION omega_corr[3]; /* global correction rotation vector */

  int vsig_mode;		/* mode of slip uncertainty 
				   determination */
  long int seed;		/* seed for random gen. */
};
