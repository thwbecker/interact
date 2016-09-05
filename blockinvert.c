/*


invert velocities for block rotations in a plane (2-D) (can run in
spherical version with _sph)

reads in GPS velocities and tries to predict velocities within blocks
as indicated by block codes for each GPS observation


procedure with respect to the inversion for velocities follows closely
the approach of Meade et al. (BSSA, 92, 208, 2002)


remember to run block_checkflt to check the fault assignment and such


$Id: blockinvert.c,v 1.28 2004/10/05 01:09:46 becker Exp $


*/
#include "interact.h"
#include "blockinvert.h"
#include "numrec_svd_routines.h"

int main(int argc, char **argv)
{
  COMP_PRECISION ldepth,sv_max,chi2[3],minchi2[3],oldchi2[3],
    *vmat,*sval,velrms,stressrms=0.0,llim,*xtry,*xfix,dxtry,*xtryvonly,
    sscale_dfac,beta,beta0,global_stress_depth,*vslip,
    fsrms,fnrms,sv_cutoff,svc_test,dummy;
  int i,j,k,l,iter_count,nlmloop,nr_random_loop,nlmres,pdim,
    siter,nr_lin_stress_iter,remove_net_trans,nsv_zero,
    nsv_zero_test,mss;
  struct prj projection;
  /* 
     
  the following are logic switches with their default
  values

  */
  int constrain_euler = 0;	/* use contrained Euler poles? 
				   1: read from file
				   2: lock last, automatically set if remove_net_trans == -4
				*/
  my_boolean rigid = FALSE, /* include interseismic 
			    deformation */
    nofiles = FALSE, 
    use_lm = FALSE, /* use LM inversion? */
    input_control = FALSE,
    verbose = FALSE,
    test_sv_cutoff = FALSE, 
    invert_for_cfac = FALSE, 	/* invert for coupling factors? */
    print_sv = FALSE,		/* print singular values? */
    load_stress_depth = FALSE,
    invert_for_ld = FALSE,	/* invert for lockign depth? */
    mod_codes_out = FALSE,
    use_numrec_svd = FALSE,	/* using numerical recipes SVD will
				   allow to extract the V matrix, for
				   LAPACK, we still have to implement
				   this */
    damp_nslip = FALSE,		/* damp slip in normal motion direction? */
    calc_ginv = FALSE,		/* calculate the general inve */
    beta_init = FALSE,
    override_locking_depth = FALSE, /* override the locking depth settings */
    fgeo_out = FALSE,
    no_rel_stress_amp_scale = FALSE,/* by default, scale all
				       stress amplitudes
				       regionally to model
				       values. this implies that
				       stress input is given
				       normalized to unity. if set
				       to TRUE, will expect
				       regional variations in
				       stress (still trace=0) and
				       only adjust the overall RMS
				    */
    use_nullspace = FALSE,	/* use the null instead of the solution space */
    override_dip = FALSE,	/* override the dips of faults */
    constrain_slip_direction = FALSE, /* constrain the slip direction of certain 
					 segments */
    compute_vel_from_omega = FALSE; /* read in omega solution from
				       previous run, and compue the
				       velocities that correspond to
				       the rigid, long term motions
				    */
  FILE *out=NULL;
  struct bmd *mod;
  /*


  initialize model structure 

  */
  init_block_mods(&mod);
  // set pointers to zero so that we can use realloc safely
  xtry = xfix = xtryvonly = vmat = sval  = vslip = NULL;
  /*

  SET DEFAULT VALUES
  
  */
  //
  // define general projection for GPS observations and 
  // such
  projection.azi = 90.0;
  projection.type = OMERC_AZI;
  //projection.type = LCONFORM;
  //
  //
  ldepth = -1;			/* locking depth, by default (if
				   negative), uses the individual 
				   fault's locking_depth. if set,
				   will override */
  mod->nslip = 2;   /* nr of possible slip modes: 
		      1: only strike 2: strike & normal or dip
		   */
  llim = 1000;			/* maximum fault half length for
				   subdivisions (in km) */
  //
  //
  nr_random_loop = 1;			/* number of random loops */
  nlmloop= 100;		        /* max nr of LM iterations */

  beta0 = 1.0;			/* scaling factor for stress fit */
                                /* SVD cutoff */
#ifdef BLOCK_SPHERICAL
  sv_cutoff = 1e-4;		/* this leads to much smaller model norms */
#else
  sv_cutoff = 1e-8;		
#endif
  nr_lin_stress_iter = 30; /* max nr of linesr stress iterations, if
			      set to zero, will not iterate */
  sscale_dfac = 1.0;		/* stress re-scaling damping factor */

  mod->nslip_damp = 0.1;		/* damping value for normal slip
				   motion, if activated
				   should be between 0 (no damping) and 
				   ~1 (full damping)
				*/

  global_stress_depth = 5;		/* stress depth */
  mod->seed = -1;			/* random nr seed */
  /* 
     0 ... nrb-1 : remove rigid best-fit rotation of a single block
     -1: determine global best-fit rotation
     -2: determine and remove global, best-fit rotation
     -3: remove rigid best-fit rotation of the last block
     -4: remove best-fit rotation from last block depending on a few sites
  */
  remove_net_trans = -3;	

  mod->xdamp = 0.05;		/* some damping of solution vector */
  nsv_zero = 0;			/* 
				   nr of singular values set to zero,
				   overrides sv_cutoff
				*/
#ifdef DEBUG
  print_sv = TRUE;		/* output of singular values */
#endif

  mod->vsig_mode = 1;		/* mode of slip uncertainty
				   determination 
				   0: by combinations of
				      omega vectors +/- uncertainty
				   1: by randomizing omega
				*/
  nlmres = 5;
  /* 


  CHECK FOR COMMAND LINE ARGUMENTS 


  */
  for(i=1;i<argc;i++){	
    if(strcmp(argv[i],"-nslip")==0){ /* strike or strike + normal/dip slip */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&mod->nslip);
    }else if(strcmp(argv[i],"-rigid")==0) /* no fault deformation */
      toggle(&rigid);
    else if(strcmp(argv[i],"-nofiles")==0) /* no output files */
      toggle(&nofiles);
    else if(strcmp(argv[i],"-verbose")==0) /* somewhat verbose */
      toggle(&verbose);
    else if(strcmp(argv[i],"-vverbose")==0) /* use input control*/
      toggle(&input_control);
    else if(strcmp(argv[i],"-uselm")==0) /* use LM algorithm */
      toggle(&use_lm);
    else if(strcmp(argv[i],"-cvo")==0) /* compute velocities from omega  */
      toggle(&compute_vel_from_omega);
    else if(strcmp(argv[i],"-cginv")==0) /* calc general inverse */
      toggle(&calc_ginv);
    else if(strcmp(argv[i],"-norel")==0) /* no rel stress amp scale*/
      toggle(&no_rel_stress_amp_scale);
    else if(strcmp(argv[i],"-fgeo_out")==0) /* output of fault
					      geometries for
					      testing purposes */
      toggle(&fgeo_out);
    else if(strcmp(argv[i],"-ver")==0) /* all dip is prescribed as 90 */
      toggle(&override_dip);
    else if(strcmp(argv[i],"-crp")==0){ /* Euler poles of some blocks are constrained */
      constrain_euler=(constrain_euler)?(0):(1);
    }else if(strcmp(argv[i],"-clrp")==0){ /* lock last block */
      constrain_euler = 2;
    }else if(strcmp(argv[i],"-csd")==0) /* constrain slip direction */
      toggle(&constrain_slip_direction);
    else if(strcmp(argv[i],"-psv")==0) /* print singular values*/
      toggle(&print_sv);
    else if(strcmp(argv[i],"-rsd")==0) /* read stress depth from file */
      toggle(&load_stress_depth);
    else if(strcmp(argv[i],"-use_nrsvd")==0) /* use the numerical recipes SVD solver */
      toggle(&use_numrec_svd);
    else if(strcmp(argv[i],"-icf")==0) /* invert for slip locking factors */
      toggle(&invert_for_cfac);
    else if(strcmp(argv[i],"-ild")==0) /* invert for locking depths*/
      toggle(&invert_for_ld);
    else if(strcmp(argv[i],"-dns")==0) /* switch damping of slip motions on */
      toggle(&damp_nslip);
    else if(strcmp(argv[i],"-gcl")==0) /* output of llim modified
					  flt.block type file
					  code.llim.dat, then exit */
      toggle(&mod_codes_out);
    else if(strcmp(argv[i],"-uns")==0) /* use the null instead of solution space*/
      toggle(&use_nullspace);
    else if(strcmp(argv[i],"-ld")==0){	/* locking depth */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,&ldepth);
    }else if(strcmp(argv[i],"-beta")==0){	/* stress weighting */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,&beta0);
    }else if(strcmp(argv[i],"-sd")==0){	/* stress evaluation depth */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,&global_stress_depth);
    }else if(strcmp(argv[i],"-nrr")==0){	/* number of random loops */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&nr_random_loop);
    }else if(strcmp(argv[i],"-nsd")==0){	/* value of normal slip damping  */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%lf",&mod->nslip_damp);
    }else if(strcmp(argv[i],"-rnt")==0){	/* remove net translation of block value 
						   from velocity field. 
						   -3: remove that of the last block
						   -2: remove the global net rotation
						   -1: compute global net rotation, but no correction
						   0...nrb-1: remove the rotation of that block

						*/
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&remove_net_trans);
    }else if(strcmp(argv[i],"-nlmloop")==0){	/* max number of LM iterations */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&nlmloop);
    }else if(strcmp(argv[i],"-nsv_zero")==0){	/* nr of small SVs set to zero */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&nsv_zero);
    }else if(strcmp(argv[i],"-nlmres")==0){	/* number of LM restarts */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&nlmres);
    }else if(strcmp(argv[i],"-si")==0){	/* nr. of stress iterations for linear 
					   inversion */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&nr_lin_stress_iter);
    }else if(strcmp(argv[i],"-llim")==0){	/* max length */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,&llim);
    }else if(strcmp(argv[i],"-xdamp")==0){	/* solution damping */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,&mod->xdamp);
    }else if(strcmp(argv[i],"-seed")==0){	/* random seed */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&j);
      mod->seed = (long int)j;
    }else if(strcmp(argv[i],"-svc")==0){	/* SVD cutoff  */
      advance_argument(&i,argc,argv);
      sscanf(argv[i],ONE_CP_FORMAT,&sv_cutoff);
    }else{
      fprintf(stderr,"%s\n\n",argv[0]);
      fprintf(stderr,"options (current setting in parentheses):\n");

      fprintf(stderr,"-rigid\t\t\tdon't take elastic interseismic deformation into account\n\n");
      fprintf(stderr,"-cvo\t\t\tread omega solution from stdin and produce rigid block velocities\n");
      fprintf(stderr,"\t\t\t(See also block_compute_vel_from_omega)\n\n");
      
      fprintf(stderr,"-nslip\t\tvalue\tnr of slip directions: 1: strike only or 2: strike + normal/dip (%i)\n\n",
	      mod->nslip);
      fprintf(stderr,"-rnt\t\tvalue\tremove net translation of block value+1 (input from 0...nrb-1) (%i)\n",
	      remove_net_trans);
      fprintf(stderr,"\t\t\tif set to -1, will determine global net rotation\n");
      fprintf(stderr,"\t\t\tif set to -2, will remove    global net rotation\n");
      fprintf(stderr,"\t\t\tif set to -3, will remove net rotation of last block (like nrb-1)\n");
      fprintf(stderr,"\t\t\tif set to -4, will remove net rotation based on a few sites, and reference to last block\n");
      fprintf(stderr,"\t\t\t              in this case, hard-coded sites will be used or read from file %s\n",
	      RIGIDBLOCKSITES_FILE);
      fprintf(stderr,"\t\t\t              in lon lat format.\n");

      fprintf(stderr,"-crp\t\t\tconstrain the rotation poles of selected blocks.\n");
      fprintf(stderr,"\t\t\tif set, will read file %s in format: block_code(1...nrb) wx wy wz\n",EULER_POLE_FILE);
      fprintf(stderr,"\t\t\tand constrain these block motions. block_code has to be at end of list. (see -clrp) \n\n");
      fprintf(stderr,"-clrp\t\t\tconstrain the rotation pole of last block in list to be zero (see -crp)\n");

      fprintf(stderr,"-dns\t\t\tdamp normal slip motion (off by default, see -nsd)\n");
      fprintf(stderr,"-nsd\t\tvalue\tvalue for normal slip damping (0: none 1: full) (%g)\n",
	      mod->nslip_damp);
      fprintf(stderr,"-xdamp\t\tvalue\tdamp the norm of the Euler vector solution (%g)\n",
	      mod->xdamp);
      fprintf(stderr,"-csd\t\t\tconstrain the slip direction of faults (off by default)\n");
      fprintf(stderr,"\t\t\treads in constraints from %s, format: code strike_cons normal_cons\n",
	      CSD_FILE);
      fprintf(stderr,"\t\t\tcode is from 1 .. norig_flts, cons can be -1, 0, or 1\n\n");
      fprintf(stderr,"-ld\t\tvalue\toverride locking depths, should be >0 and in km (%g)\n",
	      ldepth);
      fprintf(stderr,"\t\t\tif set to negative value, will read locking depth for each fault\n");
      fprintf(stderr,"-ver\t\t\toverride dip setting and prescribe only vertical faults\n");
      fprintf(stderr,"-llim\t\tvalue\tmaximum flt half length in km (%g) for subdivisions\n\n",
	      llim);
      fprintf(stderr,"-beta\t\tvalue\tweighting for stress observations (%g)\n",
	      beta0);
      fprintf(stderr,"-sd\t\tvalue\tdepth to evaluate stresses at (%g)\n",
	      global_stress_depth);
      fprintf(stderr,"-rsd\t\t\tread stress evaluation depth for each stress observations from file\n");
      fprintf(stderr,"\t\t\tinstead of assigning constant of %g, if not changed with -sd\n",
	      global_stress_depth);
      fprintf(stderr,"-norel\t\t\tdo not scale the stress amplitudes to the model regionally, only global RMS\n\n");

      fprintf(stderr,"-nrr\t\tvalue\tnr of random optimization loops (%i)\n",
	      nr_random_loop);
      fprintf(stderr,"-si\t\tvalue\tmax nr of linear stress re-solving iterations (%i)\n",
	      nr_lin_stress_iter);
      fprintf(stderr,"-seed\t\tvalue\trandom number seed (<0) (%i)\n\n",
	      (int)mod->seed);

      fprintf(stderr,"-uselm\t\t\tuse LM optimization\n");
      fprintf(stderr,"-nlmloop\tvalue\tmaximum number of LM iteration (%i)\n",
	      nlmloop);
      fprintf(stderr,"-nlmres\t\tvalue\tnumber of random LM restarts (%i)\n\n",
	      nlmres);


      fprintf(stderr,"-svc\t\tvalue\tfraction of the largest singular value set to zero (%g)\n",
	      sv_cutoff);
      fprintf(stderr,"-nsv_zero\tvalue\tnr of small singular values set to zero, will override sv (%i)\n",
	      nsv_zero);
      fprintf(stderr,"\t\t\tif nsv_zero < 0, will test successive deletion of singular values\n");
      fprintf(stderr,"-use_nrsvd\t\tuse the numerical recipes implementation of SVD (for error bars)  (%i)\n",
	      use_numrec_svd);
      fprintf(stderr,"-uns\t\tuse the nullspace, not the solution space (for testing purposes, %i)\n",
	      use_nullspace);
      fprintf(stderr,"-psv\t\t\tprint singular values to file (0/1: %i)\n\n",
	      print_sv);

 
      fprintf(stderr,"-cginv\t\t\tcalculate and print the general inverse, resolution, and data matrices (0/1: %i)\n\n",
	      calc_ginv);
     
      fprintf(stderr,"-ild\t\t\tinvert for locking depths\n");
      fprintf(stderr,"-icf\t\t\tinvert for cfac slip (locking) factors\n\n");

      fprintf(stderr,"-nofiles\t\tno output files, no stress obs input\n");
      fprintf(stderr,"-fgeo_out\t\toutput of fault surface endpoints and center for testing\n");
      fprintf(stderr,"-gcl\t\t\toutput of %s type file with modified llim: code.llim.dat, then exit\n",
	      FLTBLOCK_FILE);
      fprintf(stderr,"-verbose\t\tproduce more output to stderr\n\n");
      fprintf(stderr,"-vverbose\t\talso echo some of the input velocities, stresses, and faults\n\n");
      exit(-1);
    }
  }
  // initialize random number generator
  if(mod->seed > 0)
    mod->seed = -mod->seed;
  myrand(&mod->seed);
  /* 

  do some input checks 
  
  */
  if(input_control)
    verbose = TRUE;
  if((mod->nslip < 1 ) || (mod->nslip > 2)){
    fprintf(stderr,"%s: error, nslip should be 1 (only strike) or 2 (strike+normal or dip depending on dip)\n",
	    argv[0]);
    exit(-1);
  }
  fprintf(stderr,"%s:\n%s: model setup\n%s:\n%s: slip modes: %i llim: %g %s beta: %g stress_depth: %g\n",
	  argv[0],argv[0],argv[0],argv[0],
	  mod->nslip,llim,(rigid)?("rigid calculation"):(""),
	  beta0,global_stress_depth);
  if(ldepth < 0.0){
     fprintf(stderr,"%s: using individual fault's locking depth as read in\n",
	     argv[0]);
     override_locking_depth = FALSE;
  }else{
    fprintf(stderr,"%s: WARNING: overriding  fault locking depth with global value of %g km\n",
	    argv[0],ldepth);
    override_locking_depth = TRUE;
  }
  if((mod->nslip_damp < 0)||(mod->xdamp < 0)){
    fprintf(stderr,"%s: error: norm or x damping negative: %g %g\n",
	    argv[0],mod->nslip_damp,mod->xdamp);
    exit(-1);
  }

  if(compute_vel_from_omega){
    /* 
       solution is read in from file, not computed 
    */
    rigid = TRUE;
    use_lm = FALSE;
    beta0 = 0.0;
    remove_net_trans = -1;
  }
  /* 
     SVD stuff
  */
  if(sv_cutoff < 0){
    fprintf(stderr,"%s: error: sv_cutoff < 0 (use nsv_zero < 0 for testing)\n",
	    argv[0]);
    exit(-1);
  }
  if(use_nullspace){
    fprintf(stderr,"%s: WARNING: using null instead of solution space\n",
	    argv[0]);
    if(!use_numrec_svd){
      fprintf(stderr,"%s: only implemented for numrec SVD\n",argv[0]);
      exit(-1);
    }
  }
  /* 
     check for nsv_zero settings 
  */
  if(nsv_zero > 0){		/* positive, leave some SVs out */
    fprintf(stderr,"%s: WARNING: setting the %i smallest singular values to zero\n",
	    argv[0],nsv_zero);
  }
  if(nsv_zero < 0){		/* negative, test SV deletion */
    test_sv_cutoff = TRUE;
    nsv_zero=0;
    if(!use_numrec_svd){
      fprintf(stderr,"%s: error: need to use Numerical Recipes SVD\n",
	      argv[0]);
      exit(-1);
    }
  }
  /* 
     rigid / beta settings
  */
  if(rigid && (fabs(beta0) > EPS_COMP_PREC)){
    fprintf(stderr,"%s: rigid calculation: suggest beta = 0 (beta: %g)\n",
	    argv[0],beta0);
    exit(-1);
  }
  beta = beta0;
  /* 
     
     constraining the slip directions requires LM inversion
  
  */
  if(constrain_slip_direction && (!use_lm)){
    fprintf(stderr,"%s: constrained slip directions are only implemented for LM inversion\n",
	    argv[0]);
    exit(-1);
  }
  /*

    read GPS velocities from in. vector allocation is done
    in this subroutine. also reads in block codes, fixed codes
    and 

    generates a general projection (different from the one used in the
    geographic adaptation of the Okada routine)

    might set constrain_euler, depending on the remove_net_trans setting

  */
  read_gps_velocities(mod,&projection,&velrms,argv,FALSE,
		      &remove_net_trans,&constrain_euler);
  if(!mod->nrgp){
    fprintf(stderr,"%s: error: no GPS data points read in, exiting\n",
	    argv[0]);
    exit(-1);
  }
  //
  //
#ifdef BLOCK_SPHERICAL
  /* 
     
  determine number of contraints from velocities
  
  */
  /* spherical */
  pdim = 3;
  mod->m1 =  mod->nrgp * pdim;	  /* nr of rows of A and E */
  mod->mgd = mod->nrgp * BLOCK_DIM; /* nr of entries in GPS data
				     vector */
#else  /* cartesian */
  pdim = BLOCK_DIM;
  mod->m1 = mod->nrgp * pdim;
  mod->mgd= mod->m1;
#endif
  //
  //
  /*
    
    read in stress observations and add them to v and sigv vectors
    the last number is the minimum stress sigma uncertainty
    the stress rms will be used for scaling beta later 

  */
  if(!mod_codes_out){
    /* if not only interested in the llim modified fault codes */
    read_stress_observations(mod,&stressrms,0.001,
			     no_rel_stress_amp_scale,
			     &mod->stress_depths,global_stress_depth,
			     load_stress_depth);
  }
  /* 
    
  determine number of contraints from stress inversion
  
  */
  if(rigid)
    mod->m2 = 0;			/* don't evaluate stresses if rigid
				   calculation */
  else{
    // nr of individual new observations
    mod->m2 = mod->nrsp * 6;
  }
  if(mod->m2){
    /* 
       save the original stress and stress uncertainty amplitudes 
       (will be modified by linear inversion) 
    */
    my_vecalloc(&mod->saved_stress,mod->m2*2,"saved_stress");
    a_equals_b_vector(mod->saved_stress,(mod->v+mod->mgd),mod->m2); /* stresses */
    a_equals_b_vector((mod->saved_stress+mod->m2),
		      (mod->sigv+mod->mgd),mod->m2); /* uncertainties */
  }
  //
  //
  if(input_control){
    for(i=j=0;i < mod->nrgp;i++,j+=BLOCK_DIM)
      fprintf(stderr,"vel: %5i: x: %11g %11g px: %11g %11g v: %11g %11g sigv: %11g %11g rho: %11g bc: %3i %s\n",
	      i+1,mod->gx[j+X], mod->gx[j+Y],mod->gpx[j+X], 
	      mod->gpx[j+Y],mod->v[j+X],mod->v[j+Y],
	      mod->sigv[j+X],mod->sigv[j+Y],
	      mod->rho[i],mod->bcode[j]+1,
	      (mod->block[mod->bcode[j]].fixed)?("fixed"):(""));
    for(i=0,j=mod->mgd;i < mod->nrsp;i++,j+=6)
      fprintf(stderr,"str: %5i: x: %11g %11g sxx: %g (%g) sxy: %g (%g) sxz: %g (%g) syy: %g (%g) syz: %g (%g) szz: %g (%g)\n",
	      i+1,*(mod->sx+i*BLOCK_DIM+X),*(mod->sx+i*BLOCK_DIM+Y),
	      mod->v[j],mod->sigv[j],mod->v[j+1],mod->sigv[j+1],mod->v[j+2],mod->sigv[j+2],mod->v[j+3],
	      mod->sigv[j+3],mod->v[j+4],mod->sigv[j+4],mod->v[j+5],mod->sigv[j+5]);
  }
  /*
    
    INCLUDE FAULTS
    
    read in faults and set up normal vectors as well as projected
    displacements at the observational velocity locations gx and
    stresses at the stress locations sx
    
  */
  read_bflt(mod,ldepth,projection,input_control,llim,
	    (mod_codes_out)?(TRUE):(rigid),TRUE,FALSE,override_locking_depth,
	    override_dip,fgeo_out,damp_nslip,mod_codes_out,
	    constrain_slip_direction);
  if(mod_codes_out){
    fprintf(stderr,"%s: written to code.%g.dat, now exiting\n",argv[0],llim);
    exit(0);
  }
  if(mod->nflt){
    //
    // nr of slip directions times number of faults
    //
    mod->nsnf = mod->nslip * mod->nflt;
  }else{
    fprintf(stderr,"%s: WARNING: %s not found or no faults\n",
	    argv[0],FLTBLOCK_FILE);
    mod->nflt = mod->nsnf = 0;
  }
  if(damp_nslip){
    if((mod->nfdamp != mod->nflt) && (mod->nfdamp != mod->nsnf)){
      fprintf(stderr,"%s: logic error: damping slip motion but nfdamp (%i) not %i (nflt) and not %i (nflt*nslip)\n",
	      argv[0],mod->nfdamp,mod->nflt,mod->nsnf);
      exit(-1);
    }
  }else{
    if(mod->nfdamp != 0){
      fprintf(stderr,"%s: logic error: no damping of slip motion but nfdamp %i\n",
	      argv[0],mod->nfdamp);
      exit(-1);
    }
  }
  /* 
     
  some block motions might be constrained, else this routine set nrbc
  to zero

  */
  read_constrained_euler_poles(mod,argv,constrain_euler);
  /*
  
    START MATRIX ASSEMBLY 

    dimensions of problem NBASE is the building block of the block
    motion parameter vector omega, vx0, and vy0

  */
  /*  
      number of free parameters (besides nc)
  */
  mod->n = BLOCK_NBASE * mod->nrb;/* free block motion parameters  */
  mod->na = mod->n;		/* 
				   general number of unknowns for
				   block motion paramters
				*/
  if(invert_for_ld)		/* 
				   increase by the locking depth
				   parameters for each fault, if we
				   invert for locking depth
				*/
    mod->na += mod->nflt;
  if(invert_for_cfac)
    mod->na += mod->nflt;		/* 
				   nflt coupling factors for slip
				   which are typically unity, add
				   those if we are inverting for them
				*/
  /* 
     number of constraints (data) (rows of K)
  */
  mod->m = mod->m1 + mod->m2; /* velocity components ((BLOCK_DIM or 3) *
			      nrgp) plus stresses (6 * nrsp)
			   */
  if(mod->m < mod->n){
    fprintf(stderr,"%s: error: %i blocks (n: %i) and only %i GPS and %i stress observations (m: %i)\n",
	    argv[0],mod->nrb,mod->n,mod->nrgp,mod->nrsp,mod->m);
    exit(-1);
  }
  mss = mod->ms = mod->m;
  if(damp_nslip && (mod->nfdamp==0)){
    fprintf(stderr,"%s: logic error: slip damping but nfdamp: %i\n",
	    argv[0],damp_nslip);
    exit(-1);
  }
  if(mod->nfdamp){	
    fprintf(stderr,"%s: damping slip motions with gamma: %g\n",
	    argv[0],mod->nslip_damp);
    /* 
       normal slip damping (this gets set in read_bflt)
       add number of flt slip directions to m
    */
    mod->ms += mod->nfdamp;
    mss = mod->ms;
  }
  if((mod->xdamp > 0)||(!compute_vel_from_omega)){
    /* 

       general norm of Euler solution damping toward rigid block
       motion from input

    */
    mod->nxdamp = mod->n;		/* elements of solution x */
    fprintf(stderr,"%s: damping Euler solution toward input rigid with alpha: %g\n",
	    argv[0],mod->xdamp);
    mod->ms += mod->nxdamp;	/* only one row necessary */
  }else{
    mod->nxdamp = 0;		/* this isn't really necessary */
  }
  /* 
     generate the zeroes for the right-hand side in case of
     normal/strike slip motion and/or norm solution damping
  */
#ifdef BLOCK_SPHERICAL
  my_vecrealloc(&mod->v,mod->ms,"v");
  my_vecrealloc(&mod->sigv,mod->ms,"v");
  my_vecrealloc(&mod->vc,mod->ms,"vc");
  my_vecrealloc(&mod->sigvc,mod->ms,"sigvc");
  /* 
     normal/strike slip damping 
  */
  for(i=mod->m,j=mod->mgd+mod->m2;i < mss;i++,j++){
    mod->vc[i]    = mod->v[j] = 0.0;	/* for damping */
    mod->sigvc[i] = mod->sigv[j] = 1.0;	/* fake uncertainty */
  }
  if(mod->nxdamp){
    if(((mod->ms - mss) != mod->nrb * BLOCK_NBASE) || (BLOCK_NBASE != 3)){
      fprintf(stderr,"%s: logic error A at nxdamp assignment: ms: %i mss: %i nrb: %i bnbase: %i\n",
	      argv[0],mod->ms,mss,mod->nrb,BLOCK_NBASE);
      exit(-1);
    }
    for(l=k=0;k < mod->nrb;k++)
      for(l=0;l < 3;l++){
	/* 
	   add the solution damping 
	*/
	mod->vc[i]    = mod->v[j] = mod->block[k].xrigid[l] * mod->xdamp;
	mod->sigvc[i] = mod->sigv[j] = 1.0;	/* fake uncertainty */
	i++;j++;
      }
    if((k*l != mod->nxdamp) || (i != mod->ms)){
      fprintf(stderr,"%s: logic error B at nxdamp assignment: k %i l: %i i: %i ms: %i nxdamp: %i\n",
	      argv[0],k,l,i,mod->ms,mod->nxdamp);
      exit(-1);
    }
  }
#else
  my_vecrealloc(&mod->v,mod->ms,"v");
  my_vecrealloc(&mod->sigv,mod->ms,"sigv");
  for(i=mod->m;i < mss;i++){
    mod->v[i]    = 0.0;	/* for damping */
    mod->sigv[i] = 1.0;	/* fake uncertainty */
  }	  
  for(k=0;i<mod->ms;i++,k++){
    mod->v[i] = mod->xrigid[k] * mod->xdamp;
    mod->sigv[i] = 1.0;
  }
#endif
  /* 

     nr of degrees of free parameters, now corrected by the possibly
     constrained blocks

  */
  mod->ncov = mod->na - mod->nc;
  if(use_lm && mod->nc){
    fprintf(stderr,"%s: error, LM scheme is not prepared for constrained blocks\n",
	    argv[0]);
    exit(-1);
  }
  //
  // allocate solution related arrays
  //
  // block motion parameters + possible cfac[nflt]
  my_vecrealloc(&mod->xsol,mod->na,"blockinvert: xsol"); /* this will
							  have been
							  already
							  allocated if
							  read_constrained_euler_poles
							  is called */
  my_vecrealloc(&xtry,mod->na,"blockinvert: xtry");
  my_vecrealloc(&xtryvonly,mod->na,"blockinvert: xtryvonly");
  my_vecrealloc(&xfix,mod->na,"blockinvert: xfix");
  if(mod->nc){			/* some Euler poles have been assigned
				   as a fixed solution */
    a_equals_b_vector((xtry+mod->first_c*BLOCK_NBASE),(mod->xsol+mod->first_c*BLOCK_NBASE),mod->nc);
    a_equals_b_vector((xtryvonly+mod->first_c*BLOCK_NBASE),(mod->xsol+mod->first_c*BLOCK_NBASE),mod->nc);
    a_equals_b_vector((xfix+mod->first_c*BLOCK_NBASE),(mod->xsol+mod->first_c*BLOCK_NBASE),mod->nc);
  }
  my_vecrealloc(&mod->sigma,mod->ncov,"blockinvert: sigma");
  my_vecrealloc(&mod->cov,mod->ncov*mod->ncov,"blockinvert: cov");
  //
  // observed/predicted velocities and stresses plus damping
  //
  my_vecrealloc(&mod->vmod,mod->ms,"blockinvert: vmod");
#ifdef BLOCK_SPHERICAL
  my_vecrealloc(&mod->vmodc,mod->ms,"blockinvert: vmodc");
  fprintf(stderr,"%s: %i (m: (vel: %i * 3 =) %i + (str: %i * 6 =) %i obs. + %i slip dp. + %i x dp.) by %i (n: %i * %i unknowns - con: %i) system\n",
	  argv[0],mod->ms,mod->nrgp,mod->m1,mod->nrsp,
	  mod->m2,mod->nfdamp * ((damp_nslip)?(1):(0)),
	  mod->nxdamp,mod->ncov,mod->nrb,BLOCK_NBASE,mod->nc);
#else
  fprintf(stderr,"%s: %i (m: (vel: %i * %i =) %i + (str: %i * 6 =) %i obs. + %i slip dp. + %i x. dp.) by %i (n: %i * %i unknowns) system\n",
	  argv[0],mod->ms,mod->nrgp,BLOCK_DIM,
	  mod->m1,mod->nrsp,mod->m2,mod->nfdamp * ((damp_nslip)?(1):(0)),
	  mod->nxdamp,mod->n,mod->nrb,BLOCK_NBASE);
#endif
  if(mod->nc)
    fprintf(stderr,"%s: WARNING: %i out of %i n-parameters are constrained, %i blocks\n",
	    argv[0],mod->nc,mod->na,mod->nrbc);
  if(compute_vel_from_omega){
    /* 
       read in previous solution from stdin
    */
    fprintf(stderr,"%s: reading omega solution in w_x s(w_x) w_y s(w_y) w_z s(w_z) dummy dummy dummy format in deg/Myr\n",
	   argv[0]);
    for(i=j=0;i < mod->nrb;i++,j+=BLOCK_NBASE){
      if(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		(mod->xsol+j+X),(mod->sigma+j+X),
		(mod->xsol+j+Y),(mod->sigma+j+Y),
		(mod->xsol+j+Z),(mod->sigma+j+Z),
		&dummy,&dummy,&dummy) != 9){
	fprintf(stderr,"%s: omega input read error\n",argv[0]);
	exit(-1);
      }
      for(k=0;k<BLOCK_NBASE;k++){
	/* scale to internal normalization */
	mod->xsol[j+k] *= BLOCK_GFAC;
      	mod->sigma[j+k] *= BLOCK_GFAC;
      }
    }
  }
  /*
    
  BEGIN SOLUTION PROCEDURE
  
  */
  minchi2[0] = minchi2[1] = minchi2[2] = FLT_MAX;
  iter_count = 0;
  do{
    if(!iter_count){  
      /*
	
      put together the A matrix that links block motion parameters
      with observed velocities at the sites. pass site coordinates in
      the globally projected frame
      
      assemble the rigid block motion matrix 
      A[nrgp*dim][nbase*nrb]
      
      only once, will not change during cfac randomization loop
      since A only depends on the block (not fault) geometry as it
      only links the translation and rotation of blocks to local
      velocity vectors
      
      dim is either BLOCK_DIM or 3
      
      */
      assemble_block_a(&mod->a,mod->gpx,mod->bcode,
		       mod->nrgp,mod->nrb,mod->block,-1,mod);
    }    
    /*
      
    if there exist faults, assemble D, G, and F matrices
    regardless of rigid or non-rigid calculation

    */
    if(mod->nflt){
      if(!iter_count){	 
	/* 
	   first pass  
	*/
	/*

	  initialize the locking depth/coupling factor vector
	  during the first call to the routine

	*/
	assign_additional_sol_values(xtry,mod->n,mod->na,
				     mod->nflt,mod->fault,
				     invert_for_ld,
				     invert_for_cfac,
				     &mod->seed,INIT_ADD_SOL);
	//
	// assemble the D matrix which relates fault local
	// to global slip 
	// dimensions: (dim*nrgp,nslip*nflt)
	// D depends on fault geometry and interaction coeff
	//
	// dim is BLOCK_DIM or 3 (for spherical)
	//
	assemble_block_d(&mod->d,mod->fault,mod->nflt,mod->nrgp,
			 mod->nslip,mod);
	//
	// assemble G, which links global-frame slip to 
	// projections on fault-strike,normal,and dip vectors
	// dimensions: (nslip*nflt, dim*nflt)
	// G depends on orientation but not locking depth
	//
	assemble_block_g(&mod->g,mod->fault,mod->nflt,mod->nslip);
	//
	// assemble I, which links fault locale slip to stresses
	// at the nrsp observational points
	// dimensions: (6*nrsp,nslip*nflt)
	// I depends on the fault geoemtry
	// 
	if(mod->nrsp)
	  assemble_block_i(&mod->imat,mod->fault,mod->nflt,
			   mod->nrsp,mod->nslip);
      }else{
	/* 
	   
	second pass

	*/
	if((!invert_for_cfac)&&(!invert_for_ld)){
	  fprintf(stderr,"%s: logic error: no ld or cfac inversion but in second random pass through loop\n",
		  argv[0]);
	  exit(-1);
	}
	/* 
	   randomize the locking depth and coupling factors for the
	   second pass through and later
	*/
	assign_additional_sol_values(xtry,mod->n,mod->na,mod->nflt,mod->fault,
				     invert_for_ld,invert_for_cfac,
				     &mod->seed,RANDOM_ADD_SOL);
	if(invert_for_ld){	/* 

				   if the locking depths change,
				   have to recalculate the fault 
				   geometry and the 
				   interaction coefficients 

				*/
	  change_locking_depths(&mod->fault,mod->nflt,
				mod->nrgp,mod->nrsp,mod->nslip,
				mod->gx,mod->sx,mod->stress_depths,
				rigid,projection,&mod->d,
				&mod->imat,(xtry+mod->n),mod);
	}
      }
      /*
	
	from now on, matrices might change since we could be fiddling
	around with the locking factors or with the locking depths,
	both of which affect the F matrix

      */
      assemble_block_fltdep_matrices(&mod->f,&mod->gf,&mod->e,
				     &mod->k2mat,mod->g,mod->d,
				     mod->imat,mod->m, mod->n, 
				     mod->nflt,mod->nrb,mod->nrgp,
				     mod->nrsp,mod->nsnf,xtry,
				     mod->fault,mod->block,
				     invert_for_ld,
				     invert_for_cfac);
#ifdef DEBUG
      if(!iter_count){
	fprintf(stderr,"%s: writing G.F to gf.dat\n",argv[0]);
	out=myopen("gf.dat","w");
	print_matrix_ftrn(mod->gf,mod->nsnf,mod->n,out,FALSE);
	fclose(out);
	fprintf(stderr,"%s: writing E = D.G.F to e.dat\n",
		argv[0]);
	out=myopen("e.dat","w");
	print_matrix_ftrn(mod->e,pdim*mod->nrgp,mod->n,out,FALSE);
	fclose(out);
	if(mod->nrsp){
	  fprintf(stderr,"%s: writing K2 = I.G.F to k2.dat\n",
		  argv[0]);
	  out=myopen("k2.dat","w");
	  print_matrix_ftrn(mod->k2mat,6*mod->nrsp,mod->n,out,FALSE);
	  fclose(out);
	}
      }
#endif
    } /* end of nflt > 0 part  */
    /*
      get the complete K matrix


          A - D G F 
	    - I G F
	gamma   G F, if strike/normal slip damping activated
        alpha unity, if normal solution vector damping is 
	             activated
	
    */
    assemble_block_k(&mod->kmat,mod->ms,mod->m1,mod->m2,
		     mod->n,mod->nrgp,mod->nrsp,mod->nflt,
		     mod->a,mod->e,mod->k2mat,rigid,damp_nslip,
		     mod->gf,mod->fault,mod->nsnf,mod->nslip,
		     mod->nfdamp,mod->nxdamp,mod->xdamp);
    if(0)
      fprintf(stderr,"1: |a|: %g |e|: %g |gf|: %g |d|: %g |k2|: %g |k|: %g\n",
	      norm(mod->a,BLOCK_DIM*mod->nrgp*BLOCK_NBASE*mod->nrb),
	      norm(mod->e,pdim*mod->nrgp*BLOCK_NBASE*mod->nrb),
	      norm(mod->gf,mod->n*mod->nslip*mod->nflt), 
	      norm(mod->d,pdim * mod->nrgp *mod->nslip* mod->nflt),
	      norm(mod->k2mat,mod->nrsp*6*BLOCK_NBASE*mod->nrb),
	      norm(mod->kmat,mod->n * mod->ms));
#ifdef DEBUG
    if(!iter_count){
      fprintf(stderr,"%s: writing K matrix %s to k.dat %s\n",
	      argv[0],(rigid)?("(equal to A)"):("(K=A-E,-K2)"),
	      (damp_nslip?("(g G.F)"):("")));
      out=myopen("k.dat","w");
      print_matrix_ftrn(mod->kmat,mod->ms,mod->n,out,FALSE);
      fclose(out);
    }
#endif
    /* 

    END OF MATRIX ASSEMBLY PART, NOW SOLVE

    
    */
    if(!compute_vel_from_omega){
      /* 

      
      real solution finding branch of program
      

      */
      if(test_sv_cutoff){
	out = myopen("damp.dat","w");
	fprintf(stderr,"%s: printing chi2 as a function of sv cutoff to damp.dat\n",
		argv[0]);	
	/* output in format 
	   SV_cutoff_fraction chi2_GPS chi2_stress \
	   ... norm(x) strike_rms normal_rms nr_zero_SV
	   
	*/
	fprintf(out,"# svc vchi2 schi2 chi2 norm(xsol) fsrms rnrms nr_zero_SV norm(xsigma)\n");
	my_vecrealloc(&vslip,mod->nsnf,"vslip");
      }
      svc_test = sv_cutoff;
      for(nsv_zero_test = nsv_zero;nsv_zero_test < mod->n;
	  nsv_zero_test += (test_sv_cutoff)?(1):(mod->n+1)){
	if(!beta_init){  
	  /*
	    normalize beta by BLOCK_DIM/6 such that each observation is
	    weighted the same, regardless of the parameters it has (two
	    velocities vs. six stress components)
	  */
	  beta = beta0*(((COMP_PRECISION)BLOCK_DIM)/6.0);
	  beta_init = TRUE;
	}
	/*
	  
	solve the system for velocities only 
	
	the solution will be in xtry
	
	if nrbc != 0, only the first n-nrbc entries of xtry will
	be overwritten

	*/
#ifdef BLOCK_SPHERICAL
	solve_block(mod->kmat,xtry,mod->vc,mod->m1,mod->m2,mod->m,
		    mod->ms,FALSE,mod->n,mod->nsnf,mod->nfdamp,
		    mod->nxdamp,mod->sigvc,beta,&sv_max,&sval,&vmat,
		    mod->vmodc,&svc_test,use_numrec_svd,
		    nsv_zero_test,calc_ginv,mod->nc,use_nullspace);
	
#else
	solve_block(mod->kmat,xtry,mod->v,mod->m1,mod->m2,mod->m,mod->ms,
		    FALSE,mod->n,mod->nsnf,mod->nfdamp,
		    mod->nxdamp,mod->sigv,beta,&sv_max,&sval,&vmat,
		    mod->vmod,&svc_test,use_numrec_svd,
		    nsv_zero_test,calc_ginv,mod->nc,use_nullspace);
#endif
	/* 
	   save velocity only solution 
	*/
	a_equals_b_vector(xtryvonly,xtry,mod->na);
	if(print_sv){
	  /* output of singular value */
	  if((!test_sv_cutoff) && (!iter_count)){
	    fprintf(stderr,"%s: printing solution for sv_cutoff: %g to %s (%i dim)\n",
		    argv[0],svc_test,SV_FILE,mod->n-mod->nc);
	    out = myopen(SV_FILE,"w");
	    print_singular_values(sval,(mod->n-mod->nc),out);
	    fclose(out);
	  }
	}
	/*
	  
	calculate the predicted observables (velocities+stresses)
	K . xtry = vmod
	
	*/
	evaluate_block_solution(mod->kmat,mod->ms,mod->n,xtry,
				mod->vmod,mod->vmodc,mod);
	/* 
	   scale stresses to model predictions 
	*/
	if(mod->m2)
	  rescale_observed_stresses(mod->vmod,mod->v,mod->sigv,
				    sscale_dfac,&stressrms,
				    no_rel_stress_amp_scale,mod,
				    TRUE,TRUE);
	/*
	  calculate chi^2 misfit, weighted sum over 
	  squared diff between v_obs and v_mod
	  
	  (do not include the possible last nsnf rows which 
	  are used for damping)
	  
	*/
	chi2[0] = block_chi_square(mod->vmod,mod->v,mod->sigv,mod->mgd,
				   mod->m2,beta,(chi2+1),(chi2+2));
	fprintf(stderr,"%s: ri: %4i vel     sol:     chi2: %12g vel_c2: %12g str_c2: %12g str_rms: %12g (beta: %11g) |x|: %g\n",
		argv[0],iter_count+1,chi2[0],chi2[1],chi2[2],
		stressrms,beta,norm(xtry,mod->na));
	/* 
	   now we have a solution in xtry and stress observables scaled
	   to that solution
	*/
	if((!rigid) && (nr_lin_stress_iter) && (mod->nrsp) 
	   && (beta != 0.0)){
	  siter = 0;
	  oldchi2[0] = oldchi2[1] = oldchi2[2] = FLT_MAX;
	  /* 
	     
	  START ITERATION FOR STRESS fit
	  
	  */
	  do{
	    if(!beta_init){
	      fprintf(stderr,"%s: error, beta not initialized\n",
		      argv[0]);
	      exit(-1);
	    }
	    /*
	      recompute solution xtry with stresses included in
	      inversion
	    */
#ifdef BLOCK_SPHERICAL
	    solve_block(mod->kmat,xtry,mod->vc,mod->m1,mod->m2,
			mod->m,mod->ms,TRUE,
			mod->n,mod->nsnf,mod->nfdamp,mod->nxdamp,
			mod->sigvc,beta,
			&sv_max,&sval,&vmat,mod->vmodc,&svc_test,
			use_numrec_svd,nsv_zero_test,calc_ginv,
			mod->nc,use_nullspace);
#else
	    solve_block(mod->kmat,xtry,mod->v,mod->m1,mod->m2,
			mod->m,mod->ms,TRUE,mod->n,mod->nsnf,
			mod->nfdamp,mod->nxdamp,mod->sigv,beta,
			&sv_max,&sval,&vmat,mod->vmod,&svc_test,
			use_numrec_svd,nsv_zero_test,calc_ginv,
			mod->nc,use_nullspace);
#endif
	    evaluate_block_solution(mod->kmat,mod->ms,mod->n,xtry,
				    mod->vmod,mod->vmodc,mod);
	    /* 
	       rescale stresses
	    */
	    rescale_observed_stresses(mod->vmod,mod->v,mod->sigv,
				      sscale_dfac,&stressrms,
				      no_rel_stress_amp_scale,
				      mod,FALSE,TRUE);
	    /* calculate new chi2 */
	    chi2[0] = block_chi_square(mod->vmod,mod->v,mod->sigv,
				       mod->mgd,mod->m2,beta,
				       (chi2+1),(chi2+2));
	    if(siter == 0)
	      a_equals_b_vector(oldchi2,chi2,3);
	    siter++;
	    fprintf(stderr,"%s: ri: %4i vel+str sol: %2i: chi2: %12g vel_c2: %12g str_c2: %12g str_rms: %12g (beta: %11g) |x|: %g\n",
		    argv[0],iter_count+1,siter,chi2[0],chi2[1],
		    chi2[2],stressrms,beta,norm(xtry,mod->na));
	    if(chi2[0] <= oldchi2[0]){	/* save good solutions in xfix and
					   chi2 ins oldchi2 */
	      a_equals_b_vector(xfix,xtry,mod->na);
	      a_equals_b_vector(oldchi2,chi2,3);
	    }
	  }while((beta > 0.0)&&(chi2[0] <= oldchi2[0])&&
		 (siter < nr_lin_stress_iter));
	  if(siter == nr_lin_stress_iter)
	    fprintf(stderr,"%s: WARNING: max nr of lin. stress iterations reached\n",
		    argv[0]);
	  /* copy best solution to xtry */
	  a_equals_b_vector(xtry,xfix,mod->na);
	  a_equals_b_vector(chi2,oldchi2,3);
	}
	/* 
	   
	now we have the best solution given velocities and stresses after
	the stress iteration of this (random) step in xtry, and the best
	chi2s in chi2[]
	
	*/
	if(test_sv_cutoff) {
	  /* 
	     test SVD cutoff and write 
	     
	     cutoff vchi2 schi2 chi2 norm(x) fsrms fnrms nr_sv_zero norm(sigma)
	     
	     to damp.dat 
	  */
	  /* get fault slip values */
	  calc_Ax_ftn(mod->gf,mod->nsnf,mod->n,xtry,vslip);
	  calc_fault_sn_rms(vslip,mod->nflt,mod->nslip,&fsrms,&fnrms);
	  // calculate cov. matrix
	  assemble_cov_svd(vmat,mod->ncov,sval,mod->cov);
	  // get uncertainties
	  for(i=j=0;i < mod->ncov;i++,j += mod->ncov)
	    mod->sigma[i] = sqrt(mod->cov[j+i]);
	  /* 
	     write output
	     
	     svc vchi2 schi2 chi2 norm(xsol) fsrms fnrms nr_zero_SV norm(sigma)
	     
	  */
	  fprintf(out,"%22.10e %22.10e %22.10e %22.10e %22.10e %22.10e %22.10e %5i %22.10e\n",
		  svc_test,chi2[1],chi2[2],chi2[0],norm(xtry,mod->n),
		  fsrms,fnrms,countzero_vec(sval,mod->n),
		  norm(mod->sigma,mod->n));
	}
      }
      if(test_sv_cutoff)	/* close output of log file */
	fclose(out);
      /* 
	 
      was there a decrease in chi^2 for this (random) iteration?
      
      */
      if((chi2[0] < minchi2[0])||(!finite(chi2[0]))){		
	// save minimum chi2
	a_equals_b_vector(minchi2,chi2,3);		
	//
	// save solution
	a_equals_b_vector(mod->xsol,xtry,mod->na);
	/*
	  calculate covariance from solution of SVD 
	*/
	assemble_cov_svd(vmat,mod->ncov,sval,mod->cov);
	fprintf(stderr,"%s: lin: i:%5i improved chi2: %12g vc2: %12g sc2: %12g |x|: %12g ",
		argv[0],iter_count+1,minchi2[0],
		minchi2[1],minchi2[2],
		norm(mod->xsol,mod->na));
	i=mod->n;
	if(invert_for_ld){
	  fprintf(stderr,"<ld>: %g ",
		  mean((mod->xsol+i),1,mod->nflt));
	  i += mod->nflt;
	}
	if(invert_for_cfac)
	  fprintf(stderr,"<cf>: %g",mean((mod->xsol+i),1,mod->nflt));
	fprintf(stderr,"\n");
      }
      iter_count++;
      if((nr_random_loop > 1)&&(!invert_for_cfac)&&
	 (!invert_for_ld)){
	/* 
	   if we are not trying to invert for locking depths of slip
	   factors, we will explore the solution space by random
	   variations around the original solution to the block vector
	   without iterating the stresses. this is meant to make sure
	   that we are not drawn to a local minimum by the stress iteration
	   
	*/
	fprintf(stderr,"%s: WARNING: more than one random iteration, but no ld or cfac inversion\n",
		argv[0]);
	fprintf(stderr,"%s: commencing Monte Carlo exploration of block motion vector space: %i trials\n",
		argv[0],nr_random_loop);
	fprintf(stderr,"%s: chi2 values in rc2.log. If no further output is found to stderr, no improved solutions\n",
		argv[0]);
	if((!mod->m2) || (beta == 0)){
	  fprintf(stderr,"%s: this iteration only makes sense (?) for beta != 0 and stresses\n",
		  argv[0]);
	  exit(-1);
	}
	/* save original solution and set range of variation */
	//a_equals_b_vector(xfix,xtryvonly,mod->na);
	a_equals_b_vector(xfix,xtry,mod->na);
	a_equals_b_vector(oldchi2,chi2,3);
	dxtry = norm(xfix,mod->na)/(COMP_PRECISION)mod->na * 0.001;
	out=myopen("rc2.log","w");
	/* 
	   set up local loop 
	*/
	while(iter_count < nr_random_loop){
	  /* add a random variation to the solution */	
	  a_equals_b_vector(xtry,xfix,mod->na);
	  for(i=0;i < mod->na;i++)
	    xtry[i] += mygauss_randnr(dxtry,&mod->seed);
	  //
	  // get new solution
	  evaluate_block_solution(mod->kmat,mod->ms,mod->n,xtry,
				  mod->vmod,mod->vmodc,mod);
	  /* rescale stresses */
	  rescale_observed_stresses(mod->vmod,mod->v,mod->sigv,
				    sscale_dfac,&stressrms,
				    no_rel_stress_amp_scale,mod,
				    TRUE,TRUE);
	  /* evaluate new misfit */
	  chi2[0] = block_chi_square(mod->vmod,mod->v,mod->sigv,
				     mod->mgd,mod->m2,beta,
				     (chi2+1),(chi2+2));
	  fprintf(out,"%20.10e %20.10e\n",chi2[0],norm(xtry,mod->na));
	  if(chi2[0] < minchi2[0]){ /* save improved solutions */
	    fprintf(stderr,"%s: ri: %4i vel     sol:     chi2: %12g vel_c2: %12g str_c2: %12g str_rms: %12g (beta: %11g) |x|: %g\n",
		    argv[0],iter_count+1,chi2[0],chi2[1],chi2[2],stressrms,beta,norm(xtry,mod->na));
	    a_equals_b_vector(minchi2,chi2,3);		
	    a_equals_b_vector(mod->xsol,xtry,mod->na);
	  }
	  iter_count++;
	}
	fclose(out);
	nr_random_loop = 1;	/* this ensures exit
				   two lines down*/
      }
      /* 
	end of solution finding branch 
      */
    }else{
      /* 
	 branch for omega solution is simply read in 
      */
      iter_count++;
    }
    /* 
       end of possible random evalutaion do-loop
    */
  }while(iter_count < nr_random_loop);		
  /* 
     call the LM routine, if needed. this will return new best fit
     solution as mod->xsol, chi2, and the covariance cov
  */
  if(use_lm){
    /* 
       
      run a Levenberg Marquardt iteraton scheme to find a solution
      this gives (nicely enough) the linear solution back for beta = 0
      and no other special iteration
      
    */
    /* 
       nr of random LM restarts 
    */
    if(beta == 0)		/* one should be enough, if we are not
				   inverting for locking depths or
				   cfactors
				*/
	nlmres = 5;
    else			/* we have to figure out if this is OK
				   or not */
      nlmres = 15;			
    /* 
       for debugging purposes
    */
    fprintf(stderr,"%s: |v|: %11g |vmod|: %11g |vsig|: %11g |s|: %11g |x|: %11g |un|: %11g\n",
	    argv[0],norm(mod->v,mod->mgd),norm(mod->vmod,mod->mgd),
	    norm(mod->sigv,mod->mgd),norm((mod->vmod+mod->mgd),mod->m2),
	    norm(mod->xsol,mod->n),
	    (mod->nfdamp)?(norm((mod->vmod+mod->mgd+mod->m2),mod->nfdamp)):(0.0));
    
    run_lm(mod,&mod->seed,projection,chi2,rigid,beta_init,beta,
	   &stressrms,no_rel_stress_amp_scale,damp_nslip,
	   sscale_dfac,invert_for_ld,invert_for_cfac,
	   constrain_slip_direction,nlmloop,nlmres,argv);
  }else{			/* no LM method */
    mod->ncov = mod->n - mod->nc;
  }
  /*
    
  the best solution is now in xsol, it's covariance in cov
  
  recompute velocity and stress model values
  
  */
  if(invert_for_ld)
    block_assemble_fit_vector(&mod->kmat,mod->m,mod->ms,mod->m1,mod->m2,
			      mod->n,mod->nrgp,mod->nrsp,mod->nflt,
			      mod->nslip,mod->xsol,-1,mod->vmod,mod->vmodc,
			      TRUE,mod->a,&mod->d,mod->g,&mod->imat,mod->nrb,
			      mod->block,mod->nsnf,&mod->fault,rigid,
			      &mod->gf,invert_for_ld,
			      invert_for_cfac,TRUE,mod->gx,mod->sx,
			      mod->stress_depths,projection,
			      damp_nslip,mod->nfdamp,mod->nxdamp,
			      mod->xdamp,mod);
  else{
    evaluate_block_solution(mod->kmat,mod->ms,mod->n,mod->xsol,
			    mod->vmod,mod->vmodc,mod);
  }
  if((!rigid) && (mod->nrsp))   {
    /* scale stresses model */
    rescale_observed_stresses(mod->vmod,mod->v,mod->sigv,
			      sscale_dfac,&stressrms,
			      no_rel_stress_amp_scale,mod,TRUE,
			      TRUE);
  }
  chi2[0] = block_chi_square(mod->vmod,mod->v,mod->sigv,mod->mgd,
			     mod->m2,beta,(chi2+1),(chi2+2));
  fprintf(stderr,"%s: total:              chi2: %12g vc2: %12g sc2: %12g |x|: %12g |y|: %12g\n",
	  argv[0],chi2[0],chi2[1],chi2[2],
	  norm(mod->xsol,mod->na),norm(mod->vmod,(mod->mgd+mod->m2)));
  if(!compute_vel_from_omega){
    /*
      sigma uncertainties in model parameters are sqrt(diagonal(cov))
    */
    my_vecrealloc(&mod->sigma,mod->ncov,"sigma");
    for(i=j=0;i < mod->ncov;i++,j += mod->ncov){
      mod->sigma[i] = sqrt(mod->cov[j+i]);
    }
  }
  
  //
  // call output routine
  //
  block_output(mod,rigid,argv,projection,nofiles,
	       ((use_lm)||(iter_count != 1))?(TRUE):(FALSE),
	       beta,invert_for_ld,damp_nslip,invert_for_cfac,
	       no_rel_stress_amp_scale,use_nullspace,verbose);
  //
  // free arrays
  //
  // vectors
  free(mod->gx);free(mod->v);free(mod->xsol);free(mod->sigv);
  free(mod->sigma);free(sval);free(xfix);free(mod->vc);
  free(mod->sigvc);free(mod->pbase);free(mod->gcx);
  free(mod->gpx);free(mod->bcode);free(mod->rho);free(xtry);
  free(mod->vmod);free(mod->vmodc);
  // matrices
  free(mod->kmat);free(mod->a);free(mod->cov);
  free(mod->imat);free(mod->k2mat);
  if(mod->nflt){// if faults were read in 
    // matrices
    free(mod->gf);free(mod->d);free(mod->g);free(mod->f);free(mod->e);
    // fault structure
    free_bflt(&mod->fault,mod->nflt,mod->nrsp);
  } 
  if(mod->nrsp){// if stresses were read in
    free(mod->sx);
    if(mod->m2)			/* actually used */
      free(mod->saved_stress);
  }
  if(test_sv_cutoff){
    free(vslip);
  }
  if(mod->nrb)
    free(mod->block);
  GMT_end (1, argv);
  fprintf(stderr,"%s: done\n",argv[0]);
  return 0;
}

