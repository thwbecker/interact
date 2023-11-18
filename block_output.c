/*


output part of the blockinvert program

normally assumes that gf, k, and vmod (the solution vector, 
K . x) have been computed already. 

if new_kgfvmod is true, will recompute those

inportant input: xsol[n] the solution, and sigma[n] the 
uncertainties in the parameters


if use_nullspace is set, will print out velocities without
the shift in reference frame


$Id: block_output.c,v 1.20 2011/01/07 07:19:58 becker Exp $


*/
#include "interact.h"
#include "blockinvert.h"

void block_output(struct bmd *mod,my_boolean rigid,
		  char **argv,struct prj projection,
		  my_boolean nofiles,my_boolean new_kgfvmod,
		  COMP_PRECISION beta,my_boolean invert_for_ld,
		  my_boolean damp_nslip,my_boolean invert_for_cfac,
		  my_boolean no_rel_stress_amp_scale,
		  my_boolean use_nullspace,my_boolean verbose)
{
  FILE *out,*out2=NULL,*out3;
  int i,j,k,os1,idir,nrbf;
  int reference_block; /* 
			  if set to -1, will not pick a reference 
			  block
		       */
  COMP_PRECISION chi2,vnorm2,*vslip,*vsigma,c2a,c2b,tchi2[3],
    bchi2,bchi2sum=0.0,*disp,sign,*cfac,frms,fsrms,fnrms,*omega,
    srms,*sconst,smod[9],mean_vmisfit,gamma_rms,omega_ref[3],
    modval[3],vec[9],mod_shear,lon,lat,mag,vshift[3];
  my_boolean calc_block_chi2,calc_best_fit_euler,
    domega_out,vref_out;
#ifdef BLOCK_SPHERICAL
  COMP_PRECISION omega_block[3],wm[2],t[2],r[2];
#endif
  /* 
     end declarations 
  */
#ifdef BLOCK_SPHERICAL
  calc_best_fit_euler=FALSE;	/* for the spherical version of the code,
				   this will be identical to the x solution
				   only for a rigid computation
				*/
#else
  calc_best_fit_euler=TRUE;
#endif
  omega = NULL;
  domega_out = TRUE;		/* output of relative rotation vector
				   matrix  */
  vref_out = TRUE;		/* output of block velocities with
				   respect to reference block */

  /* 
     reference block is last block
  */
  reference_block = mod->nrb - 1;
  calc_block_chi2 = TRUE;
  fprintf(stderr,"%s: output of solution:\n",argv[0]);
  vsigma = vslip = disp = cfac = NULL;	/* init  */

  /* 
     number free blocks, ie. total - constrained blocks 
  */
  nrbf = mod->nrb - mod->nrbc;	
  if(new_kgfvmod){
    //
    // calculate K, gf (possibly anew) and vmod solution velocities
    // and stresses
    //
    block_assemble_fit_vector(&mod->kmat,mod->m,mod->ms,
			      mod->m1,mod->m2,
			      mod->n,mod->nrgp,mod->nrsp,
			      mod->nflt,mod->nslip,mod->xsol,-1,
			      mod->vmod,mod->vmodc,
			      TRUE,mod->a,&mod->d,
			      mod->g,&mod->imat,mod->nrb,
			      mod->block,mod->nsnf,&mod->fault,
			      rigid,&mod->gf,invert_for_ld,
			      invert_for_cfac,
			      invert_for_ld,mod->gx,mod->sx,
			      mod->stress_depths,
			      projection,damp_nslip,mod->nfdamp,
			      mod->nxdamp,mod->xdamp,mod);
#ifdef DEBUG
    fprintf(stderr,"block_output: writing K matrix to k.dat\n");
    out=myopen("k.dat","w");
    print_matrix_ftrn(mod->kmat,mod->ms,mod->n,out,FALSE);
    fclose(out);
    fprintf(stderr,"block_output: writing G.F to gf.dat\n");
    out=myopen("gf.dat","w");
    print_matrix_ftrn(mod->gf,mod->nsnf,mod->n,out,FALSE);
    fclose(out);
#endif
    if(mod->m2)
      rescale_observed_stresses(mod->vmod,mod->v,mod->sigv,
				1.0,&srms,
				no_rel_stress_amp_scale,
				mod,TRUE,TRUE);
    fprintf(stderr,"block_output: reevaluated K (|K|: %g) and x (|x|: %g)\n",
	    norm(mod->kmat,mod->n*mod->ms),
	    norm(mod->xsol,mod->na));
  }
  srms = rms((mod->v+mod->mgd),mod->m2);
  //
  // calculate chi2 of velocities per block
  //
  if(calc_block_chi2){
    fprintf(stderr,"%s: block sorted GPS chi2: ",argv[0]);
    bchi2sum = 0.0;
    for(i=0;i < mod->nrb;i++){
      bchi2 = 0.0;
      for(j=0;j < mod->mgd;j++)
	if(mod->bcode[j] == i)
	  bchi2 += square((mod->v[j] - mod->vmod[j])/mod->sigv[j]);
      fprintf(stderr,"%c: %12g ",bname(i),bchi2);
      bchi2sum += bchi2;
    }
    fprintf(stderr,"\n");
  }
  // the total chi^2 misfit \sum_i ((x_i - y_i)/sig_i)^2
  chi2 = block_chi_square(mod->vmod,mod->v,mod->sigv,
			  mod->mgd,mod->m2,beta,&c2a,&c2b);

  /* 
     test constant stress as chi2 
  */
  my_vecalloc(&sconst,mod->m,"sconst");
  for(i=0;i<mod->mgd;i++)
    sconst[i]=0.0;
  if(mod->m2){
    for(j=mod->mgd,i=0;i < mod->nrsp;i++,j+=6){
      /* 
	 
      this corresponds to a constant compressive stress with 
      azimuth of 0.08 degrees and e1-e2 = 1 - (-1) = 2
      
      */
      sconst[j+0]= 0.999984; sconst[j+1]=-0.00558502;sconst[j+2]=0.0;
      sconst[j+3]= -sconst[j+0];sconst[j+4]=0.0;sconst[j+5]=0.0;
      // determine stress observation shear 
      a_equals_b_vector(smod,(mod->vmod+j),6);
      expand_stress_matrix6to9(smod); 
      calc_eigensystem_sym3d(smod,modval,vec,FALSE);
      /* scaling factor is model shear stress/2 */
      mod_shear = (modval[EIGEN_E1] - modval[EIGEN_E3])/2.0; 
      for(k=0;k<6;k++)
	sconst[j+k] *= mod_shear;
    }
    tchi2[0] = block_chi_square(mod->v,sconst,mod->sigv,mod->mgd,
				mod->m2,beta,(tchi2+1),(tchi2+2));
    fprintf(stderr,"%s: stress chi2 with constant stress: %g\n",
	    argv[0],tchi2[2]);
  
    free(sconst);
  }
  if((calc_block_chi2) && (fabs(c2a - bchi2sum) > 1e-7)){
    fprintf(stderr,"%s: error, block GPS chi2 doesn't match up: %g vs. %g\n",
	    argv[0],c2a,bchi2sum);
  }
  /* 
     calculate weighted square norm of orginal v_obs for variance
     reduction
  */
  vnorm2 = block_data_norm(mod->v,mod->sigv,mod->mgd,mod->m2,beta,
			   mod);
  //
  // output of solution vector to stderr
  //
  // also, add the correction back in, if necessary
  //
  fprintf(stderr,"%s: writing w_x s(w_x) w_y s(w_y) w_z s(w_z) w^r_x w^r_y w^r_z to %s\n",
	  argv[0],OMEGA_OUT);
  out = myopen(OMEGA_OUT,"w");
  if(domega_out){		/* difference between Euler vectors
				   for debugging */
    out2 = myopen("domega.dat","w");
    fprintf(stderr,"%s: writing relative Euler vectors to domega.dat\n",
	    argv[0]);
  }
  fprintf(stderr,"%s: output of Euler vectors with uncertainties in original reference frame:\n",
	  argv[0]);
  if(norm_3d(mod->omega_corr) > EPS_COMP_PREC)
    fprintf(stderr,"%s: omega output will include the previous correction of %g, %g, %g\n",
	    argv[0],mod->omega_corr[INT_X]/BLOCK_GFAC, 
	    mod->omega_corr[INT_Y]/BLOCK_GFAC,mod->omega_corr[INT_Z]/BLOCK_GFAC);
  for(i=0;i < mod->nrb;i++){
    /* 



    THIS IS WHERE WE ADD THE CORRECTION AGAIN to the solution vector
    
    DON'T ELIMINATE THIS LOOP

    */
    add_b_to_a_vector_3d((mod->xsol+i*BLOCK_NBASE),mod->omega_corr);
    /* 
       ADD THE CORRECTION TO THE BEST FIT RIGID MOTION  
    */
    add_b_to_a_vector_3d(mod->block[i].xrigid,mod->omega_corr);
#ifndef BLOCK_SPHERICAL
    /* 
       flat 
    */
    fprintf(stderr,"%s: block %c (%2i): w: %12g (%12g) v_0: %12g (%12g) %12g (%12g) %s\n",
	    argv[0],bname(i),i+1,
	    mod->xsol[i*BLOCK_NBASE],  (i<nrbf)?(mod->sigma[i*BLOCK_NBASE]):(NAN),
	    mod->xsol[i*BLOCK_NBASE+1],(i<nrbf)?(mod->sigma[i*BLOCK_NBASE+1]):(NAN),
	    mod->xsol[i*BLOCK_NBASE+2],(i<nrbf)?(mod->sigma[i*BLOCK_NBASE+2]):(NAN),
	    (mod->changed_reference_frame)?("again in orig RF"):("in orig RF"));
#else
    /* 

    output of solution parameters to stderr, possibly reshifted into old frame
    
    */
    if(verbose)
      fprintf(stderr,"%s: block %c (%2i): wx: %12g (%12g) wy: %12g (%12g) wz: %12g (%12g) (deg/Myr) %s\n",
	      argv[0],bname(i),i+1,
	      mod->xsol[i*BLOCK_NBASE+INT_X]/BLOCK_GFAC,(i<nrbf)?(mod->sigma[i*BLOCK_NBASE+INT_X]/BLOCK_GFAC):(NAN),
	      mod->xsol[i*BLOCK_NBASE+INT_Y]/BLOCK_GFAC,(i<nrbf)?(mod->sigma[i*BLOCK_NBASE+INT_Y]/BLOCK_GFAC):(NAN),
	      mod->xsol[i*BLOCK_NBASE+INT_Z]/BLOCK_GFAC,(i<nrbf)?(mod->sigma[i*BLOCK_NBASE+INT_Z]/BLOCK_GFAC):(NAN),
	      (mod->changed_reference_frame)?("again in orig RF"):("in orig RF"));
    /* 

    OUTPUT TO FILE: format 

    w_x sig(w_x)  w_y sig(w_y) w_z sig(w_z)       w^r_x w^r_y w^r_z

    where w_i are the best-fit, long Euler vectors and w^r_i the best
    fit rigid ones from the input processing

    */
    fprintf(out,"%.7e %.7e %.7e %.7e %.7e %.7e\t\t%.7e %.7e %.7e\n",
	    mod->xsol[i*BLOCK_NBASE+INT_X]/BLOCK_GFAC,mod->sigma[i*BLOCK_NBASE+INT_X]/BLOCK_GFAC,
	    mod->xsol[i*BLOCK_NBASE+INT_Y]/BLOCK_GFAC,mod->sigma[i*BLOCK_NBASE+INT_Y]/BLOCK_GFAC,
	    mod->xsol[i*BLOCK_NBASE+INT_Z]/BLOCK_GFAC,mod->sigma[i*BLOCK_NBASE+INT_Z]/BLOCK_GFAC,
	    mod->block[i].xrigid[INT_X]/BLOCK_GFAC,mod->block[i].xrigid[INT_Y]/BLOCK_GFAC,
	    mod->block[i].xrigid[INT_Z]/BLOCK_GFAC);
    if(0){
      /* 
	 
      compute the rotational and translational component
    
      */
      /* 
	 long-term solution 
      */
      wm[0] = norm_3d((mod->xsol+i*BLOCK_NBASE));	/* magnitude */
      cross_product((mod->xsol+i*BLOCK_NBASE),mod->block[i].center,vec);
      t[0] = norm_3d(vec)/wm[0];	/* translational component */
      /* rotational component */
      r[0] = fabs(dotp_3d((mod->xsol+i*BLOCK_NBASE),mod->block[i].center))/wm[0];
      /* 
	 rigid block best-fit 
      */
      wm[1] = norm_3d(mod->block[i].xrigid);
      cross_product(mod->block[i].xrigid,mod->block[i].center,vec);
      t[1] = norm_3d(vec)/wm[1];
      r[1] = fabs(dotp_3d(mod->block[i].xrigid,mod->block[i].center))/wm[1];
      fprintf(stderr,"%s: block %c: rigid: t: %.2f r: %.2f |w|: %7.2f long-term: t: %.2f r: %.2f |w|: %7.2f\n",
	      argv[0],bname(i),t[1],r[1],wm[1],t[0],r[0],wm[0]);
    }
    if(domega_out)		/* leave redundant combinations in for
				   convenience in post-processing  */
      for(j=0;j < mod->nrb;j++)
	fprintf(out2,"%c %c %.7e %.7e %.7e\n",bname(i),bname(j),
		(mod->xsol[i*BLOCK_NBASE+INT_X]-mod->xsol[j*BLOCK_NBASE+INT_X])/BLOCK_GFAC,
		(mod->xsol[i*BLOCK_NBASE+INT_Y]-mod->xsol[j*BLOCK_NBASE+INT_Y])/BLOCK_GFAC,
		(mod->xsol[i*BLOCK_NBASE+INT_Z]-mod->xsol[j*BLOCK_NBASE+INT_Z])/BLOCK_GFAC);
#endif
  }
  fclose(out);
  if(domega_out)
    fclose(out2);
#ifdef BLOCK_SPHERICAL
  /* 

  output of rigid plate Euler motion vectors in lon lat mag format,
  possibly referenced to a certain block, and corrected back into the
  original system 

  */
  if(reference_block >= 0){	
    /* 
       reference to certain block 
    */
    if(reference_block >= mod->nrb){
      fprintf(stderr,"%s: error, reference block %i out of range, nrb: %i\n",
	      argv[0],reference_block, mod->nrb);
      exit(-1);
    }
    /* 
       assign reference solution 
    */
    
    a_equals_b_vector_3d(omega_ref,(mod->xsol+reference_block*3));
    fprintf(stderr,"%s: output of relative motion vectors with respect to block %c\n",
	    argv[0],bname(reference_block));
    for(i=0;i < mod->nrb;i++){	/* loop through blocks */
      a_equals_b_vector_3d(omega_block,(mod->xsol+i*BLOCK_NBASE));
      /* 
	 calculated reference solution with respect to omega_ref 
      */
      sub_b_from_a_vector_3d(omega_block,omega_ref);
      calc_geo_euler_pole(omega_block,&lon,&lat,&mag);
      if(norm_3d(omega_ref) < 1e-7){
	fprintf(stderr,"%s: block %c (%2i): lon: %8.3f lat: %8.3f mag: %6.3f (%6.3f %6.3f %6.3f) (deg/Myr, wrt: SCEC)\n",
		argv[0],bname(i),i+1,lon,lat,mag,omega_block[INT_X]/BLOCK_GFAC,omega_block[INT_Y]/BLOCK_GFAC,omega_block[INT_Z]/BLOCK_GFAC);
      }else{
	if(i != reference_block)
	  fprintf(stderr,"%s: block %c (%2i): lon: %8.3f lat: %8.3f mag: %6.3f (%6.3f %6.3f %6.3f) (deg/Myr, wrt: block %c)\n",
		  argv[0],bname(i),i+1,lon,lat,mag,omega_block[INT_X]/BLOCK_GFAC,omega_block[INT_Y]/BLOCK_GFAC,omega_block[INT_Z]/BLOCK_GFAC,
		  bname(reference_block));
	else
	  fprintf(stderr,"%s: block %c (%2i): lon: %8.3f lat: %8.3f mag: %6.3f (%6.3f %6.3f %6.3f) (deg/Myr, reference block)\n",
		  argv[0],bname(i),i+1,lon,lat,mag,omega_block[INT_X]/BLOCK_GFAC,omega_block[INT_Y]/BLOCK_GFAC,omega_block[INT_Z]/BLOCK_GFAC);
      }
    }
    if(vref_out){
      /* 
	 
      output of the predicted velocities (including deformation) with
      respect to the reference block, if any
      
      */
      my_vecrealloc(&omega,mod->nrb*BLOCK_NBASE,"omega");
      for(i=0;i<mod->nrb;i++){
	a_equals_b_vector_3d((omega+i*3),mod->xsol+i*3);
	sub_b_from_a_vector_3d((omega+i*3),omega_ref);
      }
      /* obtain the shifted reference frame solution */
      evaluate_block_solution(mod->kmat,mod->ms,mod->n,omega,
			      mod->vmod,mod->vmodc,mod);
      /* 
	 output of the shifted reference frame solution
      */
      fprintf(stderr,"%s: writing velocities with respect to ref block %c to vrefblock.dat\n",
	      argv[0],bname(reference_block));
      print_simple_vel(mod->gx,mod->vmod,mod->sigv,mod->rho,
		       mod->nrgp,"vrefblock.dat");
      /* re-compute the original solution */
      evaluate_block_solution(mod->kmat,mod->ms,mod->n,mod->xsol,
			      mod->vmod,mod->vmodc,mod);
      /* free memory */
      my_vecrealloc(&omega,1,"omega");
    }
  } /* end reference block output */
#endif
  fprintf(stderr,"%s: solution x norm: %g uncertainty xsigma norm: %g\n",
	  argv[0],norm(mod->xsol,mod->n),norm(mod->sigma,mod->n));
#ifdef DEBUG
  fprintf(stderr,"block_output: writing sigma from SVD to sigma.asc\n");
  out=myopen("sigma.asc","w");
  print_matrix_ftrn(mod->sigma,mod->ncov,1,out,FALSE);
  fclose(out);
#endif    
  i = mod->n;
  if(invert_for_ld){
    /* 
       locking depths, output in format =

       locking_depth lon1 lat1 lon2 lat2
    */
    fprintf(stderr,"%s: mean inverted locking depth: %g, results in %s\n",
	    argv[0],mean((mod->xsol+i),1,mod->nflt),
	    LDOUT_FILE);
    /* write to ld.fit file */
    out=myopen(LDOUT_FILE,"w");
    for(j=0;j < mod->nflt;i++,j++) /* LEAVE this loop in to inc i */
      fprintf(out,"%15.7f %15.7f %15.7f %15.7f %15.7f %11g %6i %6i\n",
	      mod->xsol[i],
	      mod->fault[j].ex[0][INT_X],mod->fault[j].ex[0][INT_Y],
	      mod->fault[j].ex[1][INT_X],mod->fault[j].ex[1][INT_Y],
	      mod->fault[j].dip,
	      mod->fault[j].block[0]+1,mod->fault[j].block[1]+1);
    fclose(out);
  }
  my_vecrealloc(&cfac, mod->nflt,"cfac");
  if(invert_for_cfac){
    /* 
       local coupling factors for each fault 
       (only if invert_for_cfac is set)
    */
    fprintf(stderr,"%s: mean locking factor: %g, results in %s\n",
	    argv[0],mean((mod->xsol+i),1,mod->nflt),
	    CFACOUT_FILE);
    out=myopen(CFACOUT_FILE,"w");
    for(j=0;j < mod->nflt;i++,j++){	
      /* WARNING:

      leave this in here, even if output is not done so that 
      the cfac's will be assigned

      */
      cfac[j] = mod->xsol[i];
      fprintf(out,"%11g\n",mod->xsol[i]);
    }
    fclose(out);
  }else{
    for(j=0;j < mod->nflt;j++)
      cfac[j] = 1.0;
  }
  if(nofiles)			/* bailout here if no file
				   output requested */
    return;
  //
  // covariance output, normalized to deg/Myr
  //
  fprintf(stderr,"%s: writing normalized (deg/Myr**2) covariance matrix to %s\n",
	  argv[0],COVOUT_FILE);
  out=myopen(COVOUT_FILE,"w");
  print_matrix_scaled_ftrn(mod->cov,mod->ncov,mod->ncov,out,FALSE,1.0/(BLOCK_GFAC*BLOCK_GFAC));
  fclose(out);
  /*
    
    output of predicted velocities and differences v_obs - v_mod
    format:
    
    lon lat ve_mod vn_mod (veo-vem) (vno-vnm) sigma(ve_obs) sigma(vn_obs) rho_obs
    
  */
  out=myopen( VELFITOUT_FILE,"w");
  mean_vmisfit = 0.0;
  if(mod->changed_reference_frame){
    if(use_nullspace){
      fprintf(stderr,"%s: block_output: WARNING: reference frame changed, NULLSPACE OUTPUT will be in changed RF\n",
	      argv[0]);
      fprintf(stderr,"%s: block_output: WARNING: look at %s for the correction which would have to be added\n",
	      argv[0],VEL_GR_FILE);
    }else{
      fprintf(stderr,"%s: block_output: reference frame for velocities was changed, output will be in original RF\n",
	      argv[0]);
    }
  }
  for(i=j=0;i < mod->nrgp;i++,j+=BLOCK_DIM){
    if((!use_nullspace)&&(mod->changed_reference_frame)){
      /* 
	 only in the case of a changed reference frame will the vcorp
	 array hold the correction. if we are interested in the nullspace,
	 we will not apply this correction
      */
      vshift[INT_X] = mod->vcorp[j+INT_X];vshift[INT_Y] = mod->vcorp[j+INT_Y];
    }else{
      vshift[INT_X] = vshift[INT_Y] = 0.0;
    }
    fprintf(out,"%12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	    mod->gx[j+INT_X],mod->gx[j+INT_Y],
	    mod->vmod[j+INT_X], 
	    mod->vmod[j+INT_Y],
	    mod->v[j+INT_X] + vshift[INT_X] - mod->vmod[j+INT_X],/* pontentially shift orig velocities back */
	    mod->v[j+INT_Y] + vshift[INT_Y] - mod->vmod[j+INT_Y],
	    mod->sigv[j+INT_X], mod->sigv[j+INT_Y],mod->rho[i]);
    mean_vmisfit += hypot(mod->v[j+INT_X] + vshift[INT_X] - mod->vmod[j+INT_X],
			  mod->v[j+INT_Y] + vshift[INT_Y] - mod->vmod[j+INT_Y]);
  }
  fclose(out);
  fprintf(stderr,"%s: written %i predicted velocities to %s, mean misfit: %g\n",
	  argv[0],mod->nrgp,VELFITOUT_FILE,
	  mean_vmisfit/(COMP_PRECISION)mod->nrgp);
  /* 
     end velocity output
  */
  if(calc_best_fit_euler){
    /* 
    output of best fit Euler rotation vectors. those are different from
    the real, rigid plate, long-term Euler vectors, which are identical
    to the solution parameters for the spherical version. 

    THIS PART OF THE CODE WILL NOT NORMALLY GET EXECUTED BY THE
    SPHERICAL VERSION, SINCE IN THIS CASE, WE CAN USE THE SOLUTION
    VECTORS
    
    
    */
    my_vecrealloc(&omega,mod->nrb*3,"omega");
    for(i=0;i<mod->nrb;i++){	/* find the best-fit rotation 
				   for block i */
      find_spherical_rotation(mod,i,(omega+i*3),TRUE,FALSE,
			      FALSE,FALSE,FALSE,mod->fblock_sites,
			      mod->fblock_nsites);
    }
    /* 
       find Euler poles with respect to reference_block
    */
    if(reference_block >= 0){
      fprintf(stderr,"%s: best fit block motion with respect to block %c\n",
	      argv[0],bname(reference_block));
      if(reference_block >= mod->nrb){
	fprintf(stderr,"%s: error, reference block out of range, nrb: %i\n",
		argv[0],mod->nrb);
	exit(-1);
      }
      /* assign relative motion vector */
      a_equals_b_vector_3d(omega_ref,(omega+reference_block*3));
    }else{
      fprintf(stderr,"%s: no reference block set\n",argv[0]);
    }
    for(i=0;i<mod->nrb;i++){
      if(reference_block >= 0)	/* reduce vector, if reference_block is set */
	sub_b_from_a_vector_3d((omega+i*3),omega_ref);
      calc_geo_euler_pole((omega+i*3),&lon,&lat,&mag);
      fprintf(stderr,"fsr: block %c (%2i): model w: %6.3f, %6.3f, %6.3f (deg/Myr) lon lat r: %8.3f %8.3f %6.3f\n",
	      bname(i),i+1,
	      omega[i*3+INT_X]/BLOCK_GFAC,omega[i*3+INT_Y]/BLOCK_GFAC,
	      omega[i*3+INT_Z]/BLOCK_GFAC,lon,lat,mag);
    }
  }
  /* 

  
  FAULT RELATED PROPERTIES


  */
  frms=0.0;
  if(mod->nflt){
    my_vecrealloc(&vslip, mod->nsnf,"vslip");
    my_vecrealloc(&vsigma, mod->nsnf,"vsigma");
    my_vecrealloc(&disp, 3*mod->nflt,"disp");
    /*
      
      output of total fault displacements at fault midpoints
      in format. total means u = u^- - u^+
      
      clon clat  s_x s_y u_s    (n_x n_y u_n     d_x d_y u_d)
      
      
      calculate the predicted slip in local system
      G . F . b = v_slip[nslip * nflt = nsnf]
      
      G.F links the solution vector b to relative velocities, and then
      to fault local slip
      
    */
    calc_Ax_ftn(mod->gf,mod->nsnf,mod->n,mod->xsol,vslip);
    frms = rms(vslip,mod->nsnf); /* total slip RMS */
    /* calculate the fault strike and normal RMSs */
    calc_fault_sn_rms(vslip,mod->nflt,mod->nslip,&fsrms,&fnrms);
    /*
      
      calculate the uncertainties in the fault local slip 
      by multiplying the solution vector uncertainties in sigma[]
      with G.F
      
      format:
      
      lon lat u^i_x u^i_y u^i sig(u^i) coupling_factor orig_code
      
      for i = strike, normal, and dip depending on nslip
      and sigma indicating uncertainties
      
      coupling_factor is the faults coupling factor from input * the
      cfac the inversion might have obtained

      orig_code is the code (nr in list, starting with zero) of the
      original fault

    */
    calculate_slipsol_sigma(mod,vslip,vsigma,mod->vsig_mode);
    //
    //
    //
    out=myopen(SLIPOUT_FILE,"w");
    for(i=0;i < mod->nflt;i++){
      fprintf(out,"%11g %11g ",mod->fault[i].x[INT_X], mod->fault[i].x[INT_Y]);
      disp[i*3+INT_X] = disp[i*3+INT_Y] = disp[i*3+INT_Z] = 0.0;
      for(j=0;j < mod->nslip;j++){	/* slip modes:
					   
					strike or normal/dip 

					*/
	idir = block_slip_direction(j,(mod->fault+i));
	/*
	  the internal convention for slip on Okada patches is like so:
	  
	  For strike: positive values of slip mean left-lateral fault,
  	  negative right-lateral; 
	  
	  For dip: positive values of slip mean thrust (up-dip motion)
  	  fault, negative normal (down-dip motion) fault;

	  For normal: positive values of slip mean explosive source,
  	  negative implosive;

	  
	  if the second fault slip mode is selected, the fault is
	  non-vertical, and has postive dip, we flip the sign such
	  that up-dip motion, indicative of thrust, gives the same,
	  negative, sign as implosive motion. if the fault has
	  negative dip, this is reversed.
	  
	  (to be checked!)

	*/
	if((j==1)&&(!mod->fault[i].vertical)&&(mod->fault[i].dip > 0))
	  sign = -1.0;
	else
	  sign = 1.0;
	fprintf(out,"%11.4e %11.4e %11.4e %11.4e ",
		vslip[i*mod->nslip+j] * mod->fault[i].evec[idir*3+INT_X],// v_e
		vslip[i*mod->nslip+j] * mod->fault[i].evec[idir*3+INT_Y],// v_u
		sign * vslip[i*mod->nslip+j],                   // u_i
		fabs(vsigma[i*mod->nslip+j]));// sig(u_i)

	//
	// assign to global slip array 
	disp[i*3+block_slip_direction(j,(mod->fault+i))] = 
	  (rigid)?(0.0):
	  (vslip[i*mod->nslip+j] * mod->fault[i].lfac * cfac[i]);
      }
      for(j=mod->nslip;j < 2;j++)
	fprintf(out,"0. 0. 0. 0. ");
      //
      // coupling factor, the one we read times the one we might
      // have inverted for. then, the code
      //
      fprintf(out,"  %g %i\n",mod->fault[i].lfac * cfac[i],
	      mod->fault[i].orig_code);
    }
    fclose(out);
    fprintf(stderr,"%s: written %i types of fault slip for %i faults %s to %s\n",
	    argv[0],mod->nslip,mod->nsnf/mod->nslip,
	    (rigid)?("(rigid blocks)"):(""),SLIPOUT_FILE );
    //
    // binary solution file with projection, solution vector
    // fault geometry and fault slip vectors
    //
    out = myopen(SOLOUT_FILE,"w");
    block_save_solution_and_faults(mod->xsol,mod->nrb,mod->nflt,
				   mod->fault,disp,&projection,
				   out,invert_for_ld,invert_for_cfac);
    fclose(out);
    fprintf(stderr,"%s: saved faults with slips (norm: %12g) to %s\n",
	    argv[0],norm(disp,3*mod->nflt),SOLOUT_FILE);
    /*
      
      
    OUTPUT OF PREDICTED STRESSES: 
    
    those are given by the background (or zero) minus the effect of
    the fault slip (deformation also enters as minus the slip times
    D) fault slip vslip is D.F.x, so we could use assemble_stress
    matrix and the slip from above like
    
    assemble_stress_matrix(smod,i,fault,vslip,nflt,nslip);
    
    for stress observation i. however, we have already calculated
    the solution vmod, which holds the stress tensor components
    at the i stress observations above mgd
    
    */
    out2 = myopen(SFITOUT_FILE,"w"); /* stress model */
    out3 = myopen(SCSTRESSOUT_FILE,"w"); /* scaled observations */
    for(os1=mod->mgd,i=0;i < mod->nrsp;i++,os1+=6){/* loop through all stress points */
      // coordinates
      fprintf(out2,"%11g %11g ",
	      mod->sx[BLOCK_DIM*i+INT_X],mod->sx[BLOCK_DIM*i+INT_Y]);
      fprintf(out3,"%11g %11g ",
	      mod->sx[BLOCK_DIM*i+INT_X],mod->sx[BLOCK_DIM*i+INT_Y]);
      // print to stress.fit file
      for(j=0;j < 6;j++){
	fprintf(out2,"%11g ",mod->vmod[os1+j]);
	fprintf(out3,"%11g ",mod->v[os1+j]);
      }
      fprintf(out2,"\n");fprintf(out3,"\n");
    }
    fclose(out2);fclose(out3);
    if(mod->nrsp){
      fprintf(stderr,"%s: written model stress           in six component format to %s, rms: %g\n",
	      argv[0],STRESSDATA_FILE,rms((mod->vmod+mod->mgd),mod->m2));
      fprintf(stderr,"%s: written scaled observed stress in six component format to %s, rms: %g\n",
	      argv[0],SCSTRESSOUT_FILE,rms((mod->v+mod->mgd),mod->m2));
    }
  }
  //
  // print variance reduction, total chi^2, vchi2, schi2, vel_rms, stress_rms, fault_rms
  //
  //
  fprintf(stderr,"%s: 1-sqrt(chi^2/|v/sigv|^2): %12.8f chi^2: %g vc2: %g sc2: %g vrms: %g srms: %g frms: %g %g %g\n",
	  argv[0],1.0 - sqrt(chi2/vnorm2),chi2,c2a,c2b,
	  rms(mod->vmod,mod->mgd),rms((mod->vmod+mod->mgd),mod->m2),
	  frms,fsrms,fnrms);
  /* 
     reduced variance reduction:

     \hat{\chi}^2  = \frac{\chi^2}{\nu_{GPS} + \nu_{\tau} - M - 1}

  */
  fprintf(stderr,"%s: reduced chi^2: %g (c:%g ng:%i mgd: %i  b: %g m2: %i ns:%i nrb:%i n: %i nc: %i ncov: %i)\n",argv[0],
	  chi2/((COMP_PRECISION)(mod->mgd + ((beta != 0)?(mod->m2):(0)) - mod->ncov - 1)),
	  chi2,mod->nrgp,mod->mgd,beta,mod->m2,mod->nrsp,mod->nrb,mod->n,mod->nc,mod->ncov);
  if(mod->nfdamp > 0){
    gamma_rms = 0.0;
    for(i=0;i < mod->nflt;i++)
      for(j=0;j < mod->nslip;j++)
	if(mod->fault[i].use_damp[j])
	  gamma_rms += square(mod->fault[i].sdamp[j]);
    gamma_rms = sqrt(gamma_rms / (COMP_PRECISION)mod->nfdamp);
    fprintf(stderr,"%s: damped slip residual rms/damping rms: %11g\n",
	    argv[0],rms((mod->vmod+mod->mgd+mod->m2),
			mod->nfdamp)/gamma_rms);
  }
  free(vslip);free(vsigma);free(disp);free(cfac);
}
/* 
   sort eigenvalues according to absolute value of eigenvalue
   which is input as f[3] vector, on return isort[3] will hold
   the indices of f such that 

*/
void sort_eigen(int *isort, COMP_PRECISION *f)
{
  int i;
  COMP_PRECISION max,tmp;
  max = 0.0;isort[0]=0;
  for(i=0;i<3;i++){	/* find largest entry */
    if((tmp = fabs(f[i])) > max){
      isort[0] = i;
      max = tmp;
    }
  }
  if(isort[0] == 0){
    isort[1] = 1;isort[2] = 2;
  }else if(isort[0] == 1){
    isort[1] = 0;isort[2] = 2;
  }else{
    isort[1] = 0;isort[2] = 1;
  }
  if(fabs(f[isort[1]]) < fabs(f[isort[2]])){
    i = isort[1];
    isort[1] = isort[2];
    isort[2] = i;
  }
}
/*

  print out horizontal component of stresses, 
  expected s matrix is in full format.

  output is e1 e2 azi(e1) (CW from north)

*/
void print_horizontal_stress(COMP_PRECISION *s,FILE *out)
{
  COMP_PRECISION loc_azi,evec[4],eval[2],sloc2d[4];
  //
  // assign horizontal components of stress tensor
  sloc2d[0] = s[INT_XX];
  sloc2d[1] = sloc2d[2] = s[INT_XY];
  sloc2d[3] = s[INT_YY];
  /* this routine returns values and vectors sorted as
     e2 < e1 
  */
  calc_eigensystem_sym2d(sloc2d,eval,evec,TRUE);
  loc_azi = vec_to_strike((evec+2));/* azimuth of e1 */
  fix_azimuth(&loc_azi);
  if(loc_azi > 180)		/* 0 .. azi .. 180 */
    loc_azi -= 180;
  fprintf(out,"%12.4e %12.4e %6.2f\n", 
	  eval[1], eval[0], loc_azi);
}
/*

  print the horizontal projection of the eigensystem of the 
  s matrix, expected and full format
  
*/
void print_projected_stress(COMP_PRECISION *s, FILE *out)
{
  COMP_PRECISION evec[9],eval[3],evalh[3],loc_azi;
  int j,isort[3];
  //
  // get the eigensystem
  calc_eigensystem_sym3d(s,eval,evec,TRUE);
  // assign hor components of eigenvectors
  for(j=0;j<3;j++)		/* horizontal component */
    evalh[j] = eval[j] * sqrt(1.0-SQUARE(evec[j*3+INT_Z]));
  sort_eigen(isort, eval);	/* sort by abs(EV) */
  for(j=0;j<3;j++){		
    loc_azi = vec_to_strike((evec+isort[j]*3));
    fix_azimuth(&loc_azi);
    //
    // output of horizontal component of eigenvalues, 
    // azimuth, and total eigenvalue 
    // sorted by total absolute eigenvalue size
    //
    fprintf(out,"%12.4e %6.2f %12.4e    ",
	    evalh[isort[j]],loc_azi,eval[isort[j]]);
  }
  fprintf(out,"\n");
}
/* 

calculate the uncertainties in fault slip values 
vslip is nsnf and input

vsigma is a nflt * nslip (nsnf) array, and will be output

mode 0: uses different combinations of the fault block vectors + 
        uncertainties

     1: uses randomized approach

*/
void calculate_slipsol_sigma(struct bmd *mod,COMP_PRECISION *vslip,
			     COMP_PRECISION *vsigma,int mode)
{
  COMP_PRECISION *test_vslip,*test_xsol,tmp_val,*vmean,*ostd,bailout;
  int i,j,k,l,m,klim,llim,lblock[2],tncount;
  /* 
     default values 
  */
  static my_boolean use_max = FALSE; /* for mode = 0, use max instead of 
					RMS for uncertainties if true
				     */
  static int nrr_max = 125000;	/* 
				   max number of random realizations
				   for mode = 1 this was determined
				   based on the mean random slip
				   compared with the mean slip of the
				   original solution. loop will be terminated 
				   earlier, if convergence criterion reached
				*/
  
  my_vecalloc(&test_vslip,mod->nsnf,"test_vslip");
  my_vecalloc(&vmean,mod->nsnf,"vmean");
  my_vecalloc(&ostd,mod->nsnf,"ostd");
  my_vecalloc(&test_xsol,mod->n,"test_xsol");
  /* 
     check if we have sigma uncertainties 
  */
  if(norm(mod->sigma,mod->n) < EPS_COMP_PREC){
    /* 
       we don't, return vsigma as zero as well
    */
    for(i=0;i<mod->nsnf;i++)
      vsigma[i] = 0.0;
    return;
  }
  switch(mode){
  case 0:
    /* 


    using a combination of the block motion vector uncertainties to 
    evaluate the slip rate deviations



    */
    fprintf(stderr,"calculate_slipsol_sigma: calculating vsigma from omega + sigma combinations\n");
    for(i=0;i<mod->nflt;i++){	/* fault loop */
      for(j=0;j<mod->nslip;j++)	/* init sigma array */
	vsigma[i*mod->nslip+j] = 0.0;
      tncount = 0;			/* init nr of contributions array */
      for(j=0;j<2;j++)		/* codes of blocks bordering fault */
	lblock[j] = mod->fault[i].block[j];
      for(j=0;j < 4;j++){		/* loop through block combinations */
	/* 
	   j=0: +0
	   j=1: 0+
	   j=2: ++ (same as --)
	   j=3: +- (same as -+)
	   
	*/
	if(j == 1)		/* loop through all directions
				   in first block only 
				   if contribution needed*/
	  klim = 1;
	else
	  klim = BLOCK_NBASE;
	if(j == 0)		/* same for second */
	  llim = 1;
	else
	  llim = BLOCK_NBASE;
	for(k=0;k < klim;k++){ /* first direction loop */
	  for(l=0;l < llim;l++){ /* second direction loop */
	    /* copy original solution vector  */
	    a_equals_b_vector(test_xsol,mod->xsol,mod->n);
	    if((j==0)||(j==2)||(j==3)){	/* left block + */
	      test_xsol[lblock[0]*BLOCK_NBASE+k] +=
		mod->sigma[lblock[0]*BLOCK_NBASE+k];
	    }
	    if((j==1)||(j==2)){	/* right block + */
	      test_xsol[lblock[1]*BLOCK_NBASE+l] +=
		mod->sigma[lblock[1]*BLOCK_NBASE+l];
	    }
	    if((j==3)){	/* right block - */
	      test_xsol[lblock[1]*BLOCK_NBASE+l] -=
		mod->sigma[lblock[1]*BLOCK_NBASE+l];
	    }
	    /* calculate a test fault slip solution for the test
	       solution vector test_xsol */
	    calc_Ax_ftn(mod->gf,mod->nsnf,mod->n,test_xsol,test_vslip);
	    //
	    tncount++;
	    for(m=0;m < mod->nslip;m++){
	      tmp_val = fabs(test_vslip[i*mod->nslip+m] - vslip[i*mod->nslip+m]);
	      if(use_max){
		if(vsigma[i*mod->nslip+m] < tmp_val) /* find max */
		  vsigma[i*mod->nslip+m] = tmp_val;
	      }else{
		vsigma[i*mod->nslip+m] += tmp_val; /* add up to mean */
	      }
	      //if(m == 0)fprintf(stderr,"flt: %3i mode: %i du_s: %11g\n",i,j,tmp_val);
	    } /* end m loop */
	  } /* end llim loop */
	} /*  end klim loop */
      }	/* end j loop */
      if(!use_max){
	for(j=0;j < mod->nslip;j++){	/* calculate mean */
	  vsigma[i*mod->nslip+j] /= (COMP_PRECISION)tncount;
	  //if(j==0)fprintf(stderr,"flt: %3i mean du_s: %11g\n",i,vsigma[i*mod->nslip+j]);
	}
      }
    } /* end the i-fault loop  */
    break;			/* end mode = 0 */
  case 1:			
    /* 

    
       USE RANDOM APPROACH

       This will assign a one standard deviation range to all the
       vsigma values based on a number of random realizations


    */
    fprintf(stderr,"calculate_slipsol_sigma: calculating vsigma from std of random realizations\n");
    for(i=0;i < mod->nsnf;i++){
      vsigma[i] = 0.0;		/* this will be the std sum */
      vmean[i] = 0.0;		/* this will be the mean sum  */
    }
    bailout = norm(mod->xsol,mod->n) * 1e-6;
    for(tncount=0;tncount < nrr_max;){		/* begin random loop */
      /* 
	 copy original solution vector  
      */
      a_equals_b_vector(test_xsol,mod->xsol,mod->n);
#ifdef DEBUG
      if(tncount != 0){
#endif
	/* add uncertainties with Gaussian distribution */
	for(j=0;j < mod->n;j++)
	  test_xsol[j] += mygauss_randnr(mod->sigma[j],&mod->seed);
#ifdef DEBUG
      }
#endif
      /* 
	 generate test slip solution 
      */
      calc_Ax_ftn(mod->gf,mod->nsnf,mod->n,test_xsol,test_vslip);
      /* 
	 sum up differences for the standard deviation computation 
      */
      for(j=0;j < mod->nsnf;j++){
	vmean[j]  += test_vslip[j];
	vsigma[j] += SQUARE(test_vslip[j]);
      }
#ifdef DEBUG
      if(tncount == 0){
	/* 
	   do a little test
	*/
	if(distance(test_vslip,vslip,mod->nsnf) > EPS_COMP_PREC){
	  fprintf(stderr,"calculate_slipsol_sigma: internal error for vslip/xsol\n");
	  exit(-1);
	}
      }
#endif
      /* increment global and local counters */
      tncount++;
      if(tncount%500==0){		/* every 500, check for
					   convergence of random
					   estimates */
	tmp_val = 0.0;
	for(j=0;j < mod->nsnf;j++)
	  tmp_val += square(ostd[j] - stddev(vmean[j],vsigma[j],tncount));
	tmp_val = sqrt(tmp_val/(COMP_PRECISION)mod->nsnf);
	if(tmp_val < bailout){
	  /* bailout */
	  nrr_max = tncount;
	}
      }
      /* 
	 compute running standard deviation
      */
      for(j=0;j < mod->nsnf;j++)
	ostd[j] = stddev(vmean[j],vsigma[j],tncount);
    } /* end random loop */
    /* calculate the std and deviation of mean from original solution */
    tmp_val = 0.0;
    for(i=0;i < mod->nflt;i++){
      for(j=0;j<mod->nslip;j++){
	k = i * mod->nslip + j;
	vsigma[k] = ostd[k];	/* has been computed already */
	vmean[k] /= (COMP_PRECISION)tncount;	/* mean */
	tmp_val += square(vslip[k] - vmean[k]); /* for RMS deviation */
      }
    }
    fprintf(stderr,"calculate_slipsol_sigma: done: %i random steps, fractional RMS of (mean sol - orig sol): %g\n",
	    tncount,sqrt(tmp_val/(COMP_PRECISION)mod->nsnf)/norm(vslip,mod->nsnf));
    break;
  default:
    fprintf(stderr,"calculate_slipsol_sigma: error: mode %i undefined\n",
	    mode);
    exit(-1);
    break;
  } /* end switch loop */

  free(vmean);free(ostd);
  free(test_xsol);free(test_vslip);
}
/* 

given a cartesian rotation vector w[3] in internal normalization,
return a geographic rotation pole with magnitude deg/Myr

*/

void calc_geo_euler_pole(COMP_PRECISION *w, 
			 COMP_PRECISION *lon,
			 COMP_PRECISION *lat, 
			 COMP_PRECISION *mag)
{
  xyz2lonlat(w,lon,lat);
  *lon=RAD2DEGF(*lon); 
  *lat=RAD2DEGF(*lat);
  *mag = norm_3d(w)/BLOCK_GFAC;
  if(*lon < 0){			/* pick the other pole */
    *lon += 180;
    *lat = - *lat;
    *mag = - (*mag);
  }

}
/* 

print a simple velocity file (different from the fit-type
file which has the deviations inlcuded, too)

input:

gx[] v[] vsig[] are nrgp * BLOCK_DIM
rho[] is nrgp dimensionalize

output format:

lon lat ve vn sigve sigvn rho

*/
void print_simple_vel(COMP_PRECISION *gx,
		      COMP_PRECISION *v,
		      COMP_PRECISION *vsig,
		      COMP_PRECISION *rho,
		      int nrgp,char *filename)
{
  int i,j;
  FILE *out;
  
  out = myopen(filename,"w");
  for(i=j=0;i < nrgp;i++,j+=BLOCK_DIM)
    fprintf(out,"%12g %12g %12g %12g %12g %12g %12g\n",
	    gx[j+INT_X],gx[j+INT_Y],
	    v[j+INT_X],v[j+INT_Y],vsig[j+INT_X],vsig[j+INT_Y],rho[i]);
  fclose(out);
}
