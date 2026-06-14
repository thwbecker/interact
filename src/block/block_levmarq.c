#include "interact.h"
#include "blockinvert.h"
/* 

part of blockinvert: handle the levenberg marquardt solution
iteraction

$Id: block_levmarq.c,v 1.9 2004/03/25 23:48:47 becker Exp $

*/
void run_lm(struct bmd *mod,long int *seed,
	    struct prj projection,COMP_PRECISION *chi2,
	    my_boolean rigid,my_boolean beta_init, COMP_PRECISION beta,
	    COMP_PRECISION *stress_rms,
	    my_boolean no_rel_stress_amp_scale,
	    my_boolean damp_nslip,COMP_PRECISION sscale_dfac,
	    my_boolean invert_for_ld, my_boolean invert_for_cfac,
	    my_boolean constrain_slip_direction,
	    int nlmloop,int nlmres, char **argv)
{

  COMP_PRECISION **nr_covar,**nr_alpha,fac,alamda,minchi2[3],
    oldchi2[3],*xfix,*xtry;
  int *ia,i,j,iter_count,ibailout,lm_restart;
  /* 

  those allocate arrays that scale with the number of parameters

  */
  my_vecalloc(&xfix,mod->na,"run_lm: xfix");
  my_vecalloc(&xtry,mod->na,"run_lm: xtry");
  my_ivecalloc(&ia,mod->na,"run_lm: ia");
  nr_covar = matrix(1,mod->na,1,mod->na);
  nr_alpha = matrix(1,mod->na,1,mod->na);
  /*
    
  begin general Levenberg Marquardt solution step
  (stresses will not change in amplitude permanently)
  
  */
  /* 
  
  copy solution from linear step to start vector, xfix

  */
  a_equals_b_vector(xfix,mod->xsol,mod->na);   
  /* 
     calculate solution 
  */
  evaluate_block_solution(mod->kmat,mod->ms,mod->n,xfix,
			  mod->vmod,mod->vmodc,mod);
  if((!rigid) && (mod->nrsp))   {
    /* 
       restore scale stresses model 
    */
    rescale_observed_stresses(mod->vmod,mod->v,mod->sigv,
			      sscale_dfac,stress_rms,
			      no_rel_stress_amp_scale,mod,TRUE,
			      FALSE);
  }
  if(!beta_init){
    fprintf(stderr,"%s: error, beta not initialized\n",argv[0]);
    exit(-1);
  }else{
    fprintf(stderr,"%s: starting LM with beta: %g, discarding linear solution\n",
	    argv[0],beta);
  }
  /* 
     initialize ia vector with 1/0 for fit/non-fit 
  */
  for(i=0;i < mod->na;i++)  /* allow all parameters to vary */
    ia[i] = 1;
  /* 
     minimum misfits 
  */
  minchi2[0] = minchi2[1] = minchi2[2] = FLT_MAX;
  /*  
      number of times we restart the LM iteration 
      
  */
  for(lm_restart=0;lm_restart < nlmres;lm_restart++){
    /*
      
    initialize test solution vector with random variations
    around solution from linear inversion
    
    */
    a_equals_b_vector(xtry,xfix,mod->na); 
    if(lm_restart == 0)
      fac = 0.0;		/* start first round
				   with the values from 
				   the linear inversion */
    else if(lm_restart < 5)
      fac = rms(xtry,mod->n);	/* these variations are 
				   for the block parameters 
				   only */
    else
      fac = rms(xtry,mod->n) * 100.0;
    if(fac != 0.0)
      for(i=0;i < mod->n;i++)
	xtry[i] += (-.5 + myrand(seed))*fac;
    /* 
       locking depths and slip locking factors
       will be randomized from the second pass on
    */
    assign_additional_sol_values(xtry,mod->n,mod->na,
				 mod->nflt,mod->fault,
				 invert_for_ld,invert_for_cfac,
				 seed,
				 (lm_restart==0)?(INIT_ADD_SOL):
				 (RANDOM_ADD_SOL));
    /* 
       initialize for each LM iteration loop
    */
    iter_count = ibailout = 0;
    /*

    initialize mrqmin

    */
    alamda = -1.0;	
    /* 
       call numrec vectors with numrec shift, others like
       normal 
    */
    mrqmin((mod->v-1),(mod->sigv-1),
	   (mod->mgd + mod->m2 + mod->nfdamp),
	   (xtry-1),(ia-1),mod->na,nr_covar,
	   nr_alpha,chi2,&alamda,mod,rigid,
	   beta,invert_for_ld,invert_for_cfac,
	   constrain_slip_direction,damp_nslip,
	   sscale_dfac,stress_rms,
	   no_rel_stress_amp_scale,projection);
    /* 
       
    start new LM iteration  
    
    */
    ibailout = 0;
    do{
      iter_count++;
      a_equals_b_vector(oldchi2,chi2,3);
      mrqmin((mod->v-1),(mod->sigv-1),
	     (mod->mgd + mod->m2 + mod->nfdamp), /* number of real data + damping */
	     (xtry-1),(ia-1),mod->na,nr_covar,nr_alpha,chi2,
	     &alamda,mod,rigid,beta,invert_for_ld,
	     invert_for_cfac,constrain_slip_direction,
	     damp_nslip,sscale_dfac,stress_rms,
	     no_rel_stress_amp_scale,projection);
      if(1){
	/* 
	   report LM iteration progress 
	*/
	print_lm_progress(chi2[0],oldchi2[0],chi2[1],chi2[2],
			  lm_restart,ibailout,xtry,alamda,
			  iter_count,argv,mod->n,mod->nflt,
			  invert_for_ld,invert_for_cfac);
      }
      if((chi2[0] - oldchi2[0] > 0) ||
	 ((chi2[0]-oldchi2[0] == 0)&&(alamda < 10000)))
	ibailout = 0; /* don't bailout when chi2 
			 increases */
      else if((oldchi2[0] - chi2[0])/oldchi2[0] < 3e-3)
	ibailout++;
    }while((ibailout < 3) && (iter_count < nlmloop) && 
	   (alamda < 1e10) && (alamda > 1e-20));
    if(iter_count == nlmloop){
      fprintf(stderr,"%s: WARNING: max LM iteration number (%i) reached without convergence\n",
	      argv[0],nlmloop);
    }
    /* 
       terminate LM and calculate covariance (we have to do this
       for each iteration for the malloc/free routines to work)
    */
    alamda = 0.0;
    mrqmin((mod->v-1),(mod->sigv-1),(mod->mgd + mod->m2 + mod->nfdamp),
	   (xtry-1),(ia-1),mod->na,nr_covar,nr_alpha,chi2,&alamda,mod,
	   rigid,beta,invert_for_ld,invert_for_cfac,
	   constrain_slip_direction,damp_nslip,sscale_dfac,
	   stress_rms,no_rel_stress_amp_scale,projection);
    fprintf(stderr,"%s: bailout: c2: %12g iter: %4i rfac: %12g vc2: %12g sc2: %12g srms: %12g\n",
	    argv[0],chi2[0],iter_count,fac,*(chi2+1),*(chi2+2),*stress_rms);
    if(chi2[0] < minchi2[0]){		/* decrease in best chi^2 */
      a_equals_b_vector(minchi2,chi2,3);
      // save solution
      a_equals_b_vector(mod->xsol,xtry,mod->na);
      /* assign local covariance to cov */
      mod->ncov = mod->na;
      my_vecrealloc(&mod->cov,mod->ncov*mod->ncov,"cov");
      for(j=0;j < mod->ncov;j++)		/* asssign to local cov matrix */
	for(i=0;i<mod->ncov;i++)
	  mod->cov[j*mod->ncov+i] = nr_covar[i+1][j+1];
    }
  } /* end randomization (nr of LM iteration restart) loop */
  free_matrix(nr_covar,1,mod->na,1,mod->na);
  free_matrix(nr_alpha,1,mod->na,1,mod->na);
  free(ia);free(xfix);free(xtry);
}
