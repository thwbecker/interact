/*

  numerical recipes in C routines for levenberg marquardt 
  nonlinear fitting 

  as always, everything in here works by addressing
  v[1...n], not v[0...n-1], as is usual 
  
  yes, this is weird
  
  see p. 685 of numerical recipes

  

  $Id: levmarq_numrec.c,v 1.8 2004/03/25 23:48:47 becker Exp $


*/
#include "interact.h"
#include "blockinvert.h"


#define NUMREC_NR_END 1
#define NUMREC_FREE_ARG char*

/* 

*/
void mrqmin(COMP_PRECISION *y,COMP_PRECISION *sig,int ndata,
	    COMP_PRECISION *a,int *ia, int ma, 
	    COMP_PRECISION **covar,
	    COMP_PRECISION **alpha,COMP_PRECISION *chisq,
	    COMP_PRECISION *alamda, 
	    /* from here, variables for fit */
	    struct bmd *mod,my_boolean fit_rigid,
	    COMP_PRECISION fit_beta, my_boolean invert_for_ld,
	    my_boolean invert_for_cfac,
	    my_boolean constrain_slip_direction,
	    my_boolean damp_nslip,
	    COMP_PRECISION fit_damp_fac,
	    COMP_PRECISION *stress_rms,
	    my_boolean no_stress_amp_scale,
	    struct prj fit_projection)
{
  int j,k,l,m;
  static int mfit;
  static COMP_PRECISION ochisq[3],
    *atry,*beta,*da,**oneda,sfac=10.;
  
  if (*alamda < 0.0) {
    /* 


       initialization step 

    */
    /* atry, beta, da allocation */
    my_vecalloc(&atry,ma+NUMREC_NR_END,"mrqmin: atry");
    my_vecalloc(&beta,ma+NUMREC_NR_END,"mrqmin: beta");
    my_vecalloc(&da,ma+NUMREC_NR_END,"mrqmin: da");
    for (mfit=0,j=1;j<=ma;j++)
      if (ia[j]) 
	mfit++;
    fprintf(stderr,"mrqmin: init: fitting %i out of %i parameters\n",
	    mfit,ma);
    /* oneade matrix allocation */
    oneda=matrix(1,mfit,1,1);
    *alamda=0.001;
    fprintf(stderr,"mrqmin: calling mrqcof for init\n");
    mrqcof(y,sig,ndata,a,ia,ma,alpha,beta,chisq,mod,
	   fit_rigid,fit_beta,invert_for_ld,
	   invert_for_cfac,constrain_slip_direction,
	   damp_nslip,fit_damp_fac,stress_rms,
	   no_stress_amp_scale,fit_projection);
    a_equals_b_vector(ochisq,chisq,3);
    for (j=1;j<=ma;j++) 
      atry[j] = a[j];
  } /* end init step */
  for (j=0,l=1;l<=ma;l++) {
    /* 

       loop through all parameters and assign to covar 
       and oneda
    */

    if (ia[l]) {
      for (j++,k=0,m=1;m<=ma;m++) {
	if (ia[m]) {
	  k++;
	  covar[j][k]=alpha[j][k];
	}
      }
      covar[j][j]=alpha[j][j]*(1.0+(*alamda));
      oneda[j][1]=beta[j];
    }
  } /* end parameters loop */
  gaussj(covar,mfit,oneda,1);
  for (j=1;j<=mfit;j++) 
    da[j]=oneda[j][1];
  if (*alamda == 0.0) {
    /* 
       
    temination step: 
    compute covariance and free memory 

    */
    covsrt(covar,ma,ia,mfit);
    fprintf(stderr,"mrqmin: freeing memory\n");
    free_matrix(oneda,1,mfit,1,1);
    free(da);
    free(beta);
    free(atry);
    return;
  }
  for (j=0,l=1;l<=ma;l++)
    if (ia[l]) 
      atry[l]=a[l]+da[++j];
  /*
    if(j==ma)
    fprintf(stderr,"|x|: %12g |dx|: %12g |beta|: %12g\n",
    norm((a+1),ma),norm((da+1),ma),
    norm((beta+1),ma));
  */
  mrqcof(y,sig,ndata,atry,ia,ma,covar,da,chisq,mod,
	 fit_rigid,fit_beta,invert_for_ld,invert_for_cfac,
	 constrain_slip_direction,damp_nslip,
	 fit_damp_fac,stress_rms,no_stress_amp_scale,
	 fit_projection);
  if (chisq[0] < ochisq[0]) {
    /* 
       
    save this solution
    
    */
    *alamda /= sfac;
    a_equals_b_vector(ochisq,chisq,3);
    for (j=0,l=1;l<=ma;l++) {
      if (ia[l]) {
	for (j++,k=0,m=1;m<=ma;m++) {
	  if (ia[m]) {
	    k++;
	    alpha[j][k] = covar[j][k];
	  }
	}
	beta[j] = da[j];
	a[l] = atry[l];
      }
    }
  } else {
    *alamda *= sfac;
     a_equals_b_vector(chisq,ochisq,3);
  }
}

void print_lm_progress(COMP_PRECISION chi2, COMP_PRECISION ochi2, 
		       COMP_PRECISION vchi2,COMP_PRECISION schi2,
		       int lm_restart, int ibailout,
		       COMP_PRECISION *x,COMP_PRECISION alamda,
		       int isc,char **argv,int n, int nflt,
		       my_boolean invert_for_ld,
		       my_boolean invert_for_cfac)
{
  COMP_PRECISION stat[4];
  fprintf(stderr,"%s:LM: restart:%i iter:%3i ib:%i chi^2:%12g rdchi^2: %12g vc2: %11g sc2: %11g al:%11g  ",
	  argv[0],lm_restart,isc,ibailout,
	  chi2,(chi2-ochi2)/ochi2,
	  vchi2,schi2,alamda);
  if(invert_for_ld){
    calc_vec_stat((x+n),nflt,stat);
    fprintf(stderr,"<ld>: %5.2f(%5.2f) ",stat[0],stat[1]);
  }
  if(invert_for_cfac)
    fprintf(stderr,"<cf>: %12g ",mean((x+n+nflt),1,nflt));
  fprintf(stderr,"\n");
}

void covsrt(COMP_PRECISION **covar,int ma,int *ia,int mfit)
{
  int i,j,k;
  COMP_PRECISION swap;
  
  for (i=mfit+1;i<=ma;i++)
    for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
  k=mfit;
  for (j=ma;j>=1;j--) {
    if (ia[j]) {
      for (i=1;i<=ma;i++) 
	NUMREC_SWAP(covar[i][k],covar[i][j]);
      for (i=1;i<=ma;i++) 
	NUMREC_SWAP(covar[k][i],covar[j][i]);
      k--;
    }
  }
}

void gaussj(COMP_PRECISION **a,int n,COMP_PRECISION **b,
	    int m)
{
  int *indxc,*indxr,*ipiv;
  int i,icol=0,irow=0,j,k,l,ll;
  COMP_PRECISION big,dum,pivinv,swap;
  my_ivecalloc(&indxc,n+NUMREC_NR_END,"gaussj: indxc");
  my_ivecalloc(&indxr,n+NUMREC_NR_END,"gaussj: indxr");
  my_ivecalloc(&ipiv,n+NUMREC_NR_END,"gaussj: ipiv");
  for (j=1;j<=n;j++) 
    ipiv[j]=0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if (ipiv[j] != 1)
	for (k=1;k<=n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) 
	    nrerror("gaussj: Singular Matrix-1");
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=1;l<=n;l++) 
	NUMREC_SWAP(a[irow][l],a[icol][l]);
      for (l=1;l<=m;l++) 
	NUMREC_SWAP(b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) 
      nrerror("gaussj: Singular Matrix-2");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=1;l<=n;l++) 
      a[icol][l] *= pivinv;
    for (l=1;l<=m;l++) 
      b[icol][l] *= pivinv;
    for (ll=1;ll<=n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=1;l<=n;l++) 
	  a[ll][l] -= a[icol][l]*dum;
	for (l=1;l<=m;l++) 
	  b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n;l>=1;l--) {
    if (indxr[l] != indxc[l])
      for (k=1;k<=n;k++)
	NUMREC_SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free(ipiv);
  free(indxr);
  free(indxc);
}
/* 


evaluate the trial solution a


*/
void mrqcof(COMP_PRECISION *y,COMP_PRECISION *sig,int ndata,
	    COMP_PRECISION *a,int *ia,int ma,
	    COMP_PRECISION **alpha,COMP_PRECISION *beta,
	    COMP_PRECISION *chisq,struct bmd *mod,
	    my_boolean fit_rigid,COMP_PRECISION fit_beta,
	    my_boolean invert_for_ld,my_boolean invert_for_cfac,
	    my_boolean constrain_slip_direction,
	    my_boolean damp_nslip,
	    COMP_PRECISION fit_damp_fac,
	    COMP_PRECISION *stressrms,
	    my_boolean no_stress_amp_scale,
	    struct prj fit_projection)
{
  int i,j,k,l,m,ilim,mfit=0;
  COMP_PRECISION *ymod,wt,sig2i,dy,*dyda,penalty,tchi2[3],
    *vslip;
#ifdef DEBUG
  fprintf(stderr,"mrqcof: starting alpha/beta assembly\n");
#endif 
  dyda = NULL;
  my_vecalloc(&ymod,ndata+ NUMREC_NR_END,"mrqcof: ymod");
  for (j=1;j<=ma;j++)
    if (ia[j]) 
      mfit++;
  for (j=1;j <= mfit;j++) {
    for (k=1;k <= j;k++) 
      alpha[j][k]=0.0;
    beta[j]=0.0;
  }
  /* 
     modify solution vector accoding to additional 
     constraints 
  */
  check_solution_vector((a+1),mod->n,mod->nflt,invert_for_cfac,
			invert_for_ld);
#ifdef DEBUG
  fprintf(stderr,"mrqcof: solution vector check OK\n");
#endif 

  /*
    
    evaluate the model ymod[1...ndata,m] at 
    a[1....ma=(fit_n+fit_nflt)]
    move all numrec type pointers up by one location

  */
  block_assemble_fit_vector(&mod->kmat,mod->m,mod->ms,
			    mod->m1,mod->m2,
			    mod->n,mod->nrgp,mod->nrsp,
			    mod->nflt,mod->nslip,
			    (a+1),-1,(ymod+1),mod->vmodc,TRUE,
			    mod->a,&mod->d,mod->g,
			    &mod->imat,mod->nrb,mod->block,
			    mod->nsnf, &mod->fault,
			    fit_rigid,&mod->gf,
			    invert_for_ld,invert_for_cfac,
			    invert_for_ld,mod->gx, mod->sx,
			    mod->stress_depths,fit_projection,
			    damp_nslip,mod->nfdamp,mod->nxdamp,
			    mod->xdamp,mod);
#ifdef DEBUG
  fprintf(stderr,"mrqcof: assemble fit vector OK\n");
#endif 

  /*
    
    evaluate the dyda matrix at a
    
  */
  block_assemble_dyda_matrix(&dyda,&mod->kmat,mod->m,mod->ms,
			     mod->m1,mod->m2,mod->n,mod->nrgp,mod->nrsp,
			     mod->nflt,mod->nslip,
			     (a+1), mod->a,&mod->d,mod->g,
			     &mod->gf,&mod->imat,mod->nrb,
			     mod->block,mod->nsnf,&mod->fault,
			     fit_rigid,invert_for_ld,
			     invert_for_cfac,mod->gx,mod->sx,
			     mod->stress_depths,fit_projection,
			     damp_nslip,mod->nfdamp,
			     mod->nxdamp,mod->xdamp,mod);
#ifdef DEBUG
  fprintf(stderr,"mrqcof: assemble dyda matrix OK\n");
#endif 
  /* 

     assemble alpha matrix and calculate the chi2 misfit

  */
  //
  // compute velocity misfit and assemble dyda
  //
  
  chisq[1] = 0.0;
  for (i=1;i <= mod->mgd;i++) {
    sig2i = 1.0/(sig[i]*sig[i]); 
    dy = y[i] - ymod[i];
    /* assembly of dyda */
    for (j=0,l=1;l <= ma;l++) {
      if (ia[l]) {
	wt = dyda[(l-1)*ndata+i-1] * sig2i;
	for (j++,k=0,m=1;m <= l;m++)
	  if (ia[m]) 
	    alpha[j][++k] += wt * dyda[(m-1)*ndata+i-1];
	beta[j] += dy * wt;
      }
    }
    chisq[1] += dy*dy*sig2i;
  }
  
  /*
    
  scale the stresses to the model
  
  when calling this function, make sure to shift 
  vectors up again since we called numrec style
  
  */
  if(mod->m2)
    rescale_observed_stresses((ymod+1),(y+1),(sig+1),
			      fit_damp_fac,stressrms,
			      no_stress_amp_scale,mod,
			      TRUE,TRUE);
  //
  // stresses
  //
  chisq[2] = 0.0;
  ilim = mod->mgd + mod->m2;
  for (i=mod->mgd + 1;i <= ilim;i++) {
    sig2i = 1.0 / (sig[i]*sig[i]); 
    dy = y[i] - ymod[i];
    for (j=0,l=1;l <= ma;l++) {
      if (ia[l]) {/* incorporate beta here */
	wt = dyda[(l-1)*ndata+i-1] * sig2i * fit_beta;
	for (j++,k=0,m=1;m <= l;m++)
	  if (ia[m]) 
	    alpha[j][++k] += wt * dyda[(m-1)*ndata+i-1];
	beta[j] += dy * wt;
      }
    }
    chisq[2] += dy*dy*sig2i;
  }
  /*
    reassign original stress amplitudes, include the shift also
    here for v and sigv
  */
  /*
    if(mod->m2)
    rescale_observed_stresses((ymod+1),(y+1),(sig+1),
			      fit_damp_fac,stressrms,
			      no_stress_amp_scale,mod,
			      TRUE,FALSE);
  */
  //
  // weighted chi2
  //
  if(mod->m2 == 0)		/* no stresses */
    chisq[0] = chisq[1];
  else
    chisq[0] = (chisq[1] + fit_beta * chisq[2])/(1.0 + fit_beta);

  if(mod->nfdamp){
    /* 
       
    damping of normal fault motion, add this to chi2
    
    */
  }
  if(0){
    /* 
       for testing purposes, compare the local chi2 and the 
       one that is normally used. tchi2 and chi2 should be the 
       same!
    
    */
    fprintf(stderr,"mrqcof: |v|: %11g |vmod|: %11g |sig|: %11g |s|: %11g |x|: %11g |un|: %11g\n",
	    norm((y+1),mod->mgd),norm((ymod+1),mod->mgd),
	    norm((sig+1),mod->mgd),norm((ymod+1+mod->mgd),mod->m2),
	    norm((a+1),mod->n),
	    (mod->nfdamp)?(norm((ymod+1+mod->mgd+mod->m2),mod->nfdamp)):(0.0));

    tchi2[0] = block_chi_square((ymod+1),(y+1),(sig+1),mod->mgd,
				mod->m2,fit_beta,(tchi2+1),(tchi2+2));
    fprintf(stderr,"mrqcof: m: %12g %12g %12g   n: %12g %12g %12g\n",
	    tchi2[0],tchi2[1],tchi2[2],chisq[0],chisq[1],chisq[2]);
  }
  //
  //
  for (j=2;j<=mfit;j++)
    for (k=1;k<j;k++) 
      alpha[k][j] = alpha[j][k];
  /* 

  penalty section if we are inverting for locking depths and such

  */
  i = mod->n+1;
  if(invert_for_ld){		/* 
				   locking depths, make sure
				   between 1.25 and 30 
				   (1.25 because of the dxdy
				   routine)
				*/
    /* 
       check the locking depths
    */
    for(penalty = 0.0,j=0;j < mod->nflt;j++,i++){	
      if(a[i] < 1.25)
	penalty += 1.25-a[i];
      if(a[i] > 30.0)
	penalty += (a[i]-30.0);
    }
    chisq[0] += penalty * 30000.0;
  }
  if(invert_for_cfac){
    /* 
       check the cfactors, if they are outside 0..1, add penalty
    */
    for(penalty = 0.0,j=0;j < mod->nflt;i++,j++){	
      if(a[i] < 0)
	penalty += -a[i];
      if(a[i] > 1.0)
	penalty += (a[i]-1.0);
    }
    chisq[0] += penalty * 10000.0;
  }
  if(constrain_slip_direction){
    /* 
       
       some faults have to slip in one direction, not the other
       add to pentalty
    */
    my_vecalloc(&vslip, mod->nsnf,"mrqcof: vslip");
    /* compute slip solution */
    calc_Ax_ftn(mod->gf,mod->nsnf,mod->n,(a+1),vslip);
    for(i=0;i<mod->nflt;i++){
      for(j=0;j<mod->nslip;j++)
	if(mod->fault[i].sc[j] != 0){
	  dy = vslip[i*mod->nslip+j] * 
	    (COMP_PRECISION)mod->fault[i].sc[j];
	  /* 
	     if the slip is in the other direction than
	     the constraint, add it to the penalty
	  */
	  if(dy < 0)
	    chisq[0] -= dy;
	}
    } /* end loop through all faults */
    free(vslip);
  }
  free(ymod);
  free(dyda);
}


void nrerror(char *error_text)
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

COMP_PRECISION **matrix(long nrl,long nrh,long ncl,long nch)
/* allocate a COMP_PRECISION matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  COMP_PRECISION **m;
#ifdef MEM_ALLOC_DEBUG
  fprintf(stderr,"matrix: newly allocating %.4f MB of memory\n",
	  (float)(sizeof(COMP_PRECISION)*(size_t)((nrow+NUMREC_NR_END)*(ncol+NUMREC_NR_END)))/
	  (float)ONE_MEGABYTE);
#endif
  /* allocate pointers to rows */
  m=(COMP_PRECISION **) malloc((unsigned int)((nrow+NUMREC_NR_END)*sizeof(COMP_PRECISION*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NUMREC_NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(COMP_PRECISION *) malloc((unsigned int)((nrow*ncol+NUMREC_NR_END)*sizeof(COMP_PRECISION)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NUMREC_NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}



void free_matrix(COMP_PRECISION **m,long nrl,long nrh,
		 long ncl,long nch)
{
  free((NUMREC_FREE_ARG) (m[nrl]+ncl-NUMREC_NR_END));
  free((NUMREC_FREE_ARG) (m+nrl-NUMREC_NR_END));
}


