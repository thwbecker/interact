/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: solve_mode_dependend.c,v 1.28 2003/01/07 02:52:29 tbecker Exp $
*/
#include "interact.h"
/*

  this source file holds all routines that depend on the way that interaction
  coefficients are calculated. to speed things up, there are always
  four different subroutines so that we can avoid numerous if-statement calls

  COMP_MODE_1: we have I matrix in memory and assemble it once
  COMP_MODE_2: we need to read the I matrix from file each time
  COMP_MODE_3: we will calculate the coefficients right now and forget them after
  COMP_MODE_4: we hold the I matrix in sparse matrix storage

*/
#ifdef USE_SLATEC_NNLS

/*

  assemble the A' matrix for mixed least squares and non-negative
  constraint solutions that depend on the slatec algorithm least
  square solutions

*/
#ifdef COMP_MODE_1
void assemble_ap_matrix_1(A_MATRIX_PREC *a,int naflt,int naflt_con,
			  my_boolean *sma,my_boolean *sma_con,
			  int nreq,int nreq_con,
			  int *nameaf,int *nameaf_con,
			  struct flt *fault,struct med *medium)
#elif defined  COMP_MODE_2
void assemble_ap_matrix_2(A_MATRIX_PREC *a,int naflt,int naflt_con,
			  my_boolean *sma,my_boolean *sma_con,
			  int nreq,int nreq_con,
			  int *nameaf,int *nameaf_con,
			  struct flt *fault,struct med *medium)
#elif defined COMP_MODE_3
void assemble_ap_matrix_3(A_MATRIX_PREC *a,int naflt,int naflt_con,
			  my_boolean *sma,my_boolean *sma_con,
			  int nreq,int nreq_con,
			  int *nameaf,int *nameaf_con,
			  struct flt *fault,struct med *medium)
#else
void assemble_ap_matrix_4(A_MATRIX_PREC *a,int naflt,int naflt_con,
			  my_boolean *sma,my_boolean *sma_con,
			  int nreq,int nreq_con,
			  int *nameaf,int *nameaf_con,
			  struct flt *fault,struct med *medium)
#endif
{
  /*
    see routine below for the logic of the stress assignments
  */
  long int i,j,k,l,eqc1,eqc2,m,namef1tmp,namef2tmp,allflt,ip1,ip2,eqc2m;
  A_MATRIX_PREC cf;
  my_boolean sma1tmp,sma2tmp;
#ifdef DEBUG
  long int n;
#endif
#ifndef COMP_MODE_1
  A_MATRIX_PREC itmp;
#endif
#ifdef COMP_MODE_3
  int iret;
  fprintf(stderr,"assemble_ap_matrix_3: assembling matrix and calculating interaction coefficients\n");
#endif
  // dimensions of A matrix (without b column)
  m = nreq + nreq_con;
#ifdef DEBUG
  n=m;
#endif
  allflt=naflt+naflt_con;
  // begin slow loop
  eqc1=0;
  for(ip1=i=0;i < allflt;i++,ip1+=3){
    for(j=0;j<3;j++){
      if(i < naflt)// get right activation mode 
	sma1tmp=sma[ip1+j];
      else
	sma1tmp=sma_con[(i-naflt)*3+j];
      if(sma1tmp){
#ifdef DEBUG
	if(eqc1 >= m){
	  fprintf(stderr,"assemble_ap_matrix: i (slow) index out of bounds %i, vs. m: %i\n",
		  (int)eqc1,(int)m);
	  exit(-1);
	}
#endif
	if(i < naflt)// get appropriate active fault name 
	  namef1tmp=nameaf[i];
	else
	  namef1tmp=nameaf_con[i-naflt];
	// assign coulomb correction for strike and dip stress directions
	// on fault namef1tmp
	cf = (j==NORMAL)?(0.0):(fault[namef1tmp].cf[j]);
	// begin fast loop
	eqc2 = eqc2m = 0;
	for(ip2=k=0;k < allflt;k++,ip2+=3){
	  for(l=0;l<3;l++){
	    if(k<naflt)
	      sma2tmp = sma[ip2+l];
	    else
	      sma2tmp = sma_con[(k-naflt)*3+l];
	    if(sma2tmp){
	      if(k<naflt)
		namef2tmp=nameaf[k];
	      else
		namef2tmp=nameaf_con[k-naflt];
#ifdef DEBUG
	      if(eqc2 >= n){
		fprintf(stderr,"assemble_ap_matrix: j (fast) index out of bounds %i, vs. n: %i\n",
			(int)eqc2,(int)n);
		exit(-1);
	      }
#endif
#ifdef COMP_MODE_1
	      *(a+eqc2m+eqc1) = (A_MATRIX_PREC)
		ICIM(medium->i,namef1tmp,namef2tmp,l,j);
	      if(cf != 0.0)
		*(a+eqc2m+eqc1) += (A_MATRIX_PREC)
		  ICIM(medium->i,namef1tmp,namef2tmp,l,NORMAL) * cf;
#elif defined COMP_MODE_2
	      *(a+eqc2m+eqc1) = (A_MATRIX_PREC)
		ic_from_file(namef1tmp,namef2tmp,l,j,medium);
	      if(cf != 0.0){
		itmp=(A_MATRIX_PREC)ic_from_file(namef1tmp,namef2tmp,l,NORMAL,medium);
		*(a+eqc2m+eqc1) +=  (A_MATRIX_PREC) itmp * cf;
	      }
#elif defined COMP_MODE_3
	      // need to calculate interaction coefficient right here
	      *(a+eqc2m+eqc1) = (A_MATRIX_PREC)
		interaction_coefficient(namef1tmp,namef2tmp,l,j,fault,&iret);
	      if(cf != 0.0){
		itmp=(A_MATRIX_PREC)interaction_coefficient(namef1tmp,namef2tmp,l,
							    NORMAL,fault,&iret);
		if(iret){
		  fprintf(stderr,"assemble_ap_matrix_3: WARNING: encountered iret: i/j/k/l: %i/%i/%i/%i\n",
			  (int)namef1tmp,(int)namef2tmp,(int)l,(int)j);
		  itmp=0.0;
		}
		*(a+eqc2m+eqc1) +=  (A_MATRIX_PREC) itmp * cf;
	      }
#else
	      // have I matrix in numerical recipes sparse matrix 
	      // storage
	      *(a+eqc2m+eqc1) = (A_MATRIX_PREC)
		get_nrs_sparse_el(POSII(namef2tmp,l),POSIJ(namef1tmp,j),
				  medium->is1,medium->val);
	      if(cf != 0.0){
		itmp=(A_MATRIX_PREC)
		  get_nrs_sparse_el(POSII(namef2tmp,l),POSIJ(namef1tmp,NORMAL),
				    medium->is1,medium->val);
		*(a+eqc2m+eqc1) +=  (A_MATRIX_PREC) itmp * cf;
	      }
#endif
	      /* reset to zero if there are no interactions between
		 faults wanted */
	      if(medium->no_interactions)
		if(fault[i].group != fault[k].group)
		  *(a+eqc2m+eqc1) = 0.0;
#ifdef SUPER_DEBUG
	      fprintf(stderr,"assemble_ap_matrix: i:%3i j:%i(%i) %e\n",
		      (int)eqc2,(int)eqc1,(int)m,*(a+eqc2m+eqc1));
#endif
	      eqc2++;
	      eqc2m += m;
	    }
	  }
	}
	eqc1++; 
      }
    }
  }
}
#endif

/*

  assemble the A matrix for plain NNLS or least square solutions

*/
#ifdef COMP_MODE_1
void assemble_a_matrix_1(A_MATRIX_PREC *a,int naflt,
			 my_boolean *sma,int nreq,int *nameaf,
			 struct flt *fault,struct med *medium)
#elif defined  COMP_MODE_2
void assemble_a_matrix_2(A_MATRIX_PREC *a,int naflt,
			 my_boolean *sma,int nreq,int *nameaf,
			 struct flt *fault,struct med *medium)
     
#elif defined COMP_MODE_3
void assemble_a_matrix_3(A_MATRIX_PREC *a,int naflt,
			 my_boolean *sma,int nreq,int *nameaf,
			 struct flt *fault,struct med *medium)
#else
void assemble_a_matrix_4(A_MATRIX_PREC *a,int naflt,
			 my_boolean *sma,int nreq,int *nameaf,
			 struct flt *fault,struct med *medium)
#endif
{
  long int i,j,k,l,eqc1,eqc2,eqc2nreq,ip1,ip2;
  A_MATRIX_PREC cf;
#ifndef COMP_MODE_1
  A_MATRIX_PREC itmp;
#endif
#ifdef COMP_MODE_3
  int iret,istep;
  istep = naflt/10+1;
#endif

  /*
    COMMENT: the flipped sign for constrained
    faults was determined in interact.c

    nameaf[i]: observing fault
    j: stress component on observing fault
    nameaf[k]: slipping fault
    l: slip mode on slipping fault

  */
  for(eqc1=ip1=i=0;i < naflt;i++,ip1+=3){
    for(j=0;j < 3;j++){
      if(sma[ip1+j]){
	// normal correction for Coulomb?
	cf = (j==NORMAL)?(0.0):fault[nameaf[i]].cf[j];
	for(eqc2=eqc2nreq=ip2=k=0;k < naflt;k++,ip2+=3){
	  for(l=0;l < 3;l++){
	    if(sma[ip2+l]){// flip around matrix ordering
#ifdef SUPER_DUPER_DEBUG
	      if(cf != 0.0)
		fprintf(stderr,"assemble_a_matrix: i: %i j: %i k: %i l: %i cf: %g\n",
			nameaf[i],nameaf[k],l,j,cf);
#endif
#ifdef COMP_MODE_1
	      // have I matrix in memory
	      *(a+eqc2nreq+eqc1) = (A_MATRIX_PREC)
		ICIM(medium->i,nameaf[i],nameaf[k],l,j);
	      if(cf != 0.0)
		*(a+eqc2nreq+eqc1) += (A_MATRIX_PREC)
		  ICIM(medium->i,nameaf[i],nameaf[k],l,NORMAL) * cf;
#elif defined COMP_MODE_2
	      // need to read it from file
	      *(a+eqc2nreq+eqc1) = (A_MATRIX_PREC)
		ic_from_file(nameaf[i],nameaf[k],l,j,medium);
	      if(cf != 0.0){// correct 
		itmp=(A_MATRIX_PREC)ic_from_file(nameaf[i],nameaf[k],l,NORMAL,medium);
		*(a+eqc2nreq+eqc1) +=  (A_MATRIX_PREC) itmp * cf;
	      }
#elif defined COMP_MODE_3
	      //
	      // calculate interaction coefficients right now
	      //
	      *(a+eqc2nreq+eqc1) = (A_MATRIX_PREC)
		interaction_coefficient(nameaf[i],nameaf[k],l,j,fault,&iret);
	       if(cf != 0.0){	/* coulomb addition */
		 itmp=(A_MATRIX_PREC)interaction_coefficient(nameaf[i],nameaf[k],l,NORMAL,fault,&iret);
		 if(iret){
		   fprintf(stderr,"assemble_a_matrix_3: WARNING: encountered iret: i/j/k/l: %i/%i/%i/%i\n",
			   nameaf[i],nameaf[k],(int)l,(int)j);
		   itmp=0.0;
		 }
		 *(a+eqc2nreq+eqc1) +=  (A_MATRIX_PREC) itmp * cf;
	       }
	       if(medium->debug)
		 if((i < 15)&&(j<15))
		   fprintf(stderr,"assemble_a_matrix_3: f1 %3i f2 %3i i1 %i i2 %i %12.3e\n",
			   (int)i,(int)k,(int)j,(int)l,*(a+eqc2nreq+eqc1));
	       
#else
	      // have I matrix in numerical recipes 
	      // sparse matrix storage
	      *(a+eqc2nreq+eqc1) = (A_MATRIX_PREC)
		get_nrs_sparse_el(POSII(nameaf[k],l),POSIJ(nameaf[i],j),
				  medium->is1,medium->val);
	      if(cf != 0.0){
		itmp=(A_MATRIX_PREC)get_nrs_sparse_el(POSII(nameaf[k],l),
						      POSIJ(nameaf[i],NORMAL),
						      medium->is1,medium->val);
		*(a+eqc2nreq+eqc1) +=  (A_MATRIX_PREC) itmp * cf;
	      }
#endif
	      if(medium->no_interactions)
		if(fault[i].group != fault[k].group)
		  *(a+eqc2nreq+eqc1) = 0.0;
	      eqc2++;
	      eqc2nreq += nreq;
	    }
	  }
	}
	eqc1++; 
      }
    }
#ifdef COMP_MODE_3
    if(!medium->debug)		/* progress report */
      if(i%istep == 0)
	fprintf(stderr,"assemble_a_matrix_3: assembling matrix and calculating interaction coefficients: %8.2f%%\r",
		(double)(i+1)/(double)naflt * 100.);
#endif
  }
#ifdef COMP_MODE_3
  fprintf(stderr,"\n");
#endif
}

/*

  add the stress change due to slip at certain fault r_flt
  to medium stresses

*/

#ifdef COMP_MODE_1
void add_quake_stress_1(my_boolean *sma,COMP_PRECISION *slip,
			int r_flt,struct flt *fault,
			struct med *medium)
#elif defined COMP_MODE_2
void add_quake_stress_2(my_boolean *sma,COMP_PRECISION *slip,
			int r_flt,struct flt *fault,
			struct med *medium)
#elif defined COMP_MODE_3
void add_quake_stress_3(my_boolean *sma,COMP_PRECISION *slip,
			int r_flt,struct flt *fault,
			struct med *medium)

#else
void add_quake_stress_4(my_boolean *sma,COMP_PRECISION *slip,
			int r_flt,struct flt *fault,
			struct med *medium)

#endif
{
  int i,j,k;
#ifdef COMP_MODE_3
  int iret;
  I_MATRIX_PREC iadbl;
#ifdef ALLOW_NON_3DQUAD_GEOM
#ifdef SUPER_DUPER_DEBUG
  
  COMP_PRECISION gstrike[3],gdip[3],gnormal[3],uglobal[3],tglobal[3];
#endif
#endif
#endif
  for(i=0;i < medium->nrflt;i++){/* loop through all flts */
   
	
#ifdef COMP_MODE_1
    /* 
       we have the interactions precomputed
    */
    for(j=0;j<3;j++){/* loop through all 
			possible slip dirs. */
      if(sma[j]){
	for(k=0;k<3;k++){/* calculate all 
			    affected components */
	  
	  fault[i].s[k] += 
	    ((COMP_PRECISION)ICIM(medium->i,i,r_flt,j,k)) * slip[j];
	  //	  fprintf(stderr,"rupnr: %5i recnr: %5i smode: %i stype: %i ic: %12g slip: %12g\n",
	  //	  r_flt,i,j,k,ICIM(medium->i,i,r_flt,j,k),slip[j]);
	}
      }
    }
#elif defined COMP_MODE_2
    /* 

       from file, really?!? 

    */
    for(j=0;j<3;j++){/* loop through all 
			possible slip dirs. */
      if(sma[j]){
	for(k=0;k<3;k++){
	  fault[i].s[k] += ((COMP_PRECISION)ic_from_file(i,r_flt,j,k,medium)) * slip[j];
	}
      }
    }
#elif defined COMP_MODE_3
    /* 
       
       calculate the effect of slip in terms of fault stress right now

    */
    if(medium->solver_mode == SPARSE_SOLVER){
      /* 

	 this does the one by one loop for some internal consistency I cannot remember

      */
      for(j=0;j<3;j++){/* loop through all 
			  possible slip dirs. */
	if(sma[j]){
	  for(k=0;k<3;k++){
	    iadbl = interaction_coefficient(i,r_flt,j,k,fault,&iret);
	    // if so, make sure that cutoff values are consistent
	    if(fabs(iadbl) < medium->i_mat_cutoff){
	      iadbl = 0.0;
	      iret = 1;// don't add small entries
	    }
	    if(!iret)
	      fault[i].s[k] += iadbl * slip[j];
	    fprintf(stderr,"sparse mode: rupnr: 1x1 %5i recnr: %5i smode: %i stype: %i ic: %12g slip: %12g iret: %i\n",r_flt,i,j,k,iadbl,slip[j],iret);
	  }
	}
      }
    }else{
      /* 
	 
	 compute stress change to due slip for all three components
	 
      */
#ifdef DEBUG
      for(j=0;j<3;j++)	/* check */
	if((fabs(slip[j]) > 0)&&(!sma[j])){
	  fprintf(stderr,"add_quake_stress_3: patch %i slip %g %g %g but sma %i %i %i\n",
		    i,slip[0],slip[1],slip[2],sma[0],sma[1],sma[2]);
	  exit(-1);
	}
#endif

#ifdef SUPER_DUPER_DEBUG
#ifdef ALLOW_NON_3DQUAD_GEOM
      if(fault[i].type == TRIANGULAR){
	fprintf(stderr,"add_quake_stress_3: rec %03i rup %03i slip %10.3e, %10.3e, %10.3e s/d %.2f %.2f",
		i,r_flt,slip[STRIKE],slip[DIP],slip[NORMAL],fault[i].strike,fault[i].dip);
	calc_global_strike_dip_from_local((fault+i),gstrike,gnormal,gdip);
	calc_global_slip_and_traction_from_local((fault+i),slip,slip,gstrike, gnormal,gdip,uglobal,tglobal,FALSE);
	fprintf(stderr," gslip %10.3e, %10.3e, %10.3e\n",uglobal[0],uglobal[2],uglobal[1]);
      }else{
	fprintf(stderr,"add_quake_stress_3: rec %03i rup %03i slip %10.3e, %10.3e, %10.3e\n",
		i,r_flt,slip[STRIKE],slip[DIP],slip[NORMAL]);
      }
#endif
#endif
      /* compute the effect of slip of r_flt on fault[i] and add to it's stress fault[i].s */
      eval_green_and_project_stress_to_fault(fault,i,r_flt,slip,fault[i].s);
    }
#else
    for(j=0;j < 3;j++){/* loop through all 
			possible slip dirs. */
      if(sma[j]){
	
	for(k=0;k<3;k++){
	  // numerical recipes sparse matrix scheme
	  fault[i].s[k] += 
	    ((COMP_PRECISION)get_nrs_sparse_el(POSII(r_flt,j),POSIJ(i,k),
					       medium->is1,medium->val))*slip[j];
	}
      }
    }
#endif
  } /* number of receiving fault loops */
}

/*
  
  calculate the actual coulomb stress changes due to 
  self-interaction and interaction and return TRUE if
  a positive feedback loop is detected

*/
#ifdef COMP_MODE_1
my_boolean check_coulomb_stress_feedback_1(int nrflt,
					int old_nrflts,
					struct flt *fault,
					struct med *medium,
					my_boolean call_out_loud,
					my_boolean bailout,
					int *evil_pair,
					COMP_PRECISION max_distance)
#elif defined COMP_MODE_2
my_boolean check_coulomb_stress_feedback_2(int nrflt,
					int old_nrflts,
					struct flt *fault,
					struct med *medium,
					my_boolean call_out_loud,
					my_boolean bailout,
					int *evil_pair,
					COMP_PRECISION max_distance)
#elif defined COMP_MODE_3
my_boolean check_coulomb_stress_feedback_3(int nrflt,
					int old_nrflts,
					struct flt *fault,
					struct med *medium,
					my_boolean call_out_loud,
					my_boolean bailout,
					int *evil_pair,
					COMP_PRECISION max_distance)
#else
my_boolean check_coulomb_stress_feedback_4(int nrflt,
					int old_nrflts,
					struct flt *fault,
					struct med *medium,
					my_boolean call_out_loud,
					my_boolean bailout,
					int *evil_pair,
					COMP_PRECISION max_distance)
#endif
{
  int i,j,k;
  MODE_TYPE mode;
  static MODE_TYPE slip_mode_arr[2]={STRIKE,DIP};
  static my_boolean first_call=TRUE;
  static COMP_PRECISION dist_sqr;
#ifdef COMP_MODE_3
  int iret;
#endif
  my_boolean hit;
  COMP_PRECISION cs_jj,cs_ij,cs_ii,cs_ji;
  if(first_call){
    fprintf(stderr,"check_coulomb_stress_feedback: using cohesion: %g, first fault has mu_s: %g\n",
	    medium->cohesion,fault[0].mu_s);
    first_call=FALSE;
    dist_sqr=SQUARE(max_distance);
  }

  evil_pair[0]=evil_pair[1]=-1;
  hit=FALSE;
  for(i=0;(i<nrflt-1)&&(!(bailout && hit));i++){
    for(j=old_nrflts;(j<nrflt)&&(!(bailout && hit));j++){
      if(i==j)continue;
      /*
	if the max distance flag is not set to -1,
	check for the fault's center point distance and
	reject those further away than max_distance
      */
      if(max_distance != -1.0)
	if(distance_squared_3d(fault[i].x,fault[j].x) > dist_sqr)
	  continue;
      for(k=0;(k<2)&&(!(bailout && hit));k++){
	mode=slip_mode_arr[k];
	/* 
	   check for Coulomb stress changes due to shear 
	   stress (strike or dip) and normal stress changes during slip of 
	   strike or dip type. the shear stress is calculated based
	   on the slipping mode, ie. either in strike or in dip
	   direction
	*/
#ifdef COMP_MODE_1 
	// in memory
	// calculate coulomb stress for self action
	cs_jj=coulomb_stress(fabs((COMP_PRECISION)ICIM(medium->i,j,j,mode,mode)),
			   fault[j].mu_s,(COMP_PRECISION)ICIM(medium->i,j,j,mode,NORMAL),medium->cohesion);
	// and effect on other fault
	cs_ij=coulomb_stress(fabs((COMP_PRECISION)ICIM(medium->i,i,j,mode,mode)),
			   fault[i].mu_s,(COMP_PRECISION)ICIM(medium->i,i,j,mode,NORMAL),medium->cohesion);
	// and turned around
	cs_ii=coulomb_stress(fabs((COMP_PRECISION)ICIM(medium->i,i,i,mode,mode)),
			   fault[i].mu_s,(COMP_PRECISION)ICIM(medium->i,i,i,mode,NORMAL),medium->cohesion);
	cs_ji=coulomb_stress(fabs((COMP_PRECISION)ICIM(medium->i,j,i,mode,mode)),
			   fault[j].mu_s,(COMP_PRECISION)ICIM(medium->i,j,i,mode,NORMAL),medium->cohesion);
#elif defined COMP_MODE_2
	// from file
	cs_jj=coulomb_stress(fabs((COMP_PRECISION)ic_from_file(j,j,mode,mode,medium)),
			   fault[j].mu_s,(COMP_PRECISION)ic_from_file(j,j,mode,NORMAL,medium),medium->cohesion);
	cs_ij=coulomb_stress(fabs((COMP_PRECISION)ic_from_file(i,j,mode,mode,medium)),
			   fault[i].mu_s,(COMP_PRECISION)ic_from_file(i,j,mode,NORMAL,medium),medium->cohesion);
	cs_ii=coulomb_stress(fabs((COMP_PRECISION)ic_from_file(i,i,mode,mode,medium)),
			   fault[i].mu_s,(COMP_PRECISION)ic_from_file(i,i,mode,NORMAL,medium),medium->cohesion);
	cs_ji=coulomb_stress(fabs((COMP_PRECISION)ic_from_file(j,i,mode,mode,medium)),
			   fault[j].mu_s,(COMP_PRECISION)ic_from_file(j,i,mode,NORMAL,medium),medium->cohesion);
#elif defined COMP_MODE_3
	// calculate now
	cs_jj=coulomb_stress(fabs((COMP_PRECISION)interaction_coefficient(j,j,mode,mode,fault,&iret)),
			   fault[j].mu_s,(COMP_PRECISION)interaction_coefficient(j,j,mode,NORMAL,fault,&iret),medium->cohesion);
	cs_ij=coulomb_stress(fabs((COMP_PRECISION)interaction_coefficient(i,j,mode,mode,fault,&iret)),
			   fault[i].mu_s,(COMP_PRECISION)interaction_coefficient(i,j,mode,NORMAL,fault,&iret),medium->cohesion);
	cs_ii=coulomb_stress(fabs((COMP_PRECISION)interaction_coefficient(i,i,mode,mode,fault,&iret)),
			   fault[i].mu_s,(COMP_PRECISION)interaction_coefficient(i,i,mode,NORMAL,fault,&iret),medium->cohesion);
	cs_ji=coulomb_stress(fabs((COMP_PRECISION)interaction_coefficient(j,i,mode,mode,fault,&iret)),
			   fault[j].mu_s,(COMP_PRECISION)interaction_coefficient(j,i,mode,NORMAL,fault,&iret),medium->cohesion);
#else
	// sparse storage
	cs_jj=coulomb_stress(fabs((COMP_PRECISION)get_nrs_sparse_el(POSII(j,mode),POSIJ(j,mode),medium->is1,medium->val)),
			     fault[j].mu_s,(COMP_PRECISION)get_nrs_sparse_el(POSII(j,mode),POSIJ(j,NORMAL),medium->is1,medium->val),
			     medium->cohesion);
	cs_ij=coulomb_stress(fabs((COMP_PRECISION)get_nrs_sparse_el(POSII(j,mode),POSIJ(i,mode),medium->is1,medium->val)),
			     fault[i].mu_s,(COMP_PRECISION)get_nrs_sparse_el(POSII(j,mode),POSIJ(i,NORMAL),medium->is1,medium->val),
			     medium->cohesion);
	cs_ii=coulomb_stress(fabs((COMP_PRECISION)get_nrs_sparse_el(POSII(i,mode),POSIJ(i,mode),medium->is1,medium->val)),
			     fault[i].mu_s,(COMP_PRECISION)get_nrs_sparse_el(POSII(i,mode),POSIJ(i,NORMAL),medium->is1,medium->val),
			   medium->cohesion);
	cs_ji=coulomb_stress(fabs((COMP_PRECISION)get_nrs_sparse_el(POSII(i,mode),POSIJ(j,mode),medium->is1,medium->val)),
			   fault[j].mu_s,(COMP_PRECISION)get_nrs_sparse_el(POSII(i,mode),POSIJ(j,NORMAL),medium->is1,medium->val),
			   medium->cohesion);
#endif
	if((cs_ij > cs_jj)&&(cs_ji > cs_ii)){// positive feedback
	  hit=TRUE;
	  evil_pair[0]=i;evil_pair[1]=j;
	  if(call_out_loud)
	    fprintf(stderr,
		    "check_coulomb_stress_feedback: two way: mode:%i, i:%5i j:%5i cii:%12.4e < cji:%12.4e cjj:%12.4e < cij:%12.4e\n",
		    k,i,j,cs_ii,cs_ji,cs_jj,cs_ij);
	}
#ifdef CHECK_CI_ONE_WAY
	//
	// self interaction smaller than effect of other fault
	// one way is enough, if patches are in different groups
	//
	if((fault[i].group != fault[j].group)&&
	   ((cs_ij > cs_jj)||(cs_ji > cs_ii))){
	  hit=TRUE;
	  evil_pair[0]=i;evil_pair[1]=j;
	  if(call_out_loud)
	    fprintf(stderr,
		    "check_coulomb_stress_feedback: one way: mode:%i, i:%5i j:%5i cii:%12.4e < cji:%12.4e cjj:%12.4e < cij:%12.4e\n",
		    k,i,j,cs_ii,cs_ji,cs_jj,cs_ij);
	}
#endif
      }
    }
  }
  return(hit);
}



