/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu


  all routines in here deal with the solution of interact 
  program specific routines


  $Id: solve.c,v 2.37 2011/01/09 02:02:43 becker Exp $

*/
#include "interact.h"

void solve(struct med *medium,struct flt *fault)
{
  A_MATRIX_PREC *a,wmax,*sv=NULL,*dummyp=NULL,wcutoff;
  char solver_name[STRLEN];
  size_t full_size,sparse_size;
  static int izero=0;
  FILE *aio;
  int i;
#ifdef USE_SLATEC_NNLS
  A_MATRIX_PREC *x,prgopt[1]={1.0},rnorm,*work,dummy;
  int n,m,*iwork,mode,l,nm,j;
  FILE *out1,*out2;
#endif
  char command_str[STRLEN];
  wcutoff = (A_MATRIX_PREC)medium->wcutoff;

  if(medium->debug){
    print_equations(medium->naflt,medium->sma,medium->nameaf,
		    medium->b,medium->nreq,"ucstr.",
		    fault);
    print_equations(medium->naflt_con,medium->sma_con,
		    medium->nameaf_con,medium->b_con,
		    medium->nreq_con," cstr.",fault);
  }
  /* 
     now target_stress_drop has the dropped stress values, 
     and sma indicates which mode should slide for every patch 
     
     assemble the A matrix to solve
     for the displacement increments x in 

     A x = b
     
     where b are the target stress drops
     A the nonzero interaction matrix with
     values for those faults and modes  that ruptured
     
     depending on the constrained (x>0) or free modes, use three different solvers
     
     solver                            free     constrained

     NNLS                                Y           Y
     LU, SVD, or other  Ax=b solver      Y           N
     Lawson-Hanson original NNLS         N           Y


  */
  if((medium->nreq) && (medium->nreq_con)){
    /*

      have to solve constrained and uncontrained equations
      mixed NNLS and least squares from NNLS package

    */
#ifndef USE_SLATEC_NNLS
    fprintf(stderr,"solve: for mixed contrained/uncontrained solution, have to use SLATEC\n");
    fprintf(stderr,"solve: library functions. Read the makefile, and compile with -lslatec_???.\n");
    exit(-1);
#else
    if(medium->debug)
      fprintf(stderr,"solve: using slatec nnls routine, uc: %i c: %i\n",
	      medium->nreq,medium->nreq_con);
    m = medium->nreq + medium->nreq_con;
    n = m;
    nm=n * m;
    // these for SLATEC routine
    //me=0;// no exact equations
    //ma=m;// least squares equations
    l=medium->nreq;// up from l+1 will be non neg
    // allocation
    a=(A_MATRIX_PREC *)calloc((n+1)*m,sizeof(A_MATRIX_PREC));
    x=(A_MATRIX_PREC *)calloc(n,sizeof(A_MATRIX_PREC));
    iwork=(int *)calloc((m+n),sizeof(int));
    iwork[0]=m+5*n;// size of work
    work=(A_MATRIX_PREC *)calloc(iwork[0],sizeof(A_MATRIX_PREC));
    if((!iwork) || (!work) || (!x) || (!a)){
      fprintf(stderr,"solve: slatec_nnls: memory allocation error, n: %i, m: %i\n",
	      n,m);
      exit(-1);
    }
    iwork[1]=m+n;// size of iwork
    /*
      assemble A part of the A' matrix 
    */
    if(medium->use_old_amat){
      fprintf(stderr,"solve: slatec_nnls: WARNING: using old A matrix!\n");
      read_a_matrix_from_file(a,m,n,A_MATRIX_FILE,"mixed A");
    }else
      assemble_ap_matrix(a,medium->naflt, medium->naflt_con,medium->sma,   
			 medium->sma_con,medium->nreq,  medium->nreq_con,
			 medium->nameaf,medium->nameaf_con,fault,medium);
    if(medium->save_amat){// write A in binary to file
      print_a_matrix_to_file(a,m,n,A_MATRIX_FILE,"mixed A");
    }
    // now b part of A'
    for(i=0;i < medium->nreq;i++)
      a[nm+i]=medium->b[i];
    for(j=0,i=medium->nreq;i < n;i++,j++)
      a[nm+i]=medium->b_con[j];
    if(medium->debug){
      fprintf(stderr,"slatec_nnls: attempting to solve %i times (%i+1) system\n",
	      m,n);
      fprintf(stderr,"slatec_nnls: writing to \"%s\" and \"%s\"\n",
	      DEBUG_A_MIX_MATRIX_ASCII_OUT,
	      DEBUG_B_MIX_VECTOR_ASCII_OUT);
      out1=myopen(DEBUG_A_MIX_MATRIX_ASCII_OUT,"w");
      out2=myopen(DEBUG_B_MIX_VECTOR_ASCII_OUT,"w");
      print_a_matrix(a,m,n,out1,&dummy,FALSE);
      print_b_vector((a+nm),m,out2,&dummy,FALSE);
      fclose(out1);fclose(out2);
#ifndef A_MATRIX_SINGLE_PREC
      if(sizeof(A_MATRIX_PREC)==4){
	fprintf(stderr,"solve: slatec_nnls: error, have to recompile with A_MATRIX_SINGLE_PREC set\n");
	exit(-1);
      }
#else
      if(sizeof(A_MATRIX_PREC)==8){
	fprintf(stderr,"solve: slatec_nnls: error, have to recompile without A_MATRIX_SINGLE_PREC set\n");
	exit(-1);
      }
#endif
    }
    //
    // solve using SLATEC routines
    //
#ifdef A_MATRIX_SINGLE_PREC
    // single prec version
    wnnls_( a, &m, &izero, &m, &n, &l, prgopt, x, &rnorm, 
	    &mode,iwork, work);
#else
    // double version
    dwnnls_(a, &m, &izero, &m, &n, &l, prgopt, x, &rnorm, 
	    &mode,iwork, work);
#endif
    free(work);free(iwork);free(a);
    //
    // assign x to xsol, this has to be a bit convoluted because of the 
    // bookkeeping scheme of the wnnls routine
    // 
    if(mode == 0){// successful return
      for(i=0;i<medium->nreq;i++)
	medium->xsol[i]=x[i];
      for(j=0,i=medium->nreq;i<n;i++,j++)
	medium->xsol_con[j]=x[i];
      free(x);
      if(medium->debug){
	fprintf(stderr,"slatec_nnls: rnorm: %e, writing to \"%s\"\n",rnorm,
		DEBUG_X_MIX_VECTOR_ASCII_OUT); 
	out1=myopen(DEBUG_X_MIX_VECTOR_ASCII_OUT,"w");
	for(i=0;i<n;i++)
	  fprintf(out1,"%g\n",(i<medium->nreq)?(medium->xsol[i]):
		  (medium->xsol_con[i-medium->nreq]));
	fclose(out1);
      }
    }else{
      fprintf(stderr,"solve: nnls_slatec: error: wnnls routine returned mode: %i\n",
	      mode);
      exit(-1);
    }
    /* end constrained solver */
#endif
  }else if(medium->nreq){
    /*
      
      SOLUTION OF UNCONSTRAINED ONLY SET OF EQUATIONS 
      
      using LU, SVD or other methods

    */
    if((a=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*SQUARE(medium->nreq)))
       ==NULL){
      fprintf(stderr,"solve: memory error for unconstrained A: n: %i by %i dsz: %i tsize: %g MB\n",
	      medium->nreq,medium->nreq,(int)sizeof(A_MATRIX_PREC),
	      (double)sizeof(A_MATRIX_PREC)
	      *SQUARE((double)medium->nreq)/(double)ONE_MEGABYTE);
	exit(-1);
    }	
    if(medium->debug)
      fprintf(stderr,"solve: unconstrained part A is %5i by %i\n",
	      medium->nreq,medium->nreq);
    if(medium->use_old_amat){
      fprintf(stderr,"solve: unconstrained: WARNING: using old A matrix!\n");
      read_a_matrix_from_file(a,medium->nreq,medium->nreq,A_MATRIX_FILE,
			      "unconstrained A");
    }else
      assemble_a_matrix(a,medium->naflt,medium->sma,medium->nreq,medium->nameaf,
			fault,medium);
    if(medium->save_amat){// write A in binary to file
      print_a_matrix_to_file(a,medium->nreq,medium->nreq,
			     A_MATRIX_FILE,"unconstrained A");
      print_vector_file(medium->b,medium->nreq,B_VECTOR_FILE,"b vector");
    }
    if(medium->debug)
      if((medium->op_mode == ONE_STEP_CALCULATION)&&(!medium->save_amat))
	// output of binary A for debugging
	print_a_matrix_to_file(a,medium->nreq,medium->nreq,A_MATRIX_FILE,
			       "unconstrained A");
    if((medium->op_mode == ONE_STEP_CALCULATION)||(medium->debug)){
      // output of solver method only for one-step calculation
      switch(medium->solver_mode){
      case LU_SOLVER:
	sprintf(solver_name,"LAPACK LU");
	break;
      case SVD_SOLVER:
#ifdef USE_NUMREC_SVD
	sprintf(solver_name,"NUMREC SVD (wmax: %.3e)",medium->wcutoff);
#else
#ifdef USE_LAPACK_DAC
	sprintf(solver_name,"LAPACK SVD (gelsd, wmax: %.3e)",medium->wcutoff);
#else
	sprintf(solver_name,"LAPACK SVD (gelss, wmax: %.3e)",medium->wcutoff);
#endif
#endif
	break;
      case SPARSE_SOLVER:
	sprintf(solver_name,"SuperLU sparse LU solver");
	break;
      default:
	fprintf(stderr,"solve: error, don't know solver mode %i\n",
		medium->solver_mode);
	exit(-1);
      }
      fprintf(stderr,"solve: unconstrained system, time: %11g, n: %5i, solver mode: %s\n",
	      medium->time,medium->nreq,solver_name);
    }
    //
    // solve A . x = b 
    //
    switch(medium->solver_mode){
    case LU_SOLVER:// should be pretty fast 
      lu_driver(a,medium->xsol,medium->b,
		medium->nreq,medium->nreq,medium);
      free(a);
      break;
    case SVD_SOLVER:// gives least-squares solution in case A is singular

#ifdef USE_NUMREC_SVD
      svd_driver_numrec(a,medium->xsol,medium->b,medium->nreq,
			medium->nreq,&wcutoff,
			medium->op_mode,&wmax,FALSE,
			&sv,&dummyp,izero,
			FALSE,0,FALSE,TRUE,medium->debug);
      free(dummyp);
#else
      svd_driver_lapack(a,medium->xsol,medium->b,medium->nreq,
			medium->nreq,&wcutoff,
			medium->op_mode,&wmax,&sv,izero,TRUE,medium->debug);
#endif
      if(medium->save_amat){
	fprintf(stderr,"solve: saving singular values to %s\n",SV_FILE);
	aio=myopen(SV_FILE,"w");
	for(i=0;i<medium->nreq;i++)
	  fprintf(aio,"%20.7e %20.7e\n",sv[i],sv[i]/sv[0]);
	fclose(aio);
      }
      free(a);free(sv);
      break;
    case SPARSE_SOLVER:
#ifndef USE_SUPERLU
      fprintf(stderr,"solve: unconstrained: the sparse SuperLU driver is only available when SuperLU\n");
      fprintf(stderr,"solve: unconstrained: support was compiled in. Recompile with -DUSE_SUPERLU set.\n");
      exit(-1);
#endif
      //
      // the sparse solver is meant only for the one-step 
      // calculation
      if(medium->op_mode != ONE_STEP_CALCULATION){
	fprintf(stderr,"solve: unconstrained: the sparse Ax=b solver is no good for loading simulations\n");
	fprintf(stderr,"solve: unconstrained: choose SVD or LU, instead\n");
	exit(-1);
      }
      // print A matrix to file if not printed already
      if(!medium->save_amat)
	print_a_matrix_to_file(a,medium->nreq,medium->nreq,
			       A_MATRIX_FILE,"unconstrained A");
      // free A and open file pointer
      free(a);aio=myopen(A_MATRIX_FILE,"r");
      full_size =  sizeof(A_MATRIX_PREC) * SQUARE(medium->nreq);
      fprintf(stderr,"solve: unconstrained converting A (full size: %g MB) to sparse, coff: %g\n",
	      (double)full_size/ONE_MEGABYTE,medium->i_mat_cutoff);
      // create a CCS sparse representation 
      sparse_size = create_ccs_sparse_from_file(medium->nreq,
						medium->i_mat_cutoff,
						&medium->is1,&medium->is2,
						&medium->val,aio);
      fclose(aio);
      if(medium->debug)
	if(!medium->save_amat)// if A is not to be saved, remove
	  sprintf(command_str,"rm %s",A_MATRIX_FILE);system(command_str);
      fprintf(stderr,"solve: unconstrained: conversion complete (sparse size: %g MB) gain: %g\n",
	      (double)sparse_size/ONE_MEGABYTE,
	      (COMP_PRECISION)full_size/(COMP_PRECISION)sparse_size);
      //
      // solve CCS sparse system
      //
      sparse_driver(medium->is1, medium->is2,medium->val, medium->xsol, medium->b,
		    medium->nreq,medium);
      // free sparse matrix pointers
      free(medium->is1);free(medium->is2);free(medium->val);
      break;
    default:
      fprintf(stderr,"solve: unconstrained: solver mode %i undefined\n",
	      medium->solver_mode);
      exit(-1);
    }
    // write x to stderr
    //print_b_vector(medium->xsol,medium->nreq,stderr,&dummy,FALSE);
  }else if(medium->nreq_con){
    /*

      SOLUTION OF CONSTRAINED SET OF EQUATIONS 

      using original lawson - hanson nnls driver

    */
    if((a=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*
				  SQUARE(medium->nreq_con)))==NULL){
      fprintf(stderr,"solve: memory error for constrained A, %g MB\n",
	      (double)sizeof(A_MATRIX_PREC)*
	      (double)SQUARE(medium->nreq_con)/ONE_MEGABYTE);
      exit(-1);
    }
    if(medium->debug)
      fprintf(stderr,"solve:   constrained part A is %5i by %i\n",
	      medium->nreq_con,medium->nreq_con);
    if(medium->use_old_amat){
      fprintf(stderr,"solve: constrained: WARNING: using old A matrix!\n");
      read_a_matrix_from_file(a,medium->nreq_con,
			      medium->nreq_con,A_MATRIX_FILE,
			      "constrained A");
    }else
      assemble_a_matrix(a,medium->naflt_con,medium->sma_con,
			medium->nreq_con,medium->nameaf_con,
			fault,medium);
    if(medium->save_amat)// write A in binary to file
      print_a_matrix_to_file(a,medium->nreq_con,medium->nreq_con,
			     A_MATRIX_FILE,"constrained A");
    if(medium->debug)
      if((medium->op_mode == ONE_STEP_CALCULATION)&&(!medium->save_amat))
	print_a_matrix_to_file(a,medium->nreq_con,medium->nreq_con,
			       A_CON_MATRIX_FILE,"  constrained A");
    if(medium->op_mode == ONE_STEP_CALCULATION)
      fprintf(stderr,"solve: solving constrained system using NNLS, n: %i\n",
	      medium->nreq_con);
    /* call NNLS */
    nnls_driver_i(a,medium->xsol_con,medium->b_con,
		  medium->nreq_con,medium->nreq_con);
    free(a);
  }
  /* 
     now the b vector, target_stress_drops, 
     was overwritten with the x vector, the slip 
     increments  
  */
  if(medium->debug)
    fprintf(stderr,"solve: solvers done\n");
}

/* 
   reset all eqn systems and counters 
*/

void init_equation_system(struct med *medium,struct flt *fault)
{
  int i;
  static my_boolean initialized=FALSE;
  // reset the active equation counters for normal and
  medium->naflt=medium->nreq=0;
  // positivity constraint equations
  medium->naflt_con=medium->nreq_con=0;
  if(!initialized){// first time allocation
    initialized=TRUE;
  }else{
    // were allocated already, need to free them first
    // why did the reallocate version cease to work at 
    // some point?
    free(medium->sma);free(medium->sma_con);
    free(medium->nameaf);free(medium->nameaf_con);
    free(medium->b);free(medium->b_con);
    free(medium->xsol);free(medium->xsol_con);
  }
  medium->sma=(my_boolean *)malloc(3*sizeof(my_boolean));
  medium->sma_con=(my_boolean *)malloc(3*sizeof(my_boolean));
  for(i=0;i<3;i++){
    medium->sma[i]=INACTIVE;
    medium->sma_con[i]=INACTIVE;
  }
  medium->nameaf=(int *)calloc(1,sizeof(int));
  medium->nameaf_con =(int *)calloc(1,sizeof(int)); 
  medium->b=(A_MATRIX_PREC *)calloc(1,sizeof(A_MATRIX_PREC));
  medium->xsol=(A_MATRIX_PREC *)calloc(1,sizeof(A_MATRIX_PREC));
  medium->b_con=(A_MATRIX_PREC *)calloc(1,sizeof(A_MATRIX_PREC));
  medium->xsol_con=(A_MATRIX_PREC *)calloc(1,sizeof(A_MATRIX_PREC));
  if(medium->debug)
    if((!medium->sma)||(!medium->sma_con)||(!medium->nameaf)||
       (!medium->nameaf_con)||(!medium->b)||(!medium->b_con)||
       (!medium->xsol)||!(medium->xsol_con))
      MEMERROR("solve: init_equation_system:");
  // reset all fault activation flags
  for(i=0;i<medium->nrflt;i++)
    fault[i].active=FALSE;
  // reset all group activation flags
  medium->nr_active_groups=0;
  for(i=0;i<medium->nrgrp;i++)
    medium->fault_group[i].active=FALSE;
}

/* 

   add fault to active fault list and increment counter 

*/
void add_to_active_fault_list(int aflt,int **al,int *naf,my_boolean **sma)
{
  int ip;
  /* assign active fault name */
  (*al)[*naf]=aflt;
#ifdef DEBUG
  // check codes which have been assigned already
  if((*(*sma+ *naf*3)>3)||(*(*sma+ *naf*3+1)>3)||(*(*sma+ *naf*3+2)>3)){
    fprintf(stderr,"add_to_active_fault_list: fault %i: code screw up: sma 1/2/3: %i %i %i\n",
	    *naf,*(*sma+ *naf*3),*(*sma+ *naf*3+1),*(*sma+ *naf*3+2));
    exit(-1);
  }
#endif
  // increment number of active faults, now *naf is new index
  *naf += 1;
  /* grow slide mode array */
  if((*sma=(my_boolean *)realloc(*sma,sizeof(my_boolean)*3*
			      (*naf+1)))==NULL)
    MEMERROR("add_to_active_fault_list: 1");
  // zero out new codes
  ip = *naf * 3;
  *(*sma + ip + STRIKE)=INACTIVE;
  *(*sma + ip + DIP)=   INACTIVE;
  *(*sma + ip + NORMAL)=INACTIVE;
  /* grow active fault array */
  if((*al=(int *)realloc(*al,sizeof(int)*(*naf+1)))==NULL)
    MEMERROR("add_to_active_fault_list: 2");
}

/* 

   add to x value to right hand side (b), grow b, grow x, 
   and increment counter 

*/
void add_to_right_hand_side(COMP_PRECISION bval,A_MATRIX_PREC **b,
			    COMP_PRECISION **xsol, int *nreq)
{
  size_t i;
  (*b)[*nreq]=(A_MATRIX_PREC)bval;
  *nreq += 1;//increment equation counter
  // resizing size
  i= (*nreq + 1)*sizeof(A_MATRIX_PREC);
  // grow b
  if((*b=(A_MATRIX_PREC *)realloc(*b,i))==NULL)
    MEMERROR("add_to_right_hand_side:");
  // grow x
  if((*xsol=(A_MATRIX_PREC *)realloc(*xsol,i))==NULL)
    MEMERROR("add_to_right_hand_side:");
#ifdef SUPER_DUPER_DEBUG
  fprintf(stderr,"add_to_right_hand_side: adding %20.10e as eq. %5i\n",
	  bval,*nreq);
#endif
}

/* 

   add the x vector times sfac to the fault system
   
   if mark quake is set, will add to plot files, total moment and the
   like
 
   if add_quake_stress is not set, will not calculate the stress
   change (since
   
   it was done before) and not add to the fault slip counter

*/
void add_solution(int naflt,my_boolean *sma,A_MATRIX_PREC *xsol,int *nameaf,
		  struct med *medium,struct flt *fault,
		  my_boolean mark_quake,my_boolean add_quake_stress,
		  COMP_PRECISION sfac)
{
  int i,j,eqc,ip;
  COMP_PRECISION slip[3];
  for(eqc=i=ip=0;i<naflt;i++,ip += 3){
    // create a slip vector  by looping through all components
    for(j=0;j<3;j++){
      slip[j]=0.0;
      // if we use all modes, change with direction j
      // else, use only mode[0]
      switch(fault[nameaf[i]].mode[(medium->nr_flt_mode==3)?(j):(0)]){
	/* give a list of constrained equations that have
	   to be taken negative */
      case COULOMB_STRIKE_SLIP_RIGHTLATERAL:
      case COULOMB_DIP_SLIP_DOWNWARD:
      case STRIKE_SLIP_RIGHTLATERAL:
      case DIP_SLIP_DOWNWARD:
      case NORMAL_SLIP_INWARD:{
	/* for these, we will have to subtract the solution
	   as non-negative should actually go in the other
	   direction */
	if(sma[ip + j]){
	  slip[j]= -((COMP_PRECISION)xsol[eqc] * sfac);
	  eqc++;
	}
	break;
      }
      default:{
	if(sma[ip + j]){
	  slip[j]=   ((COMP_PRECISION)xsol[eqc] * sfac);
	  eqc++;
	}
	break;
      }}
    }
    /* assign slip and stress increments to all faults */
    quake((sma+ip),slip,nameaf[i],fault,medium,add_quake_stress,mark_quake);
  }
}

/*
  assemble the A matrix
*/

void assemble_a_matrix(A_MATRIX_PREC *a,int naflt,my_boolean *sma,
		       int nreq,int *nameaf,
		       struct flt *fault,struct med *medium)
{
  int imatmode;
  // determine way to get interaction coefficient
  imatmode=select_i_coeff_calc_mode(medium);
  switch(imatmode){
  case I_MAT_IN_MEMORY:{
    assemble_a_matrix_1(a,naflt,sma,nreq,nameaf,fault,medium);
    break;
  }
  case I_MAT_ON_FILE:{
    assemble_a_matrix_2(a,naflt,sma,nreq,nameaf,fault,medium);
    break;
  }
  case CALC_I_COEFF_NOW:{
    assemble_a_matrix_3(a,naflt,sma,nreq,nameaf,fault,medium);
    break;
  }
  case SPARSE_I_MAT:{
    assemble_a_matrix_4(a,naflt,sma,nreq,nameaf,fault,medium);
    break;
  }
  default:{
    fprintf(stderr,"assemble_a_matrix: can not deal with I matrix code %i\n",
	    imatmode);
    exit(-1);
  }}
}


#ifdef USE_SLATEC_NNLS

// matrix for combines constrained and unconstrained 

void assemble_ap_matrix(A_MATRIX_PREC *a,int naflt,int naflt_con,
			my_boolean *sma,my_boolean *sma_con,
			int nreq,int nreq_con,int *nameaf,
			int *nameaf_con,struct flt *fault,
			struct med *medium)

{
  int imatmode;
  // determine way to get interaction coefficient
  imatmode=select_i_coeff_calc_mode(medium);
  switch(imatmode){
  case I_MAT_IN_MEMORY:{
    assemble_ap_matrix_1(a,naflt,naflt_con,sma,sma_con,nreq,nreq_con,nameaf,nameaf_con,fault,medium);
    break;
  }
  case I_MAT_ON_FILE:{
    assemble_ap_matrix_2(a,naflt,naflt_con,sma,sma_con,nreq,nreq_con,nameaf,nameaf_con,fault,medium);
    break;
  }
  case CALC_I_COEFF_NOW:{
    assemble_ap_matrix_3(a,naflt,naflt_con,sma,sma_con,nreq,nreq_con,nameaf,nameaf_con,fault,medium);
    break;
  }
  case SPARSE_I_MAT:{
    assemble_ap_matrix_4(a,naflt,naflt_con,sma,sma_con,nreq,nreq_con,nameaf,nameaf_con,fault,medium);
    break;
  }
  default:{
    fprintf(stderr,"assemble_ap_matrix: can not deal with I matrix code %i\n",
	    imatmode);
    exit(-1);
  }}
}
#endif


