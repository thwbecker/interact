/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@post.harvard.edu


  all routines in here deal with the solution of interact program
  specific routines


  parallel Petsc implementation based on an example by Dave May (UCSD)

*/
#include "interact.h"

int solve(struct med *medium,struct flt *fault)
{
  A_MATRIX_PREC *a,wmax,*sv=NULL,wcutoff;
#ifdef USE_NUMREC_SVD
  A_MATRIX_PREC *dummyp=NULL;
#endif
  char solver_name[STRLEN];
  size_t full_size,sparse_size;
  static int izero=0;
  FILE *aio;
  long int isize,index_numbers;
#ifdef USE_SLATEC_NNLS
  A_MATRIX_PREC *x,prgopt[1]={1.0},rnorm,*work,dummy;
  int *iwork,mode,l,nm;
  FILE *out1,*out2;
#endif
#ifdef USE_PETSC
#define PETSC_HELPER_STR_LEN 256  
  char mattype[PETSC_HELPER_STR_LEN],out_string[300];
  PetscBool pset = PETSC_FALSE;
  Vec         px, pr,pxout;
  KSP         pksp;
  PC          ppc;
  PetscInt    i, m, n;
  PetscScalar *values=NULL;
  PetscReal   norm;
  PetscInt lm, ln, dn, on;
  VecScatter ctx;
#else
  unsigned int i,m,n;
#endif

  
  char command_str[STRLEN];
  /* test valgrind */
  /* int test_a; */
  /* if (test_a==1) printf("random value is 1\n"); */


  
  wcutoff = (A_MATRIX_PREC)medium->wcutoff;

  HEADNODE
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
    if(medium->comm_size > 1){
      HEADNODE
	fprintf(stderr,"solve: SLATEC MIX NNLS only serial, %i cores requested\n",medium->comm_size);
      exit(-1);
    }
#ifndef USE_SLATEC_NNLS
    HEADNODE{
      fprintf(stderr,"solve: for mixed contrained/uncontrained solution, have to use SLATEC\n");
      fprintf(stderr,"solve: library functions. Read the makefile, and compile with -lslatec_???.\n");
    }
    exit(-1);
#else
    if(medium->debug&&(medium->comm_rank==0))
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
      HEADNODE
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
    }
    //
    // solve using SLATEC routines
    //
#ifndef A_MATRIX_PREC_IN_DOUBLE
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
    if(medium->force_petsc || (medium->comm_size>1)){
      /* parallel solve using PETSC */
      HEADNODE{
	sprintf(out_string,"attempting Petsc LU solve, %i core(s) requested",
	       medium->comm_size);
	time_report("solve",out_string,medium);
      }
      if(medium->use_old_amat || medium->save_amat){
	HEADNODE
	  fprintf(stderr,"solve: matrix I/O not implemented yet\n");
	exit(-1);
      }
      if(medium->solver_mode != LU_SOLVER){
	HEADNODE
	  fprintf(stderr,"solve: only LU solve availabled\n");
	exit(-1);
      }
#ifndef USE_PETSC
      HEADNODE
	fprintf(stderr,"solve: parallel solve requires Petsc compile\n");
      exit(-1);
#else
      /* 

	 Petsc parallel LU 
	 
      */
      m = n = medium->nreq;
      HEADNODE{
	index_numbers = (long int)m*(long int)m;
	if(index_numbers/medium->comm_size > INT_MAX){
	  /* PetSc is not compiled with 64 bit indexing by default,
	     this can cause a core dump 
	     
	     --with-64-bit-indices=yes
	     
	     can change that
	     
	  */
	  fprintf(stderr,"solve: WARNING: %.3f G ints on %i cores, above 32b-int limit by %.1f%%\n",
		  (double)index_numbers/1e9,medium->comm_size,
		  ((double)index_numbers/(double)medium->comm_size)/(double)INT_MAX*100);
	}
      }
      /* set up A matrix */
      PetscCall(MatCreate(PETSC_COMM_WORLD, &(medium->pA)));
      PetscCall(MatSetSizes(medium->pA, PETSC_DECIDE, PETSC_DECIDE, m, n));
      PetscCall(MatSetType(medium->pA, MATDENSE));
      /*PetscCall(MatSetFromOptions(medium->pA));*//* Always assemble into MATDENSE format */
      PetscCall(MatSetUp(medium->pA));
      /* preallocate */
      PetscCall(MatGetLocalSize(medium->pA, &lm, &ln));
      dn = ln;on = n - ln;

      PetscCall(MatSeqAIJSetPreallocation(medium->pA, n, NULL));
      PetscCall(MatMPIAIJSetPreallocation(medium->pA, dn, NULL, on, NULL));
      /* 
	 parallel assembly 
      */
      PetscCall(MatGetOwnershipRange(medium->pA, &medium->rs, &medium->re));
      medium->rn = medium->re  - medium->rs; /* number of local elements */
      
#ifdef DEBUG
      fprintf(stderr,"solve: core %i: dn %i on %i n %i rs %i re %i \n",
	      medium->comm_rank,dn,on,n,medium->rs,medium->re);
#endif
      
      /* assemble A */
      par_assemble_a_matrix(medium->naflt,medium->sma,medium->nreq,medium->nameaf,
			    fault,medium);
      /*
	Always call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY).
	MatAssemblyBegin() is collective and must be called on all ranks.
      */
      PetscCall(MatAssemblyBegin(medium->pA, MAT_FINAL_ASSEMBLY));
      PetscCall(MatAssemblyEnd(medium->pA, MAT_FINAL_ASSEMBLY));
      HEADNODE
	time_report("solve","parallel assembly done",medium);
     
      /* Convert MATDENSE to another format required by solver package */
      PetscCall(PetscOptionsGetString(NULL, NULL, "-mat_type", mattype, PETSC_HELPER_STR_LEN, &pset));
      if (pset) { /* Convert MATDENSE to desired format */
	PetscCall(MatConvert(medium->pA, mattype, MAT_INPLACE_MATRIX, &medium->pA));
      }

      //MatView(pA,PETSC_VIEWER_STDOUT_WORLD);
  
      PetscCall(MatCreateVecs(medium->pA, &medium->pb, &px)); /* For A x = b: x -> left, b -> right */
      /* 
	 insert right hand side 
      */
      for (i = medium->rs; i < medium->re; i++) {
	PetscCall(VecSetValue(medium->pb, (PetscInt)i, (PetscScalar)(medium->b[i]), INSERT_VALUES));
      }
      PetscCall(VecAssemblyBegin(medium->pb));
      PetscCall(VecAssemblyEnd(medium->pb));

      //PetscCall(VecView(medium->pb,PETSC_VIEWER_STDOUT_WORLD)); 
      /* 
	 solver
      */
      PetscCall(KSPCreate(PETSC_COMM_WORLD, &pksp));
      PetscCall(KSPSetOperators(pksp, medium->pA, medium->pA));
      PetscCall(KSPSetType(pksp, KSPPREONLY));
      PetscCall(KSPGetPC(pksp, &ppc));
      PetscCall(PCSetType(ppc, PCLU));
      /* override at run time via -pc_factor_mat_solver_type xxx */
      PetscCall(PCFactorSetMatSolverType(ppc, MATSOLVERPETSC));
      PetscCall(KSPSetFromOptions(pksp));
      /* solve step */
      PetscCall(KSPSolve(pksp, medium->pb, px));
      if(medium->debug)
	PetscCall(KSPView(pksp, PETSC_VIEWER_STDERR_WORLD));
    
      /* 
	 distribute to zero node
      */
      PetscCall(VecScatterCreateToZero(px,&ctx,&pxout));
      // scatter as many times as you need
      PetscCall(VecScatterBegin(ctx,px,pxout,INSERT_VALUES,SCATTER_FORWARD));
      PetscCall(VecScatterEnd(ctx,px,pxout,INSERT_VALUES,SCATTER_FORWARD));
      /* assign to x solution vector */
      if(medium->comm_rank == 0){
	PetscCall(VecGetArray(pxout,&values));
	for(i=0;i<m;i++)
	  medium->xsol[i] = (A_MATRIX_PREC)values[i];
	PetscCall(VecRestoreArray(pxout,&values));
      }
      HEADNODE
	time_report("solve","parallel solve done, broadcasting",medium);
      /* broadcast solution */
#ifdef A_MATRIX_PREC_IN_DOUBLE
      MPI_Bcast(medium->xsol,m,MPI_DOUBLE,0, MPI_COMM_WORLD);
#else
      MPI_Bcast(medium->xsol,m,MPI_FLOAT,0, MPI_COMM_WORLD);
#endif
      // destroy scatter context and local vector when no longer needed
      PetscCall(VecScatterDestroy(&ctx));
      PetscCall(VecDestroy(&pxout));
 
      if(medium->debug){
	PetscCall(VecDuplicate(px, &pr));
	/* check residual */
	PetscCall(MatMult(medium->pA, px, pr)); /* r = A x */
	PetscCall(VecAXPY(pr, -1.0, medium->pb)); /* r <- r - b */
	PetscCall(VecNorm(pr, NORM_2, &norm));
	PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Norm of residual %1.10e\n", (double)norm));
	PetscCall(VecDestroy(&pr));
      }
      
      PetscCall(VecDestroy(&px));
      PetscCall(VecDestroy(&medium->pb));
      PetscCall(MatDestroy(&medium->pA));
      PetscCall(KSPDestroy(&pksp));
      /* 
	 petsc done 
      */
      HEADNODE
	time_report("solve","parallel solve completed",medium);

#endif
    }else{			/* serial assembly */
      isize  = ((long int)sizeof(A_MATRIX_PREC)) * ((long int)medium->nreq*(long int)medium->nreq);
      if(medium->debug)
	fprintf(stderr,"solve: attempting to allocate %g GB for A matrix\n",(double)isize/ONE_MEGABYTE/1024.);
      if((a=(A_MATRIX_PREC *)malloc((size_t)isize)) == NULL){ /* do some testing for memory error */
	fprintf(stderr,"solve: memory error for unconstrained A: (%p) n: %i by %i dsz: %i tsize: %g MB\n",
		a,medium->nreq,medium->nreq,(int)sizeof(A_MATRIX_PREC),((double)isize)/ONE_MEGABYTE);
	/* check what would have been ok */
	isize  = (long int)((double)isize/10);
	for(i=1;i<10;i++){
	  if((a=(A_MATRIX_PREC *)malloc((size_t)(isize*(long int)i)))==NULL){
	    fprintf(stderr,"solve: failing at %i %g MB\n",i,((double)isize)/(double)(ONE_MEGABYTE));
	    exit(-1);
	  }else{
	    fprintf(stderr,"solve: OK at %i/10\n",i);
	    free(a);
	  }
	}
	exit(-1);
      }	
      if(medium->debug)
	fprintf(stderr,"solve: unconstrained part A is %5i by %i\n",
		medium->nreq,medium->nreq);
      if(medium->use_old_amat){
	fprintf(stderr,"solve: unconstrained: WARNING: using old A matrix!\n");
	read_a_matrix_from_file(a,medium->nreq,medium->nreq,A_MATRIX_FILE,
				"unconstrained A");
      }else{
	assemble_a_matrix(a,medium->naflt,medium->sma,medium->nreq,medium->nameaf,
			  fault,medium);
      }
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
	if(medium->debug){
	  if(!medium->save_amat){// if A is not to be saved, remove
	    sprintf(command_str,"rm %s",A_MATRIX_FILE);
	    system(command_str);
	  }
	}
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
    } /* end serial direct solve part */
  }else if(medium->nreq_con){
    if(medium->comm_size>1){
      HEADNODE
	fprintf(stderr,"solve: NNLS only serial, %i cores requested\n",medium->comm_size);
      exit(-1);
    }
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

#ifdef USE_PETSC


int par_assemble_a_matrix(int naflt,my_boolean *sma,int nreq,int *nameaf,
			  struct flt *fault,struct med *medium)
{
  /*
    see routine below for the logic of the stress assignments
  */
  PetscInt i,j,k,l,eqc1,eqc2,ip1,ip2;
  PetscScalar cf,*avalues=NULL,fac_itmp;
  PetscInt *col_idx=NULL;
  my_boolean in_range;
  int iret;
#ifdef DEBUG
  my_boolean *assigned;
  PetscScalar amin,amax;
  
  fprintf(stderr,"par_assemble_a_matrix: core %03i/%03i: assigning row %5i to %i\n",
	  medium->comm_rank,medium->comm_size,medium->rs,medium->re);

#endif

  PetscCall(PetscCalloc(nreq*sizeof(PetscScalar), &avalues));
  PetscCall(PetscCalloc(nreq*sizeof(PetscInt), &col_idx));
#ifdef DEBUG
  assigned = (my_boolean *)malloc(sizeof(my_boolean)*nreq);
#endif
  for (j=0; j < nreq; j++) 
    col_idx[j] = j;
  
  
  // dimensions of A matrix (without b column)
  // m = n = nreq;
  /*
    COMMENT: the flipped sign for constrained
    faults was determined in interact.c
    
    nameaf[i]: observing fault
    j: stress component on observing fault
    nameaf[k]: slipping fault
    l: slip mode on slipping fault
    
  */
  for(eqc1=ip1=i=0;i < naflt;i++,ip1+=3){ /* slow fault loop */
    for(j=0;j < 3;j++){			  /* direction */
      if(sma[ip1+j]){
	/* 
	   do we have work on this processor? 
	*/
	in_range = ((eqc1 >= medium->rs)&&(eqc1 < medium->re))?(TRUE):(FALSE);
	if(in_range){
#ifdef DEBUG
	  for(eqc2=0;eqc2<nreq;eqc2++)
	    assigned[eqc2] = FALSE;
#endif

	  /* 
	     actually compute
	  */
	  // normal correction for Coulomb?
	  cf = (j==NORMAL)?(0.0):fault[nameaf[i]].cf[j];

	  for(eqc2=ip2=k=0;k < naflt;k++,ip2+=3){ /* fast fault loop */
	    for(l=0;l < 3;l++){			  /* direction */
	      if(sma[ip2+l]){// flip around matrix ordering
		/* 
		   actually compute
		*/
#ifdef SUPER_DUPER_DEBUG
		if(cf != 0.0)
		  fprintf(stderr,"par_assemble_a_matrix: i: %i j: %i k: %i l: %i cf: %g\n",
			  (int)nameaf[i],(int)nameaf[k],(int)l,(int)j,cf);
#endif
#ifdef DEBUG
		if(eqc2 >= nreq){
		  fprintf(stderr,"par_assemble_a_matrix: core %i: eqc2 %i out of bounds, nreq %i\n",
			  medium->comm_rank,eqc2,nreq);
		  fprintf(stderr,"par_assemble_a_matrix: i: %i j: %i k: %i l: %i naflt: %i\n",
			  (int)nameaf[i],(int)nameaf[k],(int)l,(int)j,naflt);
		  exit(-1);
		}
#endif
		if(medium->no_interactions && (fault[i].group != fault[k].group)){
		  avalues[eqc2] = 0;
		}else{
		  //
		  // calculate interaction coefficients right now
		  //
		  avalues[eqc2] = (PetscScalar)
		    interaction_coefficient(nameaf[i],nameaf[k],l,j,fault,&iret);
		  if(cf != 0.0){	/* coulomb addition */
		    fac_itmp=(PetscScalar)
		      interaction_coefficient(nameaf[i],nameaf[k],l,NORMAL,fault,&iret);
		    if(iret){
		      fprintf(stderr,"par_assemble_a_matrix: WARNING: encountered iret: i/j/k/l: %i/%i/%i/%i\n",
			      (int)nameaf[i],(int)nameaf[k],(int)l,(int)j);
		      fac_itmp = 0.0;
		    }
		    avalues[eqc2] +=  fac_itmp * cf;
		  }
		}
#ifdef DEBUG
		assigned[eqc2] = TRUE;
#endif
	   	//PetscCall(MatSetValue(medium->pA, eqc1, eqc2, avalues[eqc2], ADD_VALUES));
	      	/* end value work part */
		eqc2++;
	      }
	    }
	  }
#ifdef DEBUG
	  amin = 1e20;amax=-1e20;
	  for(eqc2=0;eqc2 < nreq;eqc2++){
	    if(avalues[eqc2]< amin)
	      amin = avalues[eqc2];
	    if(avalues[eqc2]> amax)
	      amax = avalues[eqc2];
	    if(!assigned[eqc2])
	      fprintf(stderr,"%06i out of %06i unassigned?!\n",eqc2,nreq);
	    if(!finite(avalues[eqc2]))
	      fprintf(stderr,"%06i out of %06i not finite\n",eqc2,nreq);
	  }
	  fprintf(stderr,"row eqc1 %6i eqc2 %6i nreq %6i amin %12g amax %12g\n",eqc1,eqc2,nreq,amin,amax);
	  
#endif 
	  PetscCall(MatSetValues(medium->pA, 1, &eqc1, nreq, col_idx,avalues, INSERT_VALUES));
	}
	eqc1++; 
      }
    }
  }
  PetscCall(PetscFree(col_idx));
  PetscCall(PetscFree(avalues));
#ifdef DEBUG
  free(assigned);
#endif
  return(0);
}

#endif
