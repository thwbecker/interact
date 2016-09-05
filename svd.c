/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  routines in here deal with the general solution of a A x = b system
  in a least squares sense by means of singular value decomposition


  $Id: svd.c,v 2.37 2003/12/23 04:04:59 becker Exp $

*/
#include <math.h>
#include <stdio.h>
#include "interact.h"
#include "numrec_svd_routines.h"
/* 

   solve a singular linear equation system by SVD decomposition

   A[m,n] x[n] = b[m]

   where A is m by n, b is m by 1 and x has dimension n

   will return the largest singular value, wmax

   other input: 
   
   wlim:

   the singular value threshold for nearly singular matrices all
   entries smaller than this fraction times the maximum singular value
   will be set to zero

   if op_mode is set to ONE_STEP_CALCULATION (1), will print
   a verbose message about the solution

   WARNING

   if ireduce_n is set to >0 values, will only use n-ireduce_n
   parameters, ie. the solution will be for

   A[m,n-ireduce_n] x[n-ireduce_n] = b[m]

   instead

   WARNING:
   
   solver assumes that matrices A and V are stored in FORTRAN
   convention, not in C style

   A MATRIX GETS OVERWRITTEN, AND WILL SHRINK IF IREDUCE_N IS NON-ZERO

*/
/* use the divide and conquer method? */

void svd_driver_lapack(A_MATRIX_PREC *a,A_MATRIX_PREC *xsol,
		       A_MATRIX_PREC *b,int m, int norig, 
		       A_MATRIX_PREC *wlim, int op_mode,
		       A_MATRIX_PREC *wmax,A_MATRIX_PREC **sval,
		       int ireduce_n, my_boolean verbose_warn,
		       my_boolean debug)
{
  A_MATRIX_PREC *vwork,*bloc;
  int rank,info,solve_n;		
  static int nrhs=1,nwork;
  char out_string1[STRLEN],out_string2[STRLEN];
  static my_boolean query_svd = TRUE;
#ifdef USE_LAPACK_DAC
  static int smlsiz, smlsizp1p2, nlvl,oldm,oldn,liwork;
  int *iwork;
#endif
  A_MATRIX_PREC dummy;
  static int print_count=0;
  FILE *out1,*out2;
  // test size arguments
  if((m < norig)||(m < 1)||(norig < 1)||
     (ireduce_n < 0)||(ireduce_n >= norig)){
    fprintf(stderr,"svd_driver_lapack: error, m: %i n: %i ireduce_n: %i, either m < n or n,m < 1, or ireduce out of bounds\n",
	    m,norig,ireduce_n);
    exit(-1);
  }
  /* 

     real solution n may be reduced from the real n of input a

  */
  solve_n = norig - ireduce_n;
  reduce_a_matrix(&a,m,norig,solve_n); /* if solve_n = n, nothing happens */
  if(debug){
    /*
      
    debugging output
    
    */
    sprintf(out_string1,"/tmp/interact/a.%i.dat",print_count);
    sprintf(out_string2,"/tmp/interact/b.%i.dat",print_count);
    fprintf(stderr,"svd_driver_lapack: m: %i times n: %i system, writing A to %s and b to %s\n",
	    m,solve_n,out_string1,out_string2);
    if(ireduce_n)
      fprintf(stderr,"svd_driver_lapack: n was reduced by %i from %i\n",
	      ireduce_n,norig);
    out1=myopen(out_string1,"w");
    out2=myopen(out_string2,"w");
    print_a_matrix(a,m,solve_n,out1,&dummy,FALSE);
    print_b_vector(b,m,out2,&dummy,FALSE);
    fclose(out1);fclose(out2);
    if(ireduce_n)
      fprintf(stderr,"svd_driver_lapack: WARNING: reducing solution n by %i\n",
	      ireduce_n);
  }
  if((*sval=(A_MATRIX_PREC *)realloc(*sval,sizeof( A_MATRIX_PREC)*solve_n))==NULL)
    MEMERROR("svd_driver: 2a:");
  //
  // workspace: the factor in front of the parentheses is 
  // a multiple of the minimum workspace needed
  //
  /* for init */
  my_ivecalloc(&iwork,1,"iwork");
  my_vecalloc(&vwork,1,"vwork");
  my_vecalloc(&bloc,1,"vwork");
#ifdef USE_LAPACK_DAC
  if((oldn != solve_n)||(oldm != m)){		/* initialize size arrays */
    /* 
       get the blocksize from the wrapper for ILAENV 
    */
    ilwvd(&m,&solve_n,&nrhs,&smlsiz);
    smlsizp1p2 = (smlsiz+1)*(smlsiz+1);
    if(query_svd){
      /* query SVD routine */
      nwork=-1;liwork=1;
#ifdef A_MATRIX_SINGLE_PREC
      sgelsd_(&m,&solve_n,&nrhs,a,&m,bloc,&m,*sval,wlim,&rank,vwork,&nwork,iwork,&info);
#else
      dgelsd_(&m,&solve_n,&nrhs,a,&m,bloc,&m,*sval,wlim,&rank,vwork,&nwork,iwork,&info);
#endif
      fprintf(stderr,"svd_driver_lapack: determined best size from query %i\n",(int)vwork[0]);
      nwork = vwork[0];
    }else{
      if(m >= solve_n)
	nwork = 4 * (12*solve_n + 2*solve_n*smlsiz + 8*solve_n*nlvl + solve_n*nrhs + smlsizp1p2);
      else      
	nwork = 4 * (12*m + 2*m*smlsiz + 8*m*nlvl + m*nrhs + smlsizp1p2);
    }
    nlvl = (int)( log2((COMP_PRECISION)MIN(m,solve_n)/(COMP_PRECISION)(smlsiz+1))) + 1;
    fprintf(stderr,"svd_driver_lapack: determined smlsiz: %i nlvl: %i nwork: %i\n",
	    smlsiz,nlvl,nwork); 
    liwork = 2 * (3 * MIN(m,solve_n) * nlvl + 11 * MIN(m,solve_n));
    oldn = solve_n;oldm = m;

  }
  my_ivecrealloc(&iwork,liwork,"iwork");
#else
  nwork = 4 * (3*MIN(m,solve_n)  +  MAX(MAX(2*MIN(m,solve_n),MAX(m,solve_n)), nrhs));
#endif
  /* reallocate with the proper sizes */
  if((vwork=(A_MATRIX_PREC *)realloc(vwork,sizeof(A_MATRIX_PREC)*nwork))
     ==NULL)MEMERROR("svd_driver_lapack: 1:");
  if((bloc=(A_MATRIX_PREC *)realloc(bloc,sizeof(A_MATRIX_PREC)*m))==NULL)
    MEMERROR("svd_driver_lapack: 2b:");
  /* bloc = b  */
  memcpy(bloc,b,m*sizeof(A_MATRIX_PREC));
  /*
    
   solve A . x = b, bloc is b[m] on input and x[n] 
   solution on output

  */
#ifdef USE_LAPACK_DAC		/* divide and conquor algorithm */
#ifdef A_MATRIX_SINGLE_PREC
  //DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)
  sgelsd_(&m,&solve_n,&nrhs,a,&m,bloc,&m,*sval,wlim,&rank,vwork,&nwork,iwork,&info);
#else
  dgelsd_(&m,&solve_n,&nrhs,a,&m,bloc,&m,*sval,wlim,&rank,vwork,&nwork,iwork,&info);
#endif
#else
#ifdef A_MATRIX_SINGLE_PREC	/* old SVD algorithm */
  sgelss_(&m,&solve_n,&nrhs,a,&m,bloc,&m,*sval,wlim,&rank,vwork,&nwork,&info);
#else
  dgelss_(&m,&solve_n,&nrhs,a,&m,bloc,&m,*sval,wlim,&rank,vwork,&nwork,&info);
#endif
#endif
  if(info != 0){
#ifdef USE_LAPACK_DAC
    fprintf(stderr,"svd_driver_lapack: (s/d)gelsd error code: %i\n",
	    info);
#else
    fprintf(stderr,"svd_driver_lapack: (s/d)gelss error code: %i\n",
	    info);
#endif
    if(solve_n < 81){
      fprintf(stderr,"svd_driver_lapack: A after return:\n");
      print_matrix_ftrn(a,m,solve_n,stderr,FALSE);
      fprintf(stderr,"svd_driver_lapack: xsol after return:\n");
      print_vector(bloc,solve_n,stderr);
    }
    exit(-1);
  }
  // assign n elements of bloc to solution vector
  /* 
     xsol = bloc
  */
  memcpy(xsol,bloc,solve_n*sizeof(A_MATRIX_PREC));
  if((rank != solve_n)&&(verbose_warn))
    fprintf(stderr,"svd_driver_lapack: WARNING: rank: %i out of %i (nz: %i) for wlim: %g (singular?)\n",
	    rank,solve_n,solve_n - rank, *wlim);
  *wmax = *(*sval+0);// largest singular value
  // free local arrays
  free(bloc); 
  free(vwork);
#ifdef USE_LAPACK_DAC
  free(iwork);
#endif
  if(debug){
    sprintf(out_string1,"/tmp/interact/x.%i.dat",print_count);
    fprintf(stderr,"svd_driver_lapack: writing solution to %s\n",
	    out_string1);
    
    out1=myopen(out_string1,"w");
    print_b_vector(xsol,solve_n,out1,&dummy,FALSE);
    fclose(out1);
    print_count++;
  }
}

/*

  same as above but using numerical recipes with the option
  to return V and s (but not sorted!)
  
  if keep_vw is set, the singular values s and the V matrix are
  input as pointers, locally re-allocated and returned
  V hold the eigenvetors of singular values 1,2,3 ...
  in columns 1,2,3...

  that means that s and v should be initialized to NULL
  
  wlim gets overwritten!

  the s and V vector and matrix are not 
  sorted! 

  if calc_gen_inverse is set, will calculate and return the 
  general inverse matrix

  G^-g = V * diag(1/w_i) * U^T 

  g has to be initialized as NULL in this case


  if use_null_space is set, will produce the solution that 
  corresponds to the singular values which would normally 
  have been zeroed out. the null space

*/

void svd_driver_numrec(A_MATRIX_PREC *a,A_MATRIX_PREC *xsol,
		       A_MATRIX_PREC *b,int m, int norig, 
		       A_MATRIX_PREC *wlim, int op_mode,
		       A_MATRIX_PREC *wmax,
		       my_boolean keep_vw,A_MATRIX_PREC **s,
		       A_MATRIX_PREC **v,int nsv_zero,
		       my_boolean calc_ginv,
		       int ireduce_n,my_boolean use_null_space,
		       my_boolean verbose_warn,my_boolean debug)
{
  A_MATRIX_PREC *vwork,*w,wminlim,*utemp,*ginv;
  int i,j,rank,*index,solve_n;
  FILE *out;
  char out_string1[STRLEN],out_string2[STRLEN];
  A_MATRIX_PREC dummy;
  static int print_count=0;
  FILE *out1,*out2;
  if(use_null_space)
    fprintf(stderr,"\nsvd_driver_numrec: WARNING: output of nullspace, not solution\n\n");
  //
  // test size arguments
  //
  if((m < norig)||(m < 1)||(norig < 1)
     ||(ireduce_n < 0)||(ireduce_n >= norig)){
    fprintf(stderr,"svd_driver_numrec: error, m: %i n: %i ireduce_n: %i, either m < n or n,m < 1, or ireduce out of bounds\n",
	    m,norig,ireduce_n);
    exit(-1);
  } 
  /* 

     real solution n may be reduced from the real n of input a
     minus the reduced n

  */
  solve_n = norig - ireduce_n;
  reduce_a_matrix(&a,m,norig,solve_n); /* if solve_n = n, nothing happens */
  if(debug){
    //
    // debugging output
    //
    sprintf(out_string1,"/tmp/interact/a.%i.dat",print_count);
    sprintf(out_string2,"/tmp/interact/b.%i.dat",print_count);
    fprintf(stderr,"svd_driver_numrec: %i times %i system, writing A to %s and b to %s\n",
	    m,solve_n,out_string1,out_string2);
    out1=myopen(out_string1,"w");
    out2=myopen(out_string2,"w");
    print_a_matrix(a,m,solve_n,out1,&dummy,FALSE);
    print_b_vector(b,m,out2,&dummy,FALSE);
    fclose(out1);fclose(out2);
    if(ireduce_n > 0)
      fprintf(stderr,"svd_driver_numrec: WARNING: reducing solution n by %i\n",
	      ireduce_n);
  }
  // allocate space for the singular values 
  if(keep_vw){
    //
    // keep them in the s array
    //
    if((*s=(A_MATRIX_PREC *)realloc(*s,sizeof( A_MATRIX_PREC)*solve_n))==NULL)
      MEMERROR("svd_driver_numrec: 2a:");
    w= *s;
    if((*v=(A_MATRIX_PREC *)realloc(*v,sizeof(A_MATRIX_PREC)*SQUARE(solve_n)))==NULL)
      MEMERROR("svd_driver_numrec: 1i:");
    vwork = *v;
  }else{
    // only locally needed
    if((w=(A_MATRIX_PREC *)malloc(sizeof( A_MATRIX_PREC)*solve_n))==NULL)
      MEMERROR("svd_driver_numrec: 2a:");
    if((vwork=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*
				      SQUARE(solve_n)))==NULL)
      MEMERROR("svd_driver_numrec: 1ii:");
  }
  /*
    
    do SVD decomposition, on output, 
    a[1..m][1..n] will be U
    vwork will be V[1..n][1..n] (NOT the transpose)
    w[1..n] the singular values
    
  */
  svdcmp(a, &m, &solve_n, &m, &solve_n, w, vwork);
  /* 
     find the max singular value
  */
  *wmax = find_max_vec(w,solve_n);
  if(nsv_zero == 0){		
    /* 
       USE TO ELIMINATE SMALL SV 
    */
    /* determine cutoff from wlim */
    wminlim = (*wmax) * (*wlim);
    if(use_null_space){
      /* produce the null space, not the solution */
      for(rank=solve_n,
	    i=0;i < solve_n;i++){
	if(w[i] >= wminlim){
	  w[i]=0.0;
	  rank--;
	}
      }
      
    }else{
      /* 
	 normal operational mode
      */
      // set all values smaller than wlim * wmax to zero
      for(rank=solve_n,
	    i=0;i < solve_n;i++){
	if(w[i] < wminlim){
	  w[i]=0.0;
	  rank--;
	}
      }
    }
  }else{
    /* 

    ELIMINATE A SPECIFIED NUMBER OF SINGULAR VALUES FROM SYSTEM

    */
    if(nsv_zero >= solve_n){
      fprintf(stderr,"svd_driver_numrec: error: n: %i nsv_zero: %i\n",
	      solve_n,nsv_zero);
      exit(-1);
    }
    /* 
       modified rank of matrix 
    */
    rank = solve_n - nsv_zero;
    /* 
       we need to sort the singular values for elimination.

       generate ascending SV index
    */
    my_ivecalloc(&index,solve_n,"svd: index");
    indexx(solve_n,(w-1),(index-1));   /* 
					  index now holds 0 .. n-1
					  indices such that
					  w[index[i]] is ascending
				       */
    /* override wlim according to the largest SV singled out */
    *wlim = w[index[nsv_zero-1]]/(*wmax);
    if(use_null_space){
      /* 
	 
	 NULLSPACE
	 
      */
      for(i=nsv_zero;i < solve_n;i++){
	/* 
	   zero out the largest (> nsv_zero) singular values
	*/
	w[index[i]]=0.0;
      }
    }else{
      /* 
	 
      NORMAL OPERATION

      */
      for(i=0;i < nsv_zero;i++){
	/* 
	    zero out the first nsv_zero singular values 
	*/
	if((i > 1)&&(w[index[i]] < w[index[i-1]])){
	  fprintf(stderr,"svd_driver_numrec: sorting error: i: %i rank: %i n: %i %g %g\n",
		  i,rank,solve_n,w[index[i-1]],w[index[i]]);
	   exit(-1);
	}
	w[index[i]]=0.0;
       }
    }
    free(index);
  }
  if((rank != solve_n)&&(verbose_warn))
    fprintf(stderr,"svd_driver_numrec: WARNING: rank: %i out of %i (nz: %i) for wlim: %g and nsv_zero: %i\n",
	    rank,solve_n,solve_n-rank,*wlim,nsv_zero);
  /*
    
     backsubstitute and solve a . x = b  

  */
  svbksb(a,w,vwork,&m,&solve_n,&m,&solve_n,b,xsol); 
  if(calc_ginv){
    /* calculate the general inverse 

    V * diag(1/w_i) * U^T 

    */
    my_vecalloc(&ginv,m*solve_n,"svd: g"); /* allocate space for G */
    my_vecalloc(&utemp,m*solve_n,"svd: utemp"); /* workspace */
    cginv(ginv,vwork,w,a,utemp,&m,&solve_n);
    free(utemp);
    /* write general inverse to file */
    fprintf(stderr,"svd_driver_numrec: writing general inverse n: %i m: %i to ginv.dat\n",
	    solve_n,m);
    out=myopen("ginv.dat","w");
    print_matrix_ftrn(ginv,solve_n,m,out,FALSE);
    fclose(out);
    free(ginv);
    /* 

    calculate the model resolution matrix 
    V .V^T

    */
    my_vecalloc(&utemp,solve_n*solve_n,"svd: utemp"); /* workspace for V^T */
    for(i=0;i < solve_n;i++)		/* utemp = V^T */
      for(j=0;j < solve_n;j++)
	utemp[i*solve_n+j] = vwork[j*solve_n+i];
    my_vecalloc(&ginv,solve_n*solve_n,"svd: ginv");	/* for V V^T */
    calc_AB_ftn(vwork,solve_n,solve_n,utemp,solve_n,ginv);
    free(utemp);
    /* 
       write model resolution matrix to file 
    */
    fprintf(stderr,"svd_driver_numrec: writing model resolution V.V^T to vvt.dat (%i by %i)\n",
	    solve_n,solve_n);
    out=myopen("vvt.dat","w");
    print_matrix_ftrn(ginv,solve_n,solve_n,out,FALSE);
    fclose(out);
    free(ginv);
    if(0){
      /* 
	 
      calculate the data resolution matrix 
      U .U^T
      
      */
      my_vecalloc(&utemp,m*solve_n,"svd: utemp"); /* workspace for U^T */
      for(i=0;i < solve_n;i++)		/* utemp = U^T */
	for(j=0;j < m;j++)
	  utemp[i+j*solve_n] = a[j+i*m];
      fprintf(stderr,"this has to be debugged\n");
      exit(-1);
      my_vecalloc(&ginv,m*m,"svd: ginv");	/* for U U^T */
      calc_AB_ftn(a,m,solve_n,utemp,m,ginv);
      free(utemp);
      /* write model resolution matrix to file */
      fprintf(stderr,"svd_driver_numrec: writing data resolution (%i by %i) in binary U.U^T to uut.bin\n",
	      m,m);
      out=myopen("uut.bin","w");
      print_matrix_ftrn(ginv,m,m,out,TRUE);
      fclose(out);
      free(ginv);
    }
  }
  if(!keep_vw){
    free(vwork);
    free(w);
  }
  if(debug){
    /* 
       print solution to file one time
    */
    sprintf(out_string1,"/tmp/interact/x.%i.dat",print_count);
    fprintf(stderr,"svd_driver_numrec: writing solution [%i dim] to %s\n",
	    solve_n,out_string1);
    out1=myopen(out_string1,"w");
    print_b_vector(xsol,solve_n,out1,&dummy,FALSE);
    fclose(out1);
    if(keep_vw){
      // print w
      fprintf(stderr,"svd_driver_numrec: writing singular values to w.asc\n");
      print_vector_file(w,solve_n,"w.asc","w vector");
      // print V
      fprintf(stderr,"svd_driver_numrec: writing V matrix to v.asc\n");
      out1=myopen("v.asc","w");
      print_matrix_ftrn(*v,solve_n,solve_n,out1,FALSE);
      fclose(out1);
    }
    print_count++;
  }
}

/*

  assemble the covariance matrix based on the singular values
  and V matrix as returned by svd_driver_numrec if requested

  COV(a_j,a_k) = sum(i=1...n) ((Vji Vki)/wi^2)
  
  v(n,n) and w(n) are input

  cvm will be assigned on output

*/

void assemble_cov_svd(COMP_PRECISION *v,int n, 
		      COMP_PRECISION *w,
		      COMP_PRECISION *cvm)
{
  int k,j,i,o1;
  COMP_PRECISION sum,*wti;
  my_vecalloc(&wti,n,"wti");
  //
  // get scaling factors
  //
  for (i=0;i < n;i++) {
    if (w[i]) 
      wti[i]=1.0/(w[i]*w[i]);
    else
      wti[i] = 0.0;
  }
  //
  // assemble matrix
  //
  for (j=0;j < n;j++) {		/* columns of cvm */
    for (k=0;k <= j;k++) {        /* rows of cvm, do only upper 
				   right */
      sum = 0.0;
      for (i=o1=0;i < n;i++,o1+=n)	/* add up all n 
					   contributions */
	sum += v[j+o1] * v[k+o1] * wti[i];
      /* assign both C(j,k) and C(k,j) */
      cvm[k+j*n] = cvm[j+k*n] = sum;
    }
  }
  free(wti); 
}
/*

  sort singular values and print to out

*/
void print_singular_values(COMP_PRECISION *w, int n, FILE *out)
{
  COMP_PRECISION *wc;
  int i,j;
  my_vecalloc(&wc,n,"wc");
  a_equals_b_vector(wc,w,n);
  qsort(wc,n,sizeof(COMP_PRECISION),compare_flt);
  j=n-1;
  for(i=j;i >= 0;i--)
    fprintf(out,"%20.12e %20.12e\n",wc[i],wc[i]/wc[j]);
  free(wc);
}


/* 
   sort numerical recipes style arrays (1...n) 
   and gives back indices  in 0 ... n-1 fashion
   i.e. pass both arr and indx with (...-1) from
   normal c style programs
*/

#define NUMREC_MINDEXX 7
#define NUMREC_NSTACK 50

void indexx(int n,COMP_PRECISION *arr,int *indx)
{
  int i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  COMP_PRECISION a;
  
  my_ivecalloc(&istack,NUMREC_NSTACK+1,"indexx");

  for (j=1;j<=n;j++) 
    indx[j]=j;
  for (;;) {
    if (ir-l < NUMREC_MINDEXX) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j];
	a=arr[indxt];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      NUMREC_ISWAP(indx[k],indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	NUMREC_ISWAP(indx[l+1],indx[ir])
	  }
      if (arr[indx[l]] > arr[indx[ir]]) {
	NUMREC_ISWAP(indx[l],indx[ir])
	  }
      if (arr[indx[l+1]] > arr[indx[l]]) {
	NUMREC_ISWAP(indx[l+1],indx[l])
	  }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	NUMREC_ISWAP(indx[i],indx[j])
	  }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NUMREC_NSTACK){
	fprintf(stderr,"index: NUMREC_NSTACK too small in indexx.");
	exit(-1);
      }
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free(istack);
  /* go to 0 ... n-1 notation */
  for(i=1;i <= n;i++)
    indx[i]--;
}
#undef NUMREC_MINDEXX
#undef NUMREC_NSTACK
/* 

given a FORTRAN style a[m,n] matrix, modify this matrix to be
a[m,reduced_n], ie. throw away a few columns of a


*/
void reduce_a_matrix(A_MATRIX_PREC **a, int m, int n,
		     int reduced_n)
{
  if((reduced_n > n)||(reduced_n <=0)){
    fprintf(stderr,"reduce_a_matrix: error: reduced n: %i cannot be bigger than original n: %i, nor 0\n",
	    reduced_n,n);
    exit(-1);
  }
  if(reduced_n == n){		/* no change needed */
    return;
  }else{			
    /* this should do the trick for fortran sorting  */
    my_vecrealloc(a,reduced_n*m,"ram: 2");
  }
}
