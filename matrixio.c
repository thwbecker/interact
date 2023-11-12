/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  
  routine takes care of matrix input and output

  internally, we store the A matrix that is passed to the 
  solvers in FORTRAN style not in C style

  i.e. for a m (row) by n (col) matrix, row i and column j is stored
  as a[j*m+i] and not as a[i*n+j] as in C

  hence, we have to flip things on I/O since we like to store
  stuff in C style to make interfacing with other programs easier


  also contains some more harmless functions to print vectors
  

  $Id: matrixio.c,v 1.9 2011/01/07 07:19:58 becker Exp $ 

*/
#include "interact.h"


/*

  print matrix in ASCII or BINARY to stream 
  out using FORTRAN convention

*/
void print_matrix_ftrn(A_MATRIX_PREC *a,int m, int n, FILE *out,
		       my_boolean binary)
{
  int i,j,k;
  if(binary){
    for(i=0;i<m;i++)
      for(k=j=0;j<n;j++,k+=m)
	if(fwrite((a+i+k),sizeof(A_MATRIX_PREC),1,out)!=1){
	  fprintf(stderr,"print_matrix: write error: i: %i j: %i\n",
		  i,j);
	  exit(-1);
	}
  }else{
    for(i=0;i<m;i++){
      for(k=j=0;j<n;j++,k+=m)
	fprintf(out,"%20.15e ",a[i+k]);
      fprintf(out,"\n");
    }
  }
}
/* same for C style */
void print_matrix_C(A_MATRIX_PREC *a,int m, int n, FILE *out,
		    my_boolean binary)
{
  int i,j,k;
  if(binary){
    for(i=k=0;i<m;i++,k+=n)
      for(j=0;j<n;j++)
	if(fwrite((a+j+k),sizeof(A_MATRIX_PREC),1,out)!=1){
	  fprintf(stderr,"print_matrix: write error: i: %i j: %i\n",
		  i,j);
	  exit(-1);
	}
  }else{
    for(i=k=0;i<m;i++,k+=n){
      for(j=0;j<n;j++)
	fprintf(out,"%20.15e ",a[j+k]);
      fprintf(out,"\n");
    }
  }
}
/* scaled version */
void print_matrix_scaled_ftrn(A_MATRIX_PREC *a,int m, int n, FILE *out,
			      my_boolean binary,COMP_PRECISION scale)
{
  int i,nm;
  A_MATRIX_PREC *aloc;
  nm = n * m;
  aloc = (A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*nm);
  if(!aloc)MEMERROR("print_matrix_scaled_ftrn");
  for(i=0;i<nm;i++)
    aloc[i] = a[i] * scale;
  print_matrix_ftrn(aloc,m, n, out,binary);
  free(aloc);
}
/* 

print FORTRAN style matrix to a file

 */
void print_matrix_ftrn_file(A_MATRIX_PREC *a, int m, int n, 
			    char *name, 
			    my_boolean binary)
{
  FILE *out;
  fprintf(stderr,"print_matrix_file: printing %i by %i matrix to %s\n",
	  m,n,name);
  out = myopen(name,"w");
  print_matrix_ftrn(a,m,n,out,binary);
  fclose(out);
}
/*

  print A matrix (fortran style) 
  to out and determine maximum absolute value

*/
void print_a_matrix(A_MATRIX_PREC *a,int m, int n, FILE *out, 
		    A_MATRIX_PREC *amax, my_boolean binary)
{
  print_matrix_ftrn(a,m,n,out,binary);
#ifdef A_MATRIX_PREC_IN_DOUBLE
  *amax = find_max_abs_vec(a,n*m);
#else
  *amax = find_max_abs_vec_float(a,n*m);
#endif
}
/*

  read in A matrix that is stored in C format and convert to
  internal FORTRAN storage
  
*/
void read_a_matrix_from_file(A_MATRIX_PREC *a,int m, int n,
			     char *filename,char *matrixname)
{
  char tmpstr[STRLEN];
  int nr,mr,size,nread,i,j,k,up;
  FILE *in;
  size_t sizeofa;
  float tmp_flt;
  double tmp_dbl;
  sizeofa = sizeof(A_MATRIX_PREC);
  sprintf(tmpstr,"%s.hdr",filename);
  fprintf(stderr,"read_a_matrix_from_file: reading %s from \"%s\" and \"%s\"\n",
	  matrixname,tmpstr,filename);
  in=myopen(tmpstr,"r");
  size=fscanf(in,"%i %i %i\n",&mr,&size,&nr);
  if(size!=3){
    fprintf(stderr,"could not read header\n");
    exit(-1);
	   
  }
  if((n!=nr)||(m!=mr)){
    fprintf(stderr,"read_a_matrix_from_file: header file indicates m: %i by n: %i system with %i byte precision\n",
	    mr,nr,size);
    fprintf(stderr,"read_a_matrix_from_file: request for A matrix was m: %i by n: %i with %i byte precision\n",
	    m,n,(int)sizeofa);
    fprintf(stderr,"read_a_matrix_from_file: error: dimension mismatch\n");
    exit(-1);
  }
  if(size != sizeofa){// precision mismatch
    fprintf(stderr,"read_a_matrix_from_file: header file indicates %i by %i system with %i byte precision\n",
	    mr,nr,size);
    fprintf(stderr,"read_a_matrix_from_file: request for A matrix was %i by %i with %i byte precision\n",
	    m,n,(int)sizeofa);
    // see if we can read in at different precision
    if((size == sizeof(double)) && (sizeofa==sizeof(float))){
      up = -1;
      fprintf(stderr,"read_a_matrix_from_file: saved matrix is double precision, will read in to internal single precision\n");
    }else if((size == sizeof(float)) && (sizeofa==sizeof(double))){
      up = 1;
      fprintf(stderr,"read_a_matrix_from_file: saved matrix is single precision, will read in to internal double precision\n");
    }else{
      fprintf(stderr,"read_a_matrix_from_file: error: saved matrix has %i bytes, internally we have %i byte precision?!\n",
	      size,(int)sizeofa);
      exit(-1);
    }
  }else// no adjustment necessary
    up = 0;
  fclose(in);
  in=myopen(filename,"r");
  // read matrix in binary format, sorted as in print_a_matrix
  switch(up){
  case 0:// no adjustment
    for(nread=i=0;i < m;i++)
      for(k=j=0;(j < n)&&
	    (fread((a+i+k),sizeof(A_MATRIX_PREC),1,in)==1);
	  j++,k+=m)
	nread++;
    break;
  case 1:// float -> double
    for(nread=i=0;i < m;i++)
      for(j=k=0;(j < n)&&(fread(&tmp_flt,sizeof(float),1,in)==1);
	  j++,k+=m){
	a[k+i] = (A_MATRIX_PREC) tmp_flt;
	nread++;
      }
    break;
  case -1:// double -> float
    for(nread=i=0;i < m;i++)
      for(j=k=0;(j < n)&&(fread(&tmp_dbl,sizeof(double),1,in)==1);
	  j++,k+=m){
	a[k+i] = (A_MATRIX_PREC) tmp_dbl;
	nread++;
      }
    break;
  default:
    fprintf(stderr,"read_a_matrix_from_file: error: up: %i, undefined\n",up);
    exit(-1);
  }
  if(nread != n*m){
    fprintf(stderr,"read_a_matrix_from_file: read error \"%s\", file too short?\n",
	    filename);
    fprintf(stderr,"read_a_matrix_from_file: %i values could be read, expected %i\n",
	    nread,n*m);
    exit(-1);
  }
  fclose(in);
}
/*

  write the A matrix of the linear eqaution solver to files, for
  debugging or restarting purposes. will generate a straight binary
  dump and a header file

*/
void print_a_matrix_to_file(A_MATRIX_PREC *a,int m, int n,
			    char *filename,char *matrixname)
{
  char tmpstr[STRLEN];
  A_MATRIX_PREC amax,acutoff;
  int maxbw;
  FILE *out;
  my_boolean calc_stats=FALSE;
  fprintf(stderr,"print_a_matrix_to_file: %s matrix (%s precision, m: %i by n: %i) size is %g MB\n",
	  matrixname,(sizeof(A_MATRIX_PREC)==4)?"float":"double",m,n,
	  (COMP_PRECISION)(n*m)*
	  sizeof(A_MATRIX_PREC)/ONE_MEGABYTE);
  sprintf(tmpstr,"%s.hdr",filename);
  fprintf(stderr,"print_a_matrix_to_file: writing to \"%s\" and \"%s\"\n",
	  tmpstr,filename);
  out=myopen(tmpstr,"w");
  fprintf(out,"%i %i %i\n",(int)m,(int)sizeof(A_MATRIX_PREC),(int)n);
  fclose(out);
  out=myopen(filename,"w");
  //
  // write matrix (and flip to C-style binary format)
  //
  print_a_matrix(a,m,n,out,&amax,TRUE);
  fclose(out);
  if(calc_stats){
    // calculate cutoff for two different cutoff levels
    acutoff = amax * 1.0e-4;
    fprintf(stderr,"print_a_matrix_to_file:  i) maximum is %10.5e, choosing %10.5e as cutoff for bandwidth\n",
	    amax,acutoff);
    maxbw = find_lde_max(a,m,n,acutoff);
    fprintf(stderr,"print_a_matrix_to_file: maximum bandwidth is %i of n: %i, %g %%\n",
	    maxbw,n,(COMP_PRECISION)maxbw/(COMP_PRECISION)(n)*100.0);
    acutoff=amax * 1.0e-6;
    fprintf(stderr,"print_a_matrix_to_file: ii) maximum is %10.5e, choosing %10.5e as cutoff for bandwidth\n",
	    amax,acutoff);
    maxbw = find_lde_max(a,m,n,acutoff);
    fprintf(stderr,"print_a_matrix_to_file: maximum bandwidth is %i of n: %i, %g %%\n",
	    maxbw,n,(COMP_PRECISION)maxbw/(COMP_PRECISION)(n)*100.0);
  }
}
/*

  print A . x = b system of equation

*/
void print_system(A_MATRIX_PREC *a, A_MATRIX_PREC *x, 
		  A_MATRIX_PREC *b,int m, int n, FILE *out)
{
  int i,j,k;
  for(i=0;i<m;i++){
    for(j=k=0;j<n;j++,k+=m)
      fprintf(out,"%12.5e ",a[i+k]);
    fprintf(out,"\t%12.5e\t%12.5e\n",x[i],b[i]);
  }
}

/*

  print a symmetric 3x3 matrix that is given in [3][3] format 
  the output is upper right half format, see below

*/
void print_sym3x3matrix(COMP_PRECISION a[3][3], FILE *out)
{
  fprintf(out,"%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	  a[INT_X][INT_X],a[INT_X][INT_Y],a[INT_X][INT_Z],a[INT_Y][INT_Y],a[INT_Y][INT_Z],a[INT_Z][INT_Z]);
}

/*

  write the contents of the interaction matrix to file

  the interaction matrix is stored in a weird way so that we
  don't have to worry about C/FORTRAN conventions

*/
void print_interaction_matrix(struct med *medium,struct flt *fault)
{
  FILE *out;
  I_MATRIX_PREC tmpflt;
  int i,j,k,l,pid;
  // get PID for name of matrix
  pid = (int)getpid();
  if(!medium->int_mat_init){
    fprintf(stderr,"print_interaction_matrix: init interaction matrix first\n");
    exit(-1);
  }
  if(medium->read_int_mat_from_file){
    fprintf(stderr,"print_interaction_matrix: will not print interaction matrix to file\n");
    fprintf(stderr,"print_interaction_matrix: since we are reading from file already\n");
    return;
  }
  sprintf(medium->hfname,"%s.%i.hdr",INTERACTION_MATRIX_FILE,pid);
  sprintf(medium->mfname,"%s.%i.dat",INTERACTION_MATRIX_FILE,pid);

  fprintf(stderr,"print_interaction_matrix: writing to \"%s\" (header) and \"%s\" (%s precision %i by %i)\n",
	  medium->hfname,medium->mfname,(medium->i_matrix_prec_size==4)?
	  ("float"):("double"),
	  medium->nmat,medium->nmat);
  // header file
  out=myopen(medium->hfname,"w");
  fprintf(out,"%i %i %i %i\n%g %g\n",medium->nmat,
	  (int)medium->i_matrix_prec_size,
	  medium->nrflt,NRMODE_DEF,medium->imean,medium->imax);
  fclose(out);
  // matrix itself
  out=myopen(medium->mfname,"w");
  if(!medium->use_sparse_storage){
    for(j=0;j<medium->nrflt;j++)
      for(k=0;k<NRMODE_DEF;k++)
	for(i=0;i<medium->nrflt;i++)
	  for(l=0;l<3;l++){
	    fwrite(&ICIM(medium->i,i,j,k,l),medium->i_matrix_prec_size,1,out);
	  }
  }else{
    for(j=0;j<medium->nrflt;j++)
      for(k=0;k<NRMODE_DEF;k++)
	for(i=0;i<medium->nrflt;i++)
	  for(l=0;l<3;l++){
	    // numerical recipsed storage
	    tmpflt=get_nrs_sparse_el(POSII(j,k),POSIJ(i,l),
				     medium->is1,medium->val);
	    fwrite(&tmpflt,medium->i_matrix_prec_size,1,out);
	  }
  }
  fclose(out);
}
/* 
   
   print the interaction matrix such that rows are for different affected faults
   and columns are for different slipping faults


 */
void print_reduced_interaction_matrix(struct med *medium,struct flt *fault)
{
  FILE *out;
  char tmpstr[STRLEN];
  I_MATRIX_PREC tmpflt;
  int i,j,k;
  COMP_PRECISION slip[3],normal,shear_strike,shear_dip;
  if(!medium->int_mat_init){
    fprintf(stderr,"print_reduced_interaction_matrix: init interaction matrix first\n");
    exit(-1);
  }
  // header
  sprintf(tmpstr,"%s.red.hdr",INTERACTION_MATRIX_FILE);
  fprintf(stderr,"print_interaction_matrix: writing to \"%s\" (header) and ",tmpstr);
  out=myopen(tmpstr,"w");
  fprintf(out,"%i %i",
	  medium->nrflt,(int)medium->i_matrix_prec_size);
  fclose(out);
  // matrix
  sprintf(tmpstr,"%s.red",INTERACTION_MATRIX_FILE);
  fprintf(stderr,"\"%s\" (%s precision %i by %i)\n",
	  tmpstr,(medium->i_matrix_prec_size==4)?("float"):("double"),
	  medium->nrflt,medium->nrflt);

  // unit slip vector, could be more complicated
  slip[STRIKE]=1.0;slip[DIP]=0.0;slip[NORMAL]=0.0;
  //
  out=myopen(tmpstr,"w");
  if(!medium->use_sparse_storage){
    for(i=0;i < medium->nrflt;i++){
      /* evaluated on fault i */
      for(j=0;j < medium->nrflt;j++){
	/* effect of slip on fault j */
   	for(normal=shear_strike=shear_dip=0.0, /* sum up stress
						  contributions from
						  all types of slip */
	      k=0;k < 3;k++){
	  normal        += ICIM(medium->i,i,j,k,NORMAL)*slip[k];
	  shear_strike  += ICIM(medium->i,i,j,k,STRIKE)*slip[k];
	  shear_dip     += ICIM(medium->i,i,j,k,DIP)   *slip[k];
	}
	tmpflt=coulomb_stress(sqrt(SQUARE(shear_strike)+SQUARE(shear_dip)),
			      (COMP_PRECISION)fault[i].mu_s,normal,
			      medium->cohesion);
	fwrite(&tmpflt,medium->i_matrix_prec_size,1,out);
      }
    }
  }else{
    for(j=0;j<medium->nrflt;j++){
      // unit slip vector, will be more complicated
      slip[STRIKE]=1.0;slip[DIP]=0.0;slip[NORMAL]=0.0;
      for(i=0;i<medium->nrflt;i++){
	for(normal=shear_strike=shear_dip=0.0,k=0;k<3;k++){
	  normal        += get_nrs_sparse_el(POSII(j,k),POSIJ(i,NORMAL),
					     medium->is1,medium->val)*slip[k];
	  shear_strike  += get_nrs_sparse_el(POSII(j,k),POSIJ(i,STRIKE),
					     medium->is1,medium->val)*slip[k];
	  shear_dip     += get_nrs_sparse_el(POSII(j,k),POSIJ(i,DIP),
					     medium->is1,medium->val)   *slip[k];
	}
	tmpflt=coulomb_stress(sqrt(SQUARE(shear_strike)+SQUARE(shear_dip)),
			      (COMP_PRECISION)fault[i].mu_s,normal,
			      medium->cohesion);
	fwrite(&tmpflt,medium->i_matrix_prec_size,1,out);
      }
    }
  }
  fclose(out);
}
/*
  print general vector in ASCII in column format
*/
void print_vector(A_MATRIX_PREC *b,int n, FILE *out)
{
  int i;
  for(i=0;i<n;i++){
    fprintf(out,"%20.15e\n",b[i]);
  }
}
void print_vector_file(A_MATRIX_PREC *b,int n, char *name,char *message)
{
  FILE *out;
  fprintf(stderr,"printing %s to %s\n",message,name);
  out = myopen(name,"w");
  print_vector(b,n,out);
  fclose(out);
}
/*
  print general vector in ASCII in row format
*/
void print_vector_row(A_MATRIX_PREC *b,int n, FILE *out)
{
  int i;
  for(i=0;i<n;i++)
    fprintf(out,"%20.15e ",b[i]);
  fprintf(out,"\n");
}
/*
  print b vector  to out and determine max abs value
*/
void print_b_vector(A_MATRIX_PREC *b,int n, FILE *out, 
		    A_MATRIX_PREC *bmax,my_boolean binary)
{
  if(binary){
    fwrite(b,sizeof(A_MATRIX_PREC),(size_t)n,out); 
    *bmax = find_max_abs_vec(b,n);
  }else{
    print_vector(b,n,out);
    *bmax = find_max_abs_vec(b,n);
  }
}
/* print 3x3 matrix */
void print_3x3_matrix(COMP_PRECISION a[3][3],FILE *out)
{
  fprintf(out,"%11g %11g %11g %11g %11g %11g %11g %11g %11g\n",
	  a[INT_X][INT_X],a[INT_X][INT_Y],a[INT_X][INT_Z],
	  a[INT_Y][INT_X],a[INT_Y][INT_Y],a[INT_Y][INT_Z],
	  a[INT_Z][INT_X],a[INT_Z][INT_Y],a[INT_Z][INT_Z]);
}
