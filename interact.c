/*
  interact: model fault interactions using dislocations in a 
            halfspace

	    (C) Thorsten Becker, thwbecker@post.harvard.edu

  Okada and other dislocation routines based elastic half-space
  fault-patch interaction program

  calculates the stresses on fault patches due to slip on other
  faults. inversion of the interaction matrix is done using various
  matrix solves, see 'interact -h' for description



  BENCHMARK:
  produce geom.in file with:
  /usr/bin/time randomflt -f 100 -m 5 -n 5 > geom.in
  bc.in file reads:
  2
  0.01 0.05 70 
  -1 0 0 10
  

*/
#include "interact.h"
#ifdef USE_PETSC
#include <petscksp.h>
PetscErrorCode GenEntries(PetscInt , PetscInt , PetscInt ,const PetscInt *, const PetscInt *, PetscScalar *, void *);

/* 
   generate interaction matrix entries in a way suitable for petsc/htools

   sdim dimension
   M local m
   N local n 
   J[M] array with global indices for sources
   K[N] array with global indices for receivers

   this is modified from the ex82.c petsc example
 */
PetscErrorCode GenEntries(PetscInt sdim, PetscInt M, PetscInt N,
			  const PetscInt *J, const PetscInt *K, PetscScalar *ptr, void *kernel_ctx)
{
  PetscInt  j, k;
  COMP_PRECISION slip[3],disp[3],stress[3][3],trac[3],sval;
  int iret;
  struct interact_ctx *ictx;
  ictx = (struct interact_ctx *)kernel_ctx;
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionBeginUser;
#endif
  get_right_slip(slip,ictx->src_slip_mode,1.0);	/* strike motion */
  for (j = 0; j < M; j++) {
    for (k = 0; k < N; k++) {
      eval_green(ictx->fault[K[k]].x,(ictx->fault+J[j]),slip,disp,stress,&iret, GC_STRESS_ONLY,TRUE);
      if(iret != 0){
	fprintf(stderr,"get_entries: WARNING: i=%3i j=%3i singular\n",j,k);
	//s[STRIKE]=s[DIP]=s[NORMAL]=0.0;
	sval = 0.0;
      }else{
	resolve_force(ictx->fault[K[k]].normal,stress,trac);
	if(ictx->rec_stress_mode == STRIKE)
	  sval = dotp_3d(trac,ictx->fault[K[k]].t_strike);
	else if(ictx->rec_stress_mode == DIP)
	  sval = dotp_3d(trac,ictx->fault[K[k]].t_dip);
	else
	  sval = dotp_3d(trac,ictx->fault[K[k]].normal);
      }
      ptr[j + M * k] = sval;
    }
  }
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}
#endif




void calc_interaction_matrix(struct med *medium,struct flt *fault,
			     my_boolean write_to_screen)
{
  COMP_PRECISION sm[3][3],disp[3]={0.,0.,0.},trac[3],adv,u[3],
    std,tmpflt;
  I_MATRIX_PREC iflt,s[3];
  int iret,i,j,k,l,nzc,nsmall,m,
    sparse_matrix_size,evil_pair[2],pid;
  size_t dense_matrix_size,isize;
  FILE *out;  
  COMP_PRECISION shear_strike,shear_dip,normal;
#ifdef SUPER_DUPER_DEBUG
  write_to_screen=TRUE;
#endif
  // PID for filenames
  pid=(int)getpid();
  // set size variables
  // nr elements in dense matrix
  isize=imatrix_size(medium);
  // precision
  medium->i_matrix_prec_size = sizeof(I_MATRIX_PREC);
  // size in bytes
  dense_matrix_size = isize * medium->i_matrix_prec_size;
  //
  // size of square matrix, has to have three slip modes
  // so that s is [3][3], if NRMODE==3, identical to sqrt(isize)
  // (will check below)
  medium->nmat1 = medium->nrflt*3;  
  medium->nmat2 = medium->nrflt*medium->nrmode;

  // check if we want to switch to sparse storage for I matrix
  if(medium->use_sparse_storage){
    fprintf(stderr,"calc_interaction_matrix: sparse storage selected\n");
  }else{
    if(dense_matrix_size > IMAT_SPARSE_LIM*ONE_MEGABYTE){
      fprintf(stderr,"calc_interaction_matrix: switching to sparse storage, dense size (%g MB) larger than %g\n",
	      ((COMP_PRECISION)(dense_matrix_size))/(COMP_PRECISION)ONE_MEGABYTE,
	      IMAT_SPARSE_LIM);
      medium->use_sparse_storage=TRUE;
    }
  }
  // check if small enough to fit in memory
  if((!medium->use_sparse_storage)&&
     (dense_matrix_size < IMAT_SIZE_LIM*ONE_MEGABYTE)){// allocate memory
    if(!(medium->i=(struct icm **)
	 malloc(sizeof(struct icm *)*medium->nrflt)))
      MEMERROR("calc_interaction_matrix: 1:");
    for(i=0;i < medium->nrflt;i++)// init with zeroes
      if(! (*(medium->i+i)=(struct icm *)calloc(medium->nrflt,sizeof(struct icm)))){
	fprintf(stderr,"calc_interaction_matrix: 2: mem error, i: %i, inc s: %g MB ts: %g MB\n",
		i,(double)sizeof(struct icm)*(double)medium->nrflt/ONE_MEGABYTE,
		SQUARE((double)medium->nrflt)*(double)sizeof(struct icm)/ONE_MEGABYTE);
	exit(-1);
      }
  }else{
    medium->i=NULL;
  }
  if(!medium->i){
    if(!medium->use_sparse_storage){
      /* 
	 we don't have enough space to store I in memory 
      */
      fprintf(stderr,"calc_interaction_matrix: i matrix too big, %g MB, limit is %g\n",
	      (COMP_PRECISION)(dense_matrix_size)/ONE_MEGABYTE,IMAT_SIZE_LIM);
      fprintf(stderr,"calc_interaction_matrix: exiting, no file operation wanted\n");
      exit(-1);
      //      fprintf(stderr,"calc_interaction_matrix: will attempt to read from file during program run, have fun!\n");
    }else{// use sparse matrix storage which somehow wants N^2
      
      fprintf(stderr,"calc_interaction_matrix: preparing for sparse matrix storage, writing to \"%s\"\n",
	      INTERACTION_MATRIX_FILE);
      if(medium->nrmode == 2){
	fprintf(stderr,"calc_interaction_matrix: for now, we need square matrix and opening modes are blocked!\n");
	fprintf(stderr,"calc_interaction_matrix: hence compile with opening modes for sparse storage\n");
	exit(-1);
      }
    }
    medium->read_int_mat_from_file=TRUE;
    /*
      if we are not using an old I matrix,
      initialize file with matrix elements
    */
    if(!medium->use_old_imat){
      //
      // write header file for external programs, will
      // add mean and max later
      //

      sprintf(medium->hfname,"%s.%i.hdr",INTERACTION_MATRIX_FILE,pid);
      out=myopen(medium->hfname,"w");
      /* header */
      fprintf(out,"%i %i %i %i %i\n",medium->nmat1,(int)medium->i_matrix_prec_size,medium->nmat2,medium->nrflt,medium->nrmode);
      fclose(out);
      // initialize I matrix file with zeroes
      sprintf(medium->mfname,"%s.%i.dat",INTERACTION_MATRIX_FILE,pid);
      medium->i_mat_in = myopen(medium->mfname,"w");
      for(iflt=0.0,i=0;i < medium->nmat1;i++)
	for(j=0;j < medium->nmat2;j++)
	  if(fwrite(&iflt,medium->i_matrix_prec_size,1,medium->i_mat_in) != 1){
	    fprintf(stderr,"interact: write error while trying to initialize the \"%s\" file\n",
		    medium->mfname);
	    fprintf(stderr,"interact: not enough space? it's %g MB large...\n",
		    (COMP_PRECISION)(medium->nmat1*medium->nmat2*medium->i_matrix_prec_size)/
		    ONE_MEGABYTE);
	    exit(-1);
	  }
      rewind(medium->i_mat_in);
    }
  }else{
    fprintf(stderr,
	    "calc_interaction_matrix: interaction matrix in memory (%g MB)\n",
	    (COMP_PRECISION)(dense_matrix_size)/ONE_MEGABYTE);
  }
  if(!medium->use_old_imat){
    /* 
       produce new interaction matrix

       calculate the effect of unity displacements in 
       direction k at fault j on stresses of component l on fault i 
       
       the funny ordering is because of the way we store things 
       on file, don't change it!
    */
    // initialize some stat quantities
    medium->imean = medium->imax = std = 0.0;
    nzc=0;
    for(j=0;j < medium->nrflt;j++){// loop over all rupturing faults
      for(k=0;k < medium->nrmode;k++){// loop over all rupture modes
	// in the i.dat matrix, the i' coordinate is j*NRMODE+k
	// and the                  j' coordinate is i*3+l
	/* 
	   different displacement components,
	   strike, dip, and tension (if not suppressed)
	*/
	get_right_slip(disp,k,1.0);
	for(i=0;i<medium->nrflt;i++){// loop over observing faults
	  // evaluate the 'Green's function'
	  eval_green(fault[i].x,(fault+j),disp,u,sm,&iret, GC_STRESS_ONLY,TRUE);
	  if(iret != 0){
	    /* 
	       Green's function is singular at this point,
	       most likely because of overlapping fault patches
	    */
	    fprintf(stderr,"calc_interaction_matrix: WARNING: i=%3i j=%3i k=%1i singular\n",
		    i,j,k);
	    fprintf(stderr,"calc_interaction_matrix: do patches i=%3i and j=%3i overlap?\n",i,j);
	    s[STRIKE]=s[DIP]=s[NORMAL]=0.0;
	  }else if(medium->suppress_interactions &&
		   (fault[j].group != fault[i].group) && 
		   (fault[j].group != MASTER_FAULT_GROUP)){
	    /*
	      if we want to suppress interactions, simply
	      zero out here. this is SLOW, but simple
	      
	      only interaction within the patch's group 
	      and effect of magic group will be considered
	    */
	    s[STRIKE]=s[DIP]=s[NORMAL]=0.0;
	  }else{
	    /* 
	       resolve the local stress tensor in terms of 
	       stresses in normal and tangential directions 
	       POS(affected fault, affecting fault,
	       affecting component, affected component) 
	    */
	    resolve_force(fault[i].normal,sm,trac);
	    s[STRIKE]=(I_MATRIX_PREC)
	      dotp_3d(trac,fault[i].t_strike);
	    s[DIP]=(I_MATRIX_PREC)
	      dotp_3d(trac,fault[i].t_dip);
	    s[NORMAL]=(I_MATRIX_PREC)
	      dotp_3d(trac,fault[i].normal);
	  }
	  /* 
	     determine absolute maximum and add for mean
	  */
	  for(l=0;l<3;l++){
	    if((tmpflt=fabs(s[l])) >= EPS_COMP_PREC){
	      // encountered non-zero entry
	      medium->imean += tmpflt;
	      std += SQUARE(tmpflt);// for standard dev
	      if(tmpflt > medium->imax)
		medium->imax=tmpflt;
	      nzc++;
	    }else{
	      s[l]=0.0;
	    }
	    /*
	      output or assign to matrix
	    */
	    if(medium->read_int_mat_from_file){// write to file
	      if(fwrite((s+l),sizeof(I_MATRIX_PREC),1,
			medium->i_mat_in) != 1){
		fprintf(stderr,"interact: write error at i/j/k/l %i/%i/%i/%i even though initialization worked?!\n",
			i,j,k,l);
		exit(-1);
	      }
	    }else{// use i matrix in memory
	      ICIM(medium->i,i,j,k,l)=(I_MATRIX_PREC)s[l];
	    }
	  }
	}
      }
      //fprintf(stderr,"calc_interaction_matrix: working on fault patch %5i\r",j);
    }
    //fprintf(stderr,"\n");
    /*
      I matrix assembly done
    */
#define cn ((COMP_PRECISION)nzc)
    if(medium->suppress_interactions)
      fprintf(stderr,"calc_interaction_matrix: ALL INTERACTIONS SUPPRESSED (but with own group and effect of group %i)\n",
	      MASTER_FAULT_GROUP);
    if(nzc){
      // stddev
      std=sqrt((cn*std - medium->imean*medium->imean) / 
	       (cn*(cn-1)));
      medium->imean /= cn;
    }
#undef cn
    if(medium->read_int_mat_from_file){
      /*
	if we are writing to i.dat file,
	add information to header file
      */
      out=myopen(medium->hfname,"a");
      fprintf(out,"%g %g\n",medium->imean,medium->imax);
      fclose(out);
    }
    /*
      
      end of matrix calculation part
      
    */
  }else{
    std=-999.0;
    /* 

       start reading old I matrix from file

    */
    sprintf(medium->hfname,"%s.hdr",INTERACTION_MATRIX_FILE);
    sprintf(medium->mfname,"%s.dat",INTERACTION_MATRIX_FILE);
    fprintf(stderr,"calc_interaction_matrix: WARNING: using old I matrix in \"%s\" and \"%s\"\n",
	    medium->mfname,medium->hfname);
    out=myopen(medium->hfname,"r");
    if(fscanf(out,"%i %i %i %i %i\n",&i,&j,&k,&l,&m)!=5)
      READ_ERROR(medium->hfname);
    if((i != medium->nmat1)||
       (j != medium->i_matrix_prec_size)||
       (k != medium->nmat2)||
       (l != medium->nrflt)||(l != medium->nrmode)){
      fprintf(stderr,"calc_interaction_matrix: mismatch of either nmat1 (%i),prec_size (%i), namt2 (%i), nrflt (%i), or nrmode (%i)\n",
	      i,j,k,l,m);
      exit(-1);
    }
    if(fscanf(out,TWO_IP_FORMAT,&medium->imean,&medium->imax)!=2)
      READ_ERROR(medium->hfname);
    fclose(out);
    medium->i_mat_in=myopen(medium->mfname,"r");
    if(!medium->read_int_mat_from_file){
      fprintf(stderr,"calc_interaction_matrix: reading into memory\n");
      for(j=0;j < medium->nrflt;j++)
	for(k=0;k < medium->nrmode;k++)
	  for(i=0;i < medium->nrflt;i++)
	    for(l=0;l<3;l++)
	      if(fread(&ICIM(medium->i,i,j,k,l),medium->i_matrix_prec_size,1,out)!=1){
		fprintf(stderr,"calc_interaction_matrix: read error i/j/k/l %i/%i/%i/%i\n",
			i,j,k,l);
		exit(-1);
	      }
    }
  }
  i=(int)sqrt((COMP_PRECISION)isize);
  fprintf(stderr,"calc_interaction_matrix: dense I matrix size %s %g MB, %i by %i in %s precision\n",
	  medium->read_int_mat_from_file?"would be":"is",
	  ((COMP_PRECISION)dense_matrix_size)/ONE_MEGABYTE,
	  i,i,(medium->i_matrix_prec_size==8)?"double":"single");
  if(medium->read_int_mat_from_file)
    fprintf(stderr,"calc_interaction_matrix: from file: max of abs values is %g, mean of non-zero amplitudes is %g\n",
	    medium->imax,medium->imean);
  else
    fprintf(stderr,"calc_interaction_matrix: max of abs values is %g, mean/stddev  of non-zero amplitudes is %g/%g\n",
	    medium->imax,medium->imean,std);
  if(medium->read_int_mat_from_file){
    // reopen in read mode
    if(!medium->use_old_imat){
      if(!freopen(medium->mfname,"r",medium->i_mat_in)){
	fprintf(stderr,"calc_interaction_matrix: issue with %s, could not reopen\n",medium->mfname);
	exit(-1);
      }
      fprintf(stderr,"calc_interaction_matrix: written to %s\n",
	      medium->mfname);
    }
    // deal with sparse matrix mode
    if(medium->use_sparse_storage){
      // intialize the spare matrix storage arrays
      fprintf(stderr,"calc_interaction_matrix: initializing sparse storage scheme from file\n");
      fprintf(stderr,"calc_interaction_matrix: cutoff: %g, %g%% of max(abs(I)): %g\n",
	      medium->i_mat_cutoff,medium->i_mat_cutoff/medium->imax*100.0,medium->imax);
      if(medium->nmat1 != medium->nmat2){
	fprintf(stderr,"trying to call create_nrs_sparse_from_file but nmat1 %i nmat2 %i\n",
		medium->nmat1,medium->nmat2);
	exit(-1);
      }
      // use numerical recipes schemes, has to be n by n
      create_nrs_sparse_from_file(medium->nmat1,medium->i_mat_cutoff,
				  &medium->is1,&medium->val,
				  medium->i_mat_in);
      sparse_matrix_size=(int)medium->is1[medium->is1[0]-1];
      sparse_matrix_size *= (sizeof(I_MATRIX_PREC)+sizeof(unsigned int));
      fprintf(stderr,"calc_interaction_matrix: sparse storage vectors fill %g MB\n",
	      (COMP_PRECISION)sparse_matrix_size/ONE_MEGABYTE);
      adv=100.0*(COMP_PRECISION)sparse_matrix_size/(COMP_PRECISION)(dense_matrix_size);
      fprintf(stderr,"calc_interaction_matrix: this is %g%% of dense storage\n",adv);
      if(adv>80)
	fprintf(stderr,"calc_interaction_matrix: maybe sparse storage wasn't such a good idea\n");
#ifdef SUPER_DEBUG
      fprintf(stderr,"calc_interaction_matrix: extracting compressed storage to i2.dat\n");
      // check if the transformation worked
      out=myopen("i2.dat","w");
      for(i=0;i<medium->nmat1;i++)
	for(j=0;j<medium->nmat2;j++){
	  iflt = get_nrs_sparse_el(i,j,medium->is1,medium->val);
	  fwrite(&iflt,sizeof(I_MATRIX_PREC),1,out);
	}
      fclose(out);
#endif
    }
  }else{
    if(write_to_screen && medium->nrflt>100)
      write_to_screen=FALSE;
    // check for small entries
    nsmall=0;
    for(i=0;i < medium->nrflt;i++)
      for(j=0;j < medium->nrflt;j++)
	for(k=0;k < medium->nrmode;k++){
	  if(write_to_screen){
	    /* if we only have a couple of patches, print out interaction coefficients */
	    get_right_slip(disp,k,1.0);
	    // obtain shear stress and normal stress change
	    for(normal=shear_strike=shear_dip=0.0,m=0;m<3;m++){
	      normal        += ICIM(medium->i,i,j,m,NORMAL)*disp[m];
	      shear_strike  += ICIM(medium->i,i,j,m,STRIKE)*disp[m];
	      shear_dip     += ICIM(medium->i,i,j,m,DIP)   *disp[m];
	    }
#ifdef DEBUG
	    fprintf(stderr,
		    "observ: %5i slip: %5i u_s/d/n:(%2g/%2g/%2g) s_s/d/n/C:(%12.5e/%12.5e/%12.5e/%12.5e)\n",
		    i,j,disp[STRIKE],disp[DIP],disp[NORMAL],
		    ICIM(medium->i,i,j,k,STRIKE),ICIM(medium->i,i,j,k,DIP),
		    ICIM(medium->i,i,j,k,NORMAL),
		    coulomb_stress(sqrt(SQUARE(shear_strike)+SQUARE(shear_dip)),(COMP_PRECISION)fault[i].mu_s,
				   normal,medium->cohesion));
#endif
	  }
	  /* 
	     check if self-interaction is zero (shouldn't happen) 
	  */
	  if(i==j)
	    if(fabs(ICIM(medium->i,i,i,k,k)) < EPS_COMP_PREC){
	      fprintf(stderr,"calc_interaction_matrix: ERROR: slip type %i on fault %i has no effect on strike type %i (%g)\n",
		      k,i,k,ICIM(medium->i,i,i,STRIKE,STRIKE));
	      fprintf(stderr,"calc_interaction_matrix: this is no good. maybe you tried a 3-D slip mode with 2-D patches?\n");
	      exit(-1);
	    }
	  for(l=0;l<3;l++)
	    if(fabs(ICIM(medium->i,i,j,k,l)) < medium->i_mat_cutoff)
	      nsmall++;
	}
    fprintf(stderr,"calc_interaction_matrix: %i elements (%5.2f%%) are below cutoff of %g\n",
	    nsmall,((COMP_PRECISION)nsmall/
		    (COMP_PRECISION)isize)*100.0,medium->i_mat_cutoff);
    fprintf(stderr,"calc_interaction_matrix: this corresponds to a density of %g (should be < 0.33 for sparse storage)\n",
	    (COMP_PRECISION)(isize-nsmall)/(COMP_PRECISION)isize);
  }
  if(medium->check_for_interaction_feedback){
    fprintf(stderr,"calc_interaction_matrix: checking for positive Coulomb feedback loops\n");
    check_coulomb_stress_feedback(medium->nrflt,0,fault,medium,
				  TRUE,select_i_coeff_calc_mode(medium),
				  FALSE,evil_pair,-1.0);
  }
  medium->int_mat_init=TRUE;
}

/*

  calculate the stress_component_l effect on fault_i when fault_j slips with unity
  displacement in slip_mode_k

 */
COMP_PRECISION interaction_coefficient(int i, int j, int k, int l,
				       struct flt *fault,int *iret)
{
  COMP_PRECISION disp[3],fac,sm[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},u[3]={0.,0.,0.};
  static COMP_PRECISION trac[3];
  static int old_indices[4]={-1,0,0,0};
  // check if we have calculated the traction vector
  // for observing fault i, slipping fault j, and 
  // slip mode k before. if so, only need to project
  if((old_indices[0] == -1)||// first call
     (old_indices[0] != i)||// or simply different indices
     (old_indices[1] != j)||
     (old_indices[2] != k)){
    // obtain appropriate slip vector for slip mode k
    get_right_slip(disp,k,1.0);
    /* obtain the stress vector at fault i (centroid) when fault j
       slips with disp[] */
    eval_green(fault[i].x,(fault+j),disp,u,sm,iret, GC_STRESS_ONLY,TRUE);
    if(! *iret){// if not singular,
      // obtain the traction vector for i,j,k 
      resolve_force(fault[i].normal,sm,trac);
    }else{
      fprintf(stderr,"interaction_coefficient: NAN returned\n");
    }
    old_indices[0]=i;old_indices[1]=j;old_indices[2]=k;
    // use fourth index for singular code
    old_indices[3]= *iret;
  }else{// same indices, only restore iret, trac was kept
    *iret=old_indices[3];
  }
  if(! *iret){// project if not singular
    switch(l){
    case STRIKE:{
      fac=  dotp_3d(trac,fault[i].t_strike);
      break;
    }
    case DIP:{
      fac=  dotp_3d(trac,fault[i].t_dip);
      break;
    }
    case NORMAL:{
      fac=  dotp_3d(trac,fault[i].normal);
      break;
    }
    default:{
      fprintf(stderr,"interaction_coefficient: affected mode %i undefined\n",l);
      exit(-1);
    }}
  }else{// return 0 if singular
    fprintf(stderr,"interaction_coefficient: WARNING: i/j/k/l %i/%i/%i/%i is singular\n",
	    i,j,k,l);
    fac=0.0;
  }
  return(fac);
}

/* 
   this is an external routine since sign issues for 
   constrained faults were once treated here, now that's
   taken care of in assign_right_hand_side
   
   anyways, creates a unity-type slip vector in the dir direction
   using `unity' for 1, ie. one can also use some other value

*/				       
void get_right_slip(COMP_PRECISION *disp,int dir,
		    COMP_PRECISION unity)
{
  disp[NORMAL]=disp[STRIKE]=disp[DIP]=0.0;
  disp[dir] = unity;
}

/* 
   read the interaction coefficient from a binary file 
   using the normal i,j,k,l indices
*/
I_MATRIX_PREC ic_from_file(int i, int j, int k, int l,
			   struct med *medium)
{
  return(aij_from_file(POSII(j,k),POSIJ(i,l), medium->nmat2, medium->i_mat_in));
}
/* 
   read the interaction coefficient from a binary file 
   using i,j indices where i and j run from
   0 .. medium->nmat2-1
*/
I_MATRIX_PREC aij_from_file(int i, int j, int n, FILE *in)
{
  I_MATRIX_PREC tmp;
  long offset;
  int error;
  static long old_offset=-LONG_MAX;
  offset  = (long)((i*n+j) * sizeof(I_MATRIX_PREC));
  // do we really have to move the head?
  if((old_offset == -LONG_MAX)||
     (offset != old_offset+sizeof(I_MATRIX_PREC)))
    fseek(in, offset, SEEK_SET);
  // read in value
  if((error = (int)fread(&tmp,sizeof(I_MATRIX_PREC),1,in)) != 1){
    fprintf(stderr,"aij_from_file: read error %i, i/j/n %i/%i/%i file not opened?\n",
	    error,i,j,n);
    exit(-1);
  }
  old_offset=offset;
  return(tmp);
}
/* 
   select an operational mode for obtaining the interaction
   coeffficients
 */
int select_i_coeff_calc_mode(struct med *medium)
{
  int mode;
  if(medium->int_mat_init){
    if(!medium->read_int_mat_from_file)
      mode = I_MAT_IN_MEMORY;
    else{
      if(!medium->use_sparse_storage)
	mode = I_MAT_ON_FILE;
      else
	mode = SPARSE_I_MAT;
    }
  }else{
    mode = CALC_I_COEFF_NOW;
  }
  return mode;
}
// calculate the size of the I matrix
size_t imatrix_size(struct med *medium)
{
  return((size_t)(SQUARE(medium->nrflt)*medium->nrmode*3));
}
