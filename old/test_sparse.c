#include "interact.h"

// reads in i matrix from file and tests sparse matrix 
// storage 

int main(int argc, char **argv)
{
  struct med *medium;
  char tmpstr[STRLEN];
  int nrmode,isize,n,i,j,nsize;
  I_MATRIX_PREC thres,thres2,*imat,tmpflt;
  my_boolean in_memory=FALSE;
  FILE *out;
  /*
    read in header for binary i matrix file 
  */

  medium=(struct med *)calloc(1,sizeof(struct med));
  sprintf(tmpstr,"%s.hdr",INTERACTION_MATRIX_FILE);
  medium->i_mat_in=myopen(tmpstr,"r");
  fscanf(medium->i_mat_in,"%i %i %i %i\n",
	 &medium->nmat,
	 &medium->i_matrix_prec_size,
	 &medium->nrflt,&nrmode);
  fscanf(medium->i_mat_in,TWO_IP_FORMAT,
	 &medium->imean,&medium->imax);
  fclose(medium->i_mat_in);
  
  fprintf(stderr,"%s: I matrix in file \"%s\", mean: %g, max: %g\n",
	  argv[0],INTERACTION_MATRIX_FILE,medium->imean,medium->imax);
  fprintf(stderr,"%s: size of matrix is %i by %i, %s precision, %i bytes\n",
	  argv[0],medium->nmat,medium->nmat,
	  medium->i_matrix_prec_size==4?"single":"double",
	  isize=SQUARE(medium->nmat)*medium->i_matrix_prec_size);
  
  
  medium->i_mat_in=myopen(INTERACTION_MATRIX_FILE,"r");
  if(in_memory){
    /*
      hold the matrix in memory for debugging purposes
    */
    imat=(I_MATRIX_PREC *)malloc(sizeof(I_MATRIX_PREC)*SQUARE(medium->nmat));
    if(!imat)MEMERROR("main");
    fread (imat, sizeof(I_MATRIX_PREC), SQUARE(medium->nmat), medium->i_mat_in);
    fclose(medium->i_mat_in);
    fprintf(stderr,"%s: keeping matrix in memory\n",argv[0]);
  }else{
    fprintf(stderr,"%s: reading from file without storing in memory\n",
	    argv[0]);
  }
  if(in_memory){// create a log file for density vs. cutoff value
    out=myopen("shrink.log","w");
    fprintf(out,"# thres/max(abs(x)) sparse_size/full_size\n");
    fclose(out);
  }
  /*
    
    loop though different cutoff values if we 
    have the dense matrix in memory
    
   */
  for(thres=1e-10;thres<=(in_memory?1e-1:1e-10);thres *= 2.0){
    fprintf(stderr,"%s: using threshold %g%% of max abs value of dense matrix\n",
	    argv[0],thres*100.0);
    create_nrs_sparse(medium,thres*medium->imax,imat,in_memory);
    // sparse storage vector length
    n=(int)medium->is1[medium->is1[0]-1];
    nsize=n*(sizeof(I_MATRIX_PREC)+sizeof(unsigned int));
    fprintf(stderr,"%s: sparse array length is %i, bytes %i, %g%% of full\n",
	    argv[0],n,nsize,
	    (COMP_PRECISION)nsize/(COMP_PRECISION)isize*100.0);
    if(in_memory){//write to log file
      out=myopen("shrink.log","a");
      fprintf(out,"%g %g\n",thres,
	      (COMP_PRECISION)nsize/(COMP_PRECISION)isize);
      fclose(out);
    }
    if(thres == 1e-10){// write a copy to the i2.dat file to
      // check if the transformation worked
      out=myopen("i2.dat","w");
      for(i=0;i<medium->nmat;i++)
	for(j=0;j<medium->nmat;j++){
	  tmpflt=get_nrs_sparse_el(i,j,medium->is1,medium->val);
	  fwrite(&tmpflt,sizeof(I_MATRIX_PREC),1,out);
	}
      fclose(out);
      fprintf(stderr,"%s: wrote reconverted matrix to i2.dat\n",
	      argv[0]);
    }
  }
  if(in_memory)
    fprintf(stderr,"%s: shrinking log in shrink.log\n",argv[0]);	  
}
