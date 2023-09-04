#include "interact.h"

/*
  reads in geometry file and calculates the interaction matrix
  
  $Id: calc_interaction_matrix.c,v 1.9 2011/01/07 07:19:58 becker Exp $ 

*/
#define FULL_IMAT_MODE 1
#define RED_IMAT_MODE 2
void print_help_local(char *);

int main(int argc, char **argv)
{
  struct med *medium;
  struct flt *fault;
  int mode=0,i;
  medium=(struct med *)calloc(1,sizeof(struct med));

  if(argc<2)
    print_help_local(argv[0]);
  for(i=2;i<argc;i++){
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){
      print_help_local(argv[0]);
    }else if(strcmp(argv[i],"-r")==0){
      mode=RED_IMAT_MODE;
    }else{
      fprintf(stderr,"%s: can not use parameter %s, use -h for help page\n",
	      argv[0],argv[i]);
      exit(-1);
    }}
  read_geometry(argv[1],&medium,&fault,TRUE,FALSE,FALSE,FALSE);
  medium->i_mat_cutoff = I_MAT_CUTOFF_DEF;
  calc_interaction_matrix(medium,fault,TRUE);
  switch(mode){
  case FULL_IMAT_MODE:{
    if(!medium->read_int_mat_from_file)
      print_interaction_matrix(medium,fault);
    break;
  }
  case RED_IMAT_MODE:{
    print_reduced_interaction_matrix(medium,fault);
    break;
  }}
  exit(0);
}
void print_help_local(char *filename)
{ 
  fprintf(stderr,"%s file.patch\n",filename);
  fprintf(stderr,"%s: computes the interaction matrix based on fault patches in file.patch\n",
	  filename);
  fprintf(stderr,"%s: these are in the same format as in the \"%s\" file\n",
	  filename,GEOMETRY_FILE);
  fprintf(stderr,"%s: options:\n",filename);
  fprintf(stderr,"%s: -r   computes reduced interaction matrix\n",filename);
  exit(-1);
}
