#include "interact.h"

/*

  reads in geometry file and location file and calculate Green's functions 

*/
void print_help_local2(char *,int,int);

int main(int argc, char **argv)
{
  struct med *medium;
  struct flt *fault;
  FILE *in;
  int disp_dim=2;		/* 2: only horizontal displacements 3: including vertical */
  int slip_modes=2;		/* 2: strike and dip 1: strike -1: dip */
  COMP_PRECISION u[3],sm[3][3],disp[3],x[3];
  int i,j,k,l,iret,opmode;
  if(argc<3)
    print_help_local2(argv[0],disp_dim,slip_modes);
  if(argc>3)
    sscanf(argv[3],"%i",&disp_dim);

  if(argc>4)
    sscanf(argv[4],"%i",&slip_modes);
  fprintf(stderr,"%s: using %i displacement dimensions and %i slip modes\n",
	  argv[0],disp_dim,slip_modes);
  
  if((disp_dim < 1)||(disp_dim >3) ||(slip_modes <-1)||(slip_modes>3)){
    fprintf(stderr,"%s: bad parameters\n",argv[0]);
    exit(-1);
  }
  read_geometry(argv[1],&medium,&fault,TRUE,FALSE,FALSE,TRUE);
  medium->i_mat_cutoff = I_MAT_CUTOFF_DEF;
  medium->olocnr=0;
  medium->xoloc = (float *)malloc(sizeof(float)*2);
  in = myopen(argv[2],"r");
  while(fscanf(in,"%f %f",
	       (medium->xoloc+medium->olocnr*2+X),
	       (medium->xoloc+medium->olocnr*2+Y))==2){
    medium->olocnr++;
    medium->xoloc = (float *)realloc(medium->xoloc,sizeof(float)*(medium->olocnr+1)*2);
  }
  fclose(in);
  fprintf(stderr,"%s: read %i locations from %s, matrix is %i (data) by %i (model)\n",
	  argv[0],medium->olocnr,argv[2],medium->olocnr*disp_dim,
	  abs(slip_modes)*medium->nrflt);
  if(0){
    /* print relative fault location */
    for(j=0;j < medium->nrflt;j++)
      fprintf(stderr,"%i x %11g y %11g z %11g l %11g w %11g\n",
	      j+1,fault[j].x[X],fault[j].x[Y],fault[j].x[Z],fault[j].pos[X],fault[j].pos[Y]);

  }
  if(1){			/* do in memory, should be faster
				   (less evaluations of Okada */
    calc_design_matrix(medium,fault,disp_dim,slip_modes);
    print_design_matrix(medium,fault,disp_dim,slip_modes,stdout);
  }else{			/* compute directly, for debugging */
    if(slip_modes == -1){
      opmode = 2;
    }else{
      opmode = 1;
    }
    slip_modes = abs(slip_modes);
    x[Z] = 0;
    for(i=0;i < medium->olocnr;i++){
      x[X] = medium->xoloc[i*2+X];
      x[Y] = medium->xoloc[i*2+Y];
      for(l=0;l < disp_dim;l++){
	for(j=0;j < medium->nrflt;j++){
	  for(k=0;k < slip_modes;k++){
	    if(opmode == 1)
	      get_right_slip(disp,k,1.0);
	    else
	      get_right_slip(disp,DIP,1.0);
	    eval_rectangle(x,(fault+j),disp,u,sm,&iret);
	    fprintf(stdout,"%g ",u[l]);
	  }
	}
	fprintf(stdout,"\n");
      }
    }
  }
  free(medium->val);
  exit(0);
}
void print_help_local2(char *filename, int disp_dim, int slip_modes)
{ 
  fprintf(stderr,"%s file1.patch file2.xy [disp_dim, %i] [slip_modes, %i]\n",
	  filename,disp_dim,slip_modes);
  fprintf(stderr,"%s: computes Greens function matrix based on fault patches in file1.patch and observations in file2.xy\n",
	  filename);
  exit(-1);
}
void calc_design_matrix(struct med *medium,struct flt *fault,int disp_dim, int slip_modes)
{
  COMP_PRECISION u[3],sm[3][3],disp[3],x[3];
  int i,j,k,l,iret,n,nmod,ndata;
  int opmode;
  if(slip_modes == -1){
    opmode = 2;
  }else{
    opmode = 1;
  }
  slip_modes = abs(slip_modes);
  
  nmod = slip_modes * medium->nrflt;
  ndata = medium->olocnr * disp_dim;
  n = nmod * ndata;
  
  medium->val = (I_MATRIX_PREC *)calloc(sizeof(I_MATRIX_PREC),n);
  if(!medium->val){fprintf(stderr," mem error, %i\n",n);exit(-1);}

  x[Z] = 0;			/* always on surface */
  disp[NORMAL]=0;		/* never normal slip */
  
  for(i=0;i < medium->olocnr;i++){
    x[X] = medium->xoloc[i*2+X];
    x[Y] = medium->xoloc[i*2+Y];
    for(j=0;j < medium->nrflt;j++){
      for(k=0;k < slip_modes;k++){
	if(opmode == 1)
	  get_right_slip(disp,k,1.0);
	else
	  get_right_slip(disp,DIP,1.0);
	eval_rectangle(x,(fault+j),disp,u,sm,&iret);
	for(l=0;l<disp_dim;l++){
	  medium->val[((i*disp_dim)+l)*nmod + j*slip_modes +k] = u[l];
	}
      }	/* slip */
    }	/* fault */
  }	/* observation */

}
void print_design_matrix(struct med *medium,struct flt *fault,int disp_dim,
			 int slip_modes, FILE *out)
{
  int i,j,k,l,nmod,ndata,ic,jc;
  slip_modes = abs(slip_modes);

  nmod = slip_modes * medium->nrflt; /* model parameters */
  ndata = medium->olocnr * disp_dim;

  for(i=ic=0;i < medium->olocnr;i++){
    for(l=0;l<disp_dim;l++,ic++){
      for(jc=j=0;j < medium->nrflt;j++){
	for(k=0;k < slip_modes;k++,jc++){
	  fprintf(out,"%g ",medium->val[((i*disp_dim)+l)*nmod + j*slip_modes +k]);
	}
      }
      if(jc != nmod){
	fprintf(stderr,"mismatch 1 %i %i\n",jc,nmod); /* kind of a
							 silly check,
							 but
							 debugging */
	exit(-1);
      }
      fprintf(out,"\n");
    }
  }
  if(ic != ndata){
    fprintf(stderr,"mismatch 2 %i %i\n",ic,ndata);
    exit(-1);
  }
}
