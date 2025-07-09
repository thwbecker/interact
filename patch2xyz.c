#include "interact.h"
#include "properties.h"

// reads in patch format and writes xyz format for GMT plotting



int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  int i,j,comp=0;
  int opmode=PSXYZ_MODE;
  COMP_PRECISION *scalar=NULL,min,max,*tmp_scalar=NULL;
  my_boolean shrink_patches=FALSE,read_slip = FALSE,
    attempt_read_slip = FALSE,verbose = TRUE;
  int use_scalar = 0;
  FILE *in;
  medium=(struct med *)calloc(1,sizeof(struct med)); 
  if((argc > 5)||((argc>1)&&(strcmp(argv[1],"-h")==0))){
    fprintf(stderr,"%s [flt.dat] [disp_component, %i] [shrink_patches, %i] [use_scalar, %i]\n\t reads in patch format from stdin and writes GMT xyz to stdout\n",
	    argv[0],comp,(int)shrink_patches,(int)use_scalar);
    fprintf(stderr,"if flt.dat is given, will print slip or stress if file found\n");
    fprintf(stderr,"\twhere 0: strike 1: dip 2: normal 3: total shear 4: total slip\n");
    fprintf(stderr,"\twhere 10: strike 11: dip 12: normal 13: vector shear stress\n");
    fprintf(stderr,"\twhere 14: is Coulomb stress for strike\n");
    fprintf(stderr,"\twhere 15: is Coulomb stress for vector shear\n");
    fprintf(stderr,"if shrink_patches is set, will make patches smaller for plotting\n");
    fprintf(stderr,"if use_scalar > 0, will interpret the flt.dat to be a list of use_scalars to use for each patch\n");
    fprintf(stderr,"\t in this case, disp_component (0,1,...) will select column\n");
    exit(-1);
  }
  if(argc > 1){
   
    attempt_read_slip = TRUE;
  }
  if(argc > 2){
    sscanf(argv[2],"%i",&comp);
  }
  if(argc > 3){
    sscanf(argv[3],"%i",&i);
    shrink_patches = (my_boolean)i;
  }
  if(argc > 4){
    sscanf(argv[4],"%i",&use_scalar);
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing xyz to stdout. shrink: %i\n",
	  argv[0],shrink_patches);
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,verbose);
  if(use_scalar){
    if((comp >= use_scalar)||(comp<0)){
      fprintf(stderr,"%s: set to read %i scalars, but component for output is %i\n",argv[0],use_scalar,comp);
      exit(-1);
    }
  }
  if(attempt_read_slip){
    scalar = (COMP_PRECISION *)realloc(scalar,medium->nrflt*sizeof(COMP_PRECISION));
    if(use_scalar){
      /*  */
      fprintf(stderr,"%s: reading %i scalar values per row from %s\n",
	      argv[0],use_scalar,argv[1]);
      in = myopen(argv[1],"r");
      
      tmp_scalar = (COMP_PRECISION *)realloc(tmp_scalar,use_scalar*sizeof(COMP_PRECISION));
      min = 1e20; max = -1e20;
      for(i=0;i < medium->nrflt;i++){
	for(j=0;j < use_scalar;j++){
	  if(fscanf(in,ONE_CP_FORMAT,(tmp_scalar+j)) !=1 ){
	    fprintf(stderr,"%s: scalar read error at row %i col %i in file %s\n",argv[0],i,j,argv[1]);
	    exit(-1);
	  }
	  scalar[i] = tmp_scalar[comp];
	}
      }
      if(scalar[i]<min)min=scalar[i];
      if(scalar[i]>max)max=scalar[i];
      free(tmp_scalar);
      read_slip = TRUE;
      fprintf(stderr,"%s: output min: %e max: %e from scalara data component %i\n",argv[0],min,max,comp);
    }else{
      fprintf(stderr,"%s: reading flt slip from %s, printing component %i\n",
	      argv[0],argv[1],comp);
      read_slip = read_fltdat(argv[1],fault,medium,verbose);
      if(read_slip){
	min = 1e20; max = -1e20;
	for(i=0;i < medium->nrflt;i++){
	  switch(comp){
	  case 0:
	    scalar[i] = fault[i].u[STRIKE];
	    break;
	  case 1:
	    scalar[i] = fault[i].u[DIP];
	    break;
	  case 2: 
	    scalar[i] = fault[i].u[NORMAL];
	    break;
	  case 3:			/* total shear */
	    scalar[i] = hypot(fault[i].u[STRIKE],fault[i].u[DIP]);
	    break;
	  case 4:			/* total slip */
	    scalar[i] = sqrt(fault[i].u[STRIKE]*fault[i].u[STRIKE] + 
			     fault[i].u[DIP]*fault[i].u[DIP] +
			     fault[i].u[NORMAL]*fault[i].u[NORMAL]);
	    break; 
	  case 10:
	    scalar[i] = fault[i].s[STRIKE];
	    break;
	  case 11:
	    scalar[i] = fault[i].s[DIP];
	    break;
	  case 12: 
	    scalar[i] = fault[i].s[NORMAL];
	    break;
	  case 13:
	    scalar[i] = hypot(fault[i].s[STRIKE],fault[i].s[DIP]);
	    break;
	  case 14:		/* Coulomb with strike shear? */
	    scalar[i] = fault[i].s[STRIKE] - STATIC_MU * fault[i].s[NORMAL];
	    break;
	  case 15:		/* Coulomb with total shear? */
	    scalar[i] = hypot(fault[i].s[STRIKE],fault[i].s[DIP]) - STATIC_MU * fault[i].s[NORMAL];
	    break;
	    
	  default:
	    fprintf(stderr,"%s: slip component %i undefined\n",argv[0],comp);
	    exit(-1);
	    break;
	  }
	  if(scalar[i]<min)min=scalar[i];
	  if(scalar[i]>max)max=scalar[i];
	}
	fprintf(stderr,"%s: output min: %e max: %e from flt data mode %i\n",argv[0],min,max,comp);
      }
    }
  }
  if(read_slip)
    opmode = PSXYZ_SCALAR_MODE;
  fprintf(stderr,"%s: printing %i patches to stdout\n",argv[0],medium->nrflt);
  for(i=0;i < medium->nrflt;i++){
    print_patch_geometry_and_bc(0,(fault+i),opmode,0.0,TRUE,stdout,
				shrink_patches,(scalar+i));
  }
  if(read_slip)
    free(scalar);
  return 0;
}

