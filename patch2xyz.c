#include "interact.h"
#include "properties.h"

// reads in patch format and writes xyz format for GMT plotting



int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  int i,comp=0;
  int opmode=PSXYZ_MODE;
  COMP_PRECISION *scalar,min,max;
  my_boolean shrink_patches=FALSE,
    read_slip = FALSE,verbose = TRUE,use_scalar = FALSE;
  FILE *in;
  if(argc > 5){
    fprintf(stderr,"%s [flt.dat] [disp_component, %i] [shrink_patches, %i] [use_scalar, %i]\n\t reads in patch format from stdin and writes GMT xyz to stdout\n",
	    argv[0],comp,(int)shrink_patches,(int)use_scalar);
    fprintf(stderr,"if flt.dat is given, will print slip or stress\n");
    fprintf(stderr,"\twhere 0: strike 1: dip 2: normal 3: vector slip\n");
    fprintf(stderr,"\twhere 10: strike 11: dip 12: normal 13: vector stress\n");
    fprintf(stderr,"\twhere 14: is Coulomb stress for strike\n");
    fprintf(stderr,"\twhere 15: is Coulomb stress for vector shear\n");
    fprintf(stderr,"if shrink_patches is set, will make patches smaller for plotting\n");
    fprintf(stderr,"if use_scalar is set, will interpret the flt.dat to be a list of scalars to use for each patch\n");
    exit(-1);
  }

  if(argc > 1){
    opmode = PSXYZ_SCALAR_MODE;
    read_slip = TRUE;
  }
  if(argc > 2){
    sscanf(argv[2],"%i",&comp);
  }
  if(argc > 3){
    sscanf(argv[3],"%i",&i);
    shrink_patches = (my_boolean)i;
  }
  if(argc > 4){
    sscanf(argv[4],"%i",&i);
    use_scalar = (my_boolean)i;
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing xyz to stdout. shrink: %i\n",
	  argv[0],shrink_patches);
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,verbose);
  
  if(read_slip){
    if(use_scalar){
      fprintf(stderr,"%s: reading scalar values from %s\n",
	      argv[0],argv[1]);
      in = myopen(argv[1],"r");
      scalar = (COMP_PRECISION *)calloc(medium->nrflt,sizeof(COMP_PRECISION));
      min = 1e20; max = -1e20;
      for(i=0;i < medium->nrflt;i++){
	if(fscanf(in,ONE_CP_FORMAT,(scalar+i)) !=1 ){
	  fprintf(stderr,"scalar read error at %i in file %s\n",i,argv[1]);
	  exit(-1);
	}
	if(scalar[i]<min)min=scalar[i];
	if(scalar[i]>max)max=scalar[i];
      }
    }else{
      fprintf(stderr,"%s: reading flt slip from %s, printing component %i\n",
	      argv[0],argv[1],comp);
      read_fltdat(argv[1],fault,medium,verbose);
      scalar = (COMP_PRECISION *)calloc(medium->nrflt,sizeof(COMP_PRECISION));
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
	case 3:
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
	  scalar[i] = sqrt(fault[i].s[STRIKE]*fault[i].s[STRIKE] + 
			   fault[i].s[DIP]*fault[i].s[DIP] +
			   fault[i].s[NORMAL]*fault[i].s[NORMAL]);
	  break;
	case 14:
	  scalar[i] = fabs(fault[i].s[STRIKE]) - STATIC_MU * fault[i].s[NORMAL];
	  break;
	case 15:
	  scalar[i] = sqrt(fault[i].s[STRIKE]*fault[i].s[STRIKE]+
			   fault[i].s[DIP]   *fault[i].s[DIP]) -
	    STATIC_MU * fault[i].s[NORMAL];
	  break;
	  
	default:
	  fprintf(stderr,"%s: slip component %i undefined\n",argv[0],comp);
	  exit(-1);
	  break;
	}
	if(scalar[i]<min)min=scalar[i];
	if(scalar[i]>max)max=scalar[i];
      }
    }
    fprintf(stderr,"%s: min: %g max: %g from scalar mode\n",argv[0],min,max);
  }

  for(i=0;i < medium->nrflt;i++){
    print_patch_geometry_and_bc(0,(fault+i),opmode,0.0,TRUE,stdout,
				shrink_patches,(scalar+i));
  }
  if(read_slip)
    free(scalar);
  return 0;
}

