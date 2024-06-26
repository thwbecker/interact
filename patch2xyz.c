#include "interact.h"
#include "properties.h"

// reads in patch format and writes xyz format for GMT plotting



int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  int i,comp=0;
  int opmode=PSXYZ_MODE;
  COMP_PRECISION *scalar=NULL,min,max;
  my_boolean shrink_patches=FALSE,read_slip,
    attempt_read_slip = FALSE,verbose = FALSE,use_scalar = FALSE;
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
    fprintf(stderr,"if use_scalar is set, will interpret the flt.dat to be a list of scalars to use for each patch\n");
    exit(-1);
  }
  if(argc > 1){
    opmode = PSXYZ_SCALAR_MODE;
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
    sscanf(argv[4],"%i",&i);
    use_scalar = (my_boolean)i;
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing xyz to stdout. shrink: %i\n",
	  argv[0],shrink_patches);
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,verbose);

  if(attempt_read_slip)
    read_slip = read_fltdat(argv[1],fault,medium,verbose);
  else
    read_slip = FALSE;
  if(read_slip){
    fprintf(stderr,"%s: reading flt slip from %s, printing component %i\n",
	    argv[0],argv[1],comp);
    
    if(use_scalar){
      fprintf(stderr,"%s: reading scalar values from %s\n",
	      argv[0],argv[1]);
      in = myopen(argv[1],"r");
      scalar = (COMP_PRECISION *)realloc(scalar,medium->nrflt*sizeof(COMP_PRECISION));
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
    }
    fprintf(stderr,"%s: output min: %e max: %e from scalar mode %i\n",argv[0],min,max,comp);
  }

  for(i=0;i < medium->nrflt;i++){
    print_patch_geometry_and_bc(0,(fault+i),opmode,0.0,TRUE,stdout,
				shrink_patches,(scalar+i));
  }
  if(read_slip)
    free(scalar);
  return 0;
}

