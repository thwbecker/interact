#include "interact.h"
#include "properties.h"

//
// reads in patch quad format and writes patch triangle format (i.e. doubles patches)
//


int main(int argc, char **argv)
{
#ifndef ALLOW_NON_3DQUAD_GEOM
  // simply for ease of storage reasons
  fprintf(stderr,"%s needs ALLOW_NON_3DQUAD_GEOM to be set\n",argv[0]);exit(-1);
#else

  struct flt *qfault,tfault[1];
  struct med *medium;
  long int seed = -1;
  int i,j,k,l,ntflt,con[2][3];
  my_boolean left,verbose = FALSE;
  COMP_PRECISION corner[4][3],*dummy=NULL;
  int assign_mode = 0;	/* 0: left 1: right 2: randomg */
  RAND_GEN(&seed);
  medium=(struct med *)calloc(1,sizeof(struct med)); 
  if(argc > 2){
    fprintf(stderr,"%s [assign_mode, %i]\n\t reads in patch quad format from stdin and writes patch tri format to stdout\nassign_mode: 0 = left, 1 = right, 2 = random\n\n",
	    argv[0],assign_mode);
    exit(-1);
  }
  if(argc>1)
    sscanf(argv[1],"%i",&assign_mode);

  fprintf(stderr,"%s: reading patch quad format from stdin, writing patch tri  to stdout, assign_mode: %i\n",argv[0],assign_mode);
  read_geometry("stdin",&medium,&qfault,FALSE,FALSE,FALSE,verbose);
  ntflt = medium->nrflt*2;
  tfault[0].xt = (COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*9);
  tfault[0].type = TRIANGULAR;
  
  for(i=0;i < medium->nrflt;i++){
    if(qfault[i].type != RECTANGULAR_PATCH){
      fprintf(stderr,"%s: patch %i is not quad on input\n",argv[0],i+1);
      exit(-1);
    }
    switch(assign_mode){
    case 0:
      left = TRUE;
      break;
    case 1:
      left = FALSE;
      break;
    case 2:
      if(myrand(&seed) < 0.5)
	left = TRUE;
      else
	left = FALSE;
      break;
    default:
      fprintf(stderr,"%s: assign mode error %i\n",argv[0],assign_mode);
      exit(-1);
      break;
    }
    calculate_bloated_corners(corner,(qfault+i),1.0);
    if(left){
      con[0][0] = 0; con[0][1] = 1; con[0][2] = 3;
      con[1][0] = 1; con[1][1] = 2; con[1][2] = 3;
    }else{
      con[0][0] = 0; con[0][1] = 1; con[0][2] = 2;
      con[1][0] = 0; con[1][1] = 2; con[1][2] = 3;
    }
    tfault[0].group =  qfault[i].group;
    /* assign quad-like strike and dip to triangular fault patch */
    tfault[0].strike = qfault[i].strike;
    tfault[0].dip =    qfault[i].dip;
    
    for(j=0;j < 2;j++){		/* one quad = two triangle */
      for(k=0;k<3;k++){		/* node loop */
	for(l=0;l<3;l++)	/* dimension loop */
	  tfault[0].xt[3*k+l] = corner[con[j][k]][l];
      }
      calc_centroid_tri(tfault[0].xt,tfault[0].x);
      tfault[0].area = triangle_area(tfault[0].xt);
      tfault[0].l = tfault[0].w = sqrt( tfault[0].area);
      print_patch_geometry_and_bc(0,tfault,PATCH_OUT_MODE,0.0,
				  0,stdout,FALSE,dummy);
    }
  }
  fprintf(stderr,"%s: written %i triangular patches to stdout\n",argv[0],ntflt);
  free(qfault);free(medium);
  return 0;
#endif
}

