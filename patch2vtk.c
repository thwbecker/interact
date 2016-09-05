#include "interact.h"
#include "properties.h"
//
// reads in patch format and writes VTK format
//

int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  int i,j,k;
  my_boolean shrink_patches=FALSE,read_slip=FALSE,verbose=FALSE;
  COMP_PRECISION leeway,corner[4][3],alpha,sin_dip,cos_dip,u[3];
  if(argc > 1)
    read_slip = TRUE;
  if(argc > 2){
    sscanf(argv[2],"%i",&i);
    shrink_patches = (my_boolean)i;
  }
  if(argc > 4){
    fprintf(stderr,"%s [flt.dat] [shrink_patches, 0] \n\t reads in patch format from stdin and writes VTK format to stdout\n",
	    argv[0]);
    fprintf(stderr,"if an argument is given, will assume it is a flt.dat type output file and assign values for coloring\n");
    fprintf(stderr,"if shrink_patches is set, will make patches smaller for plotting\n");
    exit(-1);
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing VTK to stdout. shrink: %i \n",
	  argv[0],shrink_patches);
  
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,verbose);
  if(read_slip)
    read_fltdat(argv[1],fault,medium,verbose);
  if(shrink_patches)
    leeway = 0.9;
  else
    leeway = 1.0;

  printf("# vtk DataFile Version 2.0\n");
  printf("from patch2vtk\n");
  printf("ASCII\n");
  printf("DATASET UNSTRUCTURED_GRID\n");
  printf("POINTS %i float\n",medium->nrflt*4);
  for(i=0;i < medium->nrflt;i++){
#ifdef ALLOW_NON_3DQUAD_GEOM
    if(fault[i].type != TRIANGULAR){
      // normal (rectangular, 2-D, or point source)
      alpha=90.0-(COMP_PRECISION)fault[i].strike;
      my_sincos_deg(&fault[i].sin_alpha,&fault[i].cos_alpha,(COMP_PRECISION)alpha);
      my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[i].dip);
      calc_base_vecs(fault[i].t_strike,fault[i].normal,fault[i].t_dip,
		     fault[i].sin_alpha,fault[i].cos_alpha,sin_dip,cos_dip);
    }else{// triangular element
      get_alpha_dip_tri_gh((fault+i)->xt,&(fault+i)->sin_alpha,
			   &(fault+i)->cos_alpha,&tmpdbl,&(fault+i)->w);
      (fault+i)->dip=(float)tmpdbl;
      (fault+i)->area = (fault+i)->w;
      alpha=RAD2DEGF(asin((fault+i)->sin_alpha));
      (fault+i)->strike= 90.0 - alpha;
    }
#else
    alpha=90.0-(COMP_PRECISION)fault[i].strike;
    my_sincos_deg(&fault[i].sin_alpha,&fault[i].cos_alpha,(COMP_PRECISION)alpha);
    my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[i].dip);
    calc_base_vecs(fault[i].t_strike,fault[i].normal,fault[i].t_dip,
		   fault[i].sin_alpha,fault[i].cos_alpha,sin_dip,cos_dip);
    

#endif
    calculate_bloated_corners(corner,(fault+i),leeway);
    for(j=0;j < 4;j++){
      for(k=0;k < 3;k++)
	if(fabs(corner[j][k]/CHAR_FAULT_DIM)>EPS_COMP_PREC)
	  printf("%g ",corner[j][k]/CHAR_FAULT_DIM);
	else
	  printf("0.0 ");
      printf("\t");
    }
    printf("\n");
  }
  
  printf("CELLS %i %i\n",medium->nrflt,medium->nrflt*5);
  for(i=0;i<medium->nrflt;i++)
    printf("4 %i %i %i %i\n",i*4,i*4+1,i*4+2,i*4+3);
  printf("CELL_TYPES %i\n",medium->nrflt);
  for(i=0;i<medium->nrflt;i++)
    printf("9\n");
  if(read_slip){
    /* use slip for coloring */
    printf("CELL_DATA %i\n",medium->nrflt);
    
    printf("SCALARS sqrt(s^2+d^2) float 1\n");
    printf("LOOKUP_TABLE default\n");
    for(i=0;i<medium->nrflt;i++)
      printf("%g\n",sqrt(fault[i].u[STRIKE]*fault[i].u[STRIKE] 
			 + fault[i].u[DIP]*fault[i].u[DIP]));

    printf("SCALARS strike_slip float 1\n");
    printf("LOOKUP_TABLE default\n");
    for(i=0;i<medium->nrflt;i++)
      printf("%g\n",fault[i].u[STRIKE]);
 
    printf("SCALARS dip_slip float 1\n");
    printf("LOOKUP_TABLE default\n");
    for(i=0;i<medium->nrflt;i++)
      printf("%g\n",fault[i].u[DIP]);

    /* vectors */
    printf("VECTORS slip float\n");
    for(i=0;i < medium->nrflt;i++){
      for(j=0;j<3;j++){
	u[j]  = fault[i].t_strike[j] * fault[i].u[STRIKE];
	u[j] += fault[i].t_dip[j]    * fault[i].u[DIP];
	u[j] += fault[i].normal[j]   * fault[i].u[NORMAL];
      }
      printf("%g %g %g\n",u[X],u[Y],u[Z]);
    }

  }

  exit(0);
}

