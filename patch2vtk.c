#include "interact.h"
#include "properties.h"
//
// reads in patch format and writes VTK format
//

int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  int i,j,k,ncon,tncon;
  my_boolean shrink_patches=FALSE,
    verbose=FALSE,
    attempt_read_slip=FALSE,
    read_slip;
  COMP_PRECISION leeway,corner[4][3],u[3];
  medium=(struct med *)calloc(1,sizeof(struct med));
  
  if(argc > 1)
    attempt_read_slip = TRUE;
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
  if(attempt_read_slip)
    read_slip = read_fltdat(argv[1],fault,medium,verbose);
  else
    read_slip = FALSE;
  if(shrink_patches)		/* to make things easier to see */
    leeway = 0.9;
  else
    leeway = 1.0;
  /* count all nodes  */
  for(tncon=i=0;i < medium->nrflt;i++)
    tncon += ncon_of_patch((fault+i)) ; 
  printf("# vtk DataFile Version 2.0\n");
  printf("from patch2vtk\n");
  printf("ASCII\n");
  printf("DATASET UNSTRUCTURED_GRID\n");
  
  printf("POINTS %i float\n",tncon);
  for(i=0;i < medium->nrflt;i++){
    calculate_bloated_corners(corner,(fault+i),leeway);
    ncon = ncon_of_patch((fault+i));
    for(j=0;j < ncon;j++){
      for(k=0;k < 3;k++)
	if(fabs(corner[j][k]/CHAR_FAULT_DIM) > EPS_COMP_PREC)
	  printf("%g ",corner[j][k]/CHAR_FAULT_DIM);
	else
	  printf("0.0 ");
      printf("\t");
    }
    printf("\n");
  }
  
  printf("CELLS %i %i\n",medium->nrflt,tncon+medium->nrflt);
  for(i=k=0;i<medium->nrflt;i++){
    ncon = ncon_of_patch((fault+i));
    printf("%i ",ncon);
    for(j=0;j < ncon;j++,k++)
      printf("%i ",k);
    printf("\n");
  }
  printf("CELL_TYPES %i\n",medium->nrflt);
  for(i=0;i<medium->nrflt;i++){
    printf("%i\n",vtk_type_of_patch((fault+i)));
  }
  if(read_slip){
    /* 
       use slip for coloring 
    */
    printf("CELL_DATA %i\n",medium->nrflt);
    
    printf("SCALARS sqrt(s^2+d^2) float 1\n");
    printf("LOOKUP_TABLE default\n");
    for(i=0;i < medium->nrflt;i++)
      printf("%g\n",sqrt(fault[i].u[STRIKE]*fault[i].u[STRIKE] 
			 + fault[i].u[DIP]*fault[i].u[DIP]));

    printf("SCALARS strike_slip float 1\n");
    printf("LOOKUP_TABLE default\n");
    for(i=0;i < medium->nrflt;i++)
      printf("%g\n",fault[i].u[STRIKE]);
 
    printf("SCALARS dip_slip float 1\n");
    printf("LOOKUP_TABLE default\n");
    for(i=0;i < medium->nrflt;i++)
      printf("%g\n",fault[i].u[DIP]);

    /* 
       vectors for displacement 
    */
    printf("VECTORS slip float\n");
    for(i=0;i < medium->nrflt;i++){
      for(j=0;j<3;j++){
	u[j]  = fault[i].t_strike[j] * fault[i].u[STRIKE];
	u[j] += fault[i].t_dip[j]    * fault[i].u[DIP];
	u[j] += fault[i].normal[j]   * fault[i].u[NORMAL];
      }
      printf("%g %g %g\n",u[INT_X],u[INT_Y],u[INT_Z]);
    }

  }
  fprintf(stderr,"%s: written VTK to stdout\n",argv[0]);
  exit(0);
}

