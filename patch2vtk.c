#include "interact.h"
#include "properties.h"
//
// reads in patch format and writes VTK format
//

int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  int i,j,k,l,ncon,tncon,ielmul,nvert,nel,tnvert;
  my_boolean shrink_patches=FALSE,
    verbose=FALSE,
    attempt_read_slip=FALSE,
    remove_centroid=FALSE,
    read_slip;
  COMP_PRECISION leeway,vertex[MAX_NR_EL_VERTICES*3],u[3],t_strike[3],t_dip[3],normal[3],area,xc[3],xout;
  COMP_PRECISION sscale = 1.;
  medium=(struct med *)calloc(1,sizeof(struct med));
  
  if(argc > 1)
    attempt_read_slip = TRUE;
  if(argc > 2){
    sscanf(argv[2],"%i",&i);
    shrink_patches = (my_boolean)i;
  }
  if(argc > 3){
    sscanf(argv[3],ONE_CP_FORMAT,&sscale);
  }
  if(argc > 4){
    sscanf(argv[4],"%i",&i);
    remove_centroid = (my_boolean)i;
  }
  if((argc > 5)||((argc>1) && (strcmp(argv[1],"-h")==0))){
    fprintf(stderr,"%s [flt.dat] [shrink_patches, %i] [scale_dim, %g] [remove_centroid, %i]\n\t reads in patch format from stdin and writes VTK format to stdout\n",
	    argv[0],shrink_patches,sscale,remove_centroid);
    fprintf(stderr,"if an argument is given, will assume it is a flt.dat type output file and assign values for coloring\n");
    fprintf(stderr,"if shrink_patches is set, will make patches smaller for plotting\n");
    exit(-1);
  }
  fprintf(stderr,"%s: reading patch format from stdin, writing VTK to stdout. shrink: %i scale: %g remove_cetroid: %i\n",
	  argv[0],shrink_patches,sscale,remove_centroid);
  
  read_geometry("stdin",&medium,&fault,FALSE,FALSE,FALSE,verbose);
  if(attempt_read_slip)
    read_slip = read_fltdat(argv[1],fault,medium,verbose);
  else
    read_slip = FALSE;
  if(shrink_patches)		/* to make things easier to see */
    leeway = 0.9;
  else
    leeway = 1.0;
  
  /* 
     count all nodes  
  */
  for(xc[0]=xc[1]=xc[2]=0.,nel=tncon=tnvert=i=0;i < medium->nrflt;i++){
    nvert = nvert_of_patch((fault+i));
    tnvert += nvert;
    ielmul = number_of_subpatches((fault+i));
    nel += ielmul;
    for(l=0;l < ielmul;l++)
      tncon += ncon_of_subpatch((fault+i),l);
    calculate_bloated_vertices(vertex,(fault+i),leeway);
    for(j=0;j < nvert;j++)
      for(k=0;k < 3;k++)
	xc[k] += vertex[j*3+k]/CHAR_FAULT_DIM*sscale;
  }
  for(k=0;k<3;k++)
    xc[k]/=(COMP_PRECISION)tnvert;
  if(remove_centroid){
    fprintf(stderr,"%s: vertex centroid at: %e %e %e, removing\n",argv[0],xc[0],xc[1],xc[2]);
  }else{
    fprintf(stderr,"%s: vertex centroid at: %e %e %e\n",argv[0],xc[0],xc[1],xc[2]);
    xc[0]=xc[1]=xc[2]=0;
  }
  printf("# vtk DataFile Version 2.0\n");
  printf("from patch2vtk\n");
  printf("ASCII\n");
  printf("DATASET UNSTRUCTURED_GRID\n");

  printf("POINTS %i float\n",tnvert);
  for(i=0;i < medium->nrflt;i++){
    nvert = nvert_of_patch((fault+i));
    calculate_bloated_vertices(vertex,(fault+i),leeway);
    for(j=0;j < nvert;j++){
      for(k=0;k < 3;k++){
	xout = vertex[j*3+k]/CHAR_FAULT_DIM*sscale - xc[k];
	if(fabs(xout) > EPS_COMP_PREC)
	  printf("%20.15e ",xout);
	else
	  printf("0.0 ");
      }
      printf("\t");
    }
    printf("\n");
  }
 
 
  printf("CELLS %i %i\n",nel,tncon+nel);
  for(i=k=0;i<medium->nrflt;i++){
    ielmul = number_of_subpatches((fault+i));
    for(l=0;l < ielmul;l++){
      ncon = ncon_of_subpatch((fault+i),l);
      printf("%i ",ncon);
      for(j=0;j < ncon;j++)
	printf("%i ",k+node_number_of_subelement((fault+i),j, l));
      printf("\n");
    }
    k += nvert_of_patch((fault+i));
  }
  printf("CELL_TYPES %i\n",nel);
  for(i=j=0;i<medium->nrflt;i++){
    ielmul = number_of_subpatches((fault+i));
    for(l=0;l < ielmul;l++,j++){
      printf("%i ",vtk_type_of_patch((fault+i),l));
      if(j>40){
	printf("\n");
	j=0;
      }
    }
  }
  if(j)
    printf("\n");
  if(read_slip){
    /* 
       use slip for coloring 
    */
    printf("CELL_DATA %i\n",nel);
    
    printf("SCALARS sqrt(s^2+d^2) float 1\n");
    printf("LOOKUP_TABLE default\n");
    for(i=0;i < medium->nrflt;i++){
      ielmul = number_of_subpatches((fault+i));
      for(l=0;l < ielmul;l++){
	u[STRIKE]=u[DIP]=0.;
	for(k=0;k<3;k++){
	  u[STRIKE] += projected_slip_major_to_minor_patch((fault+i),k,STRIKE,l) * fault[i].u[k];
	  u[DIP] +=    projected_slip_major_to_minor_patch((fault+i),k,DIP,   l) * fault[i].u[k];
	}
	printf("%e\n",sqrt(u[STRIKE]*u[STRIKE] + u[DIP]*u[DIP])*sscale);
      }
    }
      

    printf("SCALARS strike_slip float 1\n");
    printf("LOOKUP_TABLE default\n");
    for(i=0;i < medium->nrflt;i++){
      ielmul = number_of_subpatches((fault+i));
      for(l=0;l < ielmul;l++){
	u[STRIKE]=0.;
	for(k=0;k<3;k++)
	  u[STRIKE] += projected_slip_major_to_minor_patch((fault+i),k, STRIKE,l) * fault[i].u[k];
	printf("%e\n",u[STRIKE]*sscale);
      }
    }
 
    printf("SCALARS dip_slip float 1\n");
    printf("LOOKUP_TABLE default\n");
    for(i=0;i < medium->nrflt;i++){
      ielmul = number_of_subpatches((fault+i));
      for(l=0;l < ielmul;l++){
	u[DIP]=0.;
	for(k=0;k<3;k++)
	  u[DIP] += projected_slip_major_to_minor_patch((fault+i),k,DIP,l) * fault[i].u[k];
	printf("%e\n",u[DIP]*sscale);
      }
    }

    /* 
       vectors for displacement 
    */
    printf("VECTORS slip float\n");
    for(i=0;i < medium->nrflt;i++){
      ielmul = number_of_subpatches((fault+i));
      for(l=0;l < ielmul;l++){
	get_sub_normal_vectors((fault+i),l,t_strike,t_dip,normal,&area);
	for(j=0;j<3;j++){
	  u[j] = 0.0;
	  for(k=0;k<3;k++){
	    u[j] += t_strike[j] * fault[i].u[k]*projected_slip_major_to_minor_patch((fault+i),k,STRIKE,l) ;
	    u[j] += t_dip[j]    * fault[i].u[k]*projected_slip_major_to_minor_patch((fault+i),k,DIP,l) ;
	    u[j] += normal[j]   * fault[i].u[k]*projected_slip_major_to_minor_patch((fault+i),k,NORMAL,l) ;
	  }
	}
	printf("%e %e %e\n",u[INT_X]*sscale,u[INT_Y]*sscale,u[INT_Z]*sscale);
      }
    }
  }
  fprintf(stderr,"%s: written VTK to stdout\n",argv[0]);
  exit(0);
}

