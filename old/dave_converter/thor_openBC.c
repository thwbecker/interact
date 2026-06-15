

void thor_openBC(FILE *in,
		 int *code, int *field_option,
		 double *xmin,double *xmax, int *n,
		 double *ymin,double *ymax, int *m,
		 double *zmin,double *zmax, int *o)
{
  
  if(fscanf(in,"%i %i",code,field_option)!=2){
    fprintf(stderr,"could not read two integers\n");
    exit(-1);
  }
  if(fscanf(in,"%lf %lf %i",xmin,xmax,n)!=3){ /* double! */
    fprintf(stderr,"could not asdfjs\n");
    exit(-1);
  }
  /*  if(fscanf(in,"%f %f %i",xmin,xmax,n)!=3){ /* float! */ */
  /*     fprintf(stderr,"could not asdfjs\n"); */
  /*     exit(-1); */
  /*   } */
  

}
	 
