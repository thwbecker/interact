#include "interact.h"
#include "properties.h"
//
// reads in sets of three points in 3-D and converts to triangle patch format
//


int main(int argc, char **argv)
{
  
#ifndef ALLOW_NON_3DQUAD_GEOM
  // simply for ease of storage reasons
  fprintf(stderr,"%s needs ALLOW_NON_3DQUAD_GEOM to be set\n",argv[0]);exit(-1);
#else
  struct flt *fault;
  struct med *medium;
  COMP_PRECISION *dummy=NULL;
  int i,j,opmode=PATCH_OUT_MODE,n,eltype=TRIANGULAR;
  COMP_PRECISION sin_dip,cos_dip;
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(argc!=2){
    fprintf(stderr,"%s mode\n\tread in three 3-D points per line from stdin\n",
	    argv[0]);
    fprintf(stderr,"\tformat:\n\tx1_x x1_y x1_z x2_x x2_y x2_z ...\n\n");
    fprintf(stderr,"\tto form a triangle, points have to be in FE (CCW) ordering\n");
    fprintf(stderr,"\twrites patch format to stdout\n");
    fprintf(stderr,"if mode=%i, will write triangular element format\n",
	    TRIANGULAR);
    fprintf(stderr,"if mode=%i, will attempt to fit the best rectangular element\n",
	    RECTANGULAR_PATCH);
    exit(-1);
  }else
    sscanf(argv[1],"%i",&eltype);
  fprintf(stderr,"%s: reading points from %s, writing patch format %i to stdout\n",
	  argv[0],"stdin",eltype);


  if((fault=malloc(sizeof(struct flt)))==NULL)MEMERROR("main");
  fault[0].xt=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*9);
  n=0;
  while(fscanf(stdin,NINE_CP_FORMAT,
	       &fault[n].xt[  INT_X],&fault[n].xt[  INT_Y],&fault[n].xt[  INT_Z],
	       &fault[n].xt[3+INT_X],&fault[n].xt[3+INT_Y],&fault[n].xt[3+INT_Z],
	       &fault[n].xt[6+INT_X],&fault[n].xt[6+INT_Y],&fault[n].xt[6+INT_Z])==9){
    for(j=1,i=2;i<9;i+= 3,j++)
      if(fault[n].xt[i] > 0.0){
	fprintf(stderr,"%s: rectangle %i: point %i: z coordinate (%g) should be <= 0\n",
		argv[0],n,j,fault[n].xt[i]);
	exit(-1);
      }
    // space for next
    if((fault=realloc(fault,sizeof(struct flt)*(n+2)))==NULL)
      MEMERROR("main");
    fault[n+1].xt=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*9);
    if(!fault[n+1].xt)MEMERROR("main");
    /*  */
    // init the triangular properties
    get_tri_prop_based_on_gh((fault+n));
    fault[n].group=0;
    
    if(eltype == TRIANGULAR){
      // simply write triangular element type to stdout
      fault[n].type=TRIANGULAR;
    }else{// convert to rectangular
      // L=W=sqrt(A/4)
      fault[n].l=sqrt(fault[n].w/4.0);
      fault[n].w=fault[n].l;
      /* recompute */
      my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[n].dip);
      calc_quad_base_vecs(fault[n].t_strike,fault[n].normal,fault[n].t_dip,
			  fault[n].sin_alpha,fault[n].cos_alpha,sin_dip,cos_dip);
      fault[n].type=RECTANGULAR_PATCH;
    }
    n++;
  }
  medium->nrflt=n;
  fault=realloc(fault,sizeof(struct flt)*medium->nrflt);
  for(i=0;i < medium->nrflt;i++)// output
    print_patch_geometry_and_bc(0,(fault+i),opmode,0.0,FALSE,stdout,FALSE,dummy);
  exit(0);
#endif

}



