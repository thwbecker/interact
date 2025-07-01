#include "interact.h"
#include "properties.h"
//
// reads in sets of three points in 3-D and converts to triangle patch format
//

#define TRI2P_IOMODE_TRI 20
#define TRI2P_IOMODE_TRI_GROUP 21
#define TRI2P_IOMODE_QUAD 0
#define TRI2P_IOMODE_QUAD_GROUP 1

int main(int argc, char **argv)
{
  
#ifndef ALLOW_NON_3DQUAD_GEOM
  // simply for ease of storage reasons
  fprintf(stderr,"%s needs ALLOW_NON_3DQUAD_GEOM to be set\n",argv[0]);exit(-1);
#else
  struct flt *fault;
  struct med *medium;
  COMP_PRECISION *dummy=NULL;
  my_boolean read_group=FALSE;
  int i,j,opmode=PATCH_OUT_MODE,n,mode,eltype,grp_min,grp_max,itmp;
  COMP_PRECISION sin_dip,cos_dip;
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(argc!=2){
    fprintf(stderr,"%s mode\n\tread in three 3-D points per line from stdin\n",
	    argv[0]);
    fprintf(stderr,"\tformat:\n\tx1_x x1_y x1_z x2_x x2_y x2_z x3_x x3_y x3_z [group]\n\n");
    fprintf(stderr,"\tto form a triangle, points have to be in FE (CCW) ordering\n");
    fprintf(stderr,"\twrites patch format to stdout\n\n");
    fprintf(stderr,"\tif mode=%i or %i, will write triangular element format\n",
	    TRI2P_IOMODE_TRI,TRI2P_IOMODE_TRI_GROUP);
    fprintf(stderr,"\tif mode=%i or %i, will attempt to fit the best rectangular element\n",
	    TRI2P_IOMODE_QUAD,TRI2P_IOMODE_QUAD_GROUP);
    fprintf(stderr,"\tif modes are %i or %i, will also expect an integer group label in the 10th column\n",
	    TRI2P_IOMODE_QUAD_GROUP,TRI2P_IOMODE_TRI_GROUP);
    exit(-1);
  }else
    sscanf(argv[1],"%i",&mode);

  switch(mode){
  case TRI2P_IOMODE_TRI:
    eltype = TRIANGULAR;
    break;
  case TRI2P_IOMODE_TRI_GROUP:
    eltype = TRIANGULAR;
    read_group = TRUE;
    break;
  case TRI2P_IOMODE_QUAD:
    eltype = OKADA_PATCH;
    break;
  case TRI2P_IOMODE_QUAD_GROUP:
    eltype = OKADA_PATCH;
    read_group = TRUE;
    break;
  default:
    fprintf(stderr,"%s: mode %i undefined\n",argv[0],mode);
    exit(-1);
  }
  fprintf(stderr,"%s: reading set of nine %s points from %s, writing patch format %i to stdout\n",
	  argv[0],"stdin",read_group?"+1 for group":"",eltype);
  

  if((fault=malloc(sizeof(struct flt)))==NULL)MEMERROR("main");
  fault[0].xn=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*9);
  n=0;
  grp_min= INT_MAX;
  grp_max=-INT_MAX;
  while(fscanf(stdin,NINE_CP_FORMAT,
	       &fault[n].xn[  INT_X],&fault[n].xn[  INT_Y],&fault[n].xn[  INT_Z],
	       &fault[n].xn[3+INT_X],&fault[n].xn[3+INT_Y],&fault[n].xn[3+INT_Z],
	       &fault[n].xn[6+INT_X],&fault[n].xn[6+INT_Y],&fault[n].xn[6+INT_Z])==9){
    if(read_group){
      if(fscanf(stdin,"%i",&itmp)!=1){
	fprintf(stderr,"%s: group read error %i\n",argv[0],n);
	  exit(-1);
      }
      if(itmp < 0){
	fprintf(stderr,"%s: group id  %i should be >=0 (%i)\n",argv[0],n,itmp);
	exit(-1);
      }
      fault[n].group = (unsigned int)itmp;
    }else{
      fault[n].group = 0;
    }
    if((int)fault[n].group > grp_max){
      grp_max = fault[n].group;
    }
    if(fault[n].group < grp_min)
      grp_min = fault[n].group;

    for(j=1,i=2;i<9;i+= 3,j++)
      if(fault[n].xn[i] > 0.0){
	fprintf(stderr,"%s: rectangle %i: point %i: z coordinate (%g) should be <= 0\n",
		argv[0],n,j,fault[n].xn[i]);
	exit(-1);
      }
    // space for next
    if((fault=realloc(fault,sizeof(struct flt)*(n+2)))==NULL)
      MEMERROR("main");
    fault[n+1].xn=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*9);
    if(!fault[n+1].xn)MEMERROR("main");
    /*  */
    // init the triangular properties
    get_tri_prop_based_on_gh((fault+n));


    
    fault[n].type = eltype;    
    if(eltype != TRIANGULAR){
      // convert to rectangular
      // L=W=sqrt(A/4)
      fault[n].l=sqrt(fault[n].w/4.0);
      fault[n].w=fault[n].l;
      /* recompute */
      my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[n].dip);
      calc_quad_base_vecs(fault[n].t_strike,fault[n].normal,fault[n].t_dip,
			  fault[n].sin_alpha,fault[n].cos_alpha,sin_dip,cos_dip);
    }
    n++;
  }
  medium->nrflt=n;
  fault=realloc(fault,sizeof(struct flt)*medium->nrflt);
  for(i=0;i < medium->nrflt;i++)// output
    print_patch_geometry_and_bc(0,(fault+i),opmode,0.0,FALSE,stdout,FALSE,dummy);
  fprintf(stderr,"%s: produced %i elements total, group code from %i to %i\n",
	  argv[0],medium->nrflt,grp_min,grp_max);
  exit(0);
#endif

}



