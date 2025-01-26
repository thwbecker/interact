/*


$Id: block_eval_geokada.c,v 1.7 2003/07/02 22:53:12 becker Exp $


contains the routine to evaluate a fault given a gegraphic location

also holds input/output routines for saveing and loading the solution
of a block velocity and stress inversion, ie. the solution vector,
global projection, fault geometry and slip

furthermore, holds the block velocity evaluation routine

block_eval_geookada:


xc: lon,lat of obsrvational point
udepth, sdepth: depths at which to evaluate displacements 
and stresses, > 0

flon,flat,fz: fault center location in degrees, fault center depth
fazi,fdip: fault azimuth and dip in degrees
fl,fw: faulthalf-length and width
flfac: fault slip factor (usually unity)
fca, fsa: cos(alpha) and sin(alpha) for fault

this returns u as displacements and s[][] as stresses only
if the eval_dip and/or eval_stress flags are set. else, both
arrays are zero. 

if we get NaN, we will also return zeroes


*/
#include "interact.h"
#include "blockinvert.h"


void block_eval_geookada(COMP_PRECISION *xl, 
			 COMP_PRECISION *disp,
			 COMP_PRECISION *u, 
			 COMP_PRECISION s[3][3],
			 COMP_PRECISION udepth,
			 COMP_PRECISION sdepth,
			 COMP_PRECISION flon,COMP_PRECISION flat,
			 COMP_PRECISION fz,COMP_PRECISION fazi,
			 COMP_PRECISION fdip,
			 COMP_PRECISION fl,COMP_PRECISION fw,
			 COMP_PRECISION fca,COMP_PRECISION fsa,
			 int *iret,
			 my_boolean eval_disp,my_boolean eval_stress)
{
  COMP_PRECISION x[3]={0,0,0},px[3],pu[3],ps[3][3],dummy=0;
  int i,j;
  //
  // project observational point in fault local system
  a_equals_b_vector(x,xl,2);
  geoproject(x, px, FLT_ROT_PROJECTION,
	     flon,flat,fazi,dummy, dummy, 
	     dummy, dummy, (int)FALSE);
  //fprintf(stderr,"p: %11g %11g %11g ",flon,flat,fazi);
  if(eval_disp){
    /* evaluate velocities at udepth */
    px[INT_Z] = -udepth;			
    //
    // evaluate displacements in rotated frame
    eval_rectangle_basic(px,fl,fw,fdip,-fz,disp,pu,ps,iret);
    /* this will also return stresses, but possibly at 
       different depth from sdepth */
    if(*iret){
      fprintf(stderr,"block_eval_geookada: WARNING: singular at obs point: %g, %g, %g\n",
	      *(xl+INT_X),*(xl+INT_Y),udepth);
      /* NaN returns zeros */
      for(i=0;i < 3;i++){
	u[i]=0.0;
	for(j=0;j < 3;j++)
	  s[i][j]=0.0;
      }
    }else{
      // rotate displacements pu back into the observational frame
      rotate_vec(pu,u,fca,-fsa);
    }
  }else{			/* don't return displacements */
    for(i=0;i<3;i++)
      u[i] = 0.0;
  }
  if(eval_stress){
    if((udepth == sdepth) && (eval_disp) && (!(*iret))){
      /* use previously determined stresses and 
	 rotate stresses back */
      rotate_mat_z(ps,s,fca,-fsa);
    }else{
      /* reevaluate stresses, since different depth */
      px[Z] = -sdepth;
      eval_rectangle_basic(px,fl,fw,fdip,-fz,disp,pu,ps,iret);
      if(*iret){
	fprintf(stderr,"block_eval_geookada: WARNING: singular at obs point: %g, %g, %g\n",
		*(xl+INT_X),*(xl+INT_Y),sdepth);
	/* return stresses as zeroes if NaN */
	for(i=0;i < 3;i++){
	  for(j=0;j < 3;j++)
	    s[i][j]=0.0;
	}
      }else{
	/* rotate into reference frame */
	rotate_mat_z(ps,s,fca, -fsa);
      }
    }
  }else{			/* don't evaluate stresses */
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	s[i][j] = 0.0;
  }
}
 
void block_save_solution_and_faults(COMP_PRECISION *x,
				    int nrb,int nrf, 
				    struct bflt *fault, 
				    COMP_PRECISION *disp,
				    struct prj *projection,
				    FILE *out,my_boolean use_ld,
				    my_boolean use_cfac)
{
  int i,izero=0;
  // master projection
  fwrite(projection,sizeof(struct prj),1,out);
  // block motion parameters
  fwrite(&nrb,sizeof(int),1,out);
  fwrite(x,sizeof(COMP_PRECISION),nrb*BLOCK_NBASE,out);
  // fault geometry
  fwrite(&nrf,sizeof(int),1,out);
  for(i=0;i < nrf;i++){
    fwrite(fault[i].x,sizeof(COMP_PRECISION),3,out);
    fwrite(&fault[i].azi,sizeof(COMP_PRECISION),1,out);
    fwrite(&fault[i].dip,sizeof(COMP_PRECISION),1,out);
    fwrite(&fault[i].ca,sizeof(COMP_PRECISION),1,out);
    fwrite(&fault[i].sa,sizeof(COMP_PRECISION),1,out);
    fwrite(&fault[i].l,sizeof(COMP_PRECISION),1,out);
    fwrite(&fault[i].w,sizeof(COMP_PRECISION),1,out);
   }
  // slip values
  fwrite(disp,sizeof(COMP_PRECISION),3*nrf,out);
  i=nrb*BLOCK_NBASE;
  if(use_ld){
    // locking depths
    fwrite(&nrf,sizeof(int),1,out);
    fwrite((x+i),sizeof(COMP_PRECISION),nrf,out);
    i += nrf;
  }else{
    fwrite(&izero,sizeof(int),1,out);
  }
  if(use_cfac){
    // coupling factors
    fwrite(&nrf,sizeof(int),1,out);
    fwrite((x+i),sizeof(COMP_PRECISION),nrf,out);
  }else{
    fwrite(&izero,sizeof(int),1,out);
  }
}
void block_load_solution_and_faults(COMP_PRECISION **x,
				    int *nrb, int *nrf, 
				    struct bflt **fault, 
				    COMP_PRECISION **disp,
				    struct prj **projection,
				    FILE *in,
				    my_boolean *use_ld,
				    my_boolean *use_cfac)
{
  int i,j,os,nf3;
  *projection=
    (struct prj *)realloc(*projection,sizeof(struct prj));
  if(!*projection)
    MEMERROR("");
  // projection
  if(fread(*projection,sizeof(struct prj),1,in)!=1)
    READ_ERROR("");
  fprintf(stderr,"block_load_solution_and_faults: projection: %i %g %g %g %g %g\n",
	  (*projection)->type,(*projection)->clon,
	  (*projection)->clat,(*projection)->azi,
	  (*projection)->lat1,(*projection)->lat2);  
  // block motion parameters
  if(fread(nrb,sizeof(int),1,in)!=1)
    READ_ERROR("");
  *x=(COMP_PRECISION *)realloc(*x,sizeof(COMP_PRECISION)
			       *(*nrb)*BLOCK_NBASE);
  if(!*x)
    MEMERROR("");
  if(fread(*x,sizeof(COMP_PRECISION),BLOCK_NBASE*(*nrb),in)!=
     BLOCK_NBASE*(*nrb))READ_ERROR("");
  fprintf(stderr,"block_load_solution_and_faults: loaded solution for %i blocks OK\n",
	  *nrb);
  //for(i=0;i<BLOCK_NBASE*(*nrb);i++)fprintf(stderr,"%g\n",*(*x+i));
  // fault parameters
  if(fread(nrf,sizeof(int),1,in)!=1)
    READ_ERROR("");
  nf3 = (*nrf) * 3;
  *fault=(struct bflt *)realloc(*fault,sizeof(struct bflt)*(*nrf));
  *disp=(COMP_PRECISION *)realloc(*disp,sizeof(COMP_PRECISION)*nf3);
  if(!*fault || !*disp)
    MEMERROR("");
  for(i=0;i<(*nrf);i++){
    j =fread((*fault+i)->x,sizeof(COMP_PRECISION),3,in);
    j+=fread(&(*fault+i)->azi,sizeof(COMP_PRECISION),1,in);
    j+=fread(&(*fault+i)->dip,sizeof(COMP_PRECISION),1,in);
    j+=fread(&(*fault+i)->ca,sizeof(COMP_PRECISION),1,in);
    j+=fread(&(*fault+i)->sa,sizeof(COMP_PRECISION),1,in);
    j+=fread(&(*fault+i)->l,sizeof(COMP_PRECISION),1,in);
    j+=fread(&(*fault+i)->w,sizeof(COMP_PRECISION),1,in);
    if(j != 9)
      READ_ERROR("");
  }
  // slip values
  if(fread(*disp,sizeof(COMP_PRECISION),nf3,in) != nf3)
    READ_ERROR("");
  //
  //
  os = (*nrb)*BLOCK_NBASE;
  /* 

     locking depths 

  */
  if(fread(&j,sizeof(int),1,in)!=1)
    READ_ERROR("");
  if(j){
    *use_ld = TRUE;
    if(j != *nrf){
      fprintf(stderr,"block_load_solution_and_faults: error: nflt for ld != nflt\n");
      exit(-1);
    }
    my_vecrealloc(x,os+(*nrf),"os1");
    if(fread((*x +os),sizeof(COMP_PRECISION),(*nrf),in) != (*nrf))
      READ_ERROR("");
    for(i=0;i < (*nrf);i++){
      if((*fault+i)->x[Z] * 2 != *(*x+os+i)){
	fprintf(stderr,"block_load_solution_and_faults: locking depth error fault %i, z: %g ld: %g\n",
		i+1,(*fault+i)->x[Z], *(*x+os+i));
	exit(-1);
      }
    }
    os += *nrf;
  }else{
    *use_ld = FALSE;
  }
  /* 
     coupling factors 
  */
  if(fread(&j,sizeof(int),1,in) != 1)
    READ_ERROR("");
  if(j){
    *use_cfac = TRUE;
    if(j != *nrf){
      fprintf(stderr,"block_load_solution_and_faults: error: nflt for cfac != nflt\n");
      exit(-1);
    }
    my_vecrealloc(x,os+(*nrf),"os2");
    if(fread((*x + os),sizeof(COMP_PRECISION),(*nrf),in) != (*nrf))
      READ_ERROR("");
  }else{
    *use_cfac = FALSE;
  }
  fprintf(stderr,
	  "block_load_solution_and_faults: loaded %i faults and slip vectors OK %s %s\n",
	  *nrf,(*use_ld)?("(using locking depth)"):(""),
	  (*use_cfac)?("(using coupling factors)"):(""));
  }
/* 
   evaluate the velocity according to the block 
   point x[2],uz is in 

*/
void block_eval_blockvec(COMP_PRECISION *xy,COMP_PRECISION z,
			 int bcode,COMP_PRECISION *lbu,
			 struct prj *projection,
			 COMP_PRECISION *xsol)
{
  COMP_PRECISION x[3],px[3],dummy=0;
  static my_boolean init=FALSE;
  if(!init){
    if((BLOCK_DIM != 2)||(BLOCK_NBASE != 3)){
      fprintf(stderr,
	      "block_eval_blockvec: not prepared for dim != 2 or block_base != 3\n");
      exit(-1);
    }
    if((projection->type != PROJECT_AZI)&&
       (projection->type != OMERC_AZI)&&
       (projection->type != LCONFORM)){
      fprintf(stderr,"block_eval_blockvec: projection type %i not implemented\n",
	      projection->type);
      exit(-1);
    }
    init = TRUE;
  }
  x[INT_X] = xy[INT_X];x[INT_Y] = xy[INT_Y];x[INT_Z]=z;
  // projection observation to global projected reference frame
  geoproject(x,px,projection->type,projection->clon,
	     projection->clat,projection->azi,dummy,dummy,
	     projection->lat1,projection->lat2,(int)FALSE);
  // solve A . x = v
  // where A is ( -o_y 1 0 )
  //            (  o_x 0 1 )
  //
  lbu[INT_X] = -px[INT_Y] * xsol[bcode*BLOCK_NBASE] + xsol[bcode*BLOCK_NBASE+1];
  lbu[INT_Y] =  px[INT_X] * xsol[bcode*BLOCK_NBASE] + xsol[bcode*BLOCK_NBASE+2];
  lbu[INT_Z] = 0.0;
}

 
