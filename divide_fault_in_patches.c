/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: divide_fault_in_patches.c,v 1.20 2011/01/09 02:02:43 becker Exp $
*/
#include "properties.h"
#include "interact.h"



void create_patches(int flt, struct flt *fault,struct flt **patch,
		    int *nrpatches,int *seg,my_boolean circular,
		    my_boolean calculate_base_vecs,COMP_PRECISION *td,
		    COMP_PRECISION srand, COMP_PRECISION drand,
		    long *seed)	/* gaussian randomization */
{
  int segn[2];
  determine_segments(seg,segn,(fault+flt),circular,td);
  divide_fault_in_patches(flt,fault,patch,nrpatches,segn,circular,
			  calculate_base_vecs,srand,drand,seed);
}
/*

  determine the number of subdivisions of a fault patch
  
*/

void determine_segments(int *seg,int *segn, struct flt *fault,
			my_boolean circular,COMP_PRECISION *td)
{
  int j;
  for(j=0;j<2;j++){
    if(td[j]>0){
      segn[j]=(int)((j==0?fault->l:fault->w)/td[j]+0.5);
      if(segn[j]<1)
	segn[j]=1;
      if(circular)
	if((segn[j] % 2) == 0)
	  segn[j]++;
    }else{
      segn[j]=seg[j];
    }
  } 
}

/* 


subdivide a rectangular surface into sub-fault patches


 */

void divide_fault_in_patches(int flt,struct flt *fault,
			     struct flt **patch,
			     int *nrpatches,int *seg,
			     my_boolean circular,
			     my_boolean calculate_base_vecs,
			     COMP_PRECISION srand, COMP_PRECISION drand,
			     long *seed)
{
  int i,j,mi,mj,irad,old_nrpatches,
    added_patches,pcnt;
  COMP_PRECISION dx[3],dy[3],sin_dip,cos_dip,
    corner[4][3],newl,neww,tmpd1,tmpd2,newarea,
    l,w;
  double alpha;
  my_boolean randomize;
  randomize = ((srand>0)||(drand>0))?(TRUE):(FALSE);
  
  added_patches=seg[0]*seg[1];
  if(!added_patches){
    fprintf(stderr,"divide_fault_in_patches: nm %i/%i error\n",
	    seg[0],seg[1]);
    exit(-1);
  }
  // increase patches structure
  old_nrpatches= *nrpatches;
  if(!*nrpatches){
    *nrpatches += added_patches;
    *patch = (struct flt *)calloc(*nrpatches,sizeof(struct flt));
    if(! *patch){
      fprintf(stderr,"divide_fault_in_patches: mem error, nr patches: %i\n",*nrpatches);
      exit(-1);
    }
  }else{
    *nrpatches += added_patches;
    *patch = (struct flt *)realloc(*patch,*nrpatches*sizeof(struct flt));
    if(! *patch){
      fprintf(stderr,"divide_fault_in_patches: mem error, nr patches: %i\n",*nrpatches);
      exit(-1);
    }
  }
  // if base vecs not yet initialized, do now
  if(calculate_base_vecs){
    alpha = 90.0 - (double)fault[flt].strike;
    my_sincos_degd(&fault[flt].sin_alpha,&fault[flt].cos_alpha,alpha);
    my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)fault[flt].dip);
    calc_base_vecs(fault[flt].t_strike,fault[flt].normal,fault[flt].t_dip,
		   fault[flt].sin_alpha,fault[flt].cos_alpha,sin_dip,cos_dip);
  }
  // obtain corner points of original patch
  calculate_corners(corner,(fault+flt),&l,&w);
  // new patch geometry
  newl = l / (COMP_PRECISION)seg[0];
  neww = w / (COMP_PRECISION)seg[1];
  newarea = newl * neww;
  for(i=0;i<3;i++){
    dx[i] = fault[flt].t_strike[i]*newl;
    dy[i] = fault[flt].t_dip[i]*neww;
  }
  // assign common specs by copying them over from the 
  // fault to the patches
  for(i=old_nrpatches;i< *nrpatches;i++){
    fltcp((fault+flt),(*patch+i));

    if(randomize){		/* randomize this patch? */
      randomize_strike_dip(srand,drand,(*patch+i),seed);
      /* need to recompute basis vectors */
      alpha = 90.0 - (double)(*patch+i)->strike;
      my_sincos_degd(&((*patch+i)->sin_alpha),&((*patch+i)->cos_alpha),alpha);
      my_sincos_deg(&sin_dip,&cos_dip,(COMP_PRECISION)(*patch+i)->dip);
      calc_base_vecs((*patch+i)->t_strike,(*patch+i)->normal,(*patch+i)->t_dip,
		     (*patch+i)->sin_alpha,(*patch+i)->cos_alpha,sin_dip,cos_dip);
    }
    (*patch+i)->area = newarea;
    (*patch+i)->w = neww/2.;
    (*patch+i)->l = newl/2.;
  }
  if(!circular){/* divide fault patch in rectangular patches */
    /* 
       
    go left right, then top down

    */

    pcnt = old_nrpatches;
    for(j=seg[1]-1;j >= 0;j--){
      for(i=0;i<seg[0];i++){
	/* calculate position of patch in fault */
	get_flt_location((*patch + pcnt),dx,dy,&corner[0][0],i,j);
	pcnt++;
      }
    }
  }else{
    /* 
       divide fault patch in circular patches 
       lk and wk are the same since n and m are equal
    */
    if((seg[0] % 2 == 0)||(seg[1] % 2 == 0)){
      fprintf(stderr,"divide_fault_in_patches: need odd numbers for circular patches\n");
      exit(-1);
    }
    mi=(seg[1]+1)/2;mj=(seg[0]+1)/2;
    irad=MIN(mi,mj);
    for(pcnt=old_nrpatches,j=0;j<seg[0];j++){
      for(i=0;i<seg[1];i++){
	tmpd1=(COMP_PRECISION)(i-mi); 
	tmpd2=(COMP_PRECISION)(j-mj);
	if(hypot(tmpd1,tmpd2) > irad)
	  continue;
	/* calculate position of patch in fault */
	get_flt_location((*patch + pcnt),dx,dy,&corner[0][0],i,j);
	pcnt++;
      }
    }
  }
  // resize patch array if we did not need all the new patches
  if(pcnt<added_patches){
    *nrpatches = old_nrpatches + pcnt;
    *patch = (struct flt *)realloc(*patch,*nrpatches * sizeof(struct flt));
  }
}

/* 

get the center point location of a sub fault patch for strike spacing
dx, dip spacing dy, and strike index i, as well as dip index j 

*/
void get_flt_location(struct flt *fault,COMP_PRECISION *dx,COMP_PRECISION *dy,
		      COMP_PRECISION *corner,int i, int j)
{
  int k;
  for(k=0;k<3;k++){
    fault->x[k]  = corner[k]+(0.5+(COMP_PRECISION)i)*dx[k];
    fault->x[k] += (0.5+(COMP_PRECISION)j)*dy[k];
  }
}  

void randomize_strike_dip(COMP_PRECISION srand,
			  COMP_PRECISION drand,struct flt *fault,
			  long *seed)
{
  COMP_PRECISION dip;
  double strike;
  /* compute Gaussian deviation with srand drand std */
  strike = fault->strike + mygauss_randnr(srand,seed); 
  dip    = fault->dip    + mygauss_randnr(drand,seed); 
  check_angles(&dip,&strike);
  fault->strike = (COMP_PRECISION)strike;
  fault->dip = dip;
}
   
		      
