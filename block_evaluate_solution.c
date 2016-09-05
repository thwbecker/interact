/*


read solution of slip and fault geometry and evaluate 
at different location

usage:


block_evaluate_solution solution.bin 



$Id: block_evaluate_solution.c,v 1.8 2011/01/07 07:19:58 becker Exp $


*/
#include "interact.h"
#include "blockinvert.h"


int main(int argc, char **argv)
{
  FILE *in;
  COMP_PRECISION *disp=NULL,x[2],*xsol=NULL,fac,px[3],dummy=0,
    s_evaldepth,u_evaldepth,fu[3],s[3][3],lfu[3],ls[3][3],bu[3];
  int i,j,n,k,nrf,iret,nrb,bcode,os1;
  struct bflt *fault=NULL;
  struct prj *projection=NULL;
  my_boolean eval_stress=TRUE,eval_disp=TRUE,rigid,
    use_cfac,warned=FALSE,use_ld,read_loc_file=TRUE,
    fault_geo_out=FALSE;
  /* 
     read in fault geometry and slip vectors 
  */
  if(argc == 1)
    in = myopen("solution.bin","r");
  else
    in = myopen(argv[1],"r");
  if(argc>2)
    for(i=2;i<argc;i++){	
      if(strcmp(argv[i],"-fout")==0){ /* output of fault 
					 geometry */
	read_loc_file=FALSE;
	fault_geo_out=TRUE;
      }
    }
  block_load_solution_and_faults(&xsol,&nrb,&nrf,&fault,&disp,
				 &projection,in,&use_ld,
				 &use_cfac);
  fclose(in);
  if(norm(disp,3*nrf) < 1e-14){
    rigid = TRUE;
    fprintf(stderr,"%s: slip vector close to zero, assuming rigid\n",
	    argv[0]);
  }else
    rigid = FALSE;
  /* 
     scale displacements with locking factors, if intialized 
  */
  if(use_cfac){
    for(i=0;i < nrf;i++){
      fac = xsol[nrb*BLOCK_NBASE + i];
      if(fabs(fac-1.0)>EPS_COMP_PREC){
	if(!warned){
	  fprintf(stderr,"%s: at least one locking factor != 1\n",
		  argv[0]);
	  warned = TRUE;
	}
	for(j=0;j<3;j++)	/* scale displacements */
	  disp[i*3+j] *= fac;
      }
    }
  }
  fprintf(stderr,"%s: total slip vector norm: %11g\n",
	  argv[0],norm(disp,3*nrf));
  if(fault_geo_out){
    /* 
       output of fault geometry in interact patch format 
       in projected coordinate system
       
    */
    for(i=0;i < nrf;i++){
      geoproject(fault[i].x,px,projection->type,
		 projection->clon,
		 projection->clat,projection->azi,dummy,dummy,
		 projection->lat1,projection->lat2,(int)FALSE);
      printf("%12g %12g %12g %12g %12g %12g %12g 0\n",
	     px[X],px[Y],px[Z],fault[i].azi,
	     fault[i].dip,fault[i].l,fault[i].w);

    }    
  }
  if(read_loc_file){
    //
    // depth at which to evaluate stuff
    u_evaldepth = 0.0;// displacements at surface
    s_evaldepth = 7.5;// stresses at depth 
    fprintf(stderr,"%s: reading lon/lat from stdin\n",argv[0]);
    n=0;
    while(fscanf(stdin,"%lf %lf %i",x,(x+1),&bcode)==3){
      n++;
      if((bcode < 1)||(bcode>nrb)){
	fprintf(stderr,"%s: error, block code: %i but nrb: %i\n",
		argv[0],bcode,nrb);
	exit(-1);
      }
      bcode--;			/* block code */
      //
      // velocity due to block motion
      block_eval_blockvec(x,u_evaldepth,bcode,bu,projection,xsol);
      /* initialize fault contributions */
      for(i=0;i<3;i++){
	fu[i]=0.0;		/* displacements */
	if(eval_stress)		/* stresses */
	  for(j=0;j<3;j++)
	    s[i][j] = 0.0;
      }
      if(!rigid){
	/* 
	   sum up fault contributions 
	*/
	for(i=os1=0;i < nrf;i++,os1+=3){
	  /* evaluate displacements and stresses due to fault i,
	     possibly at different depths. NaNs will be returned 
	     as zeroes
	  */
	  block_eval_geookada(x,(disp+os1),lfu,ls,u_evaldepth,
			      s_evaldepth,fault[i].x[X],
			      fault[i].x[Y],fault[i].x[Z],
			      fault[i].azi,fault[i].dip,
			      fault[i].l,fault[i].w,
			      fault[i].ca, fault[i].sa,&iret,
			      eval_disp,eval_stress);
	  /* add to global field */
	  for(j=0;j<3;j++){
	    /* the local displacement and stress arrays,
	       lfu and ls, will be zero if eval_disp and/or 
	       eval_stress are FALSE 
	    */
	    fu[j] += lfu[j];
	    for(k=0;k<3;k++)
	      s[j][k] += ls[j][k];
	  }
	}
      }
      /* 
	 output 
      */
      fprintf(stdout,"%12g %12g ",x[X],x[Y]);
      /* output of block and fault vel */
      fprintf(stdout,"%17.8e %17.8e %17.8e\t %17.8e %17.8e %17.8e ",
	      bu[X],bu[Y],bu[Z],fu[X],fu[Y],fu[Z]);
      fprintf(stderr,"%s: observation %5i\r",argv[0],n);
      /* output of fault-related stresses (there are no others)  */
      fprintf(stdout,"%17.8e %17.8e %17.8e %17.8e %17.8e %17.8e ",
	      s[X][X],s[X][Y],s[X][Z],
	      s[Y][Y],s[Y][Z],s[Z][Z]);
      fprintf(stdout,"\n");
    }
    fprintf(stderr,"\n");
  }
  exit(0);
}
