/*
  part of blockinvert

  
  read in stress observations

  sx[nrsp*BLOCK_DIM] will hold the lon lat of the observations
  st[6*nrsp] will hold the stress tensor at each point
  and sigv the variance likewise

  if no_stress_amp_scale is set, will only check for trace == 0,
  else shear magnitude has to be unity

  $Id: read_stress_observations.c,v 1.8 2003/07/10 22:21:22 becker Exp $
    
*/
#include "interact.h"
#include "blockinvert.h"

void read_stress_observations(struct bmd *mod,
			      COMP_PRECISION *stressrms,
			      COMP_PRECISION minssig,
			      my_boolean no_stress_amp_scale,
			      COMP_PRECISION **stress_depths,
			      COMP_PRECISION global_stress_depth,
			      my_boolean load_stress_depths)
{
  int npdim,np6,i;
  FILE *in;
  my_boolean warned=FALSE;
  COMP_PRECISION s[9],sval[3],dummy,lon,lat;
  /* set stress observation counter to zero  */
  mod->nrsp = 0;
  in = fopen(STRESSDATA_FILE,"r");
  if(!in){
    fprintf(stderr,"read_stress_observations: couldn't open %s, no stresses\n",
	    STRESSDATA_FILE);
    return;
  }
  my_vecrealloc(&mod->sx,BLOCK_DIM,"sx");
  my_vecrealloc(&mod->v,mod->mgd+6,"v");
  my_vecrealloc(&mod->sigv,mod->mgd+6,"sigv");
  /*
    
    expecting 
    
    lon lat sxx sxy sxz syy syz szz  sig(sxx) sig(sxy) sig(sxz) sig(syy) sig(syz) sig(szz)
    
    format
    
  */
  npdim = 0;
  np6 = mod->mgd;
  while(fscanf(in,"%lf %lf",(mod->sx+npdim+INT_X),
	       (mod->sx+npdim+INT_Y)) == 2){
    for(i=0;i < 6;i++)		/* read in tensor */
      if(fscanf(in,"%lf",(mod->v+np6+i)) != 1){
	fprintf(stderr,"read_stress_observations: read error\n");
	exit(-1);
      }
    for(i=0;i < 6;i++)		/* read in variance */
      if(fscanf(in,"%lf",(mod->sigv+np6+i)) != 1){
	fprintf(stderr,"read_stress_observations: read error: i: %i\n",
		mod->nrsp+1);
	exit(-1);
      }else{
	if(mod->sigv[np6+i] <= minssig){
	  fprintf(stderr,"read_stress_observations: WARNING: sigmas shouldn't be <= %g: i: %i sig: %g, adjusting to %g\n",
		  minssig,mod->nrsp+1,mod->sigv[np6+i],minssig);
	  mod->sigv[np6+i] = minssig;
	}
      }
    /* 
       check if given as trace == 0,
       and 
       max shear stress == 1 OR no_stress_amp_scale == TRUE
    */
    a_equals_b_vector(s,(mod->v+np6),6);
    expand_stress_matrix6to9(s);
    calc_eigensystem_sym3d(s,sval,&dummy,FALSE);
    if(((dummy = sval[EIGEN_E1] - sval[EIGEN_E3])>1.001)||
       (dummy < 0.999)){	/* shear stress != 1  */
      if(fabs(dummy) < EPS_COMP_PREC){ /* zero shear */
	if(!warned){
	  fprintf(stderr,"read_stress_observations: WARNING: max shear: %g at least at pt %i, can not adjust to unity\n",
		  dummy,np6/6+1);
	  exit(-1);
	}
      }else{
	if(no_stress_amp_scale){
	  if(!warned){
	    fprintf(stderr,"read_stress_observations: WARNING: max shear: %g at least at pt %i, will not adjust\n",
		    dummy,np6/6+1);
	    warned=TRUE;
	  }
	}else{			/* rescale this observation to unity */
	  if(!warned){
	    fprintf(stderr,"read_stress_observations: WARNING: max shear: %g at least at pt %i, adjusting to unity\n",
		    dummy,np6/6+1);
	    warned=TRUE;
	  }
	  // rescale
	  dummy = 1.0/dummy;
	  scale_vector((mod->v+np6),dummy,6);
	  scale_vector(s,dummy,9);
	}
      }
    }
    if(fabs(dummy = tracemat9(s)) > 0.001){ /* trace is non-zero */
      fprintf(stderr,"read_stress_observations: WARNING: trace: %g at least at pt %i, adjusting to zero\n",
	      tracemat9(s),np6/6+1);
      dummy /= 3.0;
      *(mod->v+np6)   -= dummy;
      *(mod->v+np6+3) -= dummy;
      *(mod->v+np6+5) -= dummy;
    }
    /*

      increment and allocate more space

    */
    mod->nrsp++;
    npdim += BLOCK_DIM;
    np6 += 6;
    
    my_vecrealloc(&mod->sx,BLOCK_DIM+npdim,"sx");
    my_vecrealloc(&mod->v,(6+np6),"v");
    my_vecrealloc(&mod->sigv,(6+np6),"v");
  }
  fclose(in);
  /* number of rows */
  mod->m2 = mod->nrsp * 6;
  /* determine RMS  */
  *stressrms = rms((mod->v+mod->mgd),mod->m2);
  fprintf(stderr,"read_stress_observations: read %i stress observations, rms: %g\n",
	  mod->nrsp,*stressrms);
  my_vecalloc(stress_depths,mod->nrsp,"stress_depths");
  if(load_stress_depths){	/* read stress evaluation depths
				   in lon lat depth format
				*/
    in=myopen(STRESSDEPTH_FILE,"r");
    for(i=0;i < mod->nrsp;i++){
      if(fscanf(in,"%lf %lf %lf",&lon,&lat,(*stress_depths+i))!=3){
	fprintf(stderr,"read_stress_observations: read error stress.depth file, observation %i out of %i\n",
		i+1,mod->nrsp);
	exit(-1);
      }
      if(dist_on_sphere_deg(lon,lat,
			    *(mod->sx+i*BLOCK_DIM+INT_X),
			    *(mod->sx+i*BLOCK_DIM+INT_Y))>
	 1e-7){
	fprintf(stderr,"read_stress_observations: read error stress.depth file\n");
	fprintf(stderr,"read_stress_observations: location of entry %i: %g %g\n",
		i+1,lon,lat);
	fprintf(stderr,"read_stress_observations: location of stress observation %i: %g %g\n",
		i+1,*(mod->sx+i*BLOCK_DIM+INT_X),*(mod->sx+i*BLOCK_DIM+INT_Y));
	exit(-1);
      }
    }
    fclose(in);
    fprintf(stderr,"read_stress_observations: read stress evaluation depths from stress.depth, mean: %g\n",
	    mean(*stress_depths,1,mod->nrsp));
  }else{
    for(i=0;i < mod->nrsp;i++)
      *(*stress_depths+i) = global_stress_depth;
    fprintf(stderr,"read_stress_observations: all stresses evaluated at %g km depth\n",
	    global_stress_depth);
  }
  for(i=0;i < mod->nrsp;i++)
    if(*(*stress_depths+i) < 0){
      fprintf(stderr,"read_stress_observations: error: stress observation depth %i: %g (has to be >0)\n",
	      i+1,*(*stress_depths+i));
      exit(-1);
    }
#ifdef BLOCK_SPHERICAL
  /* copy stress data and uncertainties over to vc and sigvc  */
  my_vecrealloc(&mod->vc,mod->m1+mod->m2,"vc");
  my_vecrealloc(&mod->sigvc,mod->m1+mod->m2,"vc");
  a_equals_b_vector((mod->vc+mod->m1),(mod->v+mod->mgd),mod->m2);
  a_equals_b_vector((mod->sigvc+mod->m1),(mod->sigv+mod->mgd),mod->m2);
#endif
}
