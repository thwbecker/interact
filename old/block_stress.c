#include "interact.h"
#include "blockinvert.h"
/* 
   block inversion type routines 
   that deal with stresses

   $Id: block_stress.c,v 1.3 2003/07/10 22:21:22 becker Exp $

*/


/*

  scale the observed stresses, which are given as deviatoric only
  (trace = 0) and normalized such that the max shear stress, 
  e1 - e3, = 1, to the predicted stresses from the model

  m1 is the offset in the observation vector to point to the 
  first stress observation, m2 the total vector length
  nrsp * 6

  m1 is the offset

  will also scale the uncertainties, vsig

  
  if no_rel_stress_amp_scale is TRUE, than only the overall RMS
  will be made similar, not the individual observations
  (trace will still be adapted to model)

*/
void block_scale_stresses(COMP_PRECISION *v,COMP_PRECISION *vsig,
			  COMP_PRECISION *vmod,
			  int m1, int m2, int nrsp,
			  COMP_PRECISION damp_fac,
			  COMP_PRECISION *stressrms,
			  my_boolean no_rel_stress_amp_scale)
{
  
  int i,os;
  COMP_PRECISION modval[3],smod[9],obs_trace,vec[9],obs_shear,
    traceadd,mod_shear,mod_trace,sobs[9],obsval[3],fac,modrms;
  my_boolean rescale_total_rms = FALSE, test_scaling = FALSE;
#ifdef DEBUG
  if(m2/6 != nrsp){
    fprintf(stderr,"block_scale_stresses: error: m1: %i m2: %i nrsp: %i\n",
	    m1,m2,nrsp);
    exit(-1);
  }
#endif
  for(os=m1,i=0;i < nrsp;i++,os+=6){
    /* model stresses */
    a_equals_b_vector(smod,(vmod+os),6);
    expand_stress_matrix6to9(smod);
    /* observed stresses */
    a_equals_b_vector(sobs,(v+os),6);
    expand_stress_matrix6to9(sobs);
    /* isotropic part of the model stresses  */
    mod_trace =  tracemat9(smod);
    if(!no_rel_stress_amp_scale){	/* scale rel magnitude */
      /* max shear of model stresses */
      calc_eigensystem_sym3d(smod,modval,vec,FALSE);
      mod_shear = modval[EIGEN_E1] - modval[EIGEN_E3]; 
      /*        observed stresses  */
      calc_eigensystem_sym3d(sobs,obsval,vec,FALSE);
      obs_shear = obsval[EIGEN_E1]-obsval[EIGEN_E3];
      /* 
	 determine the scaling factor for the stress observation
      */
      if(fabs(obs_shear) > EPS_COMP_PREC)
	fac = mod_shear / obs_shear * damp_fac;
      else
	fac = 1.0;
      scale_vector(sobs,fac,9);	/* scale observed stress matrix */
      // scale uncertainties, too
      scale_vector((vsig+os),fac,6);
    }
    /* 
       ensure that the observed trace corresponds to the
       predicted trace of the stress tensor
    */
    obs_trace =  tracemat9(sobs);
    traceadd = (mod_trace - obs_trace)/3.0; /* to be added */
    sobs[XX] += traceadd;sobs[YY] += traceadd;sobs[ZZ] += traceadd;
    //
    // re-assign to actual observation vector
    //
    v[os  ] = sobs[XX];v[os+1] = sobs[XY];v[os+2] = sobs[XZ];
    v[os+3] = sobs[YY];v[os+4] = sobs[YZ];v[os+5] = sobs[ZZ];
  }
  if((no_rel_stress_amp_scale)||(rescale_total_rms)){
    /* 
       rescale such that the total RMS as given by the model 
       stresses doesn't decrease 
    */
    modrms = rms((vmod+m1),m2);
    *stressrms = rms((v+m1),m2);
    fac = (fabs(*stressrms)>EPS_COMP_PREC)?
      (modrms/(*stressrms)):(1.0);
    scale_vector((v+m1),fac,m2);	/* stresses */
    scale_vector((vsig+m1),fac,m2); /* uncertainties */
  }
  *stressrms = rms((v+m1),m2);
  if((!no_rel_stress_amp_scale)&&(test_scaling)){
    /* test how the scaling worked */
    for(os=m1,i=0;i < nrsp;i++,os+=6){
      a_equals_b_vector(smod,(vmod+os),6);
      a_equals_b_vector(sobs,(v+os),6);
      expand_stress_matrix6to9(sobs);
      expand_stress_matrix6to9(smod);
      calc_eigensystem_sym3d(sobs,obsval,vec,FALSE);
      obs_shear = obsval[EIGEN_E1]-obsval[EIGEN_E3];
      calc_eigensystem_sym3d(smod,modval,vec,FALSE);
      mod_shear = modval[EIGEN_E1]-modval[EIGEN_E3];
      obs_trace =  tracemat9(sobs);
      mod_trace =  tracemat9(smod);
      fprintf(stderr,"scale_stress: obs: s: %12.3e t: %12.3e obs-mod: s: %12.3e t: %12.3e\n",
	      obs_shear,obs_trace,obs_shear-mod_shear,obs_trace-mod_trace);
    }
  }
}
/* calculate a vector that holds the azimuth of the horizontal
   stress components major eigenvector */
void calc_horizontal_stress_vec(COMP_PRECISION *s,
				COMP_PRECISION *azi,
				int n)
{
  int i;
  COMP_PRECISION e1,e2;
  for(i=0;i<n;i++)
    calc_horizontal_stress((s+6*i),&e1,&e2,(azi+i));
}
/* 
   given a matrix in 6 short sotrage, return e1 > e2 and azimuth
   of e1 of the horizontal parts 
*/
void calc_horizontal_stress(COMP_PRECISION *s6,
			    COMP_PRECISION *e1,
			    COMP_PRECISION *e2,
			    COMP_PRECISION *azi)
{
  COMP_PRECISION s[9],sloc2d[4],eval[2],evec[4];
  a_equals_b_vector(s,s6,6);
  expand_stress_matrix6to9(s);
  sloc2d[0] = s[XX];
  sloc2d[1] = sloc2d[2] = s[XY];
  sloc2d[3] = s[YY];
  /* this routine returns values and vectors sorted as
     e2 < e1 */
  calc_eigensystem_sym2d(sloc2d,eval,evec,TRUE);
  *azi = vec_to_strike((evec+2));/* azimuth of e1 */
  fix_azimuth(azi);
  if(*azi > 180)		/* 0 .. azi .. 180 */
    *azi -= 180;
  *e1 = eval[1];*e2 = eval[0];
}
/* calculate a vector with angular differences */
void calc_dir_diff_vec(COMP_PRECISION *a1,
		       COMP_PRECISION *a2, 
		       COMP_PRECISION *da,
		       int n,
		       my_boolean allow_negative)
{
  int i;
  for(i=0;i<n;i++)
    da[i] = calc_dir_diff(a1[i],a2[i],allow_negative);
}
/* 
   calculate the difference between two azimuths of orientationl
   data (180 deg periodicity), azimuths a[2] are in degrees

*/
COMP_PRECISION calc_dir_diff(COMP_PRECISION a1,
			     COMP_PRECISION a2,
			     my_boolean allow_negative)
{
  COMP_PRECISION da,a[2];
  int i;
  a[0]=a1;a[1]=a2;
  for(i=0;i<2;i++){
    while(a[i] < 0)
      a[i] += 360.0;
    while(a[i] > 360.0)
      a[i] -= 360.0;
  }
  da = a[1]-a[0];
  if(allow_negative){
    if(da > 0){
      if(da > 180)da -= 180;
      if(da > 90)da -= 180;
    }else{
      if(da < -180)da += 180;
      if(da < -90)da += 180;
    }
  }else{
    if(da < 0)da = -da;
    while(da > 180.0)da -= 180.0;
    if(da > 90.)da = 180.-da;
  }
  return da;
}
/* 
   given e1h, e2h, and azi, calculate the cartesian stress matrix
   s6[6] is short format

*/
void cart_mat_from_horsym(COMP_PRECISION e1,COMP_PRECISION e2,
			  COMP_PRECISION azi,
			  COMP_PRECISION *s6)
{
  COMP_PRECISION strike[3],dip[3],eval[3],evec[9],
    s[3][3],r[3][3],sr[3][3];
  int i,i3,j;
  /* those are the eigenvector direction */
  strike[0] = azi;strike[1]=0.0;strike[2]=azi+90.0;
  dip[0] = dip[2] = 0.0;dip[1]=90.0;
  eval[0] = e1;eval[1]=0.0;eval[2]=e2;
  /* convert to eigenvector */
  for(i=0;i<3;i++)
    angles_to_vec(dip[i],strike[i],(evec+i*3));
  for(i=i3=0;i<3;i++,i3+=3){
    // assemble rotion matrix, prepare stress matrix
    for(j=0;j<3;j++){
      r[j][i] = evec[i3+j];
      s[i][j] = 0.0;
    }
    s[i][i] = eval[i];
  }
  // get cartesian stresses by rotation
  rotate_mat(s,sr,r);
  /* asign to short format matrix */
  s6[0]=sr[X][X];s6[1]=sr[X][Y];s6[2]=sr[X][Z];
  s6[3]=sr[Y][Y];s6[4]=sr[Y][Z];s6[5]=sr[Z][Z];
}

/* 


   scale 'data' stresses contained in v 
   to those from a solution vector vmod 
   given uncertainties in sigv


*/
void rescale_observed_stresses(COMP_PRECISION *vmod,
			       COMP_PRECISION *v,
			       COMP_PRECISION *sigv,
			       COMP_PRECISION sscale_dfac,
			       COMP_PRECISION *stressrms,
			       my_boolean no_rel_stress_amp_scale,
			       struct bmd *mod,
			       my_boolean restore_amplitudes,
			       my_boolean scale_to_model)
{
  if(restore_amplitudes){
    /*  reassign stresses to original values (have been rescaled) */
    a_equals_b_vector((v+mod->mgd),mod->saved_stress,mod->m2); /* stresses */
    a_equals_b_vector((sigv+mod->mgd),(mod->saved_stress+mod->m2),mod->m2); /* uncertainties */
  }
  if(scale_to_model){
    block_scale_stresses(v,sigv,vmod,mod->mgd,mod->m2,
			 mod->nrsp,sscale_dfac,
			 stressrms,no_rel_stress_amp_scale);
#ifdef BLOCK_SPHERICAL
    /* 
       assign stresses and uncertainties 
       to cartesian solution vectors 
    */
    a_equals_b_vector((mod->vc+mod->m1),(v+mod->mgd),mod->m2);
    a_equals_b_vector((mod->sigvc+mod->m1),
		      (sigv+mod->mgd),mod->m2);
#endif
  }
}
