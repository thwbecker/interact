/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu



  $Id: fit_simple_stress_from_cart.c,v 1.3 2011/01/07 07:19:58 becker Exp $

*/
#include "interact.h"
/*

given a cartesian stress field try to fit it by a stress state in
which s1 and s3 are in the horizontal plane, and s2 strictly vertical

*/
int main(int argc, char **argv)
{
  COMP_PRECISION sr[3][3],eval[3],evec[9],sc_ref6[6],sv=0.0,
    best_par[7]={0,0,0,0,0,0,0},s1lim[3],s2lim[3],svlim[3],azilim[3],
    s[3][3],r[3][3],s1h=0.0,s2h=0.0,azi=0.0,sc_ref[3][3],chi2,
    ref_shear,ref_mean;
  int i,j,i3;
  
  /*  read in cartesian stress state in 

  sxx sxy sxz syy syz szz 

  format
  
  */
  fprintf(stderr,"%s: reading six component stress state from stdin\n",
	  argv[0]);
  i = fscanf(stdin, SIX_CP_FORMAT,
	     (sc_ref6),(sc_ref6+1),(sc_ref6+2),(sc_ref6+3),
	     (sc_ref6+4),(sc_ref6+5));
  if(i != 6){
    fprintf(stderr,"%s: read error six entries stdin\n",
	    argv[0]);
    exit(-1);
  }
  convert_6sym_to_9_matrix(sc_ref6,sc_ref);
  /*
    calculate eigenvalues of input stress state
  */
  eigensystem3d(sc_ref,eval,evec);
  /* s1 > s2 > s3 */
  fprintf(stderr,"%s: s1: %g s2: %g s3: %g\n",
	  argv[0],eval[2],eval[1],eval[0]);
  ref_mean =  (eval[0]+eval[1]+eval[2])/3.0; /* mean stress */
  ref_shear = (eval[2] - eval[0])/2.0; /* max shear */
  fprintf(stderr,"%s: mean stress: %11g max shear: %11g\n",
	  argv[0],ref_mean,ref_shear);

  best_par[0] = 1e20;
  s1lim[0] = ref_shear * 2;s1lim[1] = 0.0-1e-7;
  s1lim[2] = -ref_shear * 2 / 200.;
  s2lim[0] = 0;s2lim[1] = -ref_shear * 2 -1e-7;
  s2lim[2] = -ref_shear * 2 / 200.;
  //svlim[0] = -ref_mean;svlim[1]=ref_mean+1e-7;
  svlim[0] = -ref_shear;svlim[1]=ref_shear + 1e-7;
  svlim[2] = ref_shear;
  azilim[0] = 0.0;azilim[1]=180+1e-8;azilim[2]=1.0;
  for(s1h = s1lim[0]; s1h >= s1lim[1]; s1h += s1lim[2])
    for(s2h = s2lim[0]; s2h >= s2lim[1]; s2h += s2lim[2])
      for(sv = svlim[0]; sv <= svlim[1]; sv += svlim[2])
	for(azi = azilim[0]; azi <= azilim[1]; azi += azilim[2]){
	  if((s1h > sv)&&(s2h < sv)){
	    /* for s1h, s2h, azi compute test stress */
	    stress_vec_from_hstate(s1h,s2h,sv,azi,1,evec,eval);
	    for(i=i3=0;i<3;i++,i3+=3){
	      for(j=0;j<3;j++){
		r[j][i] = evec[i3+j];
		s[i][j] = 0.0;
	      }
	      s[i][i] = eval[i];
	    }
	    // get cartesian stresses by rotation
	    rotate_mat(s,sr,r);
	    /* compute misfit */
	    chi2 = stress_misfit(sc_ref,sr);
	    if(chi2 < best_par[0]){
	      best_par[1] = s1h;
	      best_par[2] = sv; 
	      best_par[3] = s2h;
	      best_par[4] = azi;
	      best_par[0] = chi2;
	    }
	  }
	}
  printf("s1h: %12g sv: %12g s2h: %12g as1h: %12g as2h: %12g shs: %12g shm: %12g chi2: %.6e\n",
	 best_par[1],		/* s1h */
	 best_par[2],		/* sv */
	 best_par[3],		/* s2h */
	 best_par[4],		/* azi of s1h */
	 best_par[4]-90.,		/* azi of s2h */
	 (best_par[1] - best_par[3])/2., /* shs */
	 (best_par[1] + best_par[3])/2., /* shm */
	 best_par[0]);
  return 0;
}
