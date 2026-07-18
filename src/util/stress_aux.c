/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@post.harvard.edu

  auxiliary functions for stress computations

*/
#include "interact.h"
/* 

   compute chi^2 deviation in stress state given 
   two [3][3] matrices 
*/
COMP_PRECISION stress_misfit(COMP_PRECISION sr[3][3], 
			     COMP_PRECISION st[3][3])
{
  COMP_PRECISION c2;
  int i,j;
  c2 = 0.0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c2 += square(sr[i][j] - st[i][j]);
  return c2;
}
/* 
   assign a vector based on s1h s2h and azi assuming a 
   particular stress state

   modes:

   0 plane stress: no stresses in the vertical direction
     s1 = s1h ; s2 = 0 ; s3 = s2h;

   1 s2 in vertical:
     s1 = s1h ; s2 = sv; s3 = s2h;


   azi is in degrees

 */
void stress_vec_from_hstate(COMP_PRECISION s1h, 
			    COMP_PRECISION s2h,
			    COMP_PRECISION sv,
			    COMP_PRECISION azi,
			    int mode,
			    COMP_PRECISION *evec,
			    COMP_PRECISION *eval)
{
  int i;
  COMP_PRECISION strike[3],dip[3];
  
  switch(mode){
  case 0:
    //
    // convert to plane stress state
    //
    strike[0] = azi;strike[1]=   0.0;strike[2]= azi+90.0;
    dip[0] =    0.0;dip[1]   =  90.0;dip[2]   = 0.0;
    eval[0] =   s1h;eval[1]  =   0.0;eval[2]  = s2h;
    break;
  case 1:
    //
    // convert to plane stress state
    //
    strike[0] = azi;strike[1]=   0.0;strike[2]= azi+90.0;
    dip[0] =    0.0;dip[1]   =  90.0;dip[2]   = 0.0;
    eval[0] =   s1h;eval[1]  =   sv;eval[2]  = s2h;
    break;
  default:
    fprintf(stderr,"stress_vec_from_hstate: error: mode %i undefined\n",
	    mode);
    exit(-1);
    break;
  }
  // convert from strike/dip to vector
  for(i=0;i<3;i++){
    angles_to_vec(dip[i],strike[i],(evec+i*3));
  }
}

//
// go from short to 3x3 matrix storage, assumes that s[6] hold
// the stress tensor components on input in xx xy xz yy yz zz
// format. on output, will hold the full, symmetric matrix
// in FORTRAN convention 3x3 
//
void expand_stress_matrix6to9(COMP_PRECISION *s)
{
  COMP_PRECISION sloc[9];
  sloc[INT_XX] = s[0];
  sloc[INT_YX] = sloc[INT_XY] = s[1];
  sloc[INT_ZX] = sloc[INT_XZ] = s[2];
  sloc[INT_YY] = s[3];
  sloc[INT_ZY] = sloc[INT_YZ] = s[4];
  sloc[INT_ZZ] = s[5];
  a_equals_b_vector(s,sloc,9);
}
/* 
   given a symmetric matrix in xx, xy, xz, yy, yz, zz storage 
   convert to [3][3] matrix format
*/
void convert_6sym_to_9_matrix(COMP_PRECISION *s6, 
			      COMP_PRECISION s[3][3])
{
  s[INT_X][INT_X] = s6[0];s[INT_X][INT_Y] = s6[1];s[INT_X][INT_Z] = s6[2];
  s[INT_Y][INT_X] = s6[1];s[INT_Y][INT_Y] = s6[3];s[INT_Y][INT_Z] = s6[4];
  s[INT_Z][INT_X] = s6[2];s[INT_Z][INT_Y] = s6[4];s[INT_Z][INT_Z] = s6[5];

}

/* 

   returns a single (eg., shear) 
   stress component given a matrix sm,
   normal on plane norm, and a tangential vector tang

*/
COMP_PRECISION resolved_stress(COMP_PRECISION *norm,COMP_PRECISION sm[3][3],
			       COMP_PRECISION *tang)
{
  COMP_PRECISION trac[3];
  resolve_force(norm,sm,trac);
  return(dotp_3d(trac,tang));
}
/* 
   - calculates the traction vector on a given plane as
   specified by a normal vector
   - resolves the traction vector onto three different
   directions, typically strike, dip, and normal
   to return the shear and normal stresses
   
   this routine is more efficient if three values are 
   needed as the traction vector can be reused
   
*/

void calc_three_stress_components(COMP_PRECISION sm[3][3],
				  COMP_PRECISION *normal_vec,
				  COMP_PRECISION *vec1,
				  COMP_PRECISION *vec2,
				  COMP_PRECISION *vec3,
				  COMP_PRECISION *s1,
				  COMP_PRECISION *s2,
				  COMP_PRECISION *s3)
{
  COMP_PRECISION trac[3];
  resolve_force(normal_vec,sm,trac);
  *s1 = dotp_3d(trac,vec1);
  *s2 = dotp_3d(trac,vec2);
  *s3 = dotp_3d(trac,vec3);
}
