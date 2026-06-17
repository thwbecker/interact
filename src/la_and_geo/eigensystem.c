/*
  interact: model fault interactions using dislocations in a 
  halfspace

  (C) Thorsten Becker, becker@eps.harvard.edu

  routine for calculating the eigensystem of a symmetric 3-D matrix 

  input/ouput array sizing:

  val[0...2]
  vec[0...8]

  uses EISPACK routine rs

  on output, val has the eigenvalues, vec the eigenvectors in ascending order

  $Id: eigensystem.c,v 1.5 2003-01-14 07:23:46-08 becker Exp tbecker $
*/
#include "interact.h"


void eigensystem3d(COMP_PRECISION a[3][3],
		   COMP_PRECISION *val, COMP_PRECISION *vec)
{
  COMP_PRECISION ac[9],fv1[3],fv2[3];
  int i,j;
  // resort matrix for FORTRAN storage
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      ac[j*3+i] = a[i][j];
  eispack_driver(ac,3,3,val,vec,fv1,fv2,1);
}


	    
/*
  
  given real, symmetric a 3 x 3 matrix A where only the XX, XY, XZ,
  YY, YZ, and ZZ components have to be filled (the corresponding other
  ones will be overwritten)

  calculate the eigenvalues eval[3]
  and the corresponding eigenvectors, evec is 3 x 3

  val will have the eigenvectors in ascending order, 
  e1 = val[2],e2 = val[1],e3 = val[0]
  or 
  e1 = val[E1],e2 = val[E2],e3 = val[E3]

  with e1 > e2 > e3

  vec will be each normalized to unity 

  r, theta, phi
  
  vec[6,7,8] is the vector that corresponds to the highest  eigenvalue,     val[2]
  vec[3,4,5] is the vector that corresponds to the intermediate eigenvalue, val[1]
  vec[0,1,2] is the vector that corresponds to the smallest eigenvalue,     val[0]


  WARNING:

  if icalc_vectors is not set to TRUE, 
  will only calculate the eigenvalues

  uses EISPACK routine

  we have also defined macros E1, E2, E3 that refer to the indices
  2, 1, and 0, respectively (see above)


  SINCE THIS ROUINTE ONLY WORKS FOR SYMMETRIC MATRICES, IT DOESN'T
  MATTER IF C OR FORTRAN STORAGE IF THE XX,XY,XZ,YY,YZ, AND ZZ
  ENTRIES WERE ASSIGNED


*/
void calc_eigensystem_sym3d(COMP_PRECISION *a,
			    COMP_PRECISION *eval,
			    COMP_PRECISION *evec,
			    my_boolean icalc_vectors)
{
  COMP_PRECISION fv1[3],fv2[3],loca[9];
  a_equals_b_vector(loca,a,9);// save the original A
  eispack_driver(loca,3,3,eval,evec,fv1,fv2,(icalc_vectors)?(1):(0));
}
void calc_eigensystem_sym2d(COMP_PRECISION *a,
			    COMP_PRECISION *eval,
			    COMP_PRECISION *evec,
			    my_boolean icalc_vectors)
{
  COMP_PRECISION fv1[2],fv2[2],loca[4];
  /*
    a = ( xx xy ) = (  0 1 )
        ( yx yy )   (  2 3 )

  */
  a[2] = a[1];
  a_equals_b_vector(loca,a,4);// save the original A
  eispack_driver(loca,2,2,eval,evec,fv1,fv2,(icalc_vectors)?(1):(0));
}

void eispack_driver(COMP_PRECISION *a,int m, int n,
		    COMP_PRECISION *val,COMP_PRECISION *vec,
		    COMP_PRECISION *fv1,COMP_PRECISION *fv2,
		    int matz)
{
  int ierr;
  //
  // call EISPACK routine
  //
  EISPACK_RS(&m,&n,a,val,&matz,vec,fv1,fv2,&ierr);
  if(ierr){
    fprintf(stderr,"eispack_driver: runtime error %i in EISPACK routine\n",
	    ierr);
    exit(-1);
  }
}
	
