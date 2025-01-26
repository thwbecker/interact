/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: test_stuff.c,v 2.12 2003/02/13 22:45:12 becker Exp $
*/
#include "interact.h"
#include <math.h>

#define N 5
#define M 100

#define GAUSS_POINTS 13

int main(int argc, char **argv)
{
  int i,j,ierr;
  FILE *out;
  struct med *medium;
  struct flt *fault;
  char tmpc;
  long seed=-1;

  // intialize random number generator
  RAND_GEN(&seed);

  for(i=0;i<10000;i++)
    fprintf(stdout,"%g\n",mypower_randnr(0.00001,1,-2,0,&seed));

}

