/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: myopen.c,v 2.2 2002/10/08 19:24:44 tbecker Exp $
*/
#include "interact.h"
/*
  open files safely 
*/

FILE *myopen(char *name,char *modus)
{
  FILE *tmp;
  if((tmp=(FILE *)fopen(name,modus))==NULL){	
    fprintf(stderr,"myopen: error opening \"%s\"\n",name);
    exit(-1);
  }
  return ((FILE *)tmp);
}	
