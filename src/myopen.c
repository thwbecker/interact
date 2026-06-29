/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu

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
