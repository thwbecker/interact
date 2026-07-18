/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu


*/
#include "interact.h"
#include <string.h>

my_boolean strings_match(char *a, char *b)
{
  if(strncmp(a,b,(sizeof(a)>sizeof(b))?sizeof(a):sizeof(b)) == 0)
    return TRUE;
  else
    return FALSE;
}
