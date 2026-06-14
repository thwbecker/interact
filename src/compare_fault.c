#include "interact.h"
//
// compares fault legths or widths
//
int compare_fault_length(const void *va,const void *vb)
{
  struct flt *a = (struct flt *)va,
    *b = (struct flt *)vb;
  if(a->l > b->l)
    return -1;
  else if(a->l == b->l)
    return 0;
  else 
    return 1;
}
int compare_fault_width(const void *va,const void *vb)
{
  struct flt *a = (struct flt *)va,
    *b = (struct flt *)vb;
  if(a->w > b->w)
    return -1;
  else if(a->w == b->w)
    return 0;
  else 
    return 1;
}
