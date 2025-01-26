/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: fltcopy.c,v 2.7 2002/10/08 19:24:44 tbecker Exp $
*/
#include "interact.h"

void fltswap(struct flt *a, struct flt *b)
{
  struct flt c;
  fltcp(a,&c);
  fltcp(b,a);
  fltcp(&c,b);
}
//
// copy fault a to fault b
//
// b = a
//
void fltcp(struct flt *a, struct flt *b)
{
  int i;
  a_equals_b_vector_3d(b->u,a->u);
  a_equals_b_vector_3d(b->s,a->s);
  a_equals_b_vector_3d(b->sinc,a->sinc);
  a_equals_b_vector_3d(b->x,a->x);

#ifdef ALLOW_NON_3DQUAD_GEOM
  b->type = a->type;
  if(a->type == TRIANGULAR){
    a_equals_b_vector(b->xt,a->xt,9);
  }
#endif 
  for(i=0;i<2;i++){
    b->pos[i] = a->pos[i];
    b->cf[i] = a->cf[i];
  }
  b->l = a->l;
  b->w = a->w;
  b->area = a->area;

  b->strike = a->strike;
  b->dip = a->dip;

  b->sin_alpha = a->sin_alpha;
  b->cos_alpha = a->cos_alpha;

  a_equals_b_vector_3d(b->normal,a->normal);
  a_equals_b_vector_3d(b->t_strike,a->t_strike);
  a_equals_b_vector_3d(b->t_dip,a->t_dip);
  
  b->mu_d = a->mu_d;
  b->mu_s = a->mu_s;
  b->f_initial = a->f_initial;
  b->taud = a->taud;

#ifdef LATENCY
  b->last_activation_time = a->last_activation_time;
#endif
  
  b->group = a->group;
  b->active = a->active;

  for(i=0;i<3;i++)
    b->mode[i] = a->mode[i];
}
