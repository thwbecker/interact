#include "interact.h"

/* 
   reset all eqn systems and counters 
*/

void init_equation_system(struct med *medium,struct flt *fault)
{
  int i;
  static my_boolean initialized=FALSE;
  // reset the active equation counters for normal and
  medium->naflt=medium->nreq=0;
  // positivity constraint equations
  medium->naflt_con=medium->nreq_con=0;
  if(!initialized){// first time allocation
    initialized=TRUE;
  }else{
    // were allocated already, need to free them first
    // why did the reallocate version cease to work at 
    // some point?
    free(medium->sma);free(medium->sma_con);
    free(medium->nameaf);free(medium->nameaf_con);
    free(medium->rhs_b);free(medium->rhs_b_con);
    free(medium->xsol);free(medium->xsol_con);
  }
  medium->sma=(my_boolean *)malloc(3*sizeof(my_boolean));
  medium->sma_con=(my_boolean *)malloc(3*sizeof(my_boolean));
  for(i=0;i<3;i++){
    medium->sma[i]=INACTIVE;
    medium->sma_con[i]=INACTIVE;
  }
  medium->nameaf=(int *)calloc(1,sizeof(int));
  medium->nameaf_con =(int *)calloc(1,sizeof(int)); 
  medium->rhs_b=(A_MATRIX_PREC *)calloc(1,sizeof(A_MATRIX_PREC));
  medium->xsol=(A_MATRIX_PREC *)calloc(1,sizeof(A_MATRIX_PREC));
  medium->rhs_b_con=(A_MATRIX_PREC *)calloc(1,sizeof(A_MATRIX_PREC));
  medium->xsol_con=(A_MATRIX_PREC *)calloc(1,sizeof(A_MATRIX_PREC));
  if(medium->debug)
    if((!medium->sma)||(!medium->sma_con)||(!medium->nameaf)||
       (!medium->nameaf_con)||(!medium->rhs_b)||(!medium->rhs_b_con)||
       (!medium->xsol)||!(medium->xsol_con))
      MEMERROR("solve: init_equation_system:");
  // reset all fault activation flags
  for(i=0;i<medium->nrflt;i++)
    fault[i].active=FALSE;
  // reset all group activation flags
  medium->nr_active_groups=0;
  for(i=0;i<medium->nrgrp;i++)
    medium->fault_group[i].active=FALSE;
}

/* 

   add fault to active fault list and increment counter 

*/
void add_to_active_fault_list(int aflt,int **al,int *naf,my_boolean **sma)
{
  int ip;
  /* assign active fault name */
  (*al)[*naf]=aflt;
#ifdef DEBUG
  // check codes which have been assigned already
  if((*(*sma+ *naf*3)>3)||(*(*sma+ *naf*3+1)>3)||(*(*sma+ *naf*3+2)>3)){
    fprintf(stderr,"add_to_active_fault_list: fault %i: code screw up: sma 1/2/3: %i %i %i\n",
	    *naf,*(*sma+ *naf*3),*(*sma+ *naf*3+1),*(*sma+ *naf*3+2));
    exit(-1);
  }
#endif
  // increment number of active faults, now *naf is new index
  *naf += 1;
  /* grow slide mode array */
  if((*sma=(my_boolean *)realloc(*sma,sizeof(my_boolean)*3*
			      (*naf+1)))==NULL)
    MEMERROR("add_to_active_fault_list: 1");
  // zero out new codes
  ip = *naf * 3;
  *(*sma + ip + STRIKE)=INACTIVE;
  *(*sma + ip + DIP)=   INACTIVE;
  *(*sma + ip + NORMAL)=INACTIVE;
  /* grow active fault array */
  if((*al=(int *)realloc(*al,sizeof(int)*(*naf+1)))==NULL)
    MEMERROR("add_to_active_fault_list: 2");
}

/* 

   add to x value to right hand side (b), grow b, grow x, 
   and increment counter 

*/
void add_to_right_hand_side(COMP_PRECISION bval,A_MATRIX_PREC **b,
			    A_MATRIX_PREC **xsol, int *nreq)
{
  size_t i;
  (*b)[*nreq]=(A_MATRIX_PREC)bval;
  *nreq += 1;//increment equation counter
  // resizing size
  i= (*nreq + 1)*sizeof(A_MATRIX_PREC);
  // grow b
  if((*b=(A_MATRIX_PREC *)realloc(*b,i))==NULL)
    MEMERROR("add_to_right_hand_side:");
  // grow x
  if((*xsol=(A_MATRIX_PREC *)realloc(*xsol,i))==NULL)
    MEMERROR("add_to_right_hand_side:");
#ifdef SUPER_DUPER_DEBUG
  fprintf(stderr,"add_to_right_hand_side: adding %20.10e as eq. %5i\n",
	  bval,*nreq);
#endif
}
