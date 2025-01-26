/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: quake.c,v 2.25 2003/01/07 02:52:29 tbecker Exp $
*/
#include "interact.h"
#include "properties.h"
/* 
   - add displacement to fault,
   - add coseismic stresses to other faults 
   - write event to file
   - if set, plot quake on X window

   if add_stress is not set, will not add slip and stress change
   if mark_quake is not set, will not add to plotting and moment lists

   sma: activation codes
   slip: slip vector, both are in strike, dip, normal frame
   r_flt: slipping fault number
   

*/


void quake(my_boolean *sma,COMP_PRECISION *slip,int r_flt,struct flt *fault,
	   struct med *medium,my_boolean add_stress, my_boolean mark_quake)
{
#ifdef USE_PGPLOT
  int offset,j;
#endif
  int i;
  float fslip[3];
  COMP_PRECISION cslip[3];
  COMP_PRECISION mom;
  //
  if(!medium->moment_release_init){
    /* first time around, initialize fields */
    medium->moment_release_init=TRUE;
    if(medium->events_init)
      medium->old_moment_time=medium->time;
  }
  if(add_stress)
    add_quake_stress(r_flt,sma,slip,fault,medium);
  /*

    register the slip event in event files and the like (a 'real earthquake')

  */
  if(mark_quake){
    // calculate moment of slipping patch
    mom  = norm_3d(slip) * fault[r_flt].area * SHEAR_MODULUS;
    // write single event to  event files
    if(medium->events_init){
      // single patch 
      for(i=0;i<3;i++)
	fslip[i]=(float)slip[i];
      write_patch_event_file((float)medium->time,medium->iter, r_flt,(float)mom,
			      fslip,medium->events_out);
     }
    // flush group output stack from previouse if we are at a different time
    if(medium->time != medium->old_moment_time)
      flush_moment_stack(medium);
    //
    // add new event to total moment
    medium->total_moment += mom;
    //
    // add moment to temporary group moment
    medium->fault_group[fault[r_flt].group].moment += mom;
    //
    // add weighted slip to temporary group slip
    medium->fault_group[fault[r_flt].group].sarea += 
      fault[r_flt].area;
    /* convert to cartesian */
    compute_cartesian_slip(cslip,slip,(fault+r_flt));
    for(i=0;i<3;i++){
      /* centroid location */
      medium->fault_group[fault[r_flt].group].cent[i]   += fault[r_flt].x[i]  * fault[r_flt].area;
      /* strike,dip,normal */
      medium->fault_group[fault[r_flt].group].slip[i]   += slip[i]  * fault[r_flt].area;
      /* x,y,z */
      medium->fault_group[fault[r_flt].group].slip[3+i] += cslip[i] * fault[r_flt].area;
    }
#ifdef USE_PGPLOT
    if(medium->moment_array_init){
      i=(int)((fault[r_flt].pos[1]+0.5)*
	      (COMP_PRECISION)(medium->grp0_n-1)+0.5);
      j=(int)((fault[r_flt].pos[0]+0.5)*
	      (COMP_PRECISION)(medium->grp0_m-1)+0.5);
      offset  = i*medium->grp0_m+j;
      medium->momrel[offset] += mom;
    }
    plot_quake(r_flt,medium,fault);
#endif
  }
}
/*

  deal with the stress and slip increments due to slip on a patch,
  i.e. add it's effect on stress on other faults

  
*/
void add_quake_stress(int r_flt, my_boolean *sma,
		      COMP_PRECISION *slip, 
		      struct flt *fault, struct med *medium)
{
  int i;
  /* 
     add the stress and slip to the fault patch 
     and other that are affected
     
     add slip to rupturing fault 
  */
  for(i=0;i<3;i++)
    if(sma[i])
      fault[r_flt].u[i] += slip[i];
  /* 
 
     add stress contributions of rupture at fault r_flt with slip[] to all other 
     affected faults 

  */
  switch(select_i_coeff_calc_mode(medium)){
  case I_MAT_IN_MEMORY:{
    add_quake_stress_1(sma,slip,r_flt,fault,medium);
    break;
  }
  case I_MAT_ON_FILE:{
    add_quake_stress_2(sma,slip,r_flt,fault,medium);
    break;
  }
  case CALC_I_COEFF_NOW:{
    add_quake_stress_3(sma,slip,r_flt,fault,medium);
    break;
  }
  case SPARSE_I_MAT:{
    add_quake_stress_4(sma,slip,r_flt,fault,medium);
    break;
  }
  default:{
    fprintf(stderr,"add_quake_stress: can not deal with I matrix code %i\n",
	    select_i_coeff_calc_mode(medium));
    exit(-1);
  }}

}
