/*
################################################################################
#                                                                              #
#  interact: model fault interactions using dislocations in an elastic         #
#            medium                                                            #
#	                                                                       #
#    Copyright (c) Thorsten W. Becker 2000 - 2025                              #
#                                                                              #
#    Interact uses dc3d.f as provided by Y. Okada as in Okada (BSSA, 1992) and #
#    linear algebra routines from SLATEC, LAPACK, and EISPACK.                 #
#    Contains f90 converted routines from Nikkhoo and Walter (GJI, 2015)       #
#    Might contain material copyrighted by others (e.g. Numerical Recipes etc) #
#                                                                              #
#    This program is free software; you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License as published by      #
#    the Free Software Foundation; either version 2 of the License, or         #
#    (at your option) any later version.                                       #
#                                                                              #
#    This program is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    In addition, we strongly suggest that iGMT users comply with the goals    #
#    as expressed in the Student Pugwash Pledge (www.spusa.org/pugwash/).      #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with this program; see the file COPYING.  If not, write to          #
#    the Free Software Foundation, Inc., 59 Temple Place - Suite 330,          #
#    Boston, MA 02111-1307, USA.                                               #
#                                                                              #
################################################################################


this is the main  subroutine

for documentation, type 

interact -h



*/
#include "interact.h"
#include "properties.h"



int main(int argc, char **argv)
{
  struct flt *fault;
  struct med *medium;
  // values for background stress
  my_boolean read_initial_fault_stress,faults_have_slipped;
  COMP_PRECISION a[6],b[6];
  char time_out_string[20];
  medium=(struct med *)calloc(1,sizeof(struct med)); /* init as zeros */
  /* time at start */
  clock_gettime(CLOCK_REALTIME, &(medium->init_time));
  strftime(time_out_string, 20, "%Y-%m-%d %H:%M:%S",
	   localtime(&medium->init_time.tv_sec));
#ifdef USE_PETSC
  PetscFunctionBegin;
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  if(medium->comm_size == 0)
    medium->comm_size = 1;	/* fix for non MPI call? */
#else
  medium->comm_size = 1;
#endif
  HEADNODE{
    fprintf(stderr,"%s: internal %s prec, A matrix %s prec\nmain: initializing on %s ncore: %i\n",
	    argv[0],(sizeof(COMP_PRECISION)==sizeof(double))?("double"):("single"),
	    (sizeof(A_MATRIX_PREC)==sizeof(double))?("double"):("single"),
	    time_out_string,medium->comm_size);
    fprintf(stderr,"%s: nu: %.5f mu: %.3e from properties.h, therefore lambda/mu: %.5f alpha: %.5f\n",
	    argv[0],POISSON_NU,SHEAR_MODULUS,LAMBDA_CONST/SHEAR_MODULUS, ALPHA_CONST);
  }
  check_parameters_and_init(argc,argv,&medium,&fault,&read_initial_fault_stress,a,b);
  // decide which mode we are in and execute main loops
  switch(medium->op_mode){
    /*

      LOADING SIMULATION

    */
  case SIMULATE_LOADING_AND_PLOT:
  case SIMULATE_LOADING:{
    /* 
       loading simulation is only serial 
    */
    if(medium->comm_size > 1){
      fprintf(stderr,"main: loading simulation only set up for serial, but CPU %i/%i",
	      medium->comm_rank,medium->comm_size);
      terminate(medium,fault);
    }
    if(medium->nr_flt_mode != 1){
      fprintf(stderr,"main: for simulate loading, only one fault mode can be used\n");
      exit(-1);
    }
    if(!medium->read_int_mat_from_file && medium->save_imat)
      print_interaction_matrix(medium,fault,FALSE);   
    // initialize the stress state of the faults based on background stress
    // and possible past activations if restart attempt
    initialize_stress_state(fault,medium,read_initial_fault_stress,a,b);
    print_fault_stress_stat(medium->stress_stat_out,0,medium,fault);
    /* 
       main loading loop 
    */
    fprintf(stderr,"main: running simulation\n");
    while(medium->time <= medium->stop_time){
      /* check for activations of faults and slip them */
      faults_have_slipped=activate_faults(fault,medium);
      // stress on single fault files at regular intervals
      print_fault_stress_and_slip(medium,fault,faults_have_slipped);
      /* slip along fault at constant intervals */
      print_slip_line(medium,fault);
      if(faults_have_slipped)
	// fault group related stress statistics
	print_fault_stress_stat(medium->stress_stat_out,0,medium,fault);
      /* X Window output */
#ifdef USE_PGPLOT
      update_plots(medium,fault);
#endif
      if(medium->variable_time_step)
	adjust_time_step(fault,medium);
      /* 
	 update stresses on faults due to background loading, ie. 
	 increment them by the current medium->dt
      */
      update_stress_state(fault,medium);
      // and update the time
      medium->time += medium->dt;
      medium->nr_timesteps++;
    }
    fprintf(stderr,"main: done, total timesteps: %i",medium->nr_timesteps);
    if(medium->variable_time_step)
      fprintf(stderr,", avg. dt: %g",
	      medium->stop_time/(COMP_PRECISION)medium->nr_timesteps);
    fprintf(stderr,", max iter/timestep: %i total iter: %i\n",
	    medium->max_iter_realized,medium->total_iter);
    print_fault_data(ONE_STEP_FAULT_DATA_FILE,medium,fault);
    break;
  }
  case ONE_STEP_CALCULATION:{
    /*

      SINGLE STEP CALCULATION
 

      if we had one-step boundary conditions, solve the equations and
      add solution without calling the quake output routine
       
    */
    HEADNODE{
      if(medium->int_mat_init)
	print_interaction_matrix(medium,fault,FALSE);
    }
    if((medium->naflt)||(medium->naflt_con)){
      solve(medium,fault);
#ifdef DEBUG
      HEADNODE
	fprintf(stderr,"main: now adding %i unconst. and %i const. solutions\n",
		medium->nreq,medium->nreq_con);
#endif
      // don't add to moment list, calculate the stress change, add use
      // the normal sign convention: FT 1.0
      if(medium->naflt)
	add_solution(medium->naflt,medium->sma,medium->xsol,
		     medium->nameaf,medium,fault,FALSE,TRUE,1.0);
      if(medium->naflt_con)
	add_solution(medium->naflt_con,medium->sma_con,
		     medium->xsol_con,medium->nameaf_con,
		     medium,fault,FALSE,TRUE,1.0);
    }
    HEADNODE{
      print_fault_data(ONE_STEP_FAULT_DATA_FILE,medium,fault);
    }
    if(medium->print_bulk_fields || medium->read_oloc_from_file){
      calc_fields(medium,fault,FALSE,FALSE,a,b); /* all nodes compute,
						    but only head node
						    will have all
						    entries */
      HEADNODE{
#ifdef DEBUG
	fprintf(stderr,"main: field calculation ok\n");
	fprintf(stderr,"main: output of solution\n");
#endif
	print_stress(medium,fault);
	print_displacement(medium,fault);
	/*  print_stress_on_fault(medium,fault,0);  */
      }
    }
   
    break;
  }}
  terminate(medium,fault);
}
