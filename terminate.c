/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu


  main initialization routine, sets defaults and optional parameters 
  from the command line
  
  calls initialize, which reads in geometry and boundary conditions

*/
#include "interact.h"
#include "properties.h"

/*

  close files and clean up

 */

int terminate(struct med *medium, struct flt *fault)
{
  int i;
  char tmpstr[STRLEN];
  HEADNODE{
    fprintf(stderr,"terminate: closing files and cleaning up\n");

    if(medium->flt_stress_init){
      for(i=0;i<medium->nrgrp;i++)
	fclose(medium->flt_stress_out[i]);
    }
    if(medium->events_init){
      fclose(medium->events_out); 
      /* purge the rest of the activations, if any */
      if(medium->moment_release_init)
	flush_moment_stack(medium);
      fclose(medium->cevents_out);
    }
    if(medium->slip_line_init)
      for(i=0;i<medium->nrgrp;i++)
	fclose(medium->slip_line_out[i]);
  
    if(medium->read_int_mat_from_file){// I matrix was written to file
      fclose(medium->i_mat_in);
      if(medium->save_imat)// we want to save it
	fprintf(stderr,"terminate: WARNING: leaving I matrix files \"%s\" and \"%s\" (big?)\n",
		medium->mfname,medium->hfname);
      else{// we shall remove the I matrix
	HEADNODE
	  fprintf(stderr,"terminate: removing the interaction matrix files \"%s\" and \"%s\"\n",
		  medium->mfname,medium->hfname);
	snprintf(tmpstr,STRLEN,"rm  %s %s.hdr",medium->mfname,medium->hfname);
	mysystem(tmpstr);
      }
    }else{// I matrix was kept in memory but might still be on file
      ;
    }
  }
  HEADNODE{
#ifdef USE_PGPPLOT
    close_plot_window(medium,fault);
#endif
  }
#ifdef USE_PETSC
  PetscCall(PetscFinalize());
#endif
  exit(medium->op_state);
}
