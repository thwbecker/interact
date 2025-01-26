#include "interact.h"

void adjust_medium_for_restart(struct med *medium, struct flt *fault)
{
  FILE *in;
  int i,nriter,aflt,iec=0;
  char filename[STRLEN];
  float time,slip[3],mom;
  COMP_PRECISION dslip[3];
  my_boolean sma[3];
#ifdef BINARY_PATCH_EVENT_FILE
  strcpy(filename,RESTART_EVENT_FILE_BINARY);
#else
  strcpy(filename,RESTART_EVENT_FILE_ASCII);
#endif
  if(!(in=fopen(filename,"r"))){
    fprintf(stderr,"adjust_medium_for_restart: can't open restart event file %s\n",
	    filename);
    fprintf(stderr,"adjust_medium_for_restart: WARNING: restart attempt abvorted\n");
    return;
  }
  /*
    
    read past events and apply the stresses. this closely mimicks quake.c BUT
    it is not identical for technical reasons

   */
  while(read_patch_event_file(&time,&nriter,&aflt,&mom,slip,in,medium)==7){
    // cumulative moment (not for faults)
    medium->total_moment += (COMP_PRECISION)mom;
    // individual events
    for(i=0;i<3;i++)
      if(slip[i] != 0.0){
	sma[i] = TRUE;
	dslip[i] = (COMP_PRECISION)slip[i];
      }else{
	sma[i] = FALSE;
	dslip[i] = 0.0;
      }
    add_quake_stress(aflt,sma,dslip,fault,medium);
    iec++;
  }
  fclose(in);
  fprintf(stderr,"adjust_medium_for_restart: read %i events and accounted for stresschanges \n",
	  iec);
  // update the time and the background stress
  medium->time = (COMP_PRECISION)time;
  if(medium->dt0 != 1.0){
    fprintf(stderr,"adjust_medium_for_restart: dt0 should be 1\n");
    exit(-1);
  }
  for(i=0;i<medium->nrflt;i++){// background
    fault[i].s[STRIKE] += medium->time * fault[i].sinc[STRIKE];
    fault[i].s[DIP]    += medium->time * fault[i].sinc[DIP];
    fault[i].s[NORMAL] += medium->time * fault[i].sinc[NORMAL];
  }
  fprintf(stderr,"adjust_medium_for_restart: updated time and background stress to t: %g\n",
	  medium->time);
  if(medium->time >= medium->stop_time){
    fprintf(stderr,"adjust_medium_for_restart: past activation time would be > stop time (%g)\n",
	    medium->stop_time);
    fprintf(stderr,"adjust_medium_for_restart: adjust time in \"%s\"?\n",
	    BC_FILE);
  }
}
