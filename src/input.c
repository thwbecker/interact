/*
  interact: model fault interactions using dislocations in a 
            halfspace
	    (c) Thorsten Becker, thbecker@post.harvard.edu


  
  some input routines, see also matrixio.c

*/
#include "interact.h"
#include "properties.h"
//
// read in the group moment release file as floats for plotting
// time, moment, and nr of events
//
int read_moment_file(float **tm, float **mr,float *minm, float *maxm, int *n,
		     my_boolean use_log_mom)
{
  FILE *in;
  size_t size;
  in=myopen(CEVENT_FILE,"r");
  size=sizeof(float);
  *tm=(float *)malloc(size);
  *mr=(float *)malloc(size);
  *n=0; 
  *minm = FLT_MAX; 
  *maxm = -FLT_MAX;
  // time grp mom tot_mom slip_strike slip_dip slip_normal
  while(fscanf(in,"%f %*i %f %*f %*f %*f %*f\n",(*tm+ *n),(*mr + *n)) == 2){
    if(use_log_mom)
      *(*mr+ *n) = log(*(*mr+ *n));
    if(*(*mr+ *n) > *maxm) *maxm = *(*mr+ *n);
    if(*(*mr+ *n) < *minm) *minm = *(*mr+ *n);
    *n += 1;
    size += sizeof(float);
    *tm=(float *)realloc(*tm, size);
    *mr=(float *)realloc(*mr, size);
    if(! *tm || ! *mr)
      MEMERROR("read_moment_file");
   }
  size -= sizeof(float);
  *tm=(float *)realloc(*tm, size);
  *mr=(float *)realloc(*mr, size);
  fclose(in);
  fprintf(stderr,"read_moment_file: %i events, (%sM)_min/max: %15.4e/%15.4e, t_max: %g\n",
	  *n,use_log_mom?"ln":"",*minm,*maxm,*(*tm+(*n-1)));
  return(*n);
}
/*

  read/write single patch activations to file
  for successful read, will return integer 7
*/
int read_patch_event_file(float *time,int *nriter, int *aflt,  float *mom, 
			  float *slip,FILE *in,struct med* medium)
{
  int cnt=0;
#ifdef BINARY_PATCH_EVENT_FILE
  cnt += fread(time,sizeof(float),1,in);
  cnt += fread(nriter,sizeof(int),1,in);
  cnt += fread(aflt,sizeof(int),1,in);
  cnt += fread(slip,sizeof(float),3,in);
  cnt += fread(mom,sizeof(float),1,in);
#ifdef DEBUG
  if(*aflt > medium->nrflt){
    fprintf(stderr,"read_patch_event: IO error? aflt (%i) > nrflt (%i)\n",
	    *aflt,medium->nrflt);
    exit(-1);
  }
#endif
#else
  cnt  = fscanf(in,"%f %i %i %f %f %f %f",time,nriter,aflt,
		slip,(slip+1),(slip+2),mom);
#endif
  return(cnt);

}
int  write_patch_event_file(float time, int nriter, int aflt, float mom,
			    float *slip, FILE *out)
{
  int cnt=0;
#ifdef BINARY_PATCH_EVENT_FILE
  cnt += fwrite(&time,sizeof(float),1,out);
  cnt += fwrite(&nriter,sizeof(int),1,out);
  cnt += fwrite(&aflt,sizeof(int),1,out);
  cnt += fwrite(slip,sizeof(float),3,out);
  cnt += fwrite(&mom,sizeof(float),1,out);
#else
  cnt = fprintf(out,"%20.15e %3i %6i %14.6e %14.6e %14.6e %14.6e\n",
		time,nriter,aflt,
		slip[STRIKE],slip[DIP],slip[NORMAL],mom);
#endif
  return cnt;
}
/* 
   read in rate state friction parameters 
*/
void read_rsf(char *filename, struct med *medium, struct flt *fault)
{
  int i;
  FILE *in;
  COMP_PRECISION aavg,bavg;
  in = myopen(filename,"r");
  aavg=bavg=0;
  for(i=0;i < medium->nrflt;i++){ /* read a and b */
    if(fscanf(in,TWO_CP_FORMAT,&(fault[i].mu_s),&(fault[i].mu_d))!=2){
      fprintf(stderr,"read_rsf: error reading a b for patch %i\n",i);
      exit(-1);
    }
    aavg += fault[i].mu_s;
    bavg += fault[i].mu_d;
  }
  fclose(in);
  aavg /= (COMP_PRECISION)medium->nrflt;
  bavg /= (COMP_PRECISION)medium->nrflt;
  fprintf(stderr,"read_rsf: read a and b values, on averag %g and %g\n",aavg,bavg);
}
