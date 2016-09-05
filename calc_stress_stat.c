/*
  calculates spatial and other statistics for the 
  stress and friction initialization files


  $Id: calc_stress_stat.c,v 1.10 2003/03/02 01:37:41 becker Exp $

*/
#include "interact.h"
#include "properties.h"

int main(int argc, char **argv)
{
  FILE *in,*out;
  int i,j,k,*nr,nr_spc_bins,mode,os[3];
  COMP_PRECISION *s=NULL,stat[20],*xr,*r,range_frac=0.5,
    cr=0.1,xcr,*mu=NULL;
  char tmpstr[STRLEN];
  float *w1,*w2,*w3;
  my_boolean check_stress,check_mu;
  struct flt *fault;
  struct med *medium;

  // defaults
  mode = 0;
  nr_spc_bins=40;
  switch(argc){
  case 1:
    break;
  case 2:
    sscanf(argv[1],"%i",&mode);
    break;
  case 3:
    sscanf(argv[1],"%i",&mode);
    sscanf(argv[2],"%i",&nr_spc_bins);
    break;
  default:
    fprintf(stderr,"%s [mode, %i] [nr_spatial_corr_bins, %i]\n",
	    argv[0],mode,nr_spc_bins);
    fprintf(stderr,"calculates stress and friction init file statistics\n");
    fprintf(stderr,"mode 0: only output of simple quantities\n");
    fprintf(stderr,"mode 1: additional spatial correlations\n");
    fprintf(stderr,"mode 2: output of stress/friction properties to single file\n");
    exit(-1);
    break;
  }
  // read in geometry 
  read_geometry(GEOMETRY_FILE,&medium,&fault,FALSE,FALSE,
		FALSE,FALSE);
  fprintf(stderr,"%s: trying to read stress inits from %s\n",argv[0],
	  FAULT_STRESS_INIT_FILE);
  // read in fault stress initialization file
  in=fopen(FAULT_STRESS_INIT_FILE,"r");
  if(!in){
    check_stress = FALSE;
    fprintf(stderr,"%s: file not found\n",argv[0]);
  }else{
    fprintf(stderr,"%s: ok, doing stress stats\n",argv[0]);
    s=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*3*
			       medium->nrflt);
    if(!s)
      MEMERROR("s");
    check_stress = TRUE;
    /*
      
      the input stresses will be sorted by patch number, 
      increasing for each component

    */
    for(j=medium->nrflt,k=j+medium->nrflt,i=0;
	i<medium->nrflt;
	i++,j++,k++)
      if(fscanf(in,THREE_CP_FORMAT,(s+i),(s+j),(s+k))!=3){
	fprintf(stderr,"%s: read error in %s, patch %i\n",
		argv[0],FAULT_STRESS_INIT_FILE,i+1);
	exit(-1);
      }else{
	fault[i].s[STRIKE] = s[i];
	fault[i].s[DIP] =    s[j];
	fault[i].s[NORMAL] = s[k];
      }
    fclose(in);
  }
  //
  // read the friction coefficients
  //
  fprintf(stderr,"%s: trying to read friction inits from %s\n",argv[0],
	   FAULT_PROP_FILE);
  in=fopen(FAULT_PROP_FILE,"r");
  if(!in){
    check_mu = FALSE;
    fprintf(stderr,"%s: file not found\n",argv[0]);
  }else{
    fprintf(stderr,"%s: ok, doing friction stats\n",argv[0]);
    mu=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*2*medium->nrflt);
    if(!mu)
      MEMERROR("mu");
    check_mu = TRUE;
    for(j=medium->nrflt,i=0;i<medium->nrflt;i++,j++)
      if(fscanf(in,TWO_CP_FORMAT,(mu+i),(mu+j))!=2){
	fprintf(stderr,"%s: read error in %s, patch %i\n",
		argv[0],FAULT_PROP_FILE,i+1);
	exit(-1);
      }else{
	fault[i].mu_s = (float)mu[i];
	fault[i].mu_d = (float)mu[j];
      }
    fclose(in);
  }
  switch(mode){
  case 0:
    //
    // do only simple stats
    //
    for(j=i=0;i<5;i++,j += medium->nrflt){
      if((i<3) && (check_stress)){//stress
	calc_vec_stat((s+j),medium->nrflt,(stat+i*4));
	fprintf(stderr,"%s: stress dir %i: min: %12g max: %12g mean: %12g stddev: %12g (%5.2f times SD)\n",
		argv[0],i,stat[i*4+2],stat[i*4+3],stat[i*4+0],stat[i*4+1],stat[i*4+1]/STRESS_DROP);
      }else if(check_mu){// friction
	calc_vec_stat((mu+(i-3)*medium->nrflt),medium->nrflt,(stat+i*4));
     	fprintf(stderr,"%s: fric type %i: min: %12g max: %12g mean: %12g stddev: %12g\n",
		argv[0],i-3,stat[i*4+2],stat[i*4+3],stat[i*4+0],stat[i*4+1]);
      }else{
	stat[i*4+0]=stat[i*4+1]=stat[i*4+2]=stat[i*4+3]=0.0;
      }
    }
    fprintf(stderr,"%s: writing means and stddevs to stdout\n",argv[0]);
    printf("%g %g %g %g %g %g %g %g %g %g\n",
	   stat[0],stat[1],stat[4],stat[5],stat[8],stat[9],
	   stat[12],stat[13],stat[16],stat[17]);
    break;
  case 1:
    //
    // calculate spatial statistics
    //
    if(nr_spc_bins<0){
      fprintf(stderr,"%s: error: nr of bins shouldn't be negative: %i\n",
	      argv[0],nr_spc_bins);
      exit(-1);
    }
    fprintf(stderr,"%s: using %i bins for spatial correlation, range frac: %g\n",
	    argv[0],nr_spc_bins,range_frac);
    r=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*nr_spc_bins);
    xr=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*nr_spc_bins);
    nr=(int *)malloc(sizeof(int)*nr_spc_bins);
    if(!r || !xr || !nr)
      MEMERROR("r, xr, or nr");
    //
    // assign coordinates based on location in fault
    //
    for(j=i=0;i<5;i++,j += medium->nrflt){
      if(i<3 && check_stress){//stress
	calc_vec_stat((s+j),medium->nrflt,stat);
	fprintf(stderr,"%s: stress dir %i: min: %12g max: %12g mean: %12g stddev: %12g\n",
		argv[0],i,stat[2],stat[3],stat[0],stat[1]);
      }else if(check_mu){// friction
	calc_vec_stat((mu+(i-3)*medium->nrflt),medium->nrflt,stat);
     	fprintf(stderr,"%s: fric type %i: min: %12g max: %12g mean: %12g stddev: %12g\n",
		argv[0],i-3,stat[2],stat[3],stat[0],stat[1]);
      }else{
	stat[0]=stat[1]=stat[2]=stat[3]=0.0;
      }
      if(i<3 && check_stress){
	//
	// calc spatial stress correlation for faults, component i
	//
	calc_spatial_correlation(fault,medium->nrflt,2,i,
				 nr_spc_bins,&r,&xr,&nr,
				 range_frac,cr,&xcr,&w1,&w2,&w3);
	fprintf(stderr,"%s: stress component %i: r falls under %g at %g\n",argv[0],i,cr,xcr);
	sprintf(tmpstr,"ssr.%i.dat",i);
	out=myopen(tmpstr,"w");
	fprintf(stderr,"%s: writing spatial corr in xr nr r format to %s\n",
		argv[0],tmpstr);
	// output of distance, nr of patches in distance, and correlation
	for(k=0;k<nr_spc_bins;k++)
	  fprintf(out,"%12g %6i %12g\n",xr[k],nr[k],r[k]);
	fclose(out);
      }else if(check_mu){
	//
	// calc spatial correlation for friction coefficients
	//
	calc_spatial_correlation(fault,medium->nrflt,2,i,nr_spc_bins,&r,&xr,&nr,
				 range_frac,cr,&xcr,&w1,&w2,&w3);
	fprintf(stderr,"%s: fric component %i: r falls under %g at %g\n",argv[0],i-3,cr,xcr);
	sprintf(tmpstr,"fsr.%i.dat",i-3);
	out=myopen(tmpstr,"w");
	fprintf(stderr,"%s: writing spatial corr in xr nr r format to %s\n",
		argv[0],tmpstr);
	for(k=0;k<nr_spc_bins;k++)
	  fprintf(out,"%12g %6i %12g\n",xr[k],nr[k],r[k]);
	fclose(out);
      }
    }
    free(xr);free(nr);free(r);
    break;
  case 2:
    //
    // output of stress and friction coefficients to file
    //
    fprintf(stderr,"%s: writing pos_0/1 s_STRIKE s_DIP s_NORMAL mu_s/d to sh.dat\n",argv[0]);
    out=myopen("sh.dat","w");
    os[0] = STRIKE*medium->nrflt;
    os[1] = DIP*medium->nrflt;
    os[2] = NORMAL*medium->nrflt;
    if(check_mu && check_stress)
      for(i=0;i<medium->nrflt;i++)
	fprintf(out,"%g %g %g %g %g %g %g\n",
		fault[i].pos[0],fault[i].pos[1],
		s[i+os[0]],s[i+os[1]],s[i+os[2]],
		fault[i].mu_s,fault[i].mu_d);
    else if(check_mu)
      for(i=0;i<medium->nrflt;i++)
	fprintf(out,"%g %g %g %g %g %g %g\n",
		fault[i].pos[0],fault[i].pos[1],0.,0.,0.,
		fault[i].mu_s,fault[i].mu_d);
    else if(check_stress)
      for(j=medium->nrflt,k=j+medium->nrflt,
	    i=0;i<medium->nrflt;i++,j++,k++)
	fprintf(out,"%g %g %g %g %g %g %g\n",
		fault[i].pos[0],fault[i].pos[1],
		s[i+os[0]],s[i+os[1]],s[i+os[2]],
		STATIC_MU,(STATIC_MU-DELTA_MU));
    else 
      fprintf(stderr,"%s: no irregular inits, no output\n",argv[0]);
    fclose(out);
    break;
  default:
    fprintf(stderr,"%s: error, mode %i undefined\n",argv[0],mode);
    exit(-1);
    break;
  }

  exit(0);
}





