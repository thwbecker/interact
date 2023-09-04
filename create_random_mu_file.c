#include "interact.h"
#include "properties.h"
//
// reads in patch format and creates a file for random friction initializarion
// (can also produce regular friction values)
//
// $Id: create_random_mu_file.c,v 1.10 2003/03/02 01:37:41 becker Exp $
//

int main(int argc, char **argv)
{
  struct flt *fault;
  struct med *medium;
  int i,j,os,iter[2],miter[2]={100000,1000},hit,iseed=-1;
  long seed;
  COMP_PRECISION mu0[2]={STATIC_MU,(STATIC_MU-DELTA_MU)},
    std[2]={0.0,0.0},*mu,stat[4],fac[2],mutmp[2];
  medium=(struct med *)calloc(1,sizeof(struct med));
  switch(argc){
  case 1:{
    ;
    break;
  }
  case 2:
    sscanf(argv[1],ONE_CP_FORMAT,std);
    break;
  case 3:
    sscanf(argv[1],ONE_CP_FORMAT,std);
    sscanf(argv[2],ONE_CP_FORMAT,(std+1));
    break;
  case 4:
    sscanf(argv[1],ONE_CP_FORMAT,std);
    sscanf(argv[2],ONE_CP_FORMAT,(std+1));
    sscanf(argv[3],ONE_CP_FORMAT,mu0);
    break;
  case 5:
    sscanf(argv[1],ONE_CP_FORMAT,std);
    sscanf(argv[2],ONE_CP_FORMAT,(std+1));
    sscanf(argv[3],ONE_CP_FORMAT,mu0);
    sscanf(argv[4],ONE_CP_FORMAT,(mu0+1));
    break;
  case 6:
    sscanf(argv[1],ONE_CP_FORMAT,std);
    sscanf(argv[2],ONE_CP_FORMAT,(std+1));
    sscanf(argv[3],ONE_CP_FORMAT,mu0);
    sscanf(argv[4],ONE_CP_FORMAT,(mu0+1));
    sscanf(argv[5],"%i",&iseed);
    break;
  default:
    fprintf(stderr,"%s [std mu_s, %g] [std mu_d, %g] [mu_s^0, %g] [mu_d^0, %g] [iseed, %i]\n\treads in patch format from %s and ",
	    argv[0],std[0],std[1],mu0[0],mu0[1],iseed,GEOMETRY_FILE);
    fprintf(stderr,"writes random friction value file to stdout\n");
    fprintf(stderr,"std is the standard deviation for the friction values (static and dyn.)\n");
    fprintf(stderr,"mu_i^0 are the base values\n");
    fprintf(stderr,"iseed is the random number generator seed\n");
    exit(-1);
    break;
  }
  if(iseed > 0){
    fprintf(stderr,"iseed: %i, should be negative\n",iseed);
    exit(-1);
  }else{
    seed = (long) iseed;
  }
  fprintf(stderr,"%s: reading patch format from %s, writing fp.in file to stdout\n",
	  argv[0],GEOMETRY_FILE);
  read_geometry(GEOMETRY_FILE,&medium,&fault,FALSE,FALSE,FALSE,
		FALSE);
  fprintf(stderr,"%s: starting values mu_s,d: %g %g, STD for mu_s,d: %g, %g\n",
	  argv[0],mu0[0],mu0[1],std[0],std[1]);
  mu=(COMP_PRECISION *)calloc(2*medium->nrflt,sizeof(COMP_PRECISION));
  if(!mu)MEMERROR("");
  //
  RAND_GEN(&seed);
  iter[0]=0;
  do{
    iter[0]++;
    for(j=medium->nrflt,i=0;i<medium->nrflt;i++,j++){
      iter[1] = 0;
      do{
	iter[1]++;
	mu[i] = mygauss_randnr(std[0],&seed);
	mu[j] = mygauss_randnr(std[1],&seed);
	mutmp[0] = mu0[0] + mu[i];
	mutmp[1] = mu0[1] + mu[j];
      }while(((mutmp[0] < 0) || (mutmp[1] < 0) || (mutmp[0] <= mutmp[1]))
	     &&(iter[1]<miter[1]));
      if(iter[1] == miter[1]){
	fprintf(stderr,"too many type 2 iterations\n");
	exit(-1);
      }
    }
    for(os=i=0;i<2;i++,os+=medium->nrflt){
      calc_vec_stat((mu+os),medium->nrflt,stat);
      if(stat[1] != 0.0){
	fac[i]=std[i]/stat[1];
      }else
	fac[i]=1.0;
      // remove mean
      if(stat[0] != 0.0){
	//fprintf(stderr,"correcting mean of %g\n",stat[0]);
	for(j=0;j<medium->nrflt;j++)
	  mu[os+j] -= stat[0];
      }
      //fprintf(stderr,"iter: %i fric coeff %i, stddev: theo: %g actual: %g scaling: %g\n",
      //      iter[0],i,std[i],stat[1],fac[i]);
    }
    // check for range
    for(hit=i=0,j=medium->nrflt;i<medium->nrflt;i++,j++){
      mutmp[0] = mu0[0]+mu[i]*fac[0]; 
      mutmp[1] = mu0[1]+mu[j]*fac[1];
      if(mutmp[1] > mutmp[0] || mutmp[0] < 0 || 
	 mutmp[1] < 0){
	hit=1;break;
      }
    }
  }while(hit && (iter[0] < miter[0]));
  if(iter[0] == miter[0]){
     fprintf(stderr,"too many type 1 iterations\n");
     exit(-1);
  }
  // output
  for(j=medium->nrflt,i=0;i<medium->nrflt;i++,j++){
    mutmp[0] = mu0[0]+mu[i]*fac[0]; 
    mutmp[1] = mu0[1]+mu[j]*fac[1];
    printf("%20.14e %20.14e\n",mutmp[0],mutmp[1]);
  }

  exit(0);
}

