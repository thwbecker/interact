#include "interact.h"
#include "properties.h"
//
// reads in patch format and creates a file for random stress initializarion
// $Id: create_random_stress_file.c,v 1.15 2003/03/02 01:37:41 becker Exp $
//

int main(int argc, char **argv)
{
  
  struct flt *fault;
  struct med *medium;
  int i,j,os,os1,osnormal,osdip,iter,siter,basei,n1,maxiter=5000,maxbiter=1000,
    *sindex,maxsiter=-1,mode,biter,iseed= -1;
  my_boolean reject_pos_coulomb=TRUE,hit;
  long seed;
  float *dist;
  COMP_PRECISION hrange[3],std[3],*st,cs,stat[4],fac[3],slength,ts[3],dfac,std0[3];
  medium=(struct med *)calloc(1,sizeof(struct med));
  // defaults
  std0[0]=std0[1]=std0[2] = 1;// STD in multiples of the stress drop
  mode = 0;
  slength = 0.0;
  switch(argc){
  case 2:
    sscanf(argv[1],ONE_CP_FORMAT,(std0+STRIKE));
    break;
  case 3:
    sscanf(argv[1],ONE_CP_FORMAT,(std0+STRIKE));
    sscanf(argv[2],ONE_CP_FORMAT,(std0+DIP));
    break;
  case 4:
    sscanf(argv[1],ONE_CP_FORMAT,(std0+STRIKE));
    sscanf(argv[2],ONE_CP_FORMAT,(std0+DIP));
    sscanf(argv[3],ONE_CP_FORMAT,(std0+NORMAL));
    break;
  case 5:
    sscanf(argv[1],ONE_CP_FORMAT,(std0+STRIKE));
    sscanf(argv[2],ONE_CP_FORMAT,(std0+DIP));
    sscanf(argv[3],ONE_CP_FORMAT,(std0+NORMAL));
    sscanf(argv[4],ONE_CP_FORMAT,&slength);
    break;
  case 6:
    sscanf(argv[1],ONE_CP_FORMAT,(std0+STRIKE));
    sscanf(argv[2],ONE_CP_FORMAT,(std0+DIP));
    sscanf(argv[3],ONE_CP_FORMAT,(std0+NORMAL));
    sscanf(argv[4],ONE_CP_FORMAT,&slength);
    sscanf(argv[5],"%i",&mode);
    break;
  case 7:
    sscanf(argv[1],ONE_CP_FORMAT,(std0+STRIKE));
    sscanf(argv[2],ONE_CP_FORMAT,(std0+DIP));
    sscanf(argv[3],ONE_CP_FORMAT,(std0+NORMAL));
    sscanf(argv[4],ONE_CP_FORMAT,&slength);
    sscanf(argv[5],"%i",&mode);
    sscanf(argv[6],"%i",&maxsiter);
    break;
  case 8:
    sscanf(argv[1],ONE_CP_FORMAT,(std0+STRIKE));
    sscanf(argv[2],ONE_CP_FORMAT,(std0+DIP));
    sscanf(argv[3],ONE_CP_FORMAT,(std0+NORMAL));
    sscanf(argv[4],ONE_CP_FORMAT,&slength);
    sscanf(argv[5],"%i",&mode);
    sscanf(argv[6],"%i",&maxsiter);
    sscanf(argv[7],"%i",&iseed);
    break;
  case 9:
    sscanf(argv[1],ONE_CP_FORMAT,(std0+STRIKE));
    sscanf(argv[2],ONE_CP_FORMAT,(std0+DIP));
    sscanf(argv[3],ONE_CP_FORMAT,(std0+NORMAL));
    sscanf(argv[4],ONE_CP_FORMAT,&slength);
    sscanf(argv[5],"%i",&mode);
    sscanf(argv[6],"%i",&maxsiter);
    sscanf(argv[7],"%i",&iseed);
    sscanf(argv[8],"%i",&i);
    reject_pos_coulomb=(my_boolean)i;
    break;
  default:
    fprintf(stderr,"%s [std_s, %g] [std_d, %g] [std_n, %g] [slength, %g] [mode, %i] [maxsiter, -1] [seed, -1] [rpc, %i]\n\treads in patch format from %s and ",
	    argv[0],std0[STRIKE],std0[DIP],std0[NORMAL],slength,mode,reject_pos_coulomb,GEOMETRY_FILE);
    fprintf(stderr,"writes random stress file to stdout\n");
    fprintf(stderr,"std_s,d,n are the standard deviations of the stresses AS MULTIPLES OF THE STRESS DROP in strike/dip/normal dir\n");
    fprintf(stderr,"slength is the smoothing length in units of fault half widths\n");
    fprintf(stderr,"mode is the operational mode: 0: Gauss 1: uniform\n");
    fprintf(stderr,"maxsiter is the number of Gauss additions for slength, if set to -1, will use nrpatch*2\n");
    fprintf(stderr,"seed is the random number generator seed\n");
    fprintf(stderr,"rpc: if 1, will reject positive Coulomb stress, if zero will only reject extensional patches\n");
    exit(-1);
    break;
  }
  if(iseed > 0){
    fprintf(stderr,"%s: error: iseed: %i, should be negative\n",argv[0],iseed);
    exit(-1);
  }else{
    seed = (long) iseed;
  }
  //
  // multiply input values with STRESS_DROP 
  //
  for(i=0;i<3;i++)
    std[i] = std0[i] * STRESS_DROP;
  if(mode == 0)
    fprintf(stderr,"%s: Gauss deviates, target SDEV stress_s,d,n/stress drop(%g): %g/%g/%g\n",
	    argv[0],STRESS_DROP,std[STRIKE]/STRESS_DROP,std[DIP]/STRESS_DROP,
	    std[NORMAL]/STRESS_DROP);
  else
    fprintf(stderr,"%s: uniform deviates, target SDEV stress_s,d,n/stress drop(%g): %g/%g/%g\n",
	    argv[0],STRESS_DROP,std[STRIKE]/STRESS_DROP,std[DIP]/STRESS_DROP,
	    std[NORMAL]/STRESS_DROP);
  fprintf(stderr,"%s: reading patch format from %s, writing fsi.in to stdout\n",
	  argv[0], GEOMETRY_FILE);
  if(reject_pos_coulomb)
    fprintf(stderr,"%s: rejecting positive Coulomb stress and extensional patches\n",argv[0]);
  else
    fprintf(stderr,"%s: rejecting extensional patches only\n",argv[0]);
  //
  // read fault geometry
  read_geometry(GEOMETRY_FILE,&medium,&fault,FALSE,FALSE,FALSE,
		FALSE);
  // allocate memerory
  st=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*3*medium->nrflt);
  sindex=(int *)malloc(sizeof(int)*medium->nrflt);
  if(!st || !sindex)MEMERROR("create_random_stress: st or sindex");
  if(slength!=0.0){
    if(maxsiter == -1)
      maxsiter = medium->nrflt*2;
    fprintf(stderr,"%s: smoothing length: %g, maxsiter: %i\n",
	    argv[0],slength,maxsiter);
  }
  if(mode == 0){// Gaussian
    ;
  }else{
    // for uniform deviates
    for(i=0;i<3;i++)
      hrange[i] = std[i]/2.0;
  }
  //
  // start rgen
  //
  RAND_GEN(&seed);
  //
  // assign initial stress values
  // offsets
  //
  osnormal=medium->nrflt*NORMAL;
  osdip=medium->nrflt*DIP;
  // STRIKE should be zero
  if(STRIKE != 0){
    fprintf(stderr,"create_random_stress: internal coding problem: offset error, STRIKE != 0: %i\n",
	    STRIKE);
    exit(-1);
  }
  biter=siter=0;
  do{
    biter++;
    if(slength == 0.0){
      /*
	assign random values without spatial pattern
      */
      for(i=0;i<medium->nrflt;i++){
	iter=0;
	do{
	  iter++;
	  get_random_stress((st+i),std,hrange,mode,&seed,
			    medium->nrflt);
	  if(st[i+osnormal] > PRESSURE_DEF)// reject extensional faults
	    cs=1.0;
	  else{
	    if(reject_pos_coulomb)
	      // check coulomb stress
	      cs = hypot(st[i],st[i+osdip]) - 
		STATIC_MU * (PRESSURE_DEF - st[i+osnormal]);
	    else// only extensional normal stress will be rejected
	      cs = -1.0;
	  }
	  // reject if we have positive Coulomb stress
	}while((cs > 0.0)&&(iter<maxiter));
	if(iter >= maxiter){
	  fprintf(stderr,"%s: max iterations (%i) reached\n",
		  argv[0],iter);
	  exit(-1);
	}
      }
    }else{
      /*
	
	smooth spatially

      */
      // init stresses with zeroes
      j=medium->nrflt*3;
      for(i=0;i<j;i++)
	st[i]=0.0;
      //
      // do spatial averaging
      //
      // determine distances; dist[i*medium->nrflt+j] is the 
      // distance between patch i and j
      //
      dist=(float *)malloc(sizeof(float)*SQUARE(medium->nrflt));
      if(!dist)MEMERROR("create_random_stress: dist");
      for(i=0;i<medium->nrflt;i++){
	os = medium->nrflt*i;
	for(j=i+1;j<medium->nrflt;j++){
	  // weights are obtained from exp(-(distance/slength)^2)
	  dist[os+j]  = distance_float(fault[i].pos,
				       fault[j].pos,2);
	  dist[os+j] /= (float)slength;
	  dist[os+j]  = exp(-SQUARE(dist[os+j]));
	  dist[medium->nrflt*j + i] = dist[os+j];
	}
      }
      //
      // add a bunch of Gauss shaped curves
      //
      n1=medium->nrflt-1;
      for(siter=0;siter <maxsiter;siter++){
	// pick random start location between 0 and n-1
	basei = myrandi(n1,&seed);
	os = medium->nrflt * basei;
	// initial values values for ts
	get_random_stress(ts,std,hrange,mode,&seed,1);
	//
	// add that number to all others, scaled by 
	// the distance factor
	//
	for(j=0;j<medium->nrflt;j++){
	  dfac = (COMP_PRECISION)dist[os+j];
	  st[j]          += ts[STRIKE] * dfac;
	  st[j+osdip]    += ts[DIP]    * dfac;
	  st[j+osnormal] += ts[NORMAL] * dfac;
	}
      }
    }// end of slength part
    //
    // scale so that standard deviation works out
    // also make sure mean=0
    //
    for(os=i=0;i<3;i++,os+=medium->nrflt){
      calc_vec_stat((st+os),medium->nrflt,stat);
      if(i!=2){
	// remove mean for all but normal
	for(os1=os,j=0;j<medium->nrflt;j++,os1++)
	  st[os1] -= stat[0];
      }
      if(stat[1] != 0.0){
	// standard deviation
	fac[i]=std[i]/stat[1];
      }else
	fac[i]=1.0;
      fprintf(stderr,"%s: biter: %i siter: %i stress dir %i: stddev: theo: %g actual: %g scaling: %g\n",
	      argv[0],biter,siter,i,std[i],stat[1],fac[i]);
    }
    if(!reject_pos_coulomb)
      // check for positive Coulomb and extensional normal stress
      for(i=0,hit=FALSE;i<medium->nrflt;i++)
	if((hypot(st[i]*fac[STRIKE],
		  st[i+osdip]*fac[DIP]) - 
	    STATIC_MU * (PRESSURE_DEF - 
			 st[i+osnormal]*fac[NORMAL]))>0){
	  hit=TRUE;break;
	}
	if(st[i+osnormal]*fac[NORMAL] > PRESSURE_DEF){
	  hit=TRUE;break;
	}
    else
      // check for extensional normal stress
      for(i=0,hit=FALSE;i<medium->nrflt;i++)
	if(st[i+osnormal]*fac[NORMAL] > PRESSURE_DEF){
	  hit=TRUE;break;
	}
  }while((hit) && (biter<maxbiter));
  if(biter == maxbiter){
    fprintf(stderr,"no solution found, maxbiter exceeded\n");
    exit(-1);
  }
      
  // output 
  for(i=0;i<medium->nrflt;i++)
    printf("%20.14e %20.14e %20.14e\n",
	   st[i]*fac[STRIKE],st[i+osdip]*fac[DIP],st[i+osnormal]*fac[NORMAL]);
  exit(0);
}
/*

  determine random stresses for ts[0], ts[n], and ts[2*n]

 */
void get_random_stress(COMP_PRECISION *ts, COMP_PRECISION *std,
		       COMP_PRECISION *hrange,
		       int mode, long *seed,int n)
{
  int i,os;
  if(mode == 0)// Gauss (normally) distributed with std
    for(os=0,i=0;i<3;i++,os += n)
      ts[os] = mygauss_randnr(std[i],seed);
  else if(mode==1)// uniform
    for(os=0,i=0;i<3;i++,os += n)
      ts[os] = -hrange[i] + myrandnr(std[i],seed);
  else{
    fprintf(stderr,"get_random_stress: error: mode: %i undefined\n",
	    mode);
    exit(-1);
  }
}
