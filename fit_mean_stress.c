/* 
   given two files in 

   lon lat sxx sxy sxz syy syz szz 

   format, determine a best-fit background stress
   that will reduce the misfit between the two 
   stress fields



 */
#include "interact.h"
#include "blockinvert.h"

int main(int argc,char **argv)
{
  COMP_PRECISION *s1,*s2,*ds,dsm[6],rms1,rms2,sfit[6],drms,
    lon1,lat1,lon2,lat2,*w,*azi1,*azi2,*dazi,dcorrect[6],
    s1h,s2h,hazi;
  int n,n6,i;
  /* read in two pasted stess fields */
  n=0;
  my_vecalloc(&s1,6,"s1");  
  my_vecalloc(&s2,6,"s2");
  my_vecalloc(&w,1,"w");
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &lon1,&lat1,(s1+n*6),(s1+n*6+1),(s1+n*6+2),(s1+n*6+3),(s1+n*6+4),(s1+n*6+5),
	       &lon2,&lat2,(s2+n*6),(s2+n*6+1),(s2+n*6+2),(s2+n*6+3),(s2+n*6+4),(s2+n*6+5))==16){
    if((lon1!=lon2)||(lat1!=lat2)){
      fprintf(stderr,"%s: error, location mismatch\n",argv[0]);
      exit(-1);
    }
    w[n]=cos(DEG2RAD*lat1);	/* weights */
    n++;
    my_vecrealloc(&s1,6*(n+1),"s1");  
    my_vecrealloc(&s2,6*(n+1),"s2");
    my_vecrealloc(&w,(n+1),"w");
  }

  n6 = n * 6;
  /* allocate memory */
  my_vecalloc(&ds,n6,"ds");
  my_vecalloc(&azi1,n6,"ds");
  my_vecalloc(&azi2,n6,"ds");
  my_vecalloc(&dazi,n6,"ds");
  /* calculate deviation vector */
  c_eq_a_minus_b(ds,s1,s2,n6);
  rms1=rms(s1,n6);rms2=rms(s2,n6);drms=rms(ds,n6);
  fprintf(stderr,"%s: read %i stress observations, rms: s1: %g s2: %g ds: %g\n",
	  argv[0],n,rms1,rms2,drms);
  /* calculate weighted mean deviation */
  for(i=0;i<6;i++)
    dsm[i] = wmean((ds+i),6,n,w);
  fprintf(stderr,"%s: original: mean devation tensor: %g %g %g %g %g %g\n",
	  argv[0],dsm[0],dsm[1],dsm[2],dsm[3],dsm[4],dsm[5]);
  /* calculate horizontal parts */
  calc_horizontal_stress_vec(s1,azi1,n);
  calc_horizontal_stress_vec(s2,azi2,n);
  /* deviation, allow for sign */
  calc_dir_diff_vec(azi1,azi2,dazi,n,TRUE);
  fprintf(stderr,"%s: original: mean h_azi dev: %g mean(abs(h_azi)): %g\n",
	  argv[0],wmean(dazi,1,n,w),wmean_abs(dazi,1,n,w));
  
  /* start by correcting for mean deviation */
  for(i=0;i<6;i++)
    dcorrect[i] = dsm[i];
  
  /* 
     
     begin evaluation procedure
  
  */
  eval_stress_correction(s1,s2,w,azi1,n,dcorrect,azi2,dazi,ds,&drms,dsm);
  fprintf(stderr,"%s: rms ds: %g mean devation tensor: %g %g %g %g %g %g\n",
	  argv[0],drms,dsm[0],dsm[1],dsm[2],dsm[3],dsm[4],dsm[5]);
  fprintf(stderr,"%s: test: mean h_azi dev: %g mean(abs(h_azi)): %g\n",
	  argv[0],wmean(dazi,1,n,w),wmean_abs(dazi,1,n,w));
  /* do grid search */
  for(s1h=-1.0;s1h<=1.0;s1h+=0.05)
    for(s2h=-1.0;s2h<=1.0;s2h+=0.05){
      for(hazi=0.0;hazi<=180;hazi+=0.05){
	/* get cartesian matrix that corresponds to the hor system */
	cart_mat_from_horsym(s1h,s2h,hazi,dcorrect);
	/* get new misfit */
	eval_stress_correction(s1,s2,w,azi1,n,dcorrect,azi2,dazi,ds,&drms,dsm);
	printf("%g %g %g %g %g %g\n",
	       s1h,s2h,hazi,drms,
	     wmean(dazi,1,n,w),wmean_abs(dazi,1,n,w));
      }
    }
  sfit[0]=0.0;
  return 0;
}

void eval_stress_correction(COMP_PRECISION *s1,COMP_PRECISION *s2,
			    COMP_PRECISION *w,COMP_PRECISION *azi1, int n,
			    COMP_PRECISION *dcorrect,
			    COMP_PRECISION *azi2, COMP_PRECISION *dazi,
			    COMP_PRECISION *ds,COMP_PRECISION *drms,
			    COMP_PRECISION *dsm)
{
  int i,j,os,n6;
  COMP_PRECISION *stest;
  n6 = n * 6;
  my_vecalloc(&stest,n6,"ds");
  /* adjust second set of stresses */
  for(i=os=0;i<n;i++,os+=6){
    for(j=0;j < 6;j++)
      stest[os+j] = s2[os+j] + dcorrect[j];
  }
  /* deviations */
  c_eq_a_minus_b(ds,s1,stest,n6);
  *drms = rms(ds,n6);
  /* calculate weighted mean deviation */
  for(i=0;i<6;i++)
    dsm[i] = wmean((ds+i),6,n,w);
  calc_horizontal_stress_vec(stest,azi2,n);
  calc_dir_diff_vec(azi1,azi2,dazi,n,TRUE);
  free(stest);
}
