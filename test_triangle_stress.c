#include "interact.h"
#include "properties.h"

int main(int argc,char **argv)
{
  double loc[3],Ss,Ds,Ts,stress1[6],stress2[6],strain[6];
  double mu = SHEAR_MODULUS,lambda = LAMBDA_CONST;
  double p1[3]={0,0,-1},p2[3]={1,0,-1},p3[3]={0,0,-0.01};
  int i;
  //Ss = 1;Ds=1;Ts=1;
  Ss = 1;Ds=0;Ts=0;  
  
  for(loc[1]=-0.5;loc[1]<=0.5;loc[1]+=0.1){
    for(loc[0]=-0.5;loc[0]<=1.5;loc[0]+=0.1){
      for(loc[2]=-1.5;loc[2] < 0;loc[2]+=0.1){
	
	
	tdstresshs_(loc,p1,p2,p3,&Ss,&Ds,&Ts,stress1,strain);
	hbi_tdstresshs_(loc,(loc+1),(loc+2),p1,p2,p3,&Ss,&Ds,&Ts,&mu,&lambda,
			(stress2),(stress2+1),(stress2+2),(stress2+3),(stress2+4),(stress2+5));
	for(i=0;i<6;i++){
	  if((finite(stress1[i])&&(!finite(stress2[i])))||
	     (!finite(stress1[i])&&(finite(stress2[i])))||
	     (fabs(stress1[i]-stress2[i])>1e-7)){
	    fprintf(stdout,"min) %11g %11g %11g ",loc[0],loc[1],loc[2]);
	    fprintf(stdout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e\n",stress1[0],stress1[1],stress1[2],stress1[3],stress1[4],stress1[5]);
	    
	    fprintf(stdout,"hbi) %11g %11g %11g ",loc[0],loc[1],loc[2]);
	    fprintf(stdout,"%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e\n",stress2[0],stress2[1],stress2[2],stress2[3],stress2[4],stress2[5]);
	    i=6;
	  }
	}
      }
    }
  }
  return 0;
}
