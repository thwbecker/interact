#include <math.h>
#include <stdio.h>

#define N 500
void main(void)
{
  double x[N][N][N];
  long int i,j,k,l,m;
  for(l=0;l<1000000;l++)
    for(m=0;m<100000;m++)
      for(i=0;i<N;i++)
	for(j=0;j<N;j++)
	  for(k=0;k<N;k++){
	    x[i][j][k] = cos((double)l) * sin((double)m);
	    x[i][j][k] *= (double)m;
	  }
  printf("ok\n");
}
