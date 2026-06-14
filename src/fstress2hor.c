#include "interact.h"
#include "blockinvert.h"
/*


read in full stress tensor in [6] symmetric tensor storage

xx xy xz yy yz zz

and print eigensytem of horizontal components xx xy yy

in format e1 e2 azi(e1) where e1 > e2 and azi is the azimuth 
(CW from north) to e1




*/
int main(void)
{
  COMP_PRECISION s[9],x[2];
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf",
	       x,(x+1),s,(s+1),(s+2),(s+3),(s+4),(s+5))==8){
    expand_stress_matrix6to9(s);
    fprintf(stdout,"%12g %12g ",x[0],x[1]);
    /* hor part of matrix, ie. x y e1 e2 azi(e1) */
    print_horizontal_stress(s,stdout); 
  }
  exit(0);
}
