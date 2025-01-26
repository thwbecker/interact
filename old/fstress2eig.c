#include "interact.h"
#include "blockinvert.h"
/*


read in full stress tensor and print projected eigensystem components


*/
int main(void)
{
  COMP_PRECISION s[9],x[2];
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf",
	       x,(x+1),s,(s+1),(s+2),(s+3),(s+4),(s+5))==8){
    expand_stress_matrix6to9(s);
    fprintf(stdout,"%12g %12g ",x[0],x[1]);
    print_projected_stress(s,stdout); /* hor part of matrix */
  }
  exit(0);
}
