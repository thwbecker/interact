/* standalone shim: same double-precision sincos call as interact's
   my_sincos_ftn (src/mysincos.c) so the kernels link without the full tree */
#define _GNU_SOURCE
#include <math.h>
void sincos(double, double *, double *);
void my_sincos_ftn_(double *sin_val, double *cos_val, double *alpha)
{
  sincos(*alpha, sin_val, cos_val);
}
