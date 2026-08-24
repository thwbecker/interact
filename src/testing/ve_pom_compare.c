/*

  step 7b comparison driver: approximate plate-over-Maxwell surface
  displacements for a 3-D dislocation source by the antiplane-style
  image construction, for quantitative comparison against Kaj M.
  Johnson's Plate_over_Maxwell propagator solution (pom maxwell)

  the approximation: the elastic half-space (Okada) solution of the
  source, plus images shifted DOWN by 2 n H, weighted by the
  time-domain step response of Gamma(s)^n; for equal plate and
  substrate rigidities Gamma(s) = a/(s + a) with a = 1/(2 t_M) and
  the weight is the Erlang ramp W_n(t) = P(n, a t)
  (ve_nur_mavko_weight with Gamma_0 = 0). this scalar reflection
  coefficient is exact for antiplane (SH) motion only; for general
  3-D sources the interface couples P and SV and the construction is
  an approximation whose quality this driver exists to measure

  reads a single-patch geom.in (any Okada rectangle); observation
  points on the free surface along the east (x) axis through the
  patch center

  usage: ve_pom_compare [H, 20] [t_M, 25] [t_yr, 50] [n_img, 60] [slip_mode, 0]

  output to stdout: x u_east u_north u_up (units of the slip, both
  the full field at time t and, as extra columns, the elastic t = 0
  field for transient computation): x ue un uz ue0 un0 uz0

  part of interact

*/
#include "interact.h"
#include "properties.h"
#include "prony_kernel.h"

int main(int argc, char **argv)
{
  struct flt *fault,img[1];
  struct med *medium;
  COMP_PRECISION slip[3],u[3],us[3],u0[3],sm[3][3],x[3];
  COMP_PRECISION plate_h,t_M,t,wfac,zsave;
  int n_img,i,j,n,iret,mode_i;
  MODE_TYPE slip_mode;
  plate_h = 20.0;t_M = 25.0;t = 50.0;n_img = 60;
  slip_mode = STRIKE;
  if(argc > 1)sscanf(argv[1],ONE_CP_FORMAT,&plate_h);
  if(argc > 2)sscanf(argv[2],ONE_CP_FORMAT,&t_M);
  if(argc > 3)sscanf(argv[3],ONE_CP_FORMAT,&t);
  if(argc > 4)sscanf(argv[4],"%i",&n_img);
  if(argc > 5){
    sscanf(argv[5],"%i",&mode_i);
    slip_mode = (MODE_TYPE)mode_i;
  }
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("ve_pom_compare");
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25);
  if(medium->nrflt != 1)
    fprintf(stderr,"ve_pom_compare: WARNING: %i patches, using all\n",
	    medium->nrflt);
  fprintf(stderr,"ve_pom_compare: H %g t_M %g t %g yr (a t = %g), n_img %i, mode %i\n",
	  plate_h,t_M,t,t/(2.0*t_M),n_img,(int)slip_mode);
  slip[STRIKE] = slip[DIP] = slip[NORMAL] = 0.0;
  slip[slip_mode] = 1.0;
  x[INT_Y] = 0.0;x[INT_Z] = 0.0;
  for(i=0;i < 80;i++){
    x[INT_X] = 1.0 + 79.0*(COMP_PRECISION)i/79.0;
    for(j=0;j<3;j++){u[j] = 0.0;u0[j] = 0.0;}
    for(n=0;n <= n_img;n++){
      /* Erlang weight of image family n at time t; unity for the
	 source itself */
      wfac = (n == 0)?(1.0):
	(ve_nur_mavko_weight(n,0.0,1.0/(2.0*t_M),t));
      if((n > 0)&&(wfac < 1e-14)&&(t <= 0.0))
	break;
      *img = fault[0];
      zsave = fault[0].x[INT_Z];
      img->x[INT_Z] = zsave - 2.0*(COMP_PRECISION)n*plate_h;
      eval_green(x,img,slip,us,sm,&iret,GC_DISP_ONLY,FALSE,
		 medium->full_space,medium->elastic);
      for(j=0;j<3;j++){
	u[j] += wfac * us[j];
	if(n == 0)u0[j] = us[j];	/* elastic (t = 0) field */
      }
    }
    printf("%10.4f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
	   x[INT_X],u[INT_X],u[INT_Y],u[INT_Z],
	   u0[INT_X],u0[INT_Y],u0[INT_Z]);
  }
  return 0;
}
