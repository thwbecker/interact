/*

  ve_pom_compare2: revisit of the scalar-image ("transplant")
  approximation for 3-D plate-over-Maxwell surface deformation,
  with the SECOND image family included (the 2 n H - d family that
  was missing from the layered antiplane machinery and from
  ve_pom_compare; see prony_kernel.c, ve_rybicki_antiplane_uz)

  construction: the elastic half-space (Okada) solution of the
  source, plus, per reflection order n >= 1 with Erlang weight
  W_n(t) = P(n, t/(2 t_M)) (equal rigidities):

    - the source shifted DOWN by 2 n H (the 2 n H + d family), and
    - the SURFACE-MIRRORED source shifted down by 2 n H (the
      2 n H - d family), i.e. the patch with dip -> 180 - dip at
      z = -z_source - 2 n H, with slip components sign-mapped by the
      z -> -z mirror

  every image source is evaluated with the half-space kernel, so the
  free surface is satisfied exactly at all times; the interface
  condition is approximated by the scalar SH reflection ladder, which
  this driver exists to measure against Kaj Johnson's propagator

  the mirror slip-sign map is determined and VERIFIED at startup by a
  full-space identity check, u_mirror(x,y,z) = M u_orig(x,y,-z),
  M = diag(1,1,-1), per slip mode; the program aborts if no sign
  reproduces the identity at rounding level

  usage:
    ve_pom_compare2 [H, 20] [t_M, 25] [t_yr, 50] [n_img, 200]
                    [slip_mode, 0] [use_both, 1] [sta_file, -]

  reads a single-patch geom.in; stations from sta_file (x y pairs),
  or the default east profile x = 1..80, y = 0 if sta_file is "-"

  output: x y ue un uz ue0 un0 uz0 (field at t and elastic field;
  transient = difference); with use_both 0 reproduces the old
  single-family construction of ve_pom_compare

  part of interact

*/
#include "interact.h"
#include "properties.h"
#include "prony_kernel.h"

/* fill the mirrored patch (about z = 0) of fault0 and the slip sign
   map; verify by the full-space mirror identity */
static void build_mirror(struct med *medium, struct flt *fault0,
			 struct flt *mirr, COMP_PRECISION *msign)
{
  COMP_PRECISION slip[3],u0[3],um[3],sm[3][3],x[3],xm[3];
  COMP_PRECISION res,scl,best;
  int mode,ip,is,iret,sgn,bsgn;
  static COMP_PRECISION tp[4][3] = {	/* generic test points */
    {13.0,7.0,-3.0},{-9.0,17.0,-8.0},{21.0,-11.0,-14.0},{6.0,3.0,-22.0}};
  *mirr = *fault0;
  mirr->x[INT_Z] = -fault0->x[INT_Z];
  mirr->dip = 180.0 - fault0->dip;
  for(mode=0;mode < 3;mode++){
    best = 1e30;bsgn = 0;
    for(sgn=-1;sgn <= 1;sgn += 2){
      res = scl = 0.0;
      for(ip=0;ip < 4;ip++){
	slip[STRIKE] = slip[DIP] = slip[NORMAL] = 0.0;
	slip[mode] = 1.0;
	x[INT_X] = tp[ip][0];x[INT_Y] = tp[ip][1];x[INT_Z] = tp[ip][2];
	eval_green(x,fault0,slip,u0,sm,&iret,GC_DISP_ONLY,FALSE,
		   TRUE,medium->elastic);	/* full space */
	slip[mode] = (COMP_PRECISION)sgn;
	xm[INT_X] = x[INT_X];xm[INT_Y] = x[INT_Y];xm[INT_Z] = -x[INT_Z];
	eval_green(xm,mirr,slip,um,sm,&iret,GC_DISP_ONLY,FALSE,
		   TRUE,medium->elastic);
	/* um(x,y,-z) should equal M u0(x,y,z), M = diag(1,1,-1) */
	for(is=0;is < 3;is++){
	  COMP_PRECISION tgt = (is == INT_Z)?(-u0[is]):(u0[is]);
	  res += (um[is]-tgt)*(um[is]-tgt);
	  scl += tgt*tgt;
	}
      }
      res = sqrt(res/((scl > 0.0)?scl:1.0));
      if(res < best){best = res;bsgn = sgn;}
    }
    if(best > 1e-8){
      fprintf(stderr,"build_mirror: mode %i: no sign reproduces the mirror identity (best %.1e)\n",
	      mode,best);
      exit(-1);
    }
    msign[mode] = (COMP_PRECISION)bsgn;
  }
  fprintf(stderr,"build_mirror: dip %g -> %g, slip sign map (strike, dip, normal): %g %g %g\n",
	  fault0->dip,mirr->dip,msign[STRIKE],msign[DIP],msign[NORMAL]);
}

int main(int argc, char **argv)
{
  struct flt *fault,img[1],mirr0[1],mimg[1];
  struct med *medium;
  COMP_PRECISION slip[3],mslip[3],msign[3],u[3],us[3],u0[3],sm[3][3],x[3];
  COMP_PRECISION plate_h,t_M,t,wfac,zsave,mzsave;
  COMP_PRECISION xsta[512],ysta[512];
  int n_img,i,j,n,iret,mode_i,use_both,nsta;
  MODE_TYPE slip_mode;
  FILE *fp;
  plate_h = 20.0;t_M = 25.0;t = 50.0;n_img = 200;
  slip_mode = STRIKE;use_both = 1;
  nsta = 0;
  if(argc > 1)sscanf(argv[1],ONE_CP_FORMAT,&plate_h);
  if(argc > 2)sscanf(argv[2],ONE_CP_FORMAT,&t_M);
  if(argc > 3)sscanf(argv[3],ONE_CP_FORMAT,&t);
  if(argc > 4)sscanf(argv[4],"%i",&n_img);
  if(argc > 5){
    sscanf(argv[5],"%i",&mode_i);
    slip_mode = (MODE_TYPE)mode_i;
  }
  if(argc > 6)sscanf(argv[6],"%i",&use_both);
  if((argc > 7) && (strcmp(argv[7],"-") != 0)){
    fp = myopen(argv[7],"r");
    while((nsta < 512) &&
	  (fscanf(fp,TWO_CP_FORMAT,&xsta[nsta],&ysta[nsta]) == 2))
      nsta++;
    fclose(fp);
  }else{
    for(i=0;i < 80;i++){
      xsta[nsta] = 1.0 + 79.0*(COMP_PRECISION)i/79.0;
      ysta[nsta] = 0.0;nsta++;
    }
  }
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("ve_pom_compare2");
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25);
  if(medium->nrflt != 1)
    fprintf(stderr,"ve_pom_compare2: WARNING: %i patches, using first\n",
	    medium->nrflt);
  fprintf(stderr,"ve_pom_compare2: H %g t_M %g t %g yr (a t = %g), n_img %i, mode %i, families %s, %i stations\n",
	  plate_h,t_M,t,t/(2.0*t_M),n_img,(int)slip_mode,
	  (use_both)?("BOTH"):("single (legacy)"),nsta);
  slip[STRIKE] = slip[DIP] = slip[NORMAL] = 0.0;
  slip[slip_mode] = 1.0;
  /* mirrored source and slip sign map, verified in full space */
  build_mirror(medium,fault,mirr0,msign);
  for(j=0;j < 3;j++)mslip[j] = 0.0;
  mslip[slip_mode] = msign[slip_mode];
  zsave = fault[0].x[INT_Z];
  mzsave = mirr0->x[INT_Z];
  x[INT_Z] = 0.0;
  for(i=0;i < nsta;i++){
    x[INT_X] = xsta[i];x[INT_Y] = ysta[i];
    for(j=0;j<3;j++){u[j] = 0.0;u0[j] = 0.0;}
    for(n=0;n <= n_img;n++){
      wfac = (n == 0)?(1.0):
	(ve_nur_mavko_weight(n,0.0,1.0/(2.0*t_M),t));
      /* plain shifted family, 2 n H + d */
      *img = fault[0];
      img->x[INT_Z] = zsave - 2.0*(COMP_PRECISION)n*plate_h;
      eval_green(x,img,slip,us,sm,&iret,GC_DISP_ONLY,FALSE,
		 FALSE,medium->elastic);
      for(j=0;j<3;j++){
	u[j] += wfac * us[j];
	if(n == 0)u0[j] = us[j];
      }
      if((n > 0) && use_both){
	/* mirrored family, 2 n H - d */
	*mimg = *mirr0;
	mimg->x[INT_Z] = mzsave - 2.0*(COMP_PRECISION)n*plate_h;
	eval_green(x,mimg,mslip,us,sm,&iret,GC_DISP_ONLY,FALSE,
		   FALSE,medium->elastic);
	for(j=0;j<3;j++)
	  u[j] += wfac * us[j];
      }
    }
    printf("%10.4f %10.4f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
	   x[INT_X],x[INT_Y],u[INT_X],u[INT_Y],u[INT_Z],
	   u0[INT_X],u0[INT_Y],u0[INT_Z]);
  }
  return 0;
}
