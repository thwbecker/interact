/*

  step 7a validation driver: layered antiplane plate over Maxwell

  checks, for a surface observation point and a shallow strike-slip
  source read from geom.in:

  (1) s -> inf anchor: the layered s-domain sample at large s equals
      the homogeneous elastic kernel (beta -> 0)
  (2) s -> 0 anchor: the sample at tiny s equals direct image
      summation with unit weights (free-bottom plate)
  (3) Prony fit quality versus term count np = 2..6: fit the layered
      kernel samples, report the held-out residual, THE step 7a
      measurement (term count for target accuracy)

  usage: ve_layered_check [plate_h, 2] [t_M, 1] [n_img, 60]
  part of interact

*/
#include "interact.h"
#include "properties.h"
#include "prony_kernel.h"

int main(int argc, char **argv)
{
  struct flt *fault;struct med *medium;struct prony_spec spec;
  COMP_PRECISION slip[3],u[3],ue[3],sm[3][3],x[3],plate_h,t_M;
  COMP_PRECISION s,beta,a,smp[VE_MAX_NS],amp[VE_MAX_NP+1],dev,res,scl,model;
  int n_img,np,i,k,it,iret,isrc;
  plate_h = 2.0;t_M = 1.0;n_img = 60;
  if(argc > 1)sscanf(argv[1],ONE_CP_FORMAT,&plate_h);
  if(argc > 2)sscanf(argv[2],ONE_CP_FORMAT,&t_M);
  if(argc > 3)sscanf(argv[3],"%i",&n_img);
  medium=(struct med *)calloc(1,sizeof(struct med));
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25);
  a = 1.0/(2.0*t_M);
  slip[STRIKE]=1.0;slip[DIP]=slip[NORMAL]=0.0;
  isrc = medium->nrflt/2;
  x[INT_X] = fault[isrc].x[INT_X] + 1.5;x[INT_Y]=0.0;x[INT_Z]=0.0;
  /* (1) large-s anchor vs homogeneous elastic */
  ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
			      a/(1e8+a),n_img,plate_h,GC_DISP_ONLY,u,sm);
  eval_green(x,(fault+isrc),slip,ue,sm,&iret,GC_DISP_ONLY,FALSE,FALSE,
	     medium->elastic);
  dev = fabs(u[INT_Z]-ue[INT_Z])/fabs(ue[INT_Z]);
  fprintf(stderr,"anchor s->inf (elastic): rel dev %8.1e\n",dev);
  /* (2) small-s anchor vs unit-weight images */
  ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
			      a/(1e-9+a),n_img,plate_h,GC_DISP_ONLY,u,sm);
  ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
			      1.0,n_img,plate_h,GC_DISP_ONLY,ue,sm);
  dev = fabs(u[INT_Z]-ue[INT_Z])/fabs(ue[INT_Z]);
  fprintf(stderr,"anchor s->0 (free-bottom plate): rel dev %8.1e\n",dev);
  /* (3) Prony fit quality vs np on the surface displacement kernel */
  for(np=2;np <= VE_MAX_NP;np++){
    ve_spec_layered_antiplane(&spec,medium->elastic.shear,
			      medium->elastic.poisson,t_M,np);
    scl = 0.0;
    for(k=0;k < spec.ns;k++){
      s = spec.sk[k];
      beta = a/(s+a);
      ve_layered_antiplane_sample(medium,fault,isrc,slip,x,
				  beta,n_img,plate_h,GC_DISP_ONLY,u,sm);
      smp[k] = u[INT_Z];
      if(fabs(smp[k]) > scl)scl = fabs(smp[k]);
    }
    for(it=0;it < spec.nterm;it++){
      amp[it] = 0.0;
      for(k=0;k < spec.nterm;k++)
	amp[it] += spec.W[it][k]*smp[k];
    }
    res = 0.0;
    for(k=spec.nterm;k < spec.ns;k++){
      model = 0.0;
      for(it=0;it < spec.nterm;it++)
	model += amp[it]*ve_basis(&spec,it,spec.sk[k]);
      dev = fabs(model - smp[k])/scl;
      if(dev > res)res = dev;
    }
    fprintf(stderr,"np %i: held-out fit residual %8.1e\n",np,res);
  }
  return 0;
}
