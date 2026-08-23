/*

  fault-based visco-elastic stress relaxation using the
  effective-modulus (implicit Euler) scheme by Dave May (Aug 2026),
  section 1

  Maxwell material in shear with constant elastic bulk modulus,

     sdot/(2G) + s/(2 eta) = edot,

  discretized with implicit Euler, leads to an elastic problem with an
  effective shear modulus

     G'(Dt) = G alpha(Dt),  alpha(Dt) = eta/(eta + Dt G) = t_M/(t_M + Dt),

  where t_M = eta/G is the Maxwell time.

  the stress update via Green's functions F(.) is (eq I on the note)

     sig^{k+1} = F(G'(Dt),K,S^{k+1}) - dev(F(G'(Dt),K,S^k)) + alpha(Dt) s^k

  with s^k = dev(sig^k) the deviatoric stress of the previous step; K
  is held fixed (viscous effects in shear only), so the effective
  Poisson ratio for the Green's function follows from K and G'(Dt) as

     nu'(Dt) = (3K - 2G')/(6K + 2G'),  with  K = 2G(1+nu)/(3(1-2nu))

  multi-source version: the initial slip distribution is taken from
  fault[].u[]; any number of patches may slip, all held constant after
  t = 0

  since the slip is constant and Dt IS ASSUMED FIXED, the two Green's function
  fields are each evaluated only once:

     sig_el[i]    = sum_j F(G,      K, u_j)   the elastic field (t = 0 state)
     sigma_eff[i] = sum_j F(G'(Dt), K, u_j)   the effective field

  and every time step reduces to (no kernel calls)

     sigma_i <-- vol(sigma_eff[i]) + alpha(Dt) dev(sigma_i),

  which is the constant-slip form of the update above
  
  for slip held constant this recursion has the closed form
  sigma_k = vol(sigma_eff) + alpha^k dev(sig_el), which serves as a
  check

  Note that the volumetric part is evaluated from the current-step
  effective solution only, so it depends on Dt but does not evolve in
  time

*/
#include "interact.h"
#include "properties.h"

void relax_stress_field(struct med *, struct flt *,
			COMP_PRECISION, COMP_PRECISION,
			COMP_PRECISION (*)[3][3]);
void relax_stress_step(struct med *,
			 COMP_PRECISION, COMP_PRECISION,
			 COMP_PRECISION (*)[3][3],
			 COMP_PRECISION (*)[3][3]);

int main(int argc, char **argv)
{
  struct flt *fault;
  struct med *medium;
  COMP_PRECISION (*sigma)[3][3],(*sigma_eff)[3][3],xloc;
  COMP_PRECISION t_M,Dt,time,t_max,max_slip,slip_amp,sval,lfrac,rfrac;
  MODE_TYPE slip_dir_mode,slip_dist_mode;
  int i,k,nslip,ileft,iright,irange;
  /* 
     defaults 
  */
  slip_dir_mode = STRIKE;		/* used for slip and resolved stress output */
  slip_dist_mode = 1;			/* 0: constant 1: elliptical */
  lfrac = rfrac = 0.5;		/* only center node slips by default */

  t_M = 1.0;			/* Maxwell time eta/G */
  Dt = 0.1;			/* time step, in units of t_M */
  
  t_max = 50*t_M;		/* stop time */
  
  if(argc > 1){			/* read Dt  */
    sscanf(argv[1],ONE_CP_FORMAT,&Dt);
  }
  if(argc > 2){		/* slip mode */
    sscanf(argv[2],"%i",&i);
    slip_dir_mode = (MODE_TYPE)i;
  }
  if(argc > 3){		/* left patch fraction */
    sscanf(argv[3],ONE_CP_FORMAT,&lfrac);
  }
  if(argc > 4){		/* left patch fraction */
    sscanf(argv[4],ONE_CP_FORMAT,&rfrac);
  }
  if((lfrac < 0)||(lfrac>=1)||(rfrac<0)||(rfrac>=1)){
    fprintf(stderr,"%s: lfrac %g rfrac %g out of bounds\n",argv[0],lfrac,rfrac);
    exit(-1);
  }
  /* 
     read geometry; read_geometry also sets the default elastic
     parameters in medium->elastic 
  */
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("relax_fault");
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  /* reset the elastic defaults to unity shear modulus and 1/4 Poisson */
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25); 
  /* 
     
     range of asignments
     
  */
  ileft=(int)((COMP_PRECISION)medium->nrflt*lfrac);
  iright=(int)((COMP_PRECISION)medium->nrflt*rfrac);
  if(ileft == medium->nrflt)ileft--;
  if(iright == medium->nrflt)iright--;
  if(ileft>iright){
    i = ileft;ileft = iright;iright = i;
  }
  irange = iright-ileft;
  if(irange == 0)irange++;
  /* initial slip */
  nslip = 0;max_slip = 0;
  for(i=0;i < medium->nrflt;i++){
    if((i >= ileft)&&(i <= iright)){
      if(slip_dist_mode == 0)
	sval = 1.0;
      else{			/* elliptical */
	xloc = ((COMP_PRECISION)i - (COMP_PRECISION)(ileft+iright)/2.0)/
	  ((COMP_PRECISION)irange/2.0);
	sval = 1.0 - xloc*xloc;
	sval = (sval > 0.0)?(sqrt(sval)):(0.0);
      }
      //fprintf(stderr,"%i %g\n",i,sval);
    }else{
      sval = 0.0;
    }
    get_right_slip(fault[i].u,slip_dir_mode,sval,(fault+i)); /* call for all to make sure others are zero */
    slip_amp = norm_3d(fault[i].u);
    if(slip_amp > EPS_COMP_PREC){
      nslip++;
      if(slip_amp > max_slip)
	max_slip = slip_amp;
    }
  }
  /* 
     end of slip assignment 
  */
  fprintf(stderr,"%s: %i patches, %i slipping (%i to %i), max_slip: %g, t_M: %g Dt: %g t_max: %g (G: %g nu: %g)\n",
	  argv[0],medium->nrflt,nslip,ileft,iright,max_slip,t_M,Dt,t_max,
	  medium->elastic.shear,medium->elastic.poisson);
  /* 
     stress tensor state and the effective field, per patch 
  */
  sigma=  (COMP_PRECISION (*)[3][3])calloc(medium->nrflt,9*sizeof(COMP_PRECISION));
  sigma_eff=(COMP_PRECISION (*)[3][3])calloc(medium->nrflt,9*sizeof(COMP_PRECISION));
  if((!sigma)||(!sigma_eff))MEMERROR("relax_fault");
  /* 
     the two kernel evaluation loops: elastic field into sigma (the t
     = 0 state), effective field into sigma_eff (reused every step)
  */
  relax_stress_field(medium,fault,0.0,t_M,sigma);
  relax_stress_field(medium,fault,Dt, t_M,sigma_eff);
  /* 
     time loop; step 0 is the elastic state, each further step is the
     O(nrflt) recursion 
  */
  time = 0.0;
  k = 0;
  while(time < t_max){
    if(k){			/* k = 0 is the elastic state */
      relax_stress_step(medium,Dt,t_M,sigma_eff,sigma);
      time += Dt;
    }
    printf("%11g %5i\t",time,k);
    for(i=0;i < medium->nrflt;i++){
      /* resolve stress in the slip_mode component */
      printf("%12.5e ",resolve_stress_on_fault(sigma[i],(fault+i),slip_dir_mode));
    }
    printf("\n");
    k++;
  }
  fprintf(stderr,"%s: computed %i steps\n",argv[0],k-1);
  free(sigma);free(sigma_eff);free(fault);free(medium);
  return 0;
}
/* 

   accumulate the stress tensor at all patch centroids from all
   patches with non-zero slip in fault[].u, using the effective
   elastic parameters for time step Dt; Dt = 0 gives the elastic
   field; sfield is overwritten

   sources are the outer loop so that only slipping patches are
   computed and, per source, dip and alpha are constant along the
   inner receiver loop

*/
void relax_stress_field(struct med *medium, struct flt *fault,
			COMP_PRECISION Dt, COMP_PRECISION t_M,
			COMP_PRECISION (*sfield)[3][3])
{
  struct el_par elastic_par;
  COMP_PRECISION g0,nu0,bulk0,alphat,gp,nueff;
  COMP_PRECISION u[3],sm[3][3];
  int i,j,jj,kk,iret;
  /* 
     reference elastic and derived effective parameters; fixed K 
  */
  g0  = medium->elastic.shear;
  nu0 = medium->elastic.poisson;
  bulk0 = bulk_mod_from_G_nu(g0,nu0);
  /* effective moduli */
  alphat = t_M/(t_M + Dt);	/* = eta/(eta + Dt G) */
  gp = g0 * alphat;		/* G'(Dt) */
  nueff = nu_from_G_bulk(gp,bulk0); /* nu'(Dt) at fixed K */
  /* set the elastic parameters in the format needed for the Green's functions */
  calc_medium_elastic_parameters(&elastic_par,gp,nueff);
  /* 
     zero the field, then accumulate source by source 
  */
  for(i=0;i < medium->nrflt;i++)
    for(jj=0;jj < 3;jj++)
      for(kk=0;kk < 3;kk++)
	sfield[i][jj][kk] = 0.0;
  for(j=0;j < medium->nrflt;j++){ /* sources */
    if(norm_3d(fault[j].u)<EPS_COMP_PREC)
      continue;			/* only slipping patches contribute */
    for(i=0;i < medium->nrflt;i++){ /* receivers */
      eval_green_at_receiver(fault,i,j,fault[j].u,u,sm,&iret,
			     GC_STRESS_ONLY,TRUE,medium->full_space,elastic_par);
      if(iret)
	fprintf(stderr,"relax_stress_field: WARNING: singular at patch %i from source %i\n",i,j);
      for(jj=0;jj < 3;jj++)
	for(kk=0;kk < 3;kk++)
	  sfield[i][jj][kk] += sm[jj][kk];
    }
  }
}
/* 

   one implicit Euler step of the effective-modulus visco-elastic
   stress update, in place on the stress state sigma[nrflt][3][3],
   for slip held constant since the last step:

   sigma <-- vol(sigma_eff) + alpha(Dt) dev(sigma)

   sigma_eff is the precomputed effective field of relax_stress_field
   at this Dt; no Green's function evaluations are needed here; if
   the slip distribution or Dt change, rebuild sigma_eff first

*/
void relax_stress_step(struct med *medium,
			 COMP_PRECISION Dt, COMP_PRECISION t_M,
			 COMP_PRECISION (*sigma_eff)[3][3],
			 COMP_PRECISION (*sigma)[3][3])
{
  COMP_PRECISION alphat,trace_eff,trace_sigma;
  int i,j,k;
  alphat = t_M/(t_M + Dt);	/* = eta/(eta + Dt G) */
  for(i=0;i < medium->nrflt;i++){
    trace_eff =   tracemat3x3(sigma_eff[i]);
    trace_sigma = tracemat3x3(sigma[i]);
    for(j=0;j < 3;j++)
      for(k=0;k < 3;k++){
	/* decayed deviator of the previous state plus the volumetric
	   part of the effective field */
	sigma[i][j][k] = alphat * sigma[i][j][k];
	if(j == k)
	  sigma[i][j][k] += (trace_eff - alphat * trace_sigma)/3.0;
      }
  }
}
