/*

  fault-based visco-elastic stress relaxation using the
  effective-modulus (implicit Euler) scheme by Dave May (Aug 2026),
  section 1

  a Maxwell material in shear with constant (elastic) bulk modulus,

     sdot/(2G) + s/(2 eta) = edot,

  discretized with implicit Euler, leads to an elastic problem with an
  effective shear modulus

     G'(Dt) = G alpha(Dt),  alpha(Dt) = eta/(eta + Dt G) = t_M/(t_M + Dt),

  where t_M = eta/G is the Maxwell time.

  the stress update via Green's functions F(.) reads (eq I on the note)

     sig^{k+1} = F(G'(Dt),K,S^{k+1}) - dev(F(G'(Dt),K,S^k)) + alpha(Dt) s^k

  with s^k = dev(sig^k) the deviatoric stress of the previous step; K
  is held fixed (viscous effects in shear only), so the effective
  Poisson ratio for the Green's function follows from K and G'(Dt) as

     nu'(Dt) = (3K - 2G')/(6K + 2G'),  with  K = 2G(1+nu)/(3(1-2nu))

  (in the K(E), G(E) helper relations of the source text the second
  symbol is the Poisson ratio nu)

  we read a fault geometry in interact patch format, apply a constant
  slip on one patch at t = 0, and track the full stress tensor at all
  patch centroids through time; for slip held constant, the deviatoric
  part then decays per step by alpha(Dt), the backward Euler
  approximation of exp(-Dt/t_M), and Dt -> 0 recovers the elastic
  solution

  a limitation inherited from the Green's function shortcut of the
  note: the volumetric part is evaluated from the current-step
  effective solution only, so it depends on Dt but does not evolve in
  time; the deviatoric (fault shear) relaxation is the meaningful part

*/
#include "interact.h"
#include "properties.h"

void relax_stress_update(struct med *, struct flt *, int,
			 COMP_PRECISION *, COMP_PRECISION *,
			 COMP_PRECISION, COMP_PRECISION,
			 COMP_PRECISION (*)[3][3]);

int main(int argc, char **argv)
{
  struct flt *fault;
  struct med *medium;
  COMP_PRECISION (*sigma)[3][3];
  COMP_PRECISION slip_new[3],slip_old[3];
  COMP_PRECISION t_M,Dt,slip,time,t_max;
  MODE_TYPE slip_mode;
  int islip,i,k;
  /* 
     defaults 
  */
  islip = 0;			/* patch that slips, 0...nflt-1*/
  slip_mode = STRIKE;		/* strike slip */
  slip = 1.0;			/* slip amplitude */
  t_M = 1.0;			/* Maxwell time eta/G */
  Dt = 0.1;			/* time step, in units of t_M */

  t_max = 50*t_M;			/* stop time */

  if(argc>1)			/* read Dt  */
    sscanf(argv[1],ONE_CP_FORMAT,&Dt);
  if(argc>2)			/* read in slipping patch */
    sscanf(argv[2],"%i",&islip);
  if(argc>3)			/* slip mode*/
    sscanf(argv[3],"%i",&i);
  slip_mode = (MODE_TYPE)i;
  /* 
     read geometry; read_geometry also sets the default elastic
     parameters in medium->elastic 
  */
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("relax_fault");
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  /* reset the elastic defaults ot unity shear modulus and 1/4 Poisson */
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25); 
  
  if((islip < 0)||(islip >= medium->nrflt)){
    fprintf(stderr,"%s: error: islip %i out of range 0...%i\n",
	    argv[0],islip,medium->nrflt-1);
    exit(-1);
  }
  fprintf(stderr,"%s: %i patches total, init slip: %g (mode: %i) on patch %i, t_M: %g Dt: %g t_max: %g (G: %g nu: %g)\n",
	  argv[0],medium->nrflt,slip,(int)slip_mode,islip,t_M,Dt,t_max,
	  medium->elastic.shear,medium->elastic.poisson);
  /* 
     full stress tensor state per patch, zero to start 
  */
  sigma=(COMP_PRECISION (*)[3][3])calloc(medium->nrflt,9*sizeof(COMP_PRECISION));
  if(!sigma)MEMERROR("relax_fault");
  /* 
     slip is applied at step 0 and then held constant 
  */
  slip_old[STRIKE] = slip_old[DIP] = slip_old[NORMAL] = 0.0;
  slip_new[STRIKE] = slip_new[DIP] = slip_new[NORMAL] = 0.0;
  slip_new[slip_mode] = slip;
  /* 
     time loop; step 0 uses Dt = 0 (alpha = 1, purely elastic) and
     zero old slip, which initializes the state to the elastic
     solution of the applied slip 
  */
  time = 0.0;
  k=0;
  while(time < t_max){
    relax_stress_update(medium,fault,islip,slip_new,slip_old,((k == 0)?(0.0):(Dt)),t_M,sigma);
    if(k == 0)			/* from now on, old slip = new slip */
      slip_old[slip_mode] = slip;
    else
      time += Dt;
    printf("%11g %5i\t",time,k);
    for(i=0;i < medium->nrflt;i++){
      /* resolve stress in the component of slip */
      printf("%12.5e ",resolve_stress_on_fault(sigma[i], (fault+i),slip_mode));
    }
    printf("\n");
    k++;
  }
  fprintf(stderr,"%s: computed %i steps\n",argv[0],k-1);
  free(sigma);free(fault);free(medium);
  return 0;
}
/* 

   one implicit Euler step of the effective-modulus visco-elastic
   stress update at all patch centroids, in place on the stress state
   sigma[nrflt][3][3]

     sigma <-- F(G'(Dt),K,slip_new) - dev(F(G'(Dt),K,slip_old)) + alpha(Dt) dev(sigma)

   slip slip_new/slip_old (strike, dip, normal components) is on patch
   islip only; Dt = 0 recovers the elastic evaluation; the reference
   elastic parameters are taken from medium->elastic, the bulk modulus
   is kept at its reference (elastic) value throughout

*/
void relax_stress_update(struct med *medium, struct flt *fault, int islip,
			 COMP_PRECISION *slip_new, COMP_PRECISION *slip_old,
			 COMP_PRECISION Dt, COMP_PRECISION t_M,
			 COMP_PRECISION (*sigma)[3][3])
{
  struct el_par elastic_par;
  COMP_PRECISION g0,nu0,bulk0,alphat,gp,nueff;
  COMP_PRECISION u[3],sm_new[3][3],sm_old[3][3];
  COMP_PRECISION trace_old,trace_sigma;
  int i,j,k,iret;
  /* 
     reference elastic and derived effective parameters; fixed K 
  */
  g0  = medium->elastic.shear;
  nu0 = medium->elastic.poisson;
  /* compute the rest */
  bulk0 = 2.0*g0*(1.0 + nu0)/(3.0*(1.0 - 2.0*nu0));
  /*  */
  alphat = t_M/(t_M + Dt);	/* = eta/(eta + Dt G) */
  
  gp = g0 * alphat;		/* G'(Dt) */
  nueff = (3.0*bulk0 - 2.0*gp)/(6.0*bulk0 + 2.0*gp); /* nu'(Dt) at fixed K */

  /* set the elastic parameters in the format needed for the Green's functions */
  calc_medium_elastic_parameters(&elastic_par,gp,nueff);
  /* 
     update all patches 
  */
  for(i=0;i < medium->nrflt;i++){
    /* stress at patch i centroid due to new and old slip on islip,
       both with the effective parameters of this time step */
    eval_green_at_receiver(fault,i,islip,slip_new,u,sm_new,&iret,
			   GC_STRESS_ONLY,TRUE,medium->full_space,elastic_par);
    if(iret)
      fprintf(stderr,"relax_stress_update: WARNING: singular at patch %i (new)\n",i);
    eval_green_at_receiver(fault,i,islip,slip_old,u,sm_old,&iret,
			   GC_STRESS_ONLY,TRUE,medium->full_space,elastic_par);
    if(iret)
      fprintf(stderr,"relax_stress_update: WARNING: singular at patch %i (old)\n",i);
    trace_old =   tracemat3x3(sm_old);
    trace_sigma = tracemat3x3(sigma[i]);
    for(j=0;j < 3;j++)
      for(k=0;k < 3;k++){
	/* new full stress minus deviator of old-slip stress plus
	   decayed deviator of the previous stress state */
	sigma[i][j][k] = sm_new[j][k] - sm_old[j][k] + alphat * sigma[i][j][k];
	if(j == k)
	  sigma[i][j][k] += (trace_old - alphat * trace_sigma)/3.0;
      }
  }
}
