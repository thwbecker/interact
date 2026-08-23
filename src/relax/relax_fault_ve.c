/*

  relax_fault_ve: prescribed-slip visco-elastic relaxation testbed
  using the exact Prony (correspondence principle) kernels of
  src/ve/prony_kernel.c; see rsf_ve_design.md and
  rsf_ve_implementation_plan.md, step 2

  a slip distribution is applied on fault[].u over a ramp of duration
  t_ramp (t_ramp = 0: instantaneous step at t = 0) and then held; the
  medium is Maxwell in shear with elastic bulk, Maxwell time t_M

  fault tractions are advanced by THREE routes that must agree, which
  exercises the machinery rsf_solve will use before any cycle run:

    A  closed form        tau(t) = sum_terms basis_term(t) [C_term u]
    B  ODE stepping       dh_p/dt = -h_p/tau_p + [C_p v(t)], RK4
    C  exact-step update  h_p <- e h_p + tau_p (1-e)/Dt [C_p dslip]

  route C is exact for piecewise-constant slip velocity, so with a
  ramp resolved by the time step it agrees with A to rounding; route
  B carries the RK4 truncation error, which shrinks as Dt^4

  surface displacements are tracked at nobs points along a profile at
  the free surface via the displacement Prony amplitudes; for the
  homogeneous medium the T_m amplitude vanishes and the constant
  term is the permanent (relaxed, effectively incompressible)
  deformation, so u_obs(t) evolves from the elastic to the relaxed
  field

  slip taper over a patch index fraction range as in relax_fault

  usage: relax_fault_ve [Dt, .1] [mode, 0] [lfrac, .4] [rfrac, .6] [t_M, 1] [t_ramp, 0] [t_max_fac, 10]

  reads geom.in; output files:

    rf_ve_stress.dat  per step: time, step, resolved slip-mode stress
                      on each patch (route A), same layout as
                      relax_fault output
    rf_ve_disp.dat    per step: time, step, then ux uy uz per
                      observation point (route A)
    rf_ve_route.dat   per step: time, max relative deviation of
                      routes B and C from A over all patches

  part of interact

*/
#include "interact.h"
#include "properties.h"


int main(int argc, char **argv)
{
  struct flt *fault;
  struct med *medium;
  struct prony_spec spec;
  /* per-patch, per-term stress amplitudes resolved in slip_mode
     direction: amp[iterm][ipatch]; per-obs displacement amplitudes
     damp[iterm][iobs][3] */
  COMP_PRECISION *amp[VE_MAX_NP+1],*damp[VE_MAX_NP+1];
  COMP_PRECISION *tau_a,*h_b[VE_MAX_NP],*h_c[VE_MAX_NP],*tmpv[4],*xobs;
  COMP_PRECISION C[VE_MAX_NP+1][3][3],D[VE_MAX_NP+1][3];
  COMP_PRECISION x[3],res,resmax_amp;
  COMP_PRECISION Dt,t_M,t_ramp,t_max,t_max_fac,time,sval,xloc;
  COMP_PRECISION efac,ifac,vfac,dev_b,dev_c,scl,k1,k2,k3,k4,vmid;
  COMP_PRECISION lfrac,rfrac,xmin,xmax,xlen,max_slip;
  MODE_TYPE slip_mode;
  FILE *fs,*fd,*fr;
  int i,j,p,it,k,iobs,nobs,ileft,iright,irange,nslip,nrflt,isrc;
  /* 
     defaults 
  */
  Dt = 0.1;			/* time step, units of t_M */
  slip_mode = STRIKE;
  lfrac = 0.4;rfrac = 0.6;	/* slipping patch fraction range */
  t_M = 1.0;			/* Maxwell time */
  t_ramp = 0.0;			/* slip ramp duration; 0 = step */
  t_max_fac = 50.0;		/* run to t_max_fac * t_M */
  nobs = 101;			/* surface observation points, a profile*/
  
  if(argc > 1)sscanf(argv[1],ONE_CP_FORMAT,&Dt);
  if(argc > 2){
    sscanf(argv[2],"%i",&i);
    slip_mode = (MODE_TYPE)i;
  }
  if(argc > 3)sscanf(argv[3],ONE_CP_FORMAT,&lfrac);
  if(argc > 4)sscanf(argv[4],ONE_CP_FORMAT,&rfrac);
  if(argc > 5)sscanf(argv[5],ONE_CP_FORMAT,&t_M);
  if(argc > 6)sscanf(argv[6],ONE_CP_FORMAT,&t_ramp);
  if(argc > 7)sscanf(argv[7],ONE_CP_FORMAT,&t_max_fac);
  t_max = t_max_fac * t_M;
  /* 
     geometry and reference elastic parameters 
  */
  medium=(struct med *)calloc(1,sizeof(struct med));
  if(!medium)MEMERROR("relax_fault_ve");
  read_geometry("geom.in",&medium,&fault,FALSE,FALSE,FALSE);
  calc_medium_elastic_parameters(&medium->elastic,1.0,0.25);
  nrflt = medium->nrflt;
  /* 
     slip assignment: elliptical taper over [lfrac, rfrac] of the
     patch index range, as in relax_fault 
  */
  ileft  = (int)(lfrac * nrflt);
  iright = (int)(rfrac * nrflt);
  if(ileft < 0)ileft = 0;
  if(iright > nrflt-1)iright = nrflt-1;
  irange = iright - ileft;
  if(irange == 0)irange++;
  nslip = 0;max_slip = 0.0;
  for(i=0;i < nrflt;i++){
    if((i >= ileft)&&(i <= iright)){
      xloc = ((COMP_PRECISION)i - (COMP_PRECISION)(ileft+iright)/2.0)/
	((COMP_PRECISION)irange/2.0);
      sval = 1.0 - xloc*xloc;
      sval = (sval > 0.0)?(sqrt(sval)):(0.0);
      /* assign slip */
      get_right_slip(fault[i].u,slip_mode,sval,(fault+i));
      sval = norm_3d(fault[i].u);
      if(sval > 0.0)
	nslip++;
      if(sval > max_slip)
	max_slip = sval;
    }
  }
  
  /* 
     Prony spec and amplitude assembly; stress amplitudes resolved in
     the slip_mode direction on each receiver, displacement amplitudes
     at nobs surface points along x spanning the fault twice over 
  */
  ve_spec_homogeneous(&spec,medium->elastic.shear,
		      medium->elastic.poisson,t_M);
  for(it=0;it < spec.nterm;it++){
    amp[it] =(COMP_PRECISION *)calloc(nrflt,sizeof(COMP_PRECISION));
    damp[it]=(COMP_PRECISION *)calloc(nobs*3,sizeof(COMP_PRECISION));
    if((!amp[it])||(!damp[it]))MEMERROR("relax_fault_ve");
  }
  tau_a=(COMP_PRECISION *)calloc(nrflt,sizeof(COMP_PRECISION));
  for(p=0;p < spec.np;p++){
    h_b[p]=(COMP_PRECISION *)calloc(nrflt,sizeof(COMP_PRECISION));
    h_c[p]=(COMP_PRECISION *)calloc(nrflt,sizeof(COMP_PRECISION));
    if((!h_b[p])||(!h_c[p]))MEMERROR("relax_fault_ve");
  }
  for(i=0;i < 4;i++){
    tmpv[i]=(COMP_PRECISION *)calloc(nrflt,sizeof(COMP_PRECISION));
    if(!tmpv[i])MEMERROR("relax_fault_ve");
  }
  xmin = xmax = fault[0].x[INT_X];
  for(i=1;i < nrflt;i++){
    if(fault[i].x[INT_X] < xmin)
      xmin = fault[i].x[INT_X];
    if(fault[i].x[INT_X] > xmax)
      xmax = fault[i].x[INT_X];
  }
  xlen = xmax - xmin;
  if(xlen <= 0.0)
    xlen = 1.0;

  /* profile location */
  xobs = (COMP_PRECISION *)malloc(nobs*sizeof(COMP_PRECISION));
  if(!xobs)MEMERROR("relax_fault_ve");
  for(iobs=0;iobs < nobs;iobs++){
    xobs[iobs] = -5 + 10*(COMP_PRECISION)iobs/((COMP_PRECISION)nobs - 1.0);
  }
  
  resmax_amp = 0.0;
  for(isrc=0;isrc < nrflt;isrc++){
    if(norm_3d(fault[isrc].u) < EPS_COMP_PREC)
      continue;
    /* stress amplitudes at all patch receivers */
    for(i=0;i < nrflt;i++){
      res = ve_prony_amplitudes_stress(&spec,medium,fault,i,isrc,
				       fault[isrc].u,C);
      if(res > resmax_amp)
	resmax_amp = res;
      for(it=0;it < spec.nterm;it++) /* this is the incremental stress effect */
	amp[it][i] += resolve_stress_on_fault(C[it],(fault+i),slip_mode);
    }
    /* displacement amplitudes at the surface profile */
    for(iobs=0;iobs < nobs;iobs++){
      x[INT_X] = xobs[iobs];x[INT_Y] = 0.0;x[INT_Z] = 0.0;
      /* compute the displacements */
      res = ve_prony_amplitudes_disp(&spec,medium,fault,x,isrc,fault[isrc].u,D);
      if(res > resmax_amp)resmax_amp = res;
      for(it=0;it < spec.nterm;it++)
	for(j=0;j < 3;j++)
	  damp[it][iobs*3+j] += D[it][j];
    }
  }
  fprintf(stderr,"%s: %i patches, %i slipping (%i to %i), max slip %g mode %i, t_M %g Dt %g t_ramp %g t_max %g, worst amplitude resid %8.1e\n",
	  argv[0],nrflt,nslip,ileft,iright,max_slip,slip_mode,
	  t_M,Dt,t_ramp,t_max,resmax_amp);
  if(resmax_amp > 1e-8)
    fprintf(stderr,"%s: WARNING: amplitude residual above expectation for homogeneous kernels\n",argv[0]);
  /* 
     time loop; routes B and C carry state h_p per patch, route A is
     evaluated in closed form each step; slip velocity is u/t_ramp on
     [0, t_ramp], zero after; for t_ramp = 0 routes B and C are
     initialized with the full step at t = 0 
  */
  fs = myopen("rf_ve_stress.dat","w");
  fd = myopen("rf_ve_disp.dat","w");
  fr = myopen("rf_ve_route.dat","w");
 
  if(t_ramp <= 0.0){		/* impulse initialization of the state */
    for(p=0;p < spec.np;p++)
      for(i=0;i < nrflt;i++)
	h_b[p][i] = h_c[p][i] = amp[p][i];
  }
  /* 

     time loop
  */
  time = 0.0;
  k = 0;
  while(time <= t_max + 1e-10){
    /* route A, closed form, also used for output */
    scl = 0.0;
    for(i=0;i < nrflt;i++){
      tau_a[i] = 0.0;
      for(it=0;it < spec.nterm;it++)
	tau_a[i] += ve_basis_time_ramp(&spec,it,time,t_ramp) * amp[it][i];
      if(fabs(tau_a[i]) > scl)
	scl = fabs(tau_a[i]);
    }
    /* route deviations */
    dev_b = dev_c = 0.0;
    if(scl > EPS_COMP_PREC){
      for(i=0;i < nrflt;i++){
	res = (t_ramp > 0.0)?
	  (((time <= t_ramp)?(time/t_ramp):(1.0)) * amp[spec.np][i]):
	  (amp[spec.np][i]);	/* constant-term part */
	k1 = res;k2 = res;	/* reuse as accumulators */
	for(p=0;p < spec.np;p++){
	  k1 += h_b[p][i];
	  k2 += h_c[p][i];
	}
	if(fabs(k1 - tau_a[i])/scl > dev_b)
	  dev_b = fabs(k1 - tau_a[i])/scl;
	if(fabs(k2 - tau_a[i])/scl > dev_c)
	  dev_c = fabs(k2 - tau_a[i])/scl;
      }
    }
    /* 
       output 
    */
    /* stress */
    fprintf(fs,"%11g %5i\t",time,k);
    for(i=0;i < nrflt;i++)
      fprintf(fs,"%12.5e ",tau_a[i]);
    fprintf(fs,"\n");
    /* displacement */
    fprintf(fd,"%11g\t",time);
    for(iobs=0;iobs < nobs;iobs++){
      fprintf(fd,"%11g ",xobs[iobs]); /* location */
      for(j=0;j < 3;j++){	       /* compute three displacements */
	res = 0.0;
	for(it=0;it < spec.nterm;it++) /* sum up */
	  res += ve_basis_time_ramp(&spec,it,time,t_ramp) * damp[it][iobs*3+j];
	fprintf(fd,"%12.5e ",res);
      }
    }
    fprintf(fd,"\n");
    /* relax */
    fprintf(fr,"%11g %12.5e %12.5e\n",time,dev_b,dev_c);
    if(time >= t_max)
      break;
    /* 
       advance routes B and C over [time, time + Dt]; slip velocity
       of patch i in this window (piecewise constant approximation
       used by both routes for the forcing) 
    */
    vfac = 0.0;			/* fraction of total slip per unit time */
    if(t_ramp > 0.0){
      /* overlap of [time, time+Dt] with [0, t_ramp] */
      res = ((time + Dt < t_ramp)?(time + Dt):(t_ramp)) - time;
      if(res > 0.0)
	vfac = res/(Dt * t_ramp);
    }
    /* Runge Kutta time advance for testing purposes */
    for(p=0;p < spec.np;p++){
      efac = exp(-Dt/spec.tau[p]);
      ifac = spec.tau[p]*(1.0 - efac)/Dt;
      /* route C, exact for piecewise-constant velocity */
      for(i=0;i < nrflt;i++)
	h_c[p][i] = efac * h_c[p][i] + ifac * amp[p][i] * vfac * Dt;
      /* route B, classical RK4 on dh/dt = -h/tau + amp*vfac */
      for(i=0;i < nrflt;i++){
	vmid = amp[p][i] * vfac;
	k1 = -h_b[p][i]/spec.tau[p] + vmid;
	k2 = -(h_b[p][i] + 0.5*Dt*k1)/spec.tau[p] + vmid;
	k3 = -(h_b[p][i] + 0.5*Dt*k2)/spec.tau[p] + vmid;
	k4 = -(h_b[p][i] +     Dt*k3)/spec.tau[p] + vmid;
	h_b[p][i] += Dt*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
      }
    }
    time += Dt;
    k++;
  }
  fprintf(stderr,"%s: computed %i steps\n",argv[0],k);
  fclose(fs);fclose(fd);fclose(fr);
  for(it=0;it < spec.nterm;it++){
    free(amp[it]);free(damp[it]);
  }
  for(p=0;p < spec.np;p++){
    free(h_b[p]);free(h_c[p]);
  }
  for(i=0;i < 4;i++)free(tmpv[i]);
  free(tau_a);free(fault);free(medium);
  return 0;
}
