#include "interact.h"
#include "properties.h"
/*

  visco-elastic hereditary stressing for rsf_solve: Prony
  representation of the step-response interaction kernel

    K(t) = C_inf + sum_{p=1..np} C_p exp(-t/tau_p),  K(0) = K^elastic

  carried by per-patch memory states h_p with

    dh_p/dt = -h_p/tau_p + [C_p (v - vpl)],
    tau_dot = [K^el (v - vpl)] - sum_p h_p/tau_p

  (the relaxed operator C_inf never appears explicitly: it is
  implicitly K^el - sum_p C_p, so the elastic response stays exactly
  the assembled Is matrix).  This is the generator-agnostic
  rsf_ve_implementation_plan.md machinery: the dynamics below consume
  only rates and amplitude operators, so 3-D kernels (e.g. from the
  plate-over-Maxwell propagator, via -ve_prony_file, mode 3) plug in
  without touching the time stepping.

  h treatment: the memory-light "step" mode (plan step 5): h lives
  OUTSIDE the TS state.  Between accepted steps the exact-exponential
  update

    h_p <- h_p e^(-dt/tau_p) + [C_p (vbar - vpl)] tau_p (1 - e^(-dt/tau_p))

  is applied in a TS monitor, with vbar = (slip increment)/dt the
  EXACT step-average velocity (slip is part of the TS state).  Inside
  a step the RHS sees the decay-exact, forcing-lagged states
  h_p(t) = h_p(t0) e^(-(t-t0)/tau_p), applied as a per-pole VecAXPY
  in rsf_compute_vel_and_stressing.  The lag error is second order in
  the step and negligible wherever VE matters (tau_p >> coseismic
  time scales).

  kernel generators (option -ve_mode):

    1  uniform Maxwell medium.  Implemented for 2-D antiplane
       geometries, where the kernel depends on the moduli only
       through G so the representation is a SINGLE exponential with
       tau_1 = t_M and C_1 = K^el, exactly (no sampling, no fit).
       (General 3-D geometries have the three-pole homogeneous
       representation of prony_kernel.c; not wired here yet.)

    2  elastic plate over Maxwell half space, 2-D antiplane: the
       validated layered image kernel (both image families, see
       ve_rybicki_antiplane_uz) sampled in the Laplace domain through
       the actual elements and fitted per pair on a geometric rate
       ladder [b/45, 1.5 b], b = g1/((g1+g2) t_M), samples
       [b/100, 5 b], with two held-out samples verified; the fitted
       constant is discarded (relaxed operator implicit).  This is
       the operator-level version of the ve_sp_cycle testbed.

  options (all parsed here):
    -ve_mode <1|2|3>       generator (required to switch VE on);
                           3 = external kernel from -ve_prony_file
    -ve_prony_file <f>     mode 3: shared taus + per-pair amplitude
                           matrices + a K0 block for the gate below;
                           the header's third field is the number of
                           TRACTION FAMILIES: 1 = shear only,
                           2 = shear plus fault-NORMAL traction (the
                           normal family uses the same tau ladder and
                           is activated exactly when -calc_sigma_dot
                           is on; a shear-only file with
                           -calc_sigma_dot is refused, since relaxing
                           shear with frozen normal is not a medium)
    -ve_prony_k0tol <t>    mode 3: relative tolerance of the K0 vs Is
                           consistency gate (default 0.05 warn, 5x abort)
    -ve_tmaxwell_yr <t>    Maxwell time of the relaxing medium [yr]
    -ve_plate_h <H>        elastic plate thickness [m] (mode 2)
    -ve_g2fac <f>          substrate/plate rigidity ratio (mode 2, default 1)
    -ve_np <n>             exponential terms (mode 2, default 6, max 6)
    -ve_h_stage <0|1>      1 (default) stage-consistent sink forcing;
                           0 frozen decay-only sink (see rsf_ve_apply_sink)
    -ve_h_init <0|1>       initial memory states: 1 (default) backslip
                           steady state (no secular VE spin-up),
                           0 virgin medium (zero memory)

  the sampling image depth and the held-out abort threshold are
  compile-time constants (VE_NIMG_DEF, VE_FIT_TOL_DEF in
  properties.h): good values are geometry independent, see there

  restart: the h states, slip_prev, and t0 travel in the checkpoint
  (version 1.1); a VE run cannot restart from an elastic checkpoint
  and vice versa.

  part of interact

*/
#ifdef USE_PETSC
#include "petsc_prototypes.h"
#include "rsf.h"
#include "prony_kernel.h"

/*
   parse options and, if VE is requested, build rates, amplitude
   operators, and state vectors.  Called from rsf_solve_run AFTER the
   Is matrix and the backslip products exist.  negvpl_in holds -vpl
   per patch (the slip-rate vector used for the sinc products).
   scale converts interact-internal stress to SI, as for Is.
*/
static int cmp_comp_precision(const void *, const void *);
/*
   -ve_mode 3 consistency gate: read one n x n block from the kernel
   file and compare it against an assembled interaction operator
   (Is for shear, In for normal) in the Frobenius sense, skipping
   near-diagonal entries where an external generator cannot
   reproduce the singular elastic self terms (those live in the
   assembled operator, i.e. in the implicit C_inf).  Catches unit,
   orientation and SIGN mismatches before any physics runs: a pure
   global sign flip shows up as a relative deviation of 2.
*/
static PetscErrorCode rsf_ve_prony_gate(FILE *pf, Mat A, struct med *medium,
					PetscReal *row, const char *what,
					const char *pfile, PetscReal tol)
{
  PetscInt i,j;
  PetscReal dev = 0.0,nrm = 0.0,dv,gred[2],lred[2];
  PetscFunctionBeginUser;
  for(i=0;i < medium->nrflt;i++){
    for(j=0;j < medium->nrflt;j++)
      if(fscanf(pf,"%lf",row+j) != 1)
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_READ,
		"rsf_ve_setup: bad %s-gate row %ld",what,(long)i);
    if((i >= medium->rs) && (i < medium->re)){
      PetscInt ncols;
      const PetscInt *cols;
      const PetscScalar *vals;
      PetscCall(MatGetRow(A,i,&ncols,&cols,&vals));
      for(j=0;j < ncols;j++){
	if(labs((long)i - (long)cols[j]) <= 2)
	  continue;
	dv = PetscRealPart(vals[j]) - row[cols[j]];
	dev += dv*dv;
	nrm += PetscRealPart(vals[j])*PetscRealPart(vals[j]);
      }
      PetscCall(MatRestoreRow(A,i,&ncols,&cols,&vals));
    }
  }
  lred[0] = dev; lred[1] = nrm;
  PetscCallMPI(MPI_Allreduce(lred,gred,2,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD));
  dev = sqrt(gred[0]/(gred[1] + 1e-300));
  HEADNODE
    fprintf(stderr,"rsf_ve_setup: mode 3 (%s): %s block, rel dev from assembled operator %.3e\n",
	    pfile,what,(double)dev);
  if(dev > 5.0*tol)
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
	    "rsf_ve_setup: prony %s block deviates %.3e from the assembled operator (tol %g): convention mismatch?  (a global sign flip gives exactly 2)",
	    what,(double)dev,(double)tol);
  if(dev > tol)
    HEADNODE
      fprintf(stderr,"rsf_ve_setup: WARNING: prony %s deviation %.3e exceeds %g\n",
	      what,(double)dev,(double)tol);
  PetscFunctionReturn(PETSC_SUCCESS);
}
PetscErrorCode rsf_ve_setup(struct interact_ctx *par, Vec negvpl_in,
			    PetscReal scale)
{
  struct med *medium = par->medium;
  struct rsf_vars *rsf = medium->rsf;
  struct flt *fault = par->fault;
  struct prony_spec spec;
  PetscReal tm_yr=0.0,plate_h=0.0,g2fac=1.0;
  const PetscReal fit_tol = VE_FIT_TOL_DEF;
  PetscReal gamma0,rate_b,worst_res;
  PetscInt np=6,mode=0;
  const PetscInt n_img = VE_NIMG_DEF;
  PetscInt i,j,k,n,it,p,jn;
  PetscBool flg;
  PetscFunctionBeginUser;
  rsf->ve_np = 0;rsf->ve_t0 = 0.0;rsf->ve_normal = PETSC_FALSE;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ve_mode",&mode,&flg));
  if((!flg) || (mode == 0))
    PetscFunctionReturn(PETSC_SUCCESS); /* VE off: nothing changes */
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-ve_tmaxwell_yr",&tm_yr,&flg));
  if(((!flg)||(tm_yr <= 0.0)) && (mode != 3))
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
	    "rsf_ve_setup: -ve_mode %ld requires -ve_tmaxwell_yr > 0",(long)mode);
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ve_np",&np,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-ve_plate_h",&plate_h,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-ve_g2fac",&g2fac,NULL));
  /* the two internal generators require pure 2-D antiplane geometry;
     mode 3 takes its kernel from a file and is element-agnostic */
  if(mode != 3)
    for(i=0;i < medium->nrflt;i++)
      if(fault[i].type != TWO_DIM_ANTIPLANE)
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,
		"rsf_ve_setup: -ve_mode %ld: patch %ld type %i, only 2-D antiplane elements are wired up (external kernels: -ve_mode 3 -ve_prony_file)",
		(long)mode,(long)i,(int)fault[i].type);
  if(rsf->calc_sigma_dot && (mode != 3))
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,
	    "rsf_ve_setup: -ve_mode %ld has no hereditary normal-stress coupling (the internal generators are 2-D ANTIPLANE, where slip induces no normal-stress change); use -ve_mode 3 with a kernel file carrying normal blocks, or drop -calc_sigma_dot",
	    (long)mode);
  rsf->ve_mode = mode;
  rsf->ve_hstage = 1;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ve_h_stage",&rsf->ve_hstage,NULL));
  if(mode == 1){
    /* uniform Maxwell, antiplane: exactly one pole, C_1 = K^el */
    rsf->ve_np = 1;
    rsf->ve_tau[0] = tm_yr*SEC_PER_YEAR;
    PetscCall(MatDuplicate(medium->Is,MAT_COPY_VALUES,&rsf->ve_C[0]));
    HEADNODE
      fprintf(stderr,"rsf_ve_setup: mode 1 (uniform Maxwell, antiplane): single EXACT pole, tau = %g yr, C_1 = Is\n",
	      tm_yr);
  }else if(mode == 3){
    /*
       GENERATOR-AGNOSTIC EXTERNAL KERNEL (-ve_prony_file): shared
       relaxation times and per-pair amplitude matrices fitted
       OUTSIDE this code from step responses of any layered/graded
       medium (current generator: the plane-strain plate-over-Maxwell
       + gravity prototype in test_relax/inplane_ve_proto; the same
       file contract is intended for PSGRN/PSCMP-derived 3-D layered
       viscoelastic-gravitational kernels).  ASCII format:

	 # comment lines
	 n np
	 tau_1 ... tau_np           [s]
	 K0 matrix, n rows x n cols [Pa/m]  (generator's ELASTIC kernel)
	 C_1 matrix ... C_np matrix [Pa/m]

       K(t) = C_inf + sum_p C_p exp(-t/tau_p) with C_inf implicit
       (= assembled elastic Is - sum_p C_p), exactly as modes 1/2.
       K0 is used only for a consistency gate against the assembled
       elastic operator: the generator and this code must agree on
       the elastic kernel (norm tolerance -ve_prony_k0tol, default
       0.05 warn / 5x abort), which pins discretization, units,
       orientation and sign conventions in one check.
       SHEAR ONLY: normal-stress hereditary relaxation is not
       implemented (-calc_sigma_dot is refused above).
    */
    char pfile[300];
    FILE *pf;
    PetscReal k0tol = 0.05;
    PetscInt nf,npf,nfam;
    PetscReal *row;
    PetscScalar *crow3;
    PetscInt *cidx3;
    PetscCall(PetscOptionsGetString(NULL,NULL,"-ve_prony_file",pfile,
				    sizeof(pfile),&flg));
    if(!flg)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
	      "rsf_ve_setup: -ve_mode 3 requires -ve_prony_file <file>");
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ve_prony_k0tol",&k0tol,NULL));
    pf = fopen(pfile,"r");
    if(!pf)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_OPEN,
	      "rsf_ve_setup: cannot open %s",pfile);
    {
      char lb[65536];
      /* skip comments */
      do{
	if(!fgets(lb,sizeof(lb),pf))
	  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_READ,"rsf_ve_setup: truncated prony file");
      }while(lb[0] == '#');
      nfam = 1;			/* default: shear only (format v1) */
      if(sscanf(lb,"%ld %ld %ld",(long *)&nf,(long *)&npf,(long *)&nfam) < 2)
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_READ,"rsf_ve_setup: bad prony header");
      if((nfam < 1) || (nfam > 2))
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
		"rsf_ve_setup: prony file declares %ld traction families, expect 1 (shear) or 2 (shear + normal)",
		(long)nfam);
      if(nf != medium->nrflt)
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
		"rsf_ve_setup: prony file n %ld != %ld patches",(long)nf,(long)medium->nrflt);
      if((npf < 1)||(npf > VE_MAX_NP))
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
		"rsf_ve_setup: prony file np %ld out of [1,%d]",(long)npf,VE_MAX_NP);
      rsf->ve_np = npf;
      for(p=0;p < npf;p++)
	if(fscanf(pf,"%lf",&rsf->ve_tau[p]) != 1)
	  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_READ,"rsf_ve_setup: bad tau list");
    }
    for(p=0;p < rsf->ve_np;p++){
      PetscCall(MatCreate(PETSC_COMM_WORLD,&rsf->ve_C[p]));
      PetscCall(MatSetSizes(rsf->ve_C[p],medium->rn,PETSC_DECIDE,
			    medium->nrflt,medium->nrflt));
      PetscCall(MatSetType(rsf->ve_C[p],MATDENSE));
      PetscCall(MatSetUp(rsf->ve_C[p]));
    }
    row = (PetscReal *)malloc((size_t)medium->nrflt*sizeof(PetscReal));
    crow3 = (PetscScalar *)malloc((size_t)medium->nrflt*sizeof(PetscScalar));
    cidx3 = (PetscInt *)malloc((size_t)medium->nrflt*sizeof(PetscInt));
    if((!row)||(!crow3)||(!cidx3))MEMERROR("rsf_ve_setup");
    for(j=0;j < medium->nrflt;j++)
      cidx3[j] = j;
    PetscCall(rsf_ve_prony_gate(pf,medium->Is,medium,row,"K0",pfile,k0tol));
    HEADNODE
      fprintf(stderr,"rsf_ve_setup: mode 3 (%s): np %ld, tau %g ... %g yr\n",
	      pfile,(long)rsf->ve_np,rsf->ve_tau[0]/SEC_PER_YEAR,
	      rsf->ve_tau[rsf->ve_np-1]/SEC_PER_YEAR);
    for(p=0;p < rsf->ve_np;p++){
      for(i=0;i < medium->nrflt;i++){
	for(j=0;j < medium->nrflt;j++)
	  if(fscanf(pf,"%lf",row+j) != 1)
	    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_READ,
		    "rsf_ve_setup: bad C_%ld row %ld",(long)p,(long)i);
	if((i >= medium->rs) && (i < medium->re)){
	  for(j=0;j < medium->nrflt;j++)
	    crow3[j] = row[j];
	  PetscCall(MatSetValues(rsf->ve_C[p],1,&i,medium->nrflt,cidx3,crow3,INSERT_VALUES));
	}
      }
      PetscCall(MatAssemblyBegin(rsf->ve_C[p],MAT_FINAL_ASSEMBLY));
      PetscCall(MatAssemblyEnd(rsf->ve_C[p],MAT_FINAL_ASSEMBLY));
    }
    /*
       SECOND FAMILY (nfam == 2): the same file continues with an In
       gate block and the normal amplitude matrices.  It is READ only
       when the elastic normal path is active, since otherwise sigma
       does not evolve at all and hereditary normal terms would have
       nothing to act on; conversely a run WITH -calc_sigma_dot and a
       shear-only file would be an inconsistent medium (relaxing
       shear, frozen normal) and is refused.
    */
    if(rsf->calc_sigma_dot && (nfam < 2))
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
	      "rsf_ve_setup: -calc_sigma_dot needs a kernel file with normal blocks (header 'n np 2'); %s is shear only.  Regenerate it (bp3_ve_kernels.py ... -normal) or drop -calc_sigma_dot",
	      pfile);
    if((nfam == 2) && (!rsf->calc_sigma_dot)){
      HEADNODE
	fprintf(stderr,"rsf_ve_setup: NOTE: %s carries normal-traction kernels, ignored without -calc_sigma_dot\n",
		pfile);
    }else if(nfam == 2){
      if(!medium->In)
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONGSTATE,
		"rsf_ve_setup: normal kernel family requested but the In matrix was not assembled");
      for(p=0;p < rsf->ve_np;p++){
	PetscCall(MatCreate(PETSC_COMM_WORLD,&rsf->ve_Cn[p]));
	PetscCall(MatSetSizes(rsf->ve_Cn[p],medium->rn,PETSC_DECIDE,
			      medium->nrflt,medium->nrflt));
	PetscCall(MatSetType(rsf->ve_Cn[p],MATDENSE));
	PetscCall(MatSetUp(rsf->ve_Cn[p]));
      }
      PetscCall(rsf_ve_prony_gate(pf,medium->In,medium,row,"N0",pfile,k0tol));
      for(p=0;p < rsf->ve_np;p++){
	for(i=0;i < medium->nrflt;i++){
	  for(j=0;j < medium->nrflt;j++)
	    if(fscanf(pf,"%lf",row+j) != 1)
	      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_READ,
		      "rsf_ve_setup: bad Cn_%ld row %ld",(long)p,(long)i);
	  if((i >= medium->rs) && (i < medium->re)){
	    for(j=0;j < medium->nrflt;j++)
	      crow3[j] = row[j];
	    PetscCall(MatSetValues(rsf->ve_Cn[p],1,&i,medium->nrflt,cidx3,crow3,INSERT_VALUES));
	  }
	}
	PetscCall(MatAssemblyBegin(rsf->ve_Cn[p],MAT_FINAL_ASSEMBLY));
	PetscCall(MatAssemblyEnd(rsf->ve_Cn[p],MAT_FINAL_ASSEMBLY));
      }
      rsf->ve_normal = PETSC_TRUE;
    }
    free(row);free(crow3);free(cidx3);
    fclose(pf);
  }else if(mode == 2){
    /* layered antiplane: sampled image kernel, per-pair Prony fit */
    struct flt img[1];
    COMP_PRECISION slip[3],u[3],sm[3][3],gam[VE_MAX_NS],wfac;
    COMP_PRECISION *srow,smp[VE_MAX_NS],amp[VE_MAX_NP+1],res,scl,elastic_dev,val;
    PetscScalar *crow;
    PetscInt *colidx;
    /* translational-invariance sample cache: the sampled kernel of a
       pair depends only on (x_i - x_j, z_i, z_j) (all elements are
       parallel vertical antiplane patches of equal half-length), so
       multi-fault geometries repeat the expensive image walk across
       receiver rows.  Key: unique dx values x unique z (receiver) x
       unique z (source); disabled if the elements have mixed
       half-lengths or the table would be too large */
    int cache_on;
    long ndx=0,nuz=0,ncache=0;
    COMP_PRECISION *udx=NULL,*uzz=NULL,*csamp=NULL;
    int *izof=NULL;
    char *cset=NULL;
    if((np < 2)||(np > VE_MAX_NP))
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"rsf_ve_setup: -ve_np %ld out of [2, %d]",(long)np,VE_MAX_NP);
    if(plate_h <= 0.0)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"rsf_ve_setup: -ve_mode 2 requires -ve_plate_h [m] > 0");
    for(i=0;i < medium->nrflt;i++)
      if(-(fault[i].x[INT_Y] - fault[i].l) > plate_h + 1e-3)
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,
		"rsf_ve_setup: patch %ld reaches depth %g m below the plate (H = %g m)",
		(long)i,(double)(-(fault[i].x[INT_Y]-fault[i].l)),(double)plate_h);
    ve_layered_gamma_pars(medium->elastic.shear,g2fac*medium->elastic.shear,
			  tm_yr*SEC_PER_YEAR,&gamma0,&rate_b);
    /* spec: ladder [b/45, 1.5 b], samples [b/100, 5 b], as validated
       in the ve_sp_cycle testbed */
    spec.g0 = medium->elastic.shear;spec.nu0 = medium->elastic.poisson;
    spec.t_M = tm_yr*SEC_PER_YEAR;
    spec.bulk0 = bulk_mod_from_G_nu(spec.g0,spec.nu0);
    spec.np = np;
    for(i=0;i < np;i++)
      spec.tau[i] = 1.0/(rate_b*(1.0/45.0)*
			 pow(67.5,(COMP_PRECISION)i/((COMP_PRECISION)(np-1))));
    spec.has_const = TRUE;	/* fitted, then discarded (implicit) */
    spec.nterm = np + 1;
    spec.ns = spec.nterm + VE_NHELD;
    for(i=0;i < spec.nterm;i++)
      spec.sk[i] = rate_b*(1.0/100.0)*
	pow(500.0,(COMP_PRECISION)i/((COMP_PRECISION)spec.nterm - 1.0));
    spec.sk[spec.nterm]   = rate_b/30.0;
    spec.sk[spec.nterm+1] = rate_b*0.6;
    ve_solve_weights(&spec);
    for(k=0;k < spec.ns;k++)
      gam[k] = ve_layered_gamma_of_s(gamma0,rate_b,spec.sk[k]);
    rsf->ve_np = np;
    for(p=0;p < np;p++){
      rsf->ve_tau[p] = spec.tau[p];
      PetscCall(MatCreate(PETSC_COMM_WORLD,&rsf->ve_C[p]));
      PetscCall(MatSetSizes(rsf->ve_C[p],medium->rn,PETSC_DECIDE,
			    medium->nrflt,medium->nrflt));
      PetscCall(MatSetType(rsf->ve_C[p],MATDENSE));
      PetscCall(MatSetUp(rsf->ve_C[p]));
    }
    /* sample rows for the owned receivers: srow[k*n + j] holds
       sample k of the kernel row; the image families are walked once
       per (i,j) and accumulated into all samples */
    srow = (COMP_PRECISION *)malloc((size_t)spec.ns*medium->nrflt*sizeof(COMP_PRECISION));
    colidx = (PetscInt *)malloc((size_t)medium->nrflt*sizeof(PetscInt));
    crow = (PetscScalar *)malloc((size_t)medium->nrflt*sizeof(PetscScalar));
    if((!srow)||(!colidx)||(!crow))MEMERROR("rsf_ve_setup");
    for(j=0;j < medium->nrflt;j++)
      colidx[j] = j;
    slip[STRIKE] = 1.0;slip[DIP] = slip[NORMAL] = 0.0;
    worst_res = elastic_dev = 0.0;
    /* build the sample cache index tables */
    cache_on = 1;
    for(i=1;i < medium->nrflt;i++)
      if(fabs(fault[i].l - fault[0].l) > 1e-6){cache_on = 0;break;}
    if(cache_on){
      COMP_PRECISION *tmpv;
      long ii,jj,nn;
      tmpv = (COMP_PRECISION *)malloc((size_t)medium->nrflt*sizeof(COMP_PRECISION));
      izof = (int *)malloc((size_t)medium->nrflt*sizeof(int));
      if((!tmpv)||(!izof))MEMERROR("rsf_ve_setup: cache");
      /* unique z (INT_Y) values, 1e-4 m quantization */
      for(ii=0;ii < medium->nrflt;ii++)tmpv[ii] = fault[ii].x[INT_Y];
      qsort(tmpv,(size_t)medium->nrflt,sizeof(COMP_PRECISION),cmp_comp_precision);
      uzz = (COMP_PRECISION *)malloc((size_t)medium->nrflt*sizeof(COMP_PRECISION));
      if(!uzz)MEMERROR("rsf_ve_setup: cache");
      for(nn=0,ii=0;ii < medium->nrflt;ii++)
	if((nn == 0)||(fabs(tmpv[ii]-uzz[nn-1]) > 1e-4))uzz[nn++] = tmpv[ii];
      nuz = nn;
      for(ii=0;ii < medium->nrflt;ii++){
	for(jj=0;jj < nuz;jj++)
	  if(fabs(fault[ii].x[INT_Y]-uzz[jj]) <= 1e-4)break;
	izof[ii] = (int)jj;
      }
      /* unique signed dx = x_i - x_j values over all pairs of unique x */
      {
	COMP_PRECISION *ux;
	long nux;
	for(ii=0;ii < medium->nrflt;ii++)tmpv[ii] = fault[ii].x[INT_X];
	qsort(tmpv,(size_t)medium->nrflt,sizeof(COMP_PRECISION),cmp_comp_precision);
	ux = (COMP_PRECISION *)malloc((size_t)medium->nrflt*sizeof(COMP_PRECISION));
	if(!ux)MEMERROR("rsf_ve_setup: cache");
	for(nn=0,ii=0;ii < medium->nrflt;ii++)
	  if((nn == 0)||(fabs(tmpv[ii]-ux[nn-1]) > 1e-4))ux[nn++] = tmpv[ii];
	nux = nn;
	udx = (COMP_PRECISION *)malloc((size_t)nux*nux*sizeof(COMP_PRECISION));
	if(!udx)MEMERROR("rsf_ve_setup: cache");
	for(nn=0,ii=0;ii < nux;ii++)
	  for(jj=0;jj < nux;jj++)
	    udx[nn++] = ux[ii]-ux[jj];
	qsort(udx,(size_t)(nux*nux),sizeof(COMP_PRECISION),cmp_comp_precision);
	for(ndx=0,ii=0;ii < nux*nux;ii++)
	  if((ndx == 0)||(fabs(udx[ii]-udx[ndx-1]) > 1e-4))udx[ndx++] = udx[ii];
	free(ux);
      }
      /* per-fault dx index against fault 0 is NOT enough; store each
	 fault's x for on-the-fly bsearch of dx below instead */
      ncache = ndx*nuz*nuz;
      if((ncache > 0) && (ncache*(long)spec.ns*8L < 512L*1024L*1024L)){
	csamp = (COMP_PRECISION *)calloc((size_t)(ncache*spec.ns),sizeof(COMP_PRECISION));
	cset  = (char *)calloc((size_t)ncache,sizeof(char));
	if((!csamp)||(!cset))cache_on = 0;
      }else{
	cache_on = 0;
      }
      free(tmpv);
      HEADNODE
	if(cache_on)
	  fprintf(stderr,"rsf_ve_setup: sample cache on: %ld dx x %ld^2 z = %ld entries (vs %ld pairs)\n",
		  ndx,nuz,ncache,(long)medium->nrflt*(long)medium->nrflt);
    }
    HEADNODE
      fprintf(stderr,"rsf_ve_setup: mode 2 (plate over Maxwell, antiplane): H %g m g2/g1 %g t_M %g yr, np %ld, n_img %ld, sampling...\n",
	      (double)plate_h,(double)g2fac,(double)tm_yr,(long)np,(long)n_img);
    for(i=medium->rs;i < medium->re;i++){
      int iret;
      memset(srow,0,(size_t)spec.ns*medium->nrflt*sizeof(COMP_PRECISION));
      for(j=0;j < medium->nrflt;j++){
	long ci = -1;
	if(cache_on){
	  COMP_PRECISION dx = fault[i].x[INT_X]-fault[j].x[INT_X];
	  long lo=0,hi=ndx-1,mid,kdx=-1;
	  while(lo <= hi){
	    mid = (lo+hi)/2;
	    if(fabs(udx[mid]-dx) <= 1e-4){kdx = mid;break;}
	    if(udx[mid] < dx)lo = mid+1;else hi = mid-1;
	  }
	  if(kdx >= 0){
	    ci = (kdx*nuz + (long)izof[i])*nuz + (long)izof[j];
	    if(cset[ci]){
	      for(k=0;k < spec.ns;k++)
		srow[k*medium->nrflt+j] = csamp[ci*spec.ns+k];
	      continue;
	    }
	  }
	}
	for(n=0;n <= n_img;n++){
	  /* image family n: plain shift and (n >= 1) the surface
	     mirror, each with its own auto free-surface image via
	     the half-plane kernel; see ve_rybicki_antiplane_uz */
	  for(jn=0;jn < ((n == 0)?(1):(2));jn++){
	    if(n == 0){
	      eval_green_at_receiver(fault,i,j,slip,u,sm,&iret,
				     GC_STRESS_ONLY,TRUE,FALSE,
				     medium->elastic);
	    }else{
	      *img = fault[j];
	      img->x[INT_Y] = ((jn)?(-fault[j].x[INT_Y]):(fault[j].x[INT_Y]))
		- 2.0*(COMP_PRECISION)n*plate_h;
	      eval_green(fault[i].x,img,slip,u,sm,&iret,
			 GC_STRESS_ONLY,FALSE,FALSE,medium->elastic);
	    }
	    val = resolve_stress_on_fault(sm,(fault+i),STRIKE);
	    if(n == 0){
	      for(k=0;k < spec.ns;k++)
		srow[k*medium->nrflt+j] += val;
	    }else{
	      for(k=0;k < spec.ns;k++){
		wfac = pow(gam[k],(COMP_PRECISION)n);
		srow[k*medium->nrflt+j] += wfac*val;
	      }
	    }
	  }
	  /* all samples converged geometrically? largest |gam| decides */
	  if(pow(fabs(gam[0]) > fabs(gam[spec.ns-1])?fabs(gam[0]):fabs(gam[spec.ns-1]),
		 (COMP_PRECISION)(n+1)) < 1e-16)
	    break;
	}
	if(cache_on && (ci >= 0)){
	  for(k=0;k < spec.ns;k++)
	    csamp[ci*spec.ns+k] = srow[k*medium->nrflt+j];
	  cset[ci] = 1;
	}
      }
      /* per-pair fit and fill of the np amplitude rows */
      for(p=0;p < np;p++){
	for(j=0;j < medium->nrflt;j++){
	  for(k=0;k < spec.ns;k++)
	    smp[k] = srow[k*medium->nrflt+j];
	  for(it=0;it < spec.nterm;it++){
	    amp[it] = 0.0;
	    for(k=0;k < spec.nterm;k++)
	      amp[it] += spec.W[it][k]*smp[k];
	  }
	  if(p == 0){		/* held-out and elastic checks, once per pair */
	    scl = 0.0;
	    for(k=0;k < spec.ns;k++)
	      if(fabs(smp[k]) > scl)scl = fabs(smp[k]);
	    if(scl > 0.0){
	      for(k=spec.nterm;k < spec.ns;k++){
		res = 0.0;
		for(it=0;it < spec.nterm;it++)
		  res += amp[it]*ve_basis(&spec,it,spec.sk[k]);
		res = fabs(res - smp[k])/scl;
		if(res > worst_res)worst_res = res;
	      }
	    }
	  }
	  crow[j] = (PetscScalar)(amp[p]*scale);
	}
	PetscCall(MatSetValues(rsf->ve_C[p],1,&i,medium->nrflt,colidx,crow,INSERT_VALUES));
      }
    }
    free(srow);free(colidx);free(crow);
    if(udx)free(udx);
    if(uzz)free(uzz);
    if(izof)free(izof);
    if(csamp)free(csamp);
    if(cset)free(cset);
    for(p=0;p < np;p++){
      PetscCall(MatAssemblyBegin(rsf->ve_C[p],MAT_FINAL_ASSEMBLY));
      PetscCall(MatAssemblyEnd(rsf->ve_C[p],MAT_FINAL_ASSEMBLY));
    }
    /* elastic-consistency: || (sum_p C_p) x ||: the exponential sum
       may deviate from Is by the implicit relaxed operator C_inf,
       which for a fault NOT reaching through the whole plate is
       genuinely nonzero; report the norm ratio for information */
    {
      Vec xr,y1,y2;
      PetscReal n1,n2;
      PetscCall(MatCreateVecs(medium->Is,&xr,&y1));
      PetscCall(VecDuplicate(y1,&y2));
      PetscCall(VecSet(xr,1.0));
      PetscCall(MatMult(medium->Is,xr,y1));
      PetscCall(VecSet(y2,0.0));
      for(p=0;p < np;p++){
	PetscCall(MatMultAdd(rsf->ve_C[p],xr,y2,y2));
      }
      PetscCall(VecNorm(y1,NORM_2,&n1));
      PetscCall(VecAXPY(y2,-1.0,y1)); /* -(implicit C_inf) x */
      PetscCall(VecNorm(y2,NORM_2,&n2));
      HEADNODE
	fprintf(stderr,"rsf_ve_setup: worst held-out residual %.2e (tol %.1e); ||C_inf 1||/||K^el 1|| = %.3e (relaxed operator, uniform probe)\n",
		(double)worst_res,(double)fit_tol,(double)(n2/((n1 > 0.0)?n1:1.0)));
      PetscCall(VecDestroy(&xr));PetscCall(VecDestroy(&y1));PetscCall(VecDestroy(&y2));
    }
    if(worst_res > fit_tol)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_CONV_FAILED,
	      "rsf_ve_setup: held-out residual %.2e above -ve_fit_tol %.2e",
	      (double)worst_res,(double)fit_tol);
  }else{
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"rsf_ve_setup: -ve_mode %ld undefined (1, 2 or 3)",(long)mode);
  }
  /* state vectors, Is row layout */
  for(p=0;p < rsf->ve_np;p++){
    PetscCall(MatCreateVecs(medium->Is,NULL,&rsf->ve_h[p]));
    PetscCall(VecSet(rsf->ve_h[p],0.0));
    if(rsf->ve_normal){
      PetscCall(MatCreateVecs(medium->Is,NULL,&rsf->ve_hn[p]));
      PetscCall(VecSet(rsf->ve_hn[p],0.0));
    }
  }
  PetscCall(MatCreateVecs(medium->Is,&rsf->ve_vrel,&rsf->ve_work));
  PetscCall(VecDuplicate(rsf->ve_vrel,&rsf->ve_slip_prev));
  PetscCall(VecSet(rsf->ve_slip_prev,0.0));
  PetscCall(VecDuplicate(rsf->ve_vrel,&rsf->ve_negvpl));
  PetscCall(VecCopy(negvpl_in,rsf->ve_negvpl));
  /*
     initial memory states (-ve_h_init, default 1): 1 = the backslip
     steady state h_p = C_p (-vpl) tau_p, i.e. the medium has been
     loaded by steady plate motion since forever, which removes the
     secular VE spin-up (time scale 3x the slowest ladder pole, about
     3 x 90 t_M for the layered default) from cycle studies; 0 = a
     virgin (unrelaxed) medium with zero memory, appropriate for
     relaxation-from-scratch tests such as the locked-fault
     hereditary-loading check.  NOTE for -ve_mode 1 (uniform Maxwell):
     the spun state has ZERO net fault loading (the memory sink
     exactly cancels the backslip products, the C_inf = 0
     saturation), so steady-state uniform-Maxwell backslip cycles do
     not exist; the events seen with -ve_h_init 0 are the (possibly
     very long lived, ~t_M) transient of the virgin start
  */
  {
    PetscInt hinit = 1;
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-ve_h_init",&hinit,NULL));
    if(hinit == 1){
      for(p=0;p < rsf->ve_np;p++){
	PetscCall(MatMult(rsf->ve_C[p],rsf->ve_negvpl,rsf->ve_h[p]));
	PetscCall(VecScale(rsf->ve_h[p],rsf->ve_tau[p]));
	if(rsf->ve_normal){	/* same steady state, normal family */
	  PetscCall(MatMult(rsf->ve_Cn[p],rsf->ve_negvpl,rsf->ve_hn[p]));
	  PetscCall(VecScale(rsf->ve_hn[p],rsf->ve_tau[p]));
	}
      }
      HEADNODE
	fprintf(stderr,"rsf_ve_setup: memory states initialized at the backslip steady state (-ve_h_init 1)\n");
    }else{
      HEADNODE
	fprintf(stderr,"rsf_ve_setup: memory states initialized to zero (virgin medium, -ve_h_init 0)\n");
    }
  }
  HEADNODE{
    fprintf(stderr,"rsf_ve_setup: VE ON, np %ld, rates ladder [yr]:",(long)rsf->ve_np);
    for(p=0;p < rsf->ve_np;p++)
      fprintf(stderr," %.3g",(double)(rsf->ve_tau[p]/SEC_PER_YEAR));
    fprintf(stderr,"; h: exact-exponential, %s sink forcing, %s\n",
	    (rsf->ve_hstage)?("stage-consistent"):("frozen (lagged)"),
	    (rsf->ve_normal)?("SHEAR + NORMAL tractions"):("shear traction only"));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
static int cmp_comp_precision(const void *a, const void *b)
{
  COMP_PRECISION d = *(const COMP_PRECISION *)a - *(const COMP_PRECISION *)b;
  return (d > 0.0)?(1):((d < 0.0)?(-1):(0));
}
/*
   apply the memory sink to tau_dot at stage time `time` (no-op when VE
   is off); called from rsf_compute_vel_and_stressing with the stage
   state X.

   Default (-ve_h_stage 1), stage-consistent forcing: the sink is the
   exact-exponential h propagated from the last accepted step USING THE
   STAGE STATE'S OWN SLIP INCREMENT as the within-step forcing,

     tau_dot -= sum_p [ h_p(t0) e^(-dt/tau_p) / tau_p
                        + C_p qbar (1 - e^(-dt/tau_p)) ],
     qbar = (slip(t) - slip(t0))/dt - vpl,   dt = time - t0,

   which is exact for velocity constant within the step and, because
   qbar depends on the stage state, is seen by the TS error controller
   (the frozen variant is not, and its lagged forcing accumulates a
   first-order bias during rapidly varying afterslip that can shift
   recurrence by tens of percent in relaxation-critical problems, e.g.
   the uniform-Maxwell patch of Miyake and Noda 2019; the bias was
   found by cross-checking against an exact-Erlang-chain integrator).
   Cost: np matvecs per RHS evaluation.

   -ve_h_stage 0 restores the frozen decay-only sink
   tau_dot -= sum_p h_p(t0) e^(-dt/tau_p)/tau_p (np matvecs per
   ACCEPTED STEP only).
*/
PetscErrorCode rsf_ve_apply_sink(struct interact_ctx *par, Vec X, PetscReal time)
{
  struct med *medium = par->medium;
  struct rsf_vars *rsf = medium->rsf;
  const PetscScalar *x;
  const PetscScalar *sp;
  PetscScalar *vr;
  PetscReal dt,efac;
  PetscInt p,i,j,k;
  PetscFunctionBeginUser;
  if(rsf->ve_np <= 0)
    PetscFunctionReturn(PETSC_SUCCESS);
  dt = time - rsf->ve_t0;
  if(dt < 0.0)dt = 0.0;
  if(rsf->ve_hstage && (dt > 1e-3)){
    /* stage-average relative velocity from the stage state's slip */
    PetscCall(VecGetArrayRead(X,&x));
    PetscCall(VecGetArrayRead(rsf->ve_slip_prev,&sp));
    PetscCall(VecGetArray(rsf->ve_vrel,&vr));
    for(i=medium->rs,j=0,k=0;i < medium->re;i++,j+=rsf->dim,k++)
      vr[k] = (x[j+3] - sp[k])/dt;
    PetscCall(VecRestoreArray(rsf->ve_vrel,&vr));
    PetscCall(VecRestoreArrayRead(rsf->ve_slip_prev,&sp));
    PetscCall(VecRestoreArrayRead(X,&x));
    PetscCall(VecAXPY(rsf->ve_vrel,1.0,rsf->ve_negvpl));
    for(p=0;p < rsf->ve_np;p++){
      efac = PetscExpReal(-dt/rsf->ve_tau[p]);
      PetscCall(VecAXPY(rsf->tau_dot,-efac/rsf->ve_tau[p],rsf->ve_h[p]));
      PetscCall(MatMult(rsf->ve_C[p],rsf->ve_vrel,rsf->ve_work));
      PetscCall(VecAXPY(rsf->tau_dot,-(1.0-efac),rsf->ve_work));
      if(rsf->ve_normal){	/* identical sink on the normal traction */
	PetscCall(VecAXPY(rsf->sigma_dot,-efac/rsf->ve_tau[p],rsf->ve_hn[p]));
	PetscCall(MatMult(rsf->ve_Cn[p],rsf->ve_vrel,rsf->ve_work));
	PetscCall(VecAXPY(rsf->sigma_dot,-(1.0-efac),rsf->ve_work));
      }
    }
  }else{
    for(p=0;p < rsf->ve_np;p++){
      efac = PetscExpReal(-dt/rsf->ve_tau[p]);
      PetscCall(VecAXPY(rsf->tau_dot,-efac/rsf->ve_tau[p],rsf->ve_h[p]));
      if(rsf->ve_normal)
	PetscCall(VecAXPY(rsf->sigma_dot,-efac/rsf->ve_tau[p],rsf->ve_hn[p]));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
/*
   accepted-step h update (TS monitor): exact-exponential recurrence
   with the step-average relative velocity from the slip increment
*/
PetscErrorCode rsf_ve_monitor(TS ts, PetscInt step, PetscReal time, Vec X, void *ptr)
{
  struct interact_ctx *par = (struct interact_ctx *)ptr;
  struct med *medium = par->medium;
  struct rsf_vars *rsf = medium->rsf;
  const PetscScalar *x;
  PetscScalar *vr,*sp;
  PetscReal dt,efac;
  PetscInt i,j,k,p;
  PetscFunctionBeginUser;
  if(rsf->ve_np <= 0)
    PetscFunctionReturn(PETSC_SUCCESS);
  dt = time - rsf->ve_t0;
  if(dt <= 0.0)
    PetscFunctionReturn(PETSC_SUCCESS);
  /* vbar - vpl = (slip - slip_prev)/dt + negvpl */
  PetscCall(VecGetArrayRead(X,&x));
  PetscCall(VecGetArray(rsf->ve_vrel,&vr));
  PetscCall(VecGetArray(rsf->ve_slip_prev,&sp));
  for(i=medium->rs,j=0,k=0;i < medium->re;i++,j+=rsf->dim,k++){
    vr[k] = (x[j+3] - sp[k])/dt;
    sp[k] = x[j+3];
  }
  PetscCall(VecRestoreArray(rsf->ve_slip_prev,&sp));
  PetscCall(VecRestoreArrayRead(X,&x));
  PetscCall(VecRestoreArray(rsf->ve_vrel,&vr));
  PetscCall(VecAXPY(rsf->ve_vrel,1.0,rsf->ve_negvpl));
  for(p=0;p < rsf->ve_np;p++){
    efac = PetscExpReal(-dt/rsf->ve_tau[p]);
    PetscCall(MatMult(rsf->ve_C[p],rsf->ve_vrel,rsf->ve_work));
    PetscCall(VecScale(rsf->ve_h[p],efac));
    PetscCall(VecAXPY(rsf->ve_h[p],rsf->ve_tau[p]*(1.0-efac),rsf->ve_work));
    if(rsf->ve_normal){
      PetscCall(MatMult(rsf->ve_Cn[p],rsf->ve_vrel,rsf->ve_work));
      PetscCall(VecScale(rsf->ve_hn[p],efac));
      PetscCall(VecAXPY(rsf->ve_hn[p],rsf->ve_tau[p]*(1.0-efac),rsf->ve_work));
    }
  }
  rsf->ve_t0 = time;
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
