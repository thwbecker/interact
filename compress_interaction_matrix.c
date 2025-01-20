#include "interact.h"
#ifdef USE_PETSC

#include <petscksp.h>

struct interact_ctx{
  struct flt *fault;
  struct med *medium;
};

/* 
   generate interaction matrix entries
   sdim dimension
   M local m
   N local n 
   J[M] array with global indices 
   K[N] array with global indices 

 */
static PetscErrorCode GenEntries(PetscInt sdim, PetscInt M, PetscInt N,
				 const PetscInt *J, const PetscInt *K, PetscScalar *ptr, void *kernel_ctx)
{
  PetscInt  j, k;
  struct flt *fault;
  COMP_PRECISION slip[3],disp[3],stress[3][3],trac[3],sval;
  int iret;
  struct interact_ctx *ictx;
  ictx = (struct interact_ctx *)kernel_ctx;
  fault = ictx->fault;
  
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionBeginUser;
#endif
  get_right_slip(slip,0,1.0);	/* strike motion */
  for (j = 0; j < M; j++) {
    for (k = 0; k < N; k++) {
      
      eval_green(fault[K[k]].x,(fault+J[j]),slip,disp,stress,&iret, GC_STRESS_ONLY,TRUE);
      if(iret != 0){
	fprintf(stderr,"get_entries: WARNING: i=%3i j=%3i singular\n",j,k);
	//s[STRIKE]=s[DIP]=s[NORMAL]=0.0;
	sval = 0.0;
      }else{
	resolve_force(fault[K[k]].normal,stress,trac);
	sval = dotp_3d(trac,fault[K[k]].t_strike);
	//s[STRIKE]=(I_MATRIX_PREC)dotp_3d(trac,fault[i].t_strike);
	//s[DIP]=(I_MATRIX_PREC)dotp_3d(trac,fault[i].t_dip);
	//s[NORMAL]=(I_MATRIX_PREC)dotp_3d(trac,fault[i].normal);
      }
      ptr[j + M * k] = sval;
    }
  }
#if !PetscDefined(HAVE_OPENMP)
  PetscFunctionReturn(PETSC_SUCCESS);
#else
  return 0;
#endif
}

#endif
/*
  reads in geometry file and calculates the interaction matrix, and
  the compresses it
  
*/

int main(int argc, char **argv)
{
#ifdef USE_PETSC
  struct med *medium;
  struct flt *fault;
  struct interact_ctx ictx[1];
  double *bglobal;
  Vec         x, b, u, d,bout;
  Mat         Ah;
  KSP         ksp;
  PetscReal   norm, eta, rtol,*coords,*avalues=NULL,*bvalues=NULL;
  PetscInt    basisord, leafsize, ndim, n, m, lm,ln,i,j,k, maxrank,dn,on, *col_idx=NULL;
  VecScatter ctx;
  if(argc<2){
    fprintf(stderr,"%s geom.in\n",argv[0]);
    exit(-1);
  }
  medium=(struct med *)calloc(1,sizeof(struct med));
  ictx->medium = medium;
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));


  ndim = 3;
  rtol = 1e-5;
  
  eta = 0.6;
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-eta", &eta, NULL));
  leafsize = 32;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-leafsize", &leafsize, NULL));
  basisord = 16;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-basisord", &basisord, NULL));

  read_geometry(argv[1],&medium,&fault,TRUE,FALSE,FALSE,FALSE);
  ictx->fault = fault;
  
  m = n = medium->nrflt;
  coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*n);
  bglobal = (double *)malloc(sizeof(double)*m);

  /* parallel version */
  PetscCall(MatCreate(PETSC_COMM_WORLD, &medium->pA));
  PetscCall(MatSetSizes(medium->pA, PETSC_DECIDE, PETSC_DECIDE, m, n));
  PetscCall(MatSetType(medium->pA, MATDENSE));
  PetscCall(MatSetUp(medium->pA));
  PetscCall(MatGetLocalSize(medium->pA, &lm, &ln));
  dn = ln;on = n - ln;
  PetscCall(MatSeqAIJSetPreallocation(medium->pA, n, NULL));
  PetscCall(MatMPIAIJSetPreallocation(medium->pA, dn, NULL, on, NULL));
  PetscCall(MatGetOwnershipRange(medium->pA, &medium->rs, &medium->re));
  /*  */
  medium->rn = medium->re  - medium->rs; /* number of local elements */
  fprintf(stderr,"%s: core %i: dn %i on %i n %i rs %i re %i \n",argv[0],medium->comm_rank,dn,on,n,medium->rs,medium->re);

  /* serial version */
  //PetscCall(MatCreateDense(PETSC_COMM_WORLD, n, n, PETSC_DECIDE, PETSC_DECIDE, NULL, &Ad));


  HEADNODE{
    fprintf(stderr,"%s computing %i by %i matrix\n",argv[0], medium->nrflt, medium->nrflt);
    fprintf(stderr,"%s: core %03i/%03i: assigning row %5i to %5i\n",
	    argv[0],medium->comm_rank,medium->comm_size,medium->rs,medium->re);
  }
  /*  */
  PetscCall(PetscCalloc(m*sizeof(PetscScalar), &avalues));
  PetscCall(PetscCalloc(n*sizeof(PetscInt), &col_idx));
  for (i=0; i < n; i++) 
    col_idx[i] = i;

  for(i=0;i<m;i++)		/* all sources or receiveers  */
    for(k=0;k<3;k++)
      coords[i*ndim+k] = fault[i].x[k];
  
  for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
    GenEntries(ndim,1,n,&j, col_idx, avalues,ictx);
    PetscCall(MatSetValues(medium->pA, 1, &j, n, col_idx,avalues, INSERT_VALUES));
  }
  PetscCall(PetscFree(avalues));
  PetscCall(PetscFree(col_idx));
  PetscCall(MatAssemblyBegin(medium->pA, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(medium->pA, MAT_FINAL_ASSEMBLY));

  MatView(medium->pA,PETSC_VIEWER_STDOUT_WORLD);

  /* hierarhical version */
  //PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD, m, m, m, n, ndim, coords, coords, kernel,ictx, &A));
  
  
  
  PetscCall(MatCreateVecs(medium->pA, &x, &b));/* For A x = b: x -> left, b -> right */
  for (i = medium->rs; i < medium->re; i++) {
    PetscCall(VecSetValue(x, (PetscInt)i, (PetscScalar)(1.0), INSERT_VALUES));
  }
  PetscCall(VecAssemblyBegin(x));
  PetscCall(VecAssemblyEnd(x));
  
  /* solver */
  PetscCall(MatMult(medium->pA, x, b));
  VecView(b,PETSC_VIEWER_STDOUT_WORLD);
  /* distribute */
  PetscCall(VecScatterCreateToZero(b,&ctx,&bout));
  PetscCall(VecScatterBegin(ctx,b,bout,INSERT_VALUES,SCATTER_FORWARD));
  PetscCall(VecScatterEnd(ctx,b,bout,INSERT_VALUES,SCATTER_FORWARD));
  HEADNODE{
    PetscCall(VecGetArray(bout,&bvalues));
    for(i=0;i<m;i++)
      bglobal[i] = bvalues[i];
    PetscCall(VecRestoreArray(bout,&bvalues));
  }
  MPI_Bcast(bglobal,m,MPI_DOUBLE,0, MPI_COMM_WORLD);
  PetscCall(VecScatterDestroy(&ctx));
  PetscCall(VecDestroy(&bout));
  HEADNODE{
    //  for(i=0;i<m;i++)
    //fprintf(stdout,"%20.10e\n",bglobal[i]);
  }
  
  //PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  //PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

 
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  
  free(coords);free(bglobal);
  PetscCall(MatDestroy(&medium->pA));
  //PetscCall(MatDestroy(&A));
  PetscCall(PetscFinalize());
#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);

}

