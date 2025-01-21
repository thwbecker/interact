#include "interact.h"
#ifdef USE_PETSC

#include <petscksp.h>

struct interact_ctx{
  struct flt *fault;
  struct med *medium;
};
static PetscErrorCode GenEntries(PetscInt , PetscInt , PetscInt ,
				 const PetscInt *, const PetscInt *, PetscScalar *, void *);

/* 
   generate interaction matrix entries in a way suitable for petsc/htools

   sdim dimension
   M local m
   N local n 
   J[M] array with global indices for sources
   K[N] array with global indices for receivers

   this is modified from the ex82.c petsc example
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
  Vec         x, xh, b, bh, bout,d;
  Mat         Ah;
  PetscReal   *coords,*avalues=NULL,*bvalues=NULL,norm[3];
  PetscInt    ndim, n, m, lm,ln,i,j,k,dn,on, *col_idx=NULL,rs,re;
  VecScatter ctx;
  PetscBool read_value;
  MatHtoolKernelFn *kernel = GenEntries;
  char geom_file[STRLEN]="geom.in";
  /* generate frameworks */
  medium=(struct med *)calloc(1,sizeof(struct med)); /* make one zero medium structure */
  ictx->medium = medium;
  ndim = 3;
  /* 
     start up petsc 
  */
  PetscFunctionBegin;
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  /* 
     
     read in geometry

  */
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  HEADNODE{
    if(read_value)
      fprintf(stderr,"%s: reading geometry from %s as set by -geom_file\n",argv[0],geom_file);
    else
      fprintf(stderr,"%s: reading geometry from default, %s\n",argv[0],geom_file);
  }
  read_geometry(geom_file,&medium,&fault,TRUE,FALSE,FALSE,FALSE);
  
  ictx->fault = fault;
  m = n = medium->nrflt;
  coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*n);
  bglobal = (double *)malloc(sizeof(double)*m);
  HEADNODE{
    fprintf(stderr,"%s: computing %i by %i matrix\n",argv[0], m,n);
    
  }
  for(i=0;i<m;i++)		/* all sources or receiveer coordinates  */
    for(k=0;k<3;k++)
      coords[i*ndim+k] = fault[i].x[k];
 
  /* dense matrix setup */
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
  //fprintf(stderr,"%s: core %i: dn %i on %i n %i rs %i re %i \n",argv[0],medium->comm_rank,dn,on,n,medium->rs,medium->re);
  /*  */
  PetscCall(PetscCalloc(m*sizeof(PetscScalar), &avalues));
  PetscCall(PetscCalloc(n*sizeof(PetscInt), &col_idx));
  for (i=0; i < n; i++) 
    col_idx[i] = i;

  /* assemble dense matrix */
  fprintf(stderr,"%s: core %03i/%03i: assigning dense row %5i to %5i\n",
	  argv[0],medium->comm_rank,medium->comm_size,medium->rs,medium->re);
  for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
    GenEntries(ndim,1,n,&j, col_idx, avalues,ictx);
    PetscCall(MatSetValues(medium->pA, 1, &j, n, col_idx,avalues, INSERT_VALUES));
  }
  PetscCall(PetscFree(avalues));
  PetscCall(PetscFree(col_idx));
  PetscCall(MatAssemblyBegin(medium->pA, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(medium->pA, MAT_FINAL_ASSEMBLY));
  if(n<20)
    MatView(medium->pA,PETSC_VIEWER_STDOUT_WORLD);


  /* hirarchical version */
  PetscCall(MatCreate(PETSC_COMM_WORLD, &Ah));
  PetscCall(MatSetSizes(Ah, PETSC_DECIDE, PETSC_DECIDE, m, n));  
  PetscCall(MatSetType(Ah,MATHTOOL));
  PetscCall(MatSetUp(Ah));
  PetscCall(MatGetLocalSize(Ah, &lm, &ln));
  dn = ln;on = n - ln;
  PetscCall(MatSeqAIJSetPreallocation(Ah, n, NULL));
  PetscCall(MatMPIAIJSetPreallocation(Ah, dn, NULL, on, NULL));
  PetscCall(MatGetOwnershipRange(Ah, &rs, &re));
  
  fprintf(stderr,"%s: core %03i/%03i: assigning htool row %5i to %5i\n",argv[0],medium->comm_rank,medium->comm_size,rs,re);
    
  PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE, m, n, ndim,(coords+rs), coords, kernel,ictx, &Ah));
    
  //PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD, m, n, m, n, ndim, coords, coords, kernel,ictx, &Ah));
  PetscCall(MatSetOption(Ah, MAT_SYMMETRIC, PETSC_FALSE));
  PetscCall(MatSetFromOptions(Ah));
    
  PetscCall(MatAssemblyBegin(Ah, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Ah, MAT_FINAL_ASSEMBLY));
  /* get info on H matrix */
  MatView(Ah,PETSC_VIEWER_STDOUT_WORLD);
  
  
  /* test matrix multiplication */
  PetscCall(MatCreateVecs(medium->pA, &x, &b));/* For A x = b: x -> left, b -> right */
  PetscCall(MatCreateVecs(Ah, &xh, &bh));/* For A x = b: x -> left, b -> right */
  PetscCall(VecSet(x, 1.0));
  PetscCall(VecSet(xh,1.0));
  /* do we need those? */
  PetscCall(VecAssemblyBegin(x));PetscCall(VecAssemblyEnd(x));
  PetscCall(VecAssemblyBegin(xh));PetscCall(VecAssemblyEnd(xh));

 

  /* dense solver */
  PetscCall(MatMult(medium->pA, x, b));
  if(m<20)
    VecView(b,PETSC_VIEWER_STDOUT_WORLD);
  /* H matrix solve */
  PetscCall(MatMult(Ah, xh, bh));
  if(m<20)
    VecView(bh,PETSC_VIEWER_STDOUT_WORLD);

  /* compute difference */
  PetscCall(VecDuplicate(b, &d));PetscCall(VecCopy(b, d));
  
  PetscCall(VecAXPY(d,-1.0,bh));
  PetscCall(VecNorm(b,NORM_2,norm));
  PetscCall(VecNorm(bh,NORM_2,(norm+1)));
  PetscCall(VecNorm(d,NORM_2,(norm+2)));
  fprintf(stdout,"%s: |b| = %20.10e |b_h| = %20.10e |b-b_h| = %20.10e\n",argv[0],norm[0],norm[1],norm[2]);
  
  
  
  /* get b values */
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
  
  
 
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(VecDestroy(&xh));
  PetscCall(VecDestroy(&bh));
  PetscCall(VecDestroy(&d));
  free(coords);free(bglobal);
  PetscCall(MatDestroy(&medium->pA));
  PetscCall(MatDestroy(&Ah));
  PetscCall(PetscFinalize());
#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);

}

