#include "interact.h"
#ifdef USE_PETSC
#include <petscksp.h>
PetscErrorCode GenEntries(PetscInt , PetscInt , PetscInt ,const PetscInt *, const PetscInt *, PetscScalar *, void *);
#endif

/*

 

  -geom_file geom.in - for the fault geometry
   
  
*/

int main(int argc, char **argv)
{
#ifdef USE_PETSC
  struct med *medium;
  struct flt *fault;
  struct interact_ctx ictx[1];
  KSP               ksp;
  PC                pc;
  Vec         x, b,xout;
  PetscScalar *values=NULL;
  PetscReal   *coords,*avalues=NULL;
  PetscInt    ndim, n, m, lm,ln,i,j,k,dn,on, *col_idx=NULL;
  PetscBool read_value,flg,test_forward=PETSC_TRUE,use_h=PETSC_FALSE;
  MatHtoolKernelFn *kernel = GenEntries;
  VecScatter ctx;
  char geom_file[STRLEN]="geom.in",fault_file[STRLEN]="flt.dat";
  /* generate frameworks */
  medium=(struct med *)calloc(1,sizeof(struct med)); /* make one zero medium structure */
  ictx->medium = medium;
  ictx->src_slip_mode = STRIKE;
  ictx->rec_stress_mode = STRIKE;
  ndim = 3;
  /* 
     start up petsc 
  */
  PetscFunctionBegin;
  /* set defaults, can always override */
  PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_eta","100")); /* not sure  */
  PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_epsilon","1e-3"));
  PetscCall(PetscOptionsSetValue(NULL,"-mat_htool_compressor","SVD"));
  /* start up Petsc proper */
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  /* 
     
     read in geometry

  */
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-use_h", &use_h,&read_value));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_forward", &test_forward,&read_value));
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
  HEADNODE{
    fprintf(stderr,"%s: computing %i by %i matrix\n",argv[0], m,n);
    
  }
  for(i=0;i<m;i++)		/* all sources or receiveer coordinates  */
    for(k=0;k<3;k++)
      coords[i*ndim+k] = fault[i].x[k];

  PetscCall(MatCreate(PETSC_COMM_WORLD, &medium->pA));
  PetscCall(MatSetSizes(medium->pA, PETSC_DECIDE, PETSC_DECIDE, m, n));  
  if(use_h)
    PetscCall(MatSetType(medium->pA,MATHTOOL));
  else
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

  if(use_h){
    fprintf(stderr,"%s: core %03i/%03i: assigning htool row %5i to %5i\n",argv[0],medium->comm_rank,medium->comm_size,medium->rs,medium->re);
    
    PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE, m, n, ndim,(coords+medium->rs), coords, kernel,ictx, &medium->pA));
    
    //PetscCall(MatCreateHtoolFromKernel(PETSC_COMM_WORLD, m, n, m, n, ndim, coords, coords, kernel,ictx, &medium->pA));
  
  }else{
    /* 
       assemble dense matrix 
    */
    fprintf(stderr,"%s: core %03i/%03i: assigning dense row %5i to %5i\n",
	    argv[0],medium->comm_rank,medium->comm_size,medium->rs,medium->re);
    for(j=medium->rs;j <  medium->re;j++){// rupturing faults for this CPU
      GenEntries(ndim,1,n,&j, col_idx, avalues,ictx);
      PetscCall(MatSetValues(medium->pA, 1, &j, n, col_idx,avalues, INSERT_VALUES));
    }
    PetscCall(PetscFree(avalues));
    PetscCall(PetscFree(col_idx));
  }
  PetscCall(MatAssemblyBegin(medium->pA, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(medium->pA, MAT_FINAL_ASSEMBLY));
  PetscCall(MatSetOption(medium->pA, MAT_SYMMETRIC, PETSC_FALSE));
  PetscCall(MatSetFromOptions(medium->pA));

  if(n<20)
    MatView(medium->pA,PETSC_VIEWER_STDOUT_WORLD);
  
  
  PetscCall(MatCreateVecs(medium->pA, &x, &b));/* For A x = b: x -> left, b -> right */
 
  if(test_forward){
    /* 
       test matrix multiplication 
       b = A x
    */
    PetscCall(VecSet(x, 1.0));
    PetscCall(VecAssemblyBegin(x));PetscCall(VecAssemblyEnd(x));
    PetscCall(MatMult(medium->pA, x, b));
    
    if(m<20)
      VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    PetscCall(VecScatterCreateToZero(b,&ctx,&xout));
    PetscCall(VecScatterBegin(ctx,b,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx,b,xout,INSERT_VALUES,SCATTER_FORWARD));
  }else{
    /* 
       test inverse  x = A^-1 b
    */
    /* make constext for solver */
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, medium->pA, medium->pA));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPGetPC(ksp, &pc));
    if(use_h){
      PetscCall(PetscObjectTypeCompare((PetscObject)pc, PCHPDDM, &flg));
    }else{
      PetscCall(PCSetType(pc, PCLU));
    }
    /*  */
    
    PetscCall(VecSet(b, 1.0));
    /* do we need those? */
    PetscCall(VecAssemblyBegin(b));PetscCall(VecAssemblyEnd(b));
   
    PetscCall(KSPSolve(ksp, b, x));
    if(m<20)
      VecView(x,PETSC_VIEWER_STDOUT_WORLD);
    PetscCall(VecScatterCreateToZero(x,&ctx,&xout));
    PetscCall(VecScatterBegin(ctx,x,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx,x,xout,INSERT_VALUES,SCATTER_FORWARD));
  }
  /* assign to x solution vector */
  HEADNODE{
    PetscCall(VecGetArray(xout,&values));
    if(test_forward){
      for(i=0;i<m;i++)
	fault[i].s[STRIKE] = values[i];
    }else{
      for(i=0;i<m;i++)
	fault[i].u[STRIKE] = values[i];
    }
    PetscCall(VecRestoreArray(xout,&values));
    print_fault_data(fault_file,medium,fault);
  }
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(VecDestroy(&xout));  PetscCall(VecScatterDestroy(&ctx));
  free(coords);
  PetscCall(MatDestroy(&medium->pA));
  PetscCall(PetscFinalize());
  
#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);

}

