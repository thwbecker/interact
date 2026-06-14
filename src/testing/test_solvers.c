#include "interact.h"


A_MATRIX_PREC mat_value(int ,int , int );

int main(int argc,char **argv)
{
  struct med medium[1];
  my_boolean use_petsc = 0;
  A_MATRIX_PREC *a,*b;
  PetscScalar *x=NULL;
  long int ind;
  
  Vec         px, pb, pr,pxout;
  Mat         pA;
  KSP         pksp;
  PC          ppc;
  PetscInt    i, j, m, n, rs, re;
  PetscMPIInt comm_size, comm_rank;
  PetscScalar *values=NULL;
  PetscInt    *col_idx=NULL;
  PetscReal   norm;
  PetscInt lm, ln, dn, on;
  VecScatter ctx;
  
  for(i=1;i<argc;i++)   
    if(strcmp(argv[i],"-use_petsc")==0)   
      use_petsc =  1;  
       
  PetscFunctionBegin;
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &comm_rank));

  
  medium->debug = FALSE;

  m = 5000;
  n = 5000;

  if(comm_rank == 0)
    x=(PetscScalar *)realloc(x,sizeof(PetscScalar)*m);
  
  if(use_petsc){
    /*  */
    PetscCall(MatCreate(PETSC_COMM_WORLD, &pA));
    PetscCall(MatSetSizes(pA, PETSC_DECIDE, PETSC_DECIDE, m, n));
    PetscCall(MatSetType(pA, MATDENSE));
    PetscCall(MatSetFromOptions(pA));
    PetscCall(MatSetUp(pA));
    
    /* preallocate */
    PetscCall(MatGetLocalSize(pA, &lm, &ln));
    dn = ln;
    on = n - ln;
    PetscCall(MatSeqAIJSetPreallocation(pA, n, NULL));
    PetscCall(MatMPIAIJSetPreallocation(pA, dn, NULL, on, NULL));


    //PetscCall(PetscCalloc(n*sizeof(PetscScalar), &values));
    PetscCall(PetscCalloc(n*sizeof(PetscInt), &col_idx));
    for (j=0; j<n; j++) {
      col_idx[j] = j;
    }

    /* parallel assembly */
    PetscCall(MatGetOwnershipRange(pA, &rs, &re));
    
    for (i=rs; i<re; i++) 
      for (j=0; j<n; j++)
	PetscCall(MatSetValue(pA, i, j, (PetscScalar)mat_value(i,j,m), ADD_VALUES));


    /*
      Always call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY).
      MatAssemblyBegin() is collective and must be called on all ranks.
    */
    PetscCall(MatAssemblyBegin(pA, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(pA, MAT_FINAL_ASSEMBLY));

    /* Shift diagaonl in the hope that this will always render the matrix non-singular */
    PetscCall(MatShift(pA, 1.0e2));
  
    //MatView(pA,PETSC_VIEWER_STDOUT_WORLD);
  
    PetscCall(MatCreateVecs(pA, &pb, &px)); /* For A x = b: x -> left, b -> right */
    PetscCall(VecDuplicate(px, &pr));

    PetscCall(VecSet(pb, 1.0));
  
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &pksp));
    PetscCall(KSPSetOperators(pksp, pA, pA));

    PetscCall(KSPSetType(pksp, KSPPREONLY));

    PetscCall(KSPGetPC(pksp, &ppc));
    PetscCall(PCSetType(ppc, PCLU));

    /* override at run time via -pc_factor_mat_solver_type xxx */
    PetscCall(PCFactorSetMatSolverType(ppc, MATSOLVERPETSC));
    PetscCall(KSPSetFromOptions(pksp));
 
    PetscCall(KSPSolve(pksp, pb, px));
    PetscCall(KSPView(pksp, PETSC_VIEWER_STDERR_WORLD));

    
    /* 
       distribute to zero node
    */
    PetscCall(VecScatterCreateToZero(px,&ctx,&pxout));
    // scatter as many times as you need
    PetscCall(VecScatterBegin(ctx,px,pxout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx,px,pxout,INSERT_VALUES,SCATTER_FORWARD));
    // destroy scatter context and local vector when no longer needed
    PetscCall(VecScatterDestroy(&ctx));
    /* assign to x solution vector */
    if(comm_rank == 0){
      PetscCall(VecGetArray(pxout,&values));
      memcpy(x,values,m*sizeof(PetscScalar));
      PetscCall(VecRestoreArray(pxout,&values));
    }
    PetscCall(PetscFree(col_idx));      
    PetscCall(VecDestroy(&pxout));
    
    if(0){
      /* check residual */
      PetscCall(MatMult(pA, px, pr)); /* r = A x */
      PetscCall(VecAXPY(pr, -1.0, pb)); /* r <- r - b */
      PetscCall(VecNorm(pr, NORM_2, &norm));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Norm of residual %1.10e\n", (double)norm));
    }


    PetscCall(PetscFree(values));
    PetscCall(VecDestroy(&pr));
    PetscCall(VecDestroy(&px));
    PetscCall(VecDestroy(&pb));
    PetscCall(MatDestroy(&pA));
    PetscCall(KSPDestroy(&pksp));

   
  }else{
    if (comm_rank == 0){
      a=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*n*m);
      b=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*m);
      
      for(i=0;i<m;i++){
	b[i] = 1.0;
	for(j=0;j<n;j++){
	  ind = i*n+j;
	  a[ind] = mat_value(j,i,m);
	  if(i==j)
	    a[ind] += 1e2;
	}
      }
      lu_driver(a,x,b,m,n,medium);
      free(a);free(b);
    }
  }
  if(comm_rank == 0)
    for(i=0;i<n;i++)
      printf("%g\n",x[i]);


  if(x)
    free(x);
  PetscCall(PetscFinalize());

  return(0);
}
A_MATRIX_PREC mat_value(int i,int j, int m)
{
  return (A_MATRIX_PREC)(i * m + j + 1.0);
}
