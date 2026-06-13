#include "../interact.h"
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
   KSP               ksp;
  PC                pc;
  Vec         x, b,xout;
  PetscScalar *values=NULL;
  PetscInt    n, m, i;
  
  PetscBool read_value,flg,test_forward=PETSC_TRUE,use_h=PETSC_FALSE;
  VecScatter ctx;
  char geom_file[STRLEN]="geom.in",fault_file[STRLEN]="flt.dat";
  /* generate frameworks */
  medium=(struct med *)calloc(1,sizeof(struct med)); /* make one zero medium structure */
  
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

  m = n = medium->nrflt;

  HEADNODE{
    fprintf(stderr,"%s: computing %i by %i matrix\n",argv[0], m,n);
    
  }
  
  calc_petsc_Isn_matrices(medium, fault,use_h);
 
  if(n<20){
    MatView(medium->Is,PETSC_VIEWER_STDOUT_WORLD);
    MatView(medium->In,PETSC_VIEWER_STDOUT_WORLD);
  }
  
  PetscCall(MatCreateVecs(medium->Is, &x, &b));/* For A x = b: x -> left, b -> right */
 
  if(test_forward){
    /* 
       test matrix multiplication 
       b = A x
    */
    PetscCall(VecSet(x, 1.0));
    PetscCall(VecAssemblyBegin(x));PetscCall(VecAssemblyEnd(x));
    PetscCall(MatMult(medium->Is, x, b));
    if(m<20){
      VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    }
    PetscCall(VecScatterCreateToZero(b,&ctx,&xout));
    PetscCall(VecScatterBegin(ctx,b,xout,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx,b,xout,INSERT_VALUES,SCATTER_FORWARD));
    
  }else{
    /* 
       test inverse  x = A^-1 b
    */
    /* make constext for solver */
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, medium->Is, medium->Is));
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

  PetscCall(MatDestroy(&medium->Is));
  PetscCall(MatDestroy(&medium->In));
  PetscCall(PetscFinalize());
  
#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);

}
