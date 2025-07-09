#include "interact.h"
#include <libgen.h>
#ifdef USE_PETSC_HMAT
#include "petsc_prototypes.h"


#endif

/*
  simple program to test dense and H matrix type interaction matrix
  type computations using Petsc for forward and inverse solves

  there is also compress_interaction_matrix 

  which can be used for speeed testing, and has a different
  implemtation, comparing dense and H matrix

  -geom_file geom.in - for the fault geometry
  -test_forward true/false
  -use_h true/false
  
*/

int main(int argc, char **argv)
{
#ifdef USE_PETSC_HMAT
  struct med *medium;
  struct flt *fault;
  KSP               ksp;
  PC                pc;
  Vec         x, b,xout;
  PetscScalar *values=NULL;
  PetscInt    n, m, i,modes=1;
  Mat Idense;

 
  
  PetscBool read_value,flg,test_forward=PETSC_TRUE,use_h=PETSC_FALSE;
  VecScatter ctx;
  /* input file name */
  char geom_file[STRLEN]="geom.in",fault_file[STRLEN]="flt.dat";
  /* default petsc settings */
  //char* this_source_file = __BASE_FILE__;char* compile_directory = dirname(this_source_file);
  char *home_dir = getenv("HOME");char par_file[STRLEN];
  sprintf(par_file,"%s/progs/src/interact/petsc_settings.yaml",home_dir);
  /* generate frameworks */
  medium=(struct med *)calloc(1,sizeof(struct med)); /* make one zero medium structure */
  
  /* 
     start up petsc 
  */
  PetscFunctionBegin;
  PetscCall(PetscInitialize(&argc, &argv, par_file, NULL));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &medium->comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &medium->comm_rank));
  /* 
     
     options for this code

  */
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-use_h", &use_h,&read_value)); /* H matrix or dense */
  PetscCall(PetscOptionsGetString(NULL, NULL, "-geom_file", geom_file, STRLEN,&read_value)); /* geometry file */
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_forward", &test_forward,&read_value)); /* read in true/false */

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
  /* 
     compute interaction matrices 
  */
  calc_petsc_Isn_matrices(medium, fault,use_h,1.0,0,&medium->Is);
  calc_petsc_Isn_matrices(medium, fault,use_h,1.0,1,&medium->In);
  /* 
     display 
  */
  if((n<20)||(medium->use_h)){
    /*  */
    PetscCall(MatView(medium->Is,PETSC_VIEWER_STDOUT_WORLD));
    if(medium->use_h && (n<20)){
      PetscCall(MatConvert(medium->Is, MATDENSE, MAT_INITIAL_MATRIX, &Idense));
      PetscCall(MatView(Idense,PETSC_VIEWER_STDOUT_WORLD));
      PetscCall(MatDestroy(&Idense));
    }
  }
  
  PetscCall(MatCreateVecs(medium->Is, &b, &x));/* For A x = b: x -> left, b -> right */
 
  if(test_forward){
    HEADNODE{
      fprintf(stderr,"%s: testing forward solve\n",argv[0]);
    }
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
    HEADNODE{
      fprintf(stderr,"%s: testing inverse solve\n",argv[0]);
    }
    /* 
       test INVERSE  x = A^-1 b
    */
    

    
    /* make constext for solver */
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    if(medium->use_h)
      PetscCall(KSPSetOptionsPrefix(ksp,"htool_"));
    PetscCall(KSPSetOperators(ksp, medium->Is, medium->Is));
    if(medium->use_h){
      PetscCall(KSPSetFromOptions(ksp));
      PetscCall(KSPGetPC(ksp, &pc));
      PetscCall(PetscObjectTypeCompare((PetscObject)pc, PCHPDDM, &flg));
    }else{
      PetscCall(KSPGetPC(ksp, &pc));
      PetscCall(PCSetType(pc, PCLU));
      PetscCall(KSPSetFromOptions(ksp)); 
    }

    PetscCall(VecSet(b, 1.0));
    PetscCall(VecAssemblyBegin(b));PetscCall(VecAssemblyEnd(b));
    
    /*  */
    //PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERPETSC));

    
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
  if(modes>1)
    PetscCall(MatDestroy(&medium->In));
  PetscCall(PetscFinalize());
  
#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);

}

