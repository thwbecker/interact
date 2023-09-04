/*
 Timings to assemble a 5000 x 5000 matrix
 
 ${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec -n 4 ./ex_dense -assemble_type 0  -log_view
 
 -assemble_type 0   -->
 MatAssemblyBegin       1 1.0 9.3047e+00
 MatAssemblyEnd         1 1.0 4.0266e+00
 
 -assemble_type 1   -->
 MatAssemblyBegin       1 1.0 3.5948e+00
 MatAssemblyEnd         1 1.0 4.0195e+00

 -assemble_type 2   -->
 MatAssemblyBegin       5 1.0 1.4970e+01
 MatAssemblyEnd         5 1.0 4.0264e+00

 -assemble_type 3   -->
 MatAssemblyBegin       1 1.0 1.8465e-02
 MatAssemblyEnd         1 1.0 3.0000e-05

 -assemble_type 4   -->
 MatAssemblyBegin       1 1.0 4.3400e-04
 MatAssemblyEnd         1 1.0 1.4000e-05
 
 Note:
   -assemble_type 2 does not behave as expected.
   This may be due to the fact I was running on my laptop, i.e.
   not a distributed memory machine.
*/


static char help[] = "Create a dense matrix, solve with LU.\n\n";

#include <petscksp.h>

int main(int argc, char **args)
{
  Vec         x, b, r;
  Mat         A;
  KSP         ksp;
  PC          pc;
  PetscInt    i, j, m, n, rs, re, rank, assemble_type = 0;
  PetscMPIInt comm_size, comm_rank;
  PetscScalar *values=NULL, val=0;
  PetscInt    *col_idx=NULL;
  PetscReal   norm;
  const PetscInt *ranges;

  
  PetscFunctionBegin;
  PetscCall(PetscInitialize(&argc, &args, (char *)0, help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &comm_size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &comm_rank));

  m = 50;
  n = 50;
  
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-assemble_type", &assemble_type, NULL));
  
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m, n));
  PetscCall(MatSetType(A, MATDENSE));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));

  {
    PetscInt lm, ln, dn, on;
    
    PetscCall(MatGetLocalSize(A, &lm, &ln));
    dn = ln;
    on = n - ln;
    PetscCall(MatSeqAIJSetPreallocation(A, n, NULL));
    PetscCall(MatMPIAIJSetPreallocation(A, dn, NULL, on, NULL));
  }

  PetscCall(MatGetOwnershipRange(A, &rs, &re));
  
  PetscCall(PetscCalloc(n*sizeof(PetscScalar), &values));
  PetscCall(PetscCalloc(n*sizeof(PetscInt), &col_idx));
  for (j=0; j<n; j++) {
    col_idx[j] = j;
  }

  /* Assemble matrix */
  switch (assemble_type) {
    case 0: /* sequential */
      /* [rank 0] basic - insert one (i,j) pair at a time */
      PetscPrintf(PETSC_COMM_WORLD,"[assem type 0] Sequential assemble one (i,j) at a time\n");
      if (comm_rank == 0) {
        for (i=0; i<m; i++) {
          for (j=0; j<n; j++) {
            val = i * m + j + 1.0;
            PetscCall(MatSetValue(A, i, j, val, ADD_VALUES));
          }
        }
      }
      break;
  
    case 1: /* sequential */
      /* [rank 0] intermediate - insert one entire row of values at a time */
      PetscPrintf(PETSC_COMM_WORLD,"[assem type 1] Sequential assemble one row at a time\n");

      if (comm_rank == 0) {
        for (i=0; i<m; i++) {
          for (j=0; j<n; j++) {
            values[j] = i * m + j + 1.0;
          }
          PetscCall(MatSetValues(A, 1, &i, n, col_idx, values, INSERT_VALUES));
        }
      }
      break;

    case 2: /* sequential */
      /* [rank 0] advanced - insert one entire row of values at a time, flush periodically to reduce memory usage on rank 0 */
      PetscPrintf(PETSC_COMM_WORLD,"[assem type 2] Sequential assemble one row at a time, flush chunks\n");

      PetscCall(MatGetOwnershipRanges(A, &ranges));
      for (rank=0; rank<comm_size; rank++) {
        if (comm_rank == 0) {
          PetscInt rs_r = ranges[rank];   /* start row index associated with rank r */
          PetscInt re_r = ranges[rank+1]; /* start row index associated with rank r */

          for (i=rs_r; i<re_r; i++) {
            for (j=0; j<n; j++) {
              values[j] = i * m + j + 1.0;
            }
            PetscCall(MatSetValues(A, 1, &i, n, col_idx, values, INSERT_VALUES));
          }
        }
        /*
         This is collective so much be called by all ranks.
         This call will flush any off-rank matrix, e.g. they will
         be migrated to the correct MPI rank. It is better to do this
         periodically when assemblying from rank 0, otherwise you have to
         store the entire matrix in memory on rank 0 before the call to
         MatAssemblyBegin(MAT_FINAL_ASSEMBLY) is executed.
        */
        PetscCall(MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY));
        PetscCall(MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY));
      }
      break;

    case 3: /* parallel */
      /* Assemble matrix one (i,j) at a time - slow */
      PetscPrintf(PETSC_COMM_WORLD,"[assem type 3] Parallel assemble one (i,j) at a time\n");
      for (i=rs; i<re; i++) {
        for (j=0; j<n; j++) {
          val = i * m + j + 1.0;
          PetscCall(MatSetValue(A, i, j, val, INSERT_VALUES));
        }
      }
      break;

    case 4: /* parallel */
      /* Assemble matrix one row at a time */
      PetscPrintf(PETSC_COMM_WORLD,"[assem type 4] Parallel assemble one row at a time\n");
      for (i=rs; i<re; i++) {
        for (j=0; j<n; j++) {
          values[j] = i * m + j + 1.0;
        }
        PetscCall(MatSetValues(A, 1, &i, n, col_idx, values, INSERT_VALUES));
      }
      break;

    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "assem type unknown");
      break;
  }
  
  /*
   Always call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY).
   MatAssemblyBegin() is collective and must be called on all ranks.
  */
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

  /* Shift diagaonl in the hope that this will always render the matrix non-singular */
  PetscCall(MatShift(A, 1.0e2));
  
  //MatView(A,PETSC_VIEWER_STDOUT_WORLD);
  
  PetscCall(MatCreateVecs(A, &b, &x)); /* For A x = b: x -> left, b -> right */
  PetscCall(VecDuplicate(x, &r));

  PetscCall(VecSet(b, 1.0));
  
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));

  PetscCall(KSPSetType(ksp, KSPPREONLY));

  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, PCLU));

  /* override at run time via -pc_factor_mat_solver_type xxx */
  PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERPETSC));
  //PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK));

  /* #ifdef PETSC_HAVE_ELEMENTAL */
  /* PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERELEMENTAL)); */
  /* #endif */

  /* #ifdef PETSC_HAVE_MUMPS */
  /* PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS)); */
  /* #endif */

  /* #ifdef PETSC_HAVE_SUPERLU_DIST */
  /* PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST)); */
  /* #endif */

  
  PetscCall(KSPSetFromOptions(ksp));

  PetscCall(KSPSolve(ksp, b, x));

  PetscCall(KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD));

  /* check residual */
  PetscCall(MatMult(A, x, r)); /* r = A x */
  PetscCall(VecAXPY(r, -1.0, b)); /* r <- r - b */
  PetscCall(VecNorm(r, NORM_2, &norm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Norm of residual %1.10e\n", (double)norm));

  PetscCall(PetscFree(col_idx));
  PetscCall(PetscFree(values));
  PetscCall(VecDestroy(&r));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(KSPDestroy(&ksp));

  PetscCall(PetscFinalize());
  return 0;
}

