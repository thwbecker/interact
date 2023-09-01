#include "interact.h"
#include <petscksp.h>
static char petsc_help[] = "Solves a tridiagonal linear system with KSP.\n\n";
int main(int argc,char **args)
{
  int m,j;
  struct med medium[1];

  A_MATRIX_PREC *a,*x,*b;

  Vec            px, pb, pu;      /* approx solution, RHS, exact solution */
  Mat            pA;            /* linear system matrix */
  KSP            ksp;          /* linear solver context */
  PC             pc;           /* preconditioner context */
  PetscReal      norm;         /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       i,n,col[3],its;
  PetscMPIInt    size;
  PetscScalar    value[3];


  
  ierr = PetscInitialize(&argc,&args,(char*)0,petsc_help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
    
  medium->debug = FALSE;

  n = m = 20;
  a=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*n*m);
  b=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*m);
  x=(A_MATRIX_PREC *)malloc(sizeof(A_MATRIX_PREC)*n);

  for(i=0;i<m;i++){
    b[i] = rand();
    for(j=0;j<n;j++)
      a[i*n+j]=rand();
  }
  lu_driver(a,x,b,m,n,medium);
  for(i=0;i<n;i++)
    printf("%g\n",x[i]);


  free(a);free(x);free(b);
  return(0);
}
