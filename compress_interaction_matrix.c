#include "interact.h"
#ifdef USE_PETCS

#include <petscksp.h>
static char help[] = "compress an interaction matrix with PCH2OPUS.\n\n";
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
  int iret;
  COMP_PRECISION slip[3],disp[3],stress[3][3],trac[3],sval;

  Vec         x, b, u, d;
  Mat         A, Ad = NULL;
  KSP         ksp;
  PC          pc;
  PetscReal   norm, eta, rtol,*coords;
  PetscInt    basisord, leafsize, ndim, n, i,j, maxrank;
  PetscMPIInt size;
  if(argc<2){
    fprintf(stderr,"%s geom.in\n",argv[0]);
    exit(-1);
  }
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "This is a uniprocessor example only!");

  ndim = 3;
  rtol = 1e-5;
  
  eta = 0.6;
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-eta", &eta, NULL));
  leafsize = 32;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-leafsize", &leafsize, NULL));
  basisord = 8;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-basisord", &basisord, NULL));


  medium=(struct med *)calloc(1,sizeof(struct med));
  read_geometry(argv[1],&medium,&fault,TRUE,FALSE,FALSE,FALSE);
  
  n = medium->nrflt;
  coords = (PetscReal *)malloc(sizeof(PetscReal)*ndim*n);

  
  PetscCall(MatCreateDense(PETSC_COMM_WORLD, n, n, PETSC_DECIDE, PETSC_DECIDE, NULL, &Ad));
  get_right_slip(slip,0,1.0);	/* strike motion */
  fprintf(stderr,"%s computing %i by %i matrix\n",argv[0], medium->nrflt, medium->nrflt);
  /* this is modified from calc interaction to make it simpler */
  for(j=0;j < medium->nrflt;j++){// loop over all rupturing faults
    for(i=0;i<ndim;i++)
      coords[j*ndim+i] = fault[j].x[i];
    for(i=0;i<medium->nrflt;i++){// loop over observing faults
      eval_green(fault[i].x,(fault+j),slip,disp,stress,&iret, GC_STRESS_ONLY,TRUE);
      if(iret != 0){
	fprintf(stderr,"calc_interaction_matrix: WARNING: i=%3i j=%3i singular\n",i,j);
	//s[STRIKE]=s[DIP]=s[NORMAL]=0.0;
	sval = 0.0;
      }else{
	resolve_force(fault[i].normal,stress,trac);
	sval = dotp_3d(trac,fault[i].t_strike);
	//s[STRIKE]=(I_MATRIX_PREC)dotp_3d(trac,fault[i].t_strike);
	//s[DIP]=(I_MATRIX_PREC)dotp_3d(trac,fault[i].t_dip);
	//s[NORMAL]=(I_MATRIX_PREC)dotp_3d(trac,fault[i].normal);
      }
      PetscCall(MatSetValue(Ad, i, j, sval, INSERT_VALUES));
    }
  }
  PetscCall(MatAssemblyBegin(Ad, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Ad, MAT_FINAL_ASSEMBLY));
  //MatView(Ad,PETSC_VIEWER_STDOUT_WORLD);
  
  PetscCall(MatCreateVecs(Ad, &x, &b));
  PetscCall(VecSet(x, 1.0));
  PetscCall(MatMult(Ad, x, b));
  

  
  if(0){
    /* compress */
    //PetscCall(MatCreateH2OpusFromKernel(PETSC_COMM_WORLD, n, n, PETSC_DECIDE, PETSC_DECIDE,ndim, coords, PETSC_FALSE, RBF, &fctx, eta, leafsize, basisord, &A));
    maxrank = 64;
    PetscCall(MatCreateH2OpusFromMat(Ad, ndim, coords, PETSC_TRUE, eta, leafsize, maxrank, basisord, rtol,&A));
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  }

  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  free(coords);
  PetscCall(MatDestroy(&Ad));
  //PetscCall(MatDestroy(&A));
#else
  fprintf(stderr,"%s only petsc version implemented, but not compiled as such\n",argv[0]);
#endif
  exit(0);

}
