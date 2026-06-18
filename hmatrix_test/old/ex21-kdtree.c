static char help[] = "Solves a RBF kernel matrix with KSP and PCH2OPUS.\n\n";
/*
 Dave May modified ex21 using KDTREE
 
*/
#include <petscksp.h>

/* ------ */

#define KDTR_MAX_DIM 3

typedef struct _p_kd_node_t* kd_node;
struct _p_kd_node_t {
  double  x[KDTR_MAX_DIM];
  kd_node left,right;
  int     index;
};

typedef struct _p_KDTree *KDTree;
struct _p_KDTree {
  int     npoints,cnt;
  kd_node root,point;
  int     dim;
  int     visited;
  int     setup;
};

/* prototypes */
void kdtr_node_init(kd_node n);
void KDTreeCreate(int dim,KDTree *_k);
void KDTreeDestroy(KDTree *_k);
void KDTreeReset(KDTree kt);
void KDTreeView(KDTree kt);
void KDTreeSetPoints(KDTree k,int np);
void KDTreeGetPoints(KDTree k,int *n,kd_node *nodes);
void KDTreeInsertPoint(KDTree k,double coor[]);
void KDTreeSetup(KDTree kt);
void KDTreeFindNearest(KDTree k,double coor[],kd_node *nearest,double *sep);

/* private prototypes */
static inline double kdtr_dist(kd_node a, kd_node b, int dim);
static inline void   kdtr_swap(kd_node x, kd_node y, int dim);
static kd_node       kdtr_find_median(kd_node start, kd_node end, int idx, int dim);
static kd_node       kdtr_make_tree(kd_node t, int len, int i, int dim);
static void          kdtr_nearest(kd_node root, kd_node nd, int i, int dim, kd_node *best, double *best_dist);

/* global variable, so sue me */
int kdtr_visited;

static inline double kdtr_dist(kd_node a, kd_node b, int dim)
{
  double t, d = 0;
  while (dim--) {
    t = a->x[dim] - b->x[dim];
    d += t * t;
  }
  return d;
}

static inline void kdtr_swap(kd_node x, kd_node y, int dim)
{
  double tmp[] = {0, 0, 0};
  int    tmp_i;
  memcpy(tmp,  x->x, sizeof(double)*dim);
  memcpy(x->x, y->x, sizeof(double)*dim);
  memcpy(y->x, tmp,  sizeof(double)*dim);
  tmp_i = x->index;
  x->index = y->index;
  y->index = tmp_i;
}


/* see quickselect method */
static kd_node kdtr_find_median(kd_node start, kd_node end, int idx, int dim)
{
  kd_node p,store,md ;
  double  pivot;

  if (end <= start) return NULL;
  if (end == start + 1)
  return start;

  md = start + (end - start) / 2;

  while (1) {
    pivot = md->x[idx];

    kdtr_swap(md, end - 1, dim);
    for (store = p = start; p < end; p++) {
      if (p->x[idx] < pivot) {
        if (p != store)
        kdtr_swap(p, store, dim);
        store++;
      }
    }
    kdtr_swap(store, end - 1, dim);

    /* median has duplicate values */
    if (store->x[idx] == md->x[idx])
    return md;

    if (store > md) end = store;
    else            start = store;
  }
}

static kd_node kdtr_make_tree(kd_node t, int len, int i, int dim)
{
  kd_node n;
  
  if (!len) return 0;
  
  if ((n = kdtr_find_median(t, t + len, i, dim))) {
    i = (i + 1) % dim;
    n->left  = kdtr_make_tree(t, n - t, i, dim);
    n->right = kdtr_make_tree(n + 1, t + len - (n + 1), i, dim);
  }
  return n;
}


static void kdtr_nearest(kd_node root, kd_node nd, int i, int dim,
                  kd_node *best, double *best_dist)
{
  double d, dx, dx2;
  
  if (!root) return;
  d = kdtr_dist(root, nd, dim);
  dx = root->x[i] - nd->x[i];
  dx2 = dx * dx;
  
  kdtr_visited ++;
  
  if (!*best || d < *best_dist) {
    *best_dist = d;
    *best = root;
  }
  
  /* if chance of exact match is high */
  if (!*best_dist) return;
  
  if (++i >= dim) i = 0;
  
  kdtr_nearest(dx > 0 ? root->left : root->right, nd, i, dim, best, best_dist);
  if (dx2 >= *best_dist) return;
  kdtr_nearest(dx > 0 ? root->right : root->left, nd, i, dim, best, best_dist);
}


/* Extensions */
void kdtr_node_init(kd_node n)
{
  memset(n->x,0,sizeof(double)*KDTR_MAX_DIM);
  n->left  = NULL;
  n->right = NULL;
  n->index = 0;
}

void KDTreeCreate(int dim,KDTree *_kt)
{
  KDTree kt;

  if (dim > KDTR_MAX_DIM) {
    printf("[kdtree error] KDTree cannot be created. dim must be <= %d\n",KDTR_MAX_DIM);
    *_kt = NULL;
  }
  kt = (KDTree)malloc(sizeof(struct _p_KDTree));
  memset(kt,0,sizeof(struct _p_KDTree));
  kt->root = NULL;
  kt->point = NULL;
  kt->npoints = 0;
  kt->dim = dim;
  kt->visited = 0;
  kt->setup = 0;
  kt->point = (kd_node)calloc(1, sizeof(struct _p_kd_node_t));
  kt->cnt = 0;
  *_kt = kt;
}

void KDTreeDestroy(KDTree *_kt)
{
  KDTree kt;

  if (!*_kt) return;
  kt = *_kt;
  if (!kt) return;
  if (kt->point) {
    free(kt->point);
  }
  kt->point = NULL;
  kt->root = NULL;
  free(kt);
  *_kt = NULL;
}

void KDTreeReset(KDTree kt)
{
  kt->npoints = 0;
  kt->cnt = 0;
  kt->root = NULL;
  kt->setup = 0;
}

void KDTreeView(KDTree kt)
{
  double size;
  printf("KDTree\n");
  printf(" npoints: %d\n",kt->npoints);
  printf(" dim:     %d\n",kt->dim);
  size = ((double)sizeof(struct _p_kd_node_t)) * ((double)kt->npoints);
  printf(" memory:  %1.4e (MB)\n",size/1.0e6);
}

void KDTreeSetPoints(KDTree kt,int np)
{
  int p;
  if (kt->setup == 1) {
    printf("[kdtree error] KDTree already setup. Cannot call KDTreeSetPoints() after KDTreeSetup() has been called.\n");
    return;
  }
  if (np != kt->npoints) {
    kd_node tmp;

    kt->npoints = np;
    tmp = (kd_node)realloc(kt->point,kt->npoints * sizeof(struct _p_kd_node_t));
    kt->point = tmp;
  }
  memset(kt->point,0,kt->npoints * sizeof(struct _p_kd_node_t));
  for (p=0; p<kt->npoints; p++) {
    kt->point[p].index = p;
  }
  kt->cnt = 0;
}

void KDTreeGetPoints(KDTree kt,int *n,kd_node *nodes)
{
  if (n) { *n = kt->npoints; }
  if (nodes) { *nodes = kt->point; }
}

void KDTreeInsertPoint(KDTree kt,double coor[])
{
  if (kt->setup == 1) {
    printf("[kdtree error] KDTree already setup. Cannot call KDTreeInsertPoint() after KDTreeSetup() has been called.\n");
    return;
  }
  if (kt->cnt >= kt->npoints) {
    printf("[kdtree error] Cannot insert into slot %d. Max. index = %d\n",kt->cnt,kt->npoints);
    return;
  }
  memcpy(kt->point[kt->cnt].x, coor, sizeof(double)*kt->dim);
  kt->cnt++;
}

void KDTreeSetup(KDTree kt)
{
  if (kt->setup == 1) return;
  kt->root = NULL;
  // make_tree
  kt->root = kdtr_make_tree(&kt->point[0], kt->npoints, 0, kt->dim);
  kt->setup = 1;
}

void KDTreeFindNearest(KDTree kt,double coor[],kd_node *nearest,double *sep)
{
  struct _p_kd_node_t test_node;
  kd_node             found = NULL;
  double              best_dist = 1.0e32;

  if (kt->setup == 0) {
    printf("[kdtree error] KDTree not setup. Must call KDTreeSetup() before KDTreeFindNearest().\n");
    *nearest = NULL;
    if (sep) { *sep = 1.0e32; }
    return;
  }

  kdtr_node_init(&test_node);
  memcpy(test_node.x, coor, sizeof(double)*kt->dim);

  kdtr_visited = 0;
  kdtr_nearest(kt->root, &test_node, 0, kt->dim, &found, &best_dist);
  kt->visited = kdtr_visited;

  *nearest = found;
  if (sep) { *sep = sqrt(best_dist); }
}


/* ------ */

typedef struct {
  PetscReal  sigma;
  PetscReal *l;
  PetscReal  lambda;
  KDTree     kdtree;
  PetscReal *_coords;
} RBFCtx;

static PetscScalar RBF(PetscInt sdim, PetscReal x[], PetscReal y[], void *ctx)
{
  RBFCtx    *rbfctx = (RBFCtx *)ctx;
  PetscInt   d;
  PetscReal  diff   = 0.0;
  PetscReal  s      = rbfctx->sigma;
  PetscReal *l      = rbfctx->l;
  PetscReal  lambda = rbfctx->lambda;
  kd_node    nearest_x,nearest_y;
  double     target_x[] = {0,0,0},target_y[] = {0,0,0}, sep_x,sep_y;

  for (d=0; d<sdim; d++) {
    target_x[d] = x[d];
    target_y[d] = y[d];
  }
  KDTreeFindNearest(rbfctx->kdtree,target_x,&nearest_x,&sep_x);
  KDTreeFindNearest(rbfctx->kdtree,target_y,&nearest_y,&sep_y);

/*
  if (sdim == 2) {
      printf(">> KDTree \nsearching for x[] (%g, %g)\n"
         "found [%d](%g, %g) dist %g seen %d nodes\n"
	 "--> true coord at index (%g, %g)\n",
         target_x[0], target_x[1],
         nearest_x->index, nearest_x->x[0], nearest_x->x[1], sep_x, kdtr_visited,
	 rbfctx->_coords[sdim*nearest_x->index+0],
	 rbfctx->_coords[sdim*nearest_x->index+1]);

      printf(">> KDTree \nsearching for y[] (%g, %g)\n"
         "found [%d](%g, %g) dist %g seen %d nodes\n"
         "--> true coord at index (%g, %g)\n",
         target_y[0], target_y[1],
         nearest_y->index, nearest_y->x[0], nearest_y->x[1], sep_y, kdtr_visited,
         rbfctx->_coords[sdim*nearest_y->index+0],
         rbfctx->_coords[sdim*nearest_y->index+1]);
  }
*/ 

  for (d = 0; d < sdim; d++) diff += (x[d] - y[d]) * (x[d] - y[d]) / (l[d] * l[d]);
  return s * s * PetscExpReal(-0.5 * diff) + (diff != 0.0 ? 0.0 : lambda);
}

int main(int argc, char **args)
{
  Vec         x, b, u, d;
  Mat         A, Ae = NULL, Ad = NULL;
  KSP         ksp;
  PetscRandom r;
  PC          pc;
  PetscReal   norm, *coords, eta, scale = 0.5;
  PetscInt    basisord, leafsize, sdim, di, n, its, i;
  PetscMPIInt size;
  RBFCtx      fctx;
  KDTree      kd;
  PetscLogDouble t0,t1;


  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &args, NULL, help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "This is a uniprocessor example only!");
  PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &r));
  PetscCall(PetscRandomSetFromOptions(r));

  sdim = 2;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-sdim", &sdim, NULL));
  n = 32;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL));
  eta = 0.6;
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-eta", &eta, NULL));
  leafsize = 32;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-leafsize", &leafsize, NULL));
  basisord = 8;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-basisord", &basisord, NULL));

  /* Create random points */
  PetscCall(PetscMalloc1(sdim * n, &coords));
  PetscCall(PetscRandomGetValuesReal(r, sdim * n, coords));

 
  PetscTime(&t0); 
  KDTreeCreate(sdim,&kd);
  KDTreeSetPoints(kd,n);
  {
    int     p,npoints;
    kd_node nodes;
    
    KDTreeGetPoints(kd,&npoints,&nodes);
    for (p=0; p<npoints; p++) {
      nodes[p].index = p;
      for (di=0; di<sdim; di++) {
        nodes[p].x[di] = (double)coords[sdim*p+di];
      }
    }
  }
  KDTreeSetup(kd);
  PetscTime(&t1);
  PetscPrintf(PETSC_COMM_WORLD,"[kdtree setup] %1.4e sec\n",t1-t0);

  {
    kd_node nearest;
    double  target[] = {0,0,0},sep;
    
    for (di=0; di<sdim; di++) {
      target[di] = coords[sdim*0+di];
    }
    KDTreeFindNearest(kd,target,&nearest,&sep);

    if (sdim == 2) {
      printf(">> KDTree-2 \nsearching for (%g, %g)\n"
         "found [%d](%g, %g) dist %g seen %d nodes\n",
         target[0], target[1],
         nearest->index,nearest->x[0], nearest->x[1], sep, kdtr_visited);
    }
    if (sdim == 3) {
      printf(">> KDTree-3 \nsearching for (%g, %g, %g)\n"
         "found [%d](%g, %g, %g) dist %g seen %d nodes\n",
         target[0], target[1], target[2],
         nearest->index,nearest->x[0], nearest->x[1], nearest->x[2], sep, kdtr_visited);
    }
  } 

  fctx.kdtree = kd;
  fctx._coords = coords;

  fctx.lambda = 0.01;
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-lambda", &fctx.lambda, NULL));
  PetscCall(PetscRandomGetValueReal(r, &fctx.sigma));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-sigma", &fctx.sigma, NULL));
  PetscCall(PetscMalloc1(sdim, &fctx.l));
  PetscCall(PetscRandomGetValuesReal(r, sdim, fctx.l));
  PetscCall(PetscOptionsGetRealArray(NULL, NULL, "-l", fctx.l, (i = sdim, &i), NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-scale", &scale, NULL));

  /* Populate dense matrix for comparisons */
  {
    PetscInt i, j;

    PetscCall(MatCreateDense(PETSC_COMM_WORLD, n, n, PETSC_DECIDE, PETSC_DECIDE, NULL, &Ad));
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) PetscCall(MatSetValue(Ad, i, j, RBF(sdim, coords + i * sdim, coords + j * sdim, &fctx), INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(Ad, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Ad, MAT_FINAL_ASSEMBLY));
  }

  /* Create and assemble the matrix */
  PetscCall(MatCreateH2OpusFromKernel(PETSC_COMM_WORLD, n, n, PETSC_DECIDE, PETSC_DECIDE, sdim, coords, PETSC_FALSE, RBF, &fctx, eta, leafsize, basisord, &A));
  PetscCall(MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE));
  PetscCall(MatSetOption(A, MAT_SYMMETRY_ETERNAL, PETSC_TRUE));
  MatSetFromOptions(A);
  PetscCall(PetscObjectSetName((PetscObject)A, "RBF"));
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatViewFromOptions(A, NULL, "-rbf_view"));

  PetscCall(MatCreateVecs(A, &x, &b));
  PetscCall(VecDuplicate(x, &u));
  PetscCall(VecDuplicate(x, &d));

  {
    PetscReal norm;
    PetscCall(MatComputeOperator(A, MATDENSE, &Ae));
    PetscCall(MatAXPY(Ae, -1.0, Ad, SAME_NONZERO_PATTERN));
    PetscCall(MatGetDiagonal(Ae, d));
    PetscCall(MatViewFromOptions(Ae, NULL, "-A_view"));
    PetscCall(MatViewFromOptions(Ae, NULL, "-D_view"));
    PetscCall(MatNorm(Ae, NORM_FROBENIUS, &norm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Approx err %g\n", norm));
    PetscCall(VecNorm(d, NORM_2, &norm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Approx err (diag) %g\n", norm));
    PetscCall(MatDestroy(&Ae));
  }

  {
    PetscLogDouble delta = 0;
    PetscInt its,maxits=100;

    delta = 0;
    for (its=0; its<maxits; its++) {
      PetscCall(VecSet(u, 1.0));
      PetscTime(&t0);
      PetscCall(MatMult(Ad, u, b));
      PetscTime(&t1); delta += (t1 - t0);
    }
    PetscPrintf(PETSC_COMM_WORLD,"[matmult-dense] x %d -> %1.4e sec\n",maxits,delta);

    delta = 0;
    for (its=0; its<maxits; its++) {
      PetscCall(VecSet(u, 1.0));
      PetscTime(&t0);
      PetscCall(MatMult(A, u, b));
      PetscTime(&t1); delta += (t1 - t0);
    }
    PetscPrintf(PETSC_COMM_WORLD,"[matmult-hmat]  x %d -> %1.4e sec\n",maxits,delta);
  }


  PetscCall(VecSet(u, 1.0));
  PetscCall(MatMult(Ad, u, b));
  PetscCall(MatViewFromOptions(Ad, NULL, "-Ad_view"));
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, Ad, A));
  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, PCH2OPUS));
  PetscCall(KSPSetFromOptions(ksp));
  /* we can also pass the points coordinates
     In this case it is not needed, since the preconditioning
     matrix is of type H2OPUS */
  PetscCall(PCSetCoordinates(pc, sdim, n, coords));

  PetscCall(KSPSolve(ksp, b, x));
  PetscCall(VecAXPY(x, -1.0, u));
  PetscCall(VecNorm(x, NORM_2, &norm));
  PetscCall(KSPGetIterationNumber(ksp, &its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Norm of error %g, Iterations %" PetscInt_FMT "\n", (double)norm, its));

  /* change lambda and reassemble */
  PetscCall(VecSet(x, (scale - 1.) * fctx.lambda));
  PetscCall(MatDiagonalSet(Ad, x, ADD_VALUES));
  fctx.lambda *= scale;
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  {
    PetscReal norm;
    PetscCall(MatComputeOperator(A, MATDENSE, &Ae));
    PetscCall(MatAXPY(Ae, -1.0, Ad, SAME_NONZERO_PATTERN));
    PetscCall(MatGetDiagonal(Ae, d));
    PetscCall(MatViewFromOptions(Ae, NULL, "-A_view"));
    PetscCall(MatViewFromOptions(Ae, NULL, "-D_view"));
    PetscCall(MatNorm(Ae, NORM_FROBENIUS, &norm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Approx err %g\n", norm));
    PetscCall(VecNorm(d, NORM_2, &norm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Approx err (diag) %g\n", norm));
    PetscCall(MatDestroy(&Ae));
  }
  PetscCall(KSPSetOperators(ksp, Ad, A));
  PetscCall(MatMult(Ad, u, b));
  PetscCall(KSPSolve(ksp, b, x));
  PetscCall(MatMult(Ad, x, u));
  PetscCall(VecAXPY(u, -1.0, b));
  PetscCall(VecNorm(u, NORM_2, &norm));
  PetscCall(KSPGetIterationNumber(ksp, &its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Residual norm error %g, Iterations %" PetscInt_FMT "\n", (double)norm, its));

  PetscCall(PetscFree(coords));
  PetscCall(PetscFree(fctx.l));
  PetscCall(PetscRandomDestroy(&r));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&d));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&Ad));
  PetscCall(KSPDestroy(&ksp));
  PetscCall(PetscFinalize());
  return 0;
}

/*TEST

  build:
    requires: h2opus

  test:
    requires: h2opus !single
    suffix: 1
    args: -ksp_error_if_not_converged -pc_h2opus_monitor

  test:
    requires: h2opus !single
    suffix: 1_ns
    output_file: output/ex21_1.out
    args: -ksp_error_if_not_converged -pc_h2opus_monitor -pc_h2opus_hyperorder 2

  test:
    requires: h2opus !single
    suffix: 2
    args: -ksp_error_if_not_converged -pc_h2opus_monitor -pc_h2opus_hyperorder 4

TEST*/
