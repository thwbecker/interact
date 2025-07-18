#ifdef USE_PETSC
#include "petsc_prototypes.h"
#ifdef USE_PETSC_HMAT

#include "kdtree.h"



/* private prototypes */
static inline double kdtr_dist(kd_node a, kd_node b, int dim);
static inline void   kdtr_swap(kd_node x, kd_node y, int dim);
static kd_node       kdtr_find_median(kd_node start, kd_node end, int idx, int dim);
static kd_node       kdtr_make_tree(kd_node t, int len, int i, int dim);
static void          kdtr_nearest(kd_node root, kd_node nd, int i, int dim, kd_node *best, double *best_dist);

/* ew */
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


#endif
#endif
