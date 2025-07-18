#ifndef __KDTREE_HEADER_READ

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
void kdtr_node_init(kd_node );
void KDTreeCreate(int ,KDTree *);
void KDTreeDestroy(KDTree *);
void KDTreeReset(KDTree );
void KDTreeView(KDTree );
void KDTreeSetPoints(KDTree ,int );
void KDTreeGetPoints(KDTree ,int *,kd_node *);
void KDTreeInsertPoint(KDTree ,double []);
void KDTreeSetup(KDTree );
void KDTreeFindNearest(KDTree ,double [],kd_node *,double *);



#define __KDTREE_HEADER_READ

#endif
