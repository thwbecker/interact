#include "interact.h"

/*


  $Id: randomize_list.c,v 1.4 2003/02/13 22:45:12 becker Exp $

  generate a random list of integers list(0..N-1) with 
  that contains all integers from 0 to N - 1

  list will be initialized locally

  if randomize is false, list will be 0, 1, 2 ... N-1
  
*/
void randomize_list(int **list, int n, my_boolean randomize)
{
  struct slist *lists;
  int i;
  long int seed=-1;
  if(n<1){
    fprintf(stderr,"randomize_list: error: called with n: %i\n",n);
    exit(-1);
  }
  *list = (int *)malloc(sizeof(int)*n);
  if(! *list)
    MEMERROR("randomize_list");
  if(!randomize){
    for(i=0;i<n;i++)
      *(*list+i) = i;
    return;
  }else{
    lists = (struct slist *)malloc(sizeof(struct slist)*n);
    /*

      generate random permutation of list by assigning random numbers
      and sorting

    */
    if(!lists)
      MEMERROR("randomize_list");
    for(i=0;i<n;i++){
      lists[i].x = myrand(&seed);
      lists[i].i = i;
    }
    qsort(lists,n,sizeof(struct slist),slist_sort);
    for(i=0;i<n;i++)
      *(*list+i) = lists[i].i;
    free(lists);
  }
}

int slist_sort(const void *va, const void *vb)
{
  struct slist *a = (struct slist *) va, 
    *b = (struct slist *) vb;
  if(a->x < b->x)
    return -1;
  else if(a->x == b->x)
    return 0;
  else
    return 1;
}
