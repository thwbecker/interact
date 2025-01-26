/*

block model solution solver related routines


$Id: block_solve.c,v 1.6 2011/01/07 07:19:58 becker Exp $


*/
#include "interact.h"
#include "blockinvert.h"

/*

  solve matrix equation K(m,n) x(n) = b(m)

  with bsig(m) uncertainties in b
  
  beta scaling factor are applied for i <= m1 and 
  m1 < i < m, respectively
  
  solvem can be either m1 or m
  
*/

/* if this flag is set, the velocities [m1] and uncertainties [m1]
will be written to ysig.dat, and the unscaled A matrix [m1][n] with n 
number of parameters and m1 number of velocity components will be written
to aunscaled.dat
*/
//#define DEBUG_SOLVER


void solve_block(COMP_PRECISION *kmat, COMP_PRECISION *x,
		 COMP_PRECISION *b, int m1, int m2,
		 int m, int ms,my_boolean solve_for_stress,
		 int n, int nsnf, int nfdamp,int nxdamp,
		 COMP_PRECISION *sigb,COMP_PRECISION beta, 
		 COMP_PRECISION *wmax,COMP_PRECISION **sval,
		 COMP_PRECISION **vmat,COMP_PRECISION *vmod,
		 COMP_PRECISION *sv_cutoff,
		 my_boolean use_numrec_svd, int nsv_zero,
		 my_boolean calc_ginv,int nrbc,
		 my_boolean use_nullspace)
{
  int i,j,k,o1,o2,n2,solvem,mss;
  COMP_PRECISION *kmatloc;
  my_boolean debug;
#ifdef DEBUG_SOLVER
  FILE *out;
  debug = TRUE;
#else
  debug = FALSE;
#endif
  /*  
      assemble model vector 
      vmod_i = v_i/sig_i where the observables
      v_i are the velocities and stresses and sig_i their
      undcertainties
  */
  mss = m1 + m2 + nfdamp;	/* without nxdamp */
  /* 
     velocities 
  */
  solvem = m1;			/* has to be at least m1 */
  for(j=i=0;i < m1;i++,j++) /* we checked earlier for 
			       sigv[]=0 */
    vmod[j] = b[i]/sigb[i]; 
  /* 
     stresses, which we may or may not solve for
  */
  if(solve_for_stress){
    for(i=m1;i < m;i++,j++)
      vmod[j] = beta * b[i]/sigb[i]; 
    solvem += m2;
  }
  /* damping of normal slip motion */
  if(nfdamp){
    for(i=m;i < mss;i++,j++)	
      vmod[j] = b[i]/sigb[i]; 
    solvem += nfdamp;
  }
  if(nxdamp){
    for(i=mss;i < ms;i++,j++)
      vmod[j] = b[i]/sigb[i]; 
    solvem += nxdamp;
  }
  /* do some stupidity tests */
  if(j != solvem){fprintf(stderr,"solve_block: assignment error 1a: j: %i solvem: %i m1: %i m2: %i m: %i ms: %i nsnf: %i nfdamp: %i nxdamp: %i\n",j,solvem,m1,m2,m,ms,nsnf,nfdamp,nxdamp);exit(-1);}
  if(m1 + m2 != m){fprintf(stderr,"solve_block: assignment error 1b: j: %i solvem: %i m1: %i m2: %i m: %i ms: %i nsnf: %i nfdamp: %i nxdamp: %i\n",j,solvem,m1,m2,m,ms,nsnf,nfdamp,nxdamp);exit(-1);}
  if(i > mss + nxdamp){fprintf(stderr,"solve_block: assignment error 1c: j: %i solvem: %i m1: %i m2: %i m: %i ms: %i nsnf: %i nfdamp: %i nxdamp: %i\n",j,solvem,m1,m2,m,ms,nsnf,nfdamp,nxdamp);exit(-1);}
  if(solvem != m1 + ((solve_for_stress)?(m2):(0)) + nfdamp + nxdamp){fprintf(stderr,"solve_block: assignment error 1d: j: %i solvem: %i m1: %i m2: %i m: %i ms: %i nsnf: %i nfdamp: %i nxdamp: %i\n",j,solvem,m1,m2,m,ms,nsnf,nfdamp,nxdamp);exit(-1);}
  if((nxdamp !=0)&&(nxdamp!=n)){fprintf(stderr,"solve_block: assignment error 1d: j: %i solvem: %i m1: %i m2: %i m: %i ms: %i nsnf: %i nfdamp: %i nxdamp: %i\n",j,solvem,m1,m2,m,ms,nsnf,nfdamp,nxdamp);exit(-1);}
  if((nfdamp) && ((ms-m-nxdamp) != nfdamp)){fprintf(stderr,"solve_block: assignment error 1e: j: %i solvem: %i m1: %i m2: %i m: %i ms: %i nsnf: %i nfdamp: %i nxdamp: %i\n",j,solvem,m1,m2,m,ms,nsnf,nfdamp,nxdamp);exit(-1);}
  if((nxdamp) && ((ms - mss) != n)){fprintf(stderr,"solve_block: assignment error 1f: j: %i solvem: %i m1: %i m2: %i m: %i ms: %i nsnf: %i nfdamp: %i nxdamp: %i\n",j,solvem,m1,m2,m,ms,nsnf,nfdamp,nxdamp);exit(-1);}
  /*
    
    save K(solvem,n)  matrix, else overwritten

  */
  my_vecalloc(&kmatloc,n*solvem,"K");
  /*
    
    assign and scale columns of K matrix with uncertainties

  */
  for(i=o1=o2=0;i < n;i++,o1 += solvem,o2 += ms){
    for(k=j=0;j < m1;j++,k++)
      kmatloc[o1+k] = kmat[o2+j] / sigb[j]; /* velocities */
    if(solve_for_stress)
      for(j=m1;j < m;j++,k++)
	kmatloc[o1+k] = beta  * kmat[o2+j] / sigb[j]; /* stresses */
    if(nfdamp)  
      for(j=m;j < mss;j++,k++){
	kmatloc[o1+k] = kmat[o2+j] / sigb[j]; /* slip damping */
      }
    if(nxdamp)
      for(j=mss;j<ms;j++,k++)
	kmatloc[o1+k] = kmat[o2+j] / sigb[j]; /* norm(xsol) damping */
#ifdef DEBUG
    if(k != solvem){
      fprintf(stderr,"solve_block: assignment error 2: k: %i solvem: %i\n",
	      k,solvem);
      exit(-1);
    }
#endif
  }
  if(use_numrec_svd){
    /*
      
    solve K/sig . xtry = v_obs/sig by SVD, solution is stored in xtry
    singular values in sval, and the V matrix in vmat
    
    */
    svd_driver_numrec(kmatloc,x,vmod,solvem,n,sv_cutoff,0,
		      wmax,TRUE,sval,vmat,nsv_zero,
		      calc_ginv,nrbc,use_nullspace,
		      (solve_for_stress)?(FALSE):(TRUE),debug);
  }else{
    if(calc_ginv){
      fprintf(stderr,"solve_block: error: calculation of general inverse not implemented for LAPACK\n");
      exit(-1);
    }
    if(nsv_zero > 0){
      fprintf(stderr,"solve_block: error: nsv_zero != 0 not implemented for LAPACK\n");
      exit(-1);
    }
    svd_driver_lapack(kmatloc,x,vmod,solvem,n,sv_cutoff,0,
		      wmax,sval,nrbc,
		      (solve_for_stress)?(FALSE):(TRUE),debug);
    /* LAPACK version is not set up to return the V matrix yet
       . instead, return all zeroes
    */
    n2 = (n-nrbc) * (n-nrbc);
    my_vecrealloc(vmat,n2,"vmat");
    for(i=0;i < n2;i++)
      *(*vmat+i) = 0.0;
  }
  free(kmatloc);		/* trash kmatloc */
}

/* 
   evaluate the solution given K and x  

*/
void evaluate_block_solution(COMP_PRECISION *kmat,int ms,
			     int n,COMP_PRECISION *x, 
			     COMP_PRECISION *vmod,
			     COMP_PRECISION *vmodc,
			     struct bmd *mod)
{
#ifdef BLOCK_SPHERICAL
  /* 
     obtain the solution for cartesian velocities and stresses 
     by computing vmodc = K * xsol
  */
  calc_Ax_ftn(kmat,ms,n,x,vmodc);
  /* 
     convert the cartesian velocities to spherical 
  */
  convert_cart_sol(vmodc,vmod,mod); 
#else
  /* 
     simply evaluate the solution for cartesian
  */
  calc_Ax_ftn(kmat,ms,n,x,vmod);
#endif

}

