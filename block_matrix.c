/*

assemble matrices for the block inversion routine

all matrices are stored FORTRAN style



$Id: block_matrix.c,v 1.16 2004/10/05 01:09:46 becker Exp $

*/

#include "interact.h"
#include "blockinvert.h"

/*


  assemble the matrices which change with the fault coupling 
  factors: F, GF, E, K2, and K

  some input: G, D, I

  the funny memory handling was added to track down a memory leak problem
  

*/

void assemble_block_fltdep_matrices(COMP_PRECISION **f,
				    COMP_PRECISION **gf,
				    COMP_PRECISION **e,
				    COMP_PRECISION **k2mat,
				    COMP_PRECISION *g,
				    COMP_PRECISION *d,
				    COMP_PRECISION *imat,
				    int m,int n,int nflt,int nrb,
				    int nrgp,int nrsp,int nsnf,
				    COMP_PRECISION *xtry,
				    struct bflt *fault,
				    struct bck *block,
				    my_boolean invert_for_ld,
				    my_boolean invert_for_cfac)
{
  int pdim,tsize;
  static int gf_size=0;
  /*
    
    assemble F which relates the block motion parameters
    to differences in fault location velocities, 
    ie. global slip. dimensions: (dim*nflt, nbase*nrb)
    F uses the globally-projected fault coordinates
    
    we are passing the locking factor vector, it is a [nflt] and has
    potentially different coupling factors for each fault, normally
    should be all unity. if we are inverting for both cfacs and lds,
    the vector also holds the locking depths, hence we have to shift
    this vector gets only reference if invert_for_cfac is set

  */
  assemble_block_f(f,fault,nflt,nrb,block,
		   (xtry+n+((invert_for_ld)?(nflt):(0))),
		   invert_for_cfac);
#ifdef BLOCK_SPHERICAL
  pdim = 3;
#else
  pdim = BLOCK_DIM;
#endif
  //
  // 
  // calculate GF = G(nslip*nflt,dim*nflt) * 
  //                            F(block_dim*nflt,nbase*nrb)
  //
  // dimensions: (nslip*nflt, nbase*nblock)
  //
  /* allocate space for GF  */
  tsize = nsnf *n;
  if(gf_size != tsize){
    my_vecrealloc(gf,tsize,"assemble_block_fltdep_matrices: gf");
    gf_size = tsize;
  }
  calc_AB_ftn(g,nsnf,BLOCK_DIM*nflt,*f,n,*gf);
  //
  // calculate E=D(pdim*nrgp,nslip*nflt).GF(nslip*nflt, nbase*nrb)
  // dimensions: E(pdim*nrgp, n = nbase*nrb)
  //
  /* allocate space for E */
  tsize = pdim * nrgp *n;
  my_vecrealloc(e,tsize,"assemble_block_fltdep_matrices: e");
  calc_AB_ftn(d,pdim * nrgp, nsnf, *gf, n, *e);
  if(nrsp){
    // 
    // assemble K2(nrsp*6,nbase*nrb)  = I G F
    //
    /* allocate space for K2 */
    tsize = nrsp * 6 * n;
    my_vecrealloc(k2mat,tsize,"assemble_block_fltdep_matrices: k2");
    calc_AB_ftn(imat,6*nrsp, nsnf, *gf, n, *k2mat);
  }
}

/*

  assemble A matrix that links block motion parameters, b, to 
  velocities at the nrp velocity observational points 
  b is NBASE * nrb (nrb: nr of blocks)

  A is DIM * nrp x NBASE * nrb
 
  DIM is either BLOCK_DIM or 3 (for spherical)

  gpx[nrp*BLOCK_DIM] are the projected 
               locations of the observational points

  bcode[nrp*BLOCK_DIM] the block code of each point

  block[nrb].fixed the no-rigid motion flag

  assemble A matrix in FORTRAN style (for solvers)

  A is only dependent on the block geometry, not the faults

  if single_block is >= 0, will assume that all points belong
  to single_block.
  
*/
void assemble_block_a(COMP_PRECISION **a,COMP_PRECISION *gpx,
		      int *bcode, int nrp,int nrb,struct bck *block,
		      int single_block,struct bmd *mod)
{
  int i,j,k,icol,n,m,offset,*loc_bcode,inode,pdim,irow;
  /* 
     normally, we will use the fixed flag for blocks whose rotation
     is set to zero. however, if we are only looking for one block's
     rotation, we will set it to FALSE (else all zero entries in A)
  */
  my_boolean use_fixed_block_flag = TRUE;
#ifdef SUPER_DEBUG
  static my_boolean pinit = FALSE;
#endif
  if((BLOCK_NBASE != 3)||(BLOCK_DIM != 2)){
    fprintf(stderr,"assemble_block_a: dim: %i or nbase: %i not implemented\n",
	    BLOCK_DIM,BLOCK_NBASE);
    exit(-1);
  }
  /* 

  number of rows of A  

  */
#ifdef BLOCK_SPHERICAL
  pdim = 3;			/* x,y,z cartesian components */
#else
  pdim = BLOCK_DIM;		/* only x and y in projected frame */
#endif
  m = pdim * nrp;
  /* 
     generate a copy of the boundary code vector, but make it only nrp
     long. (bcode is nrp * BLOCK_DIM)
  */
  my_ivecalloc(&loc_bcode,nrp,"assmble_block_a: loc_bcode");
  if(single_block >= 0){
    if(block[single_block].fixed){
      fprintf(stderr,"assmble_block_a: WARNING: assuming all points belong to block %i\n",
	      single_block+1);
      fprintf(stderr,"assmble_block_a: WARNING: disregarding block fixedness for all blocks (nbase: %i)\n",
	      BLOCK_NBASE);
      use_fixed_block_flag = FALSE;
    }
    for(i=0;i < nrp;i++)	/* all points are in single_block */
      loc_bcode[i] = single_block;
  }else{
    /* original */
    for(i=0;i < nrp;i++)
      loc_bcode[i] = bcode[i*BLOCK_DIM];
  }
  /* 
     number of columns of A 
  */
  n = BLOCK_NBASE * nrb;	
  /* 
     re-allocate A array
  */
  my_vecrealloc(a,n*m,"assemble_block_a: a");
  if((BLOCK_DIM != 2)||(BLOCK_NBASE != 3)){
    fprintf(stderr,"assemble_block_a: not prepared for dim != 2 or block_base != 3\n");
    exit(-1);
  }
#ifdef BLOCK_SPHERICAL
  if(!mod->pbase_init)		/* init basis vectors */
    init_gps_pbase(mod);
#endif
  for(k=0;k < nrb;k++){		/* block loop */
    for(j=0;j < BLOCK_NBASE;j++){
      icol = k * BLOCK_NBASE + j;		/* column counter 
						   0 ... n-1 */
      offset = icol * m;
      for(i=inode=irow=0;i < m;i++){/* 
				  
				    row counter, loops through 
				    all data points and pdim dimensions 
			       
				    */
	if(loc_bcode[inode] == k){      
	  /* 
	     point i is in block k
	  */
	  /*
	    
	    rigid body translation and rotation part
	    
	  */
	  if((!block[loc_bcode[inode]].fixed)||
	     (!use_fixed_block_flag)){ /* this block is not fixed with
					  respect to rigid body
					  motions and we are using the
					  fixed block flags
				       */

#ifdef BLOCK_SPHERICAL
	    /* 

	    spherical setup:
	    
	           |    0  r_z  -r_y  |
	    A' w = | -r_z    0   r_x  | w = v
	           |  r_y -r_x     0  |

	    where r_x, r_y, r_z are the cartesian coordinates and the
	    solution vector will be omega_x, omega_y, omega_z

	    NOTE THAT THIS ASSIGNMENT WORKS COLUMN BY COLUMN!
	    
	    */
	    if(j == 0){	/* first sub-col */
	      if(irow == 0)	                        /*     0 */
		*(*a+offset+i)= 0.0; 
	      else if(irow == 1) 
		*(*a+offset+i)=  -mod->gcx[inode*3+INT_Z];  /*  -r_z */
	      else
		*(*a+offset+i)=   mod->gcx[inode*3+INT_Y];  /*   r_y */

	    }else if(j == 1){	/* second sub-col */
	      if(irow == 0)
		*(*a+offset+i)=   mod->gcx[inode*3+INT_Z];	/*   r_z */
	      else if(irow == 1)
		*(*a+offset+i)=  0.0;                   /*     0 */
	      else
		*(*a+offset+i)=  -mod->gcx[inode*3+INT_X];   /* -r_x */
	    }else{  /* third sub-col  */
	      if(irow == 0)
		*(*a+offset+i)= -mod->gcx[inode*3+INT_Y];   /* -r_y */
	      else if(irow == 1)
		*(*a+offset+i)=  mod->gcx[inode*3+INT_X];   /*  r_x */
	      else
		*(*a+offset+i)=  0.0;                   /*  0  */
	    }
#else
	    /* 

	    cartesian setup:
	    
	    | -x_y    1    0   |
	    |  x_x    0    1   |

	    where x_x, x_y are the mercator projected coordinates  
	    and the solution vector will be omega, v_x, v_y

	    */
	    if(j == 0){	/* first sub-col: omega */
	      if(irow==0)
		*(*a+offset+i)= -gpx[inode*BLOCK_DIM+INT_Y]; /* - x_data_y */
	      else 
		*(*a+offset+i)=  gpx[inode*BLOCK_DIM+INT_X]; /*   x_data_x */
	    }else if(j == 1){	/* second sub-col: vx_0 */
	      if(irow==0)
		*(*a+offset+i)=1.0;
	      else 
		*(*a+offset+i)=0.0;
	    }else{  /* third sub-col: vy_0  */
	      if(irow==0)
		*(*a+offset+i)=0.0;
	      else 
		*(*a+offset+i)=1.0;
	    }
#endif
	  }else{
	    /*
	      block is fixed with respect to rigid body motions, SVD
	      should find the minimum norm solution that corresponds
	      to vx_0,vy_0,omega = 0 for this point
	    */
	    *(*a+offset+i)=0.0;
	  }
	}else{			
	  /* 
	     point is not in block 
	  */
	  *(*a+offset+i)=0.0;
	}
	irow++;/* row counter */
	if(irow == pdim){	/* increment node counter */
	  irow=0;
	  inode++;
	}
      }
    }
  }
  free(loc_bcode);
#ifdef SUPER_DEBUG
  if(!pinit){
    print_matrix_ftrn_file(*a,m,n,"a.dat",FALSE);
    pinit=TRUE;
  }
#endif
}
/* 

   assemble the D matrix that links fault-local slip, 
   \hat{vec{u}}, and global displacement rates 
   (velocities), \vec{v}

   C is the coseismic displacement part as in 

   v = R - C

   D has dimensions DIM * nrp  times nr_of_slip_dir * nrf

   C = \vec{v} = D \hat\vec{u}

   DIM is either BLOCK_DIM or 3 (for spherical) 

   D will be assembled in FORTRAN convention

*/
void assemble_block_d(COMP_PRECISION **d,struct bflt *fault, 
		      int nflt,int nrp,int nslip,
		      struct bmd *mod)
{
  int i,j,k,m,n,irow,icol,offset,idir,tsize;
  static int d_size=0;
#ifdef SUPER_DEBUG
  static my_boolean pinit=FALSE;
#endif
#ifdef BLOCK_SPHERICAL
  COMP_PRECISION xp[3];
  m = 3 * nrp;
  if(!mod->cvel_init)		/* initialize cartesian GPS velocities */
    init_cart_gps_velsig(mod,FALSE);
#else
  int l;
  m = BLOCK_DIM * nrp;
#endif
  n = nslip * nflt;		/* nslip is nr of types of slip  */
  /* allocate space for D */
  tsize = n * m;
  if(d_size != tsize){
    my_vecrealloc(d,tsize,"assemble_block_d: d");
    d_size = tsize;
  }
#ifdef BLOCK_SPHERICAL
  xp[INT_R] = 0.0;
#endif
  for(i=0;i < nflt;i++){		/* i: fault loop */
    if((!fault[i].vertical)&&(nslip < 2)){
      fprintf(stderr,"assemble_block_d: fault %i has dip %g but nslip: %i\n",
	      i+1,fault[i].dip,nslip);
      exit(-1);
    }
    for(j=0;j < nslip;j++){		
      /* j: fault slip direction loop */
      /* fault slip directions have been sorted: 
	 strike, dip, normal in read_bflt 
	 if j=0: strike direction
	    j=1: normal or dip direction, depending on the fault
       	         being vertical or not

      */
      icol = nslip * i + j;
      offset = icol * m;
      idir = block_slip_direction(j,(fault+i));
      for(k=0;k < nrp;k++){	/* k: loop through observational points */
#ifdef BLOCK_SPHERICAL
	// convert polar, rotated displacements to cartesian
	xp[INT_THETA]= -(fault+i)->v[k].vc[idir][INT_Y];
	xp[INT_PHI] =   (fault+i)->v[k].vc[idir][INT_X];
	irow = k * 3;
	pv2cv(xp,(*d + offset + irow),(mod->pbase + k*9));
#else
	for(l=0;l < BLOCK_DIM;l++){	/* l: global directions  */
	  irow = k * BLOCK_DIM + l;
	  *(*d + offset + irow) = (fault+i)->v[k].vc[idir][l];
	}
#endif
      }
    }
  }
#ifdef SUPER_DEBUG
  if(!pinit){
    print_matrix_ftrn_file(*d,m,n,"d.dat",FALSE);
    pinit=TRUE;
  }
#endif
}
/* 

   assemble the I matrix that links fault-local slip, \hat{vec{u}},
   and global stresses

   I has dimensions 6 * nrp x nr_of_slip_dir * nrf

   I will be assembled in FORTRAN convention

*/
void assemble_block_i(COMP_PRECISION **imat,struct bflt *fault, 
		      int nflt,int nrp,int nslip)
{
  int i,j,k,l,m,n,irow,icol,offset,idir,tsize;
  static int i_size=0;
#ifdef DEBUG
  static my_boolean pinit=FALSE;
#endif
  m = 6 * nrp;
  n = nslip * nflt;		/* nslip is nr of types of slip  */
  /* 
     realloc space for I
  */
  tsize = n * m;
  if(i_size != tsize){
    my_vecrealloc(imat,tsize,"assemble_block_i: I");
    i_size = tsize;
  }
  for(i=0;i < nflt;i++){		/* i: fault loop */
    if((!fault[i].vertical)&&(nslip < 2)){
      fprintf(stderr,"assemble_block_i: fault %i has dip %g but nslip: %i\n",
	      i+1,fault[i].dip,nslip);
      exit(-1);
    }
    for(j=0;j < nslip;j++){		/* j: fault slip direction loop */
      icol = nslip * i + j;
      offset = icol * m;
      idir = block_slip_direction(j,(fault+i));
      for(k=irow=0;k < nrp;k++,irow += 6){/* k: loop through observational points */
	for(l=0;l < 6;l++)	/* l: stress components  */
	  *(*imat + offset + irow + l) = 
	    (fault+i)->s[k].sc[idir][l];
      }
    }
  }
#ifdef DEBUG
  if(!pinit){
    print_matrix_ftrn_file(*imat,m,n,"i.dat",FALSE);
    pinit=TRUE;
  }
#endif
}
/* 

   assemble the G matrix that links global slip, \vec{u}, in a
   cartesian lon/lat system to fault-local strike, normal (and dip)
   components,\hat{\vec{u}}

   G has dimensions 
   nslip  * nrf  times BLOCK_DIM * nrf
   
   \hat{u} = G . u

   G will be assembled in FORTRAN convention

   G depends on the fault geometry orientation, but not the 
   locking depth

*/
void assemble_block_g(COMP_PRECISION **g,struct bflt *fault, 
		      int nflt, int nslip)
{
  int i,j,k,l,m,n,idir,irow,icol,offset,tsize;
  static int g_size=0;
#ifdef SUPER_DEBUG
  static my_boolean pinit=FALSE;
#endif
  m = nslip * nflt;	    /* nslip is nr of types of slip  */
  n = BLOCK_DIM * nflt;
  if(BLOCK_DIM != 2){
    fprintf(stderr,"assemble_block_g should be adpted for DIM != 2\n");
    exit(-1);
  }
  /* 
     
  realloc space for G
  
  */
  tsize = n * m;
  if(g_size != tsize){
    my_vecrealloc(g,tsize,"assemble_block_g: g");
    g_size = tsize;
  }
  for(i=0;i < nflt;i++){		/* i: fault loop */
    for(j=0;j < BLOCK_DIM;j++){	/* j: global direction loop */
      icol = BLOCK_DIM * i + j;
      offset = icol * m;
      for(k=0;k < nflt;k++)	/* k: second (row) fault loop */
	for(l=0;l < nslip;l++){/* l: slip directions  */
	  irow = k * nslip + l;
	  if(i == k){		/* 
				   this matrix corresponds to t
				   the unity vectors in 
				   strike and normal/dip 
				   direction and is zero else 
				*/
	    idir = block_slip_direction(l,(fault+i));
	    *(*g + offset + irow) = (fault+i)->evec[idir*3+j];	
	  }else{
	    *(*g + offset + irow) = 0.0;
	  }
	}
    }
  }
#ifdef SUPER_DEBUG
  if(!pinit){
    print_matrix_ftrn_file(*g,m,n,"g.dat",FALSE);
    pinit++;
  }
#endif
}
/*

  assemble F matrix that links block motion parameters, b, to 
  global velocity differences (v_left - v_right) at the 
  fault midpoints, interpreted as slip, in a global, lon/lat
  cartesian system

  b is NBASE * nrb (nrb: nr of blocks)

  F is DIM * nflt x NBASE * nrb
 
  assemble F matrix in FORTRAN style (for solvers)


  block[nrb].fixed is TRUE if the block is fixed with respect to 
  rigid body motions

  cfac[nflt] is a coupling factor vector (different from the 
  constant factor we have for each fault during input), that should
  normally be initialized to all unity

  cfac only gets referenced if invert_for_cfac is set

  for a spherical calculation, we determine F from 

  B . S 

  where S links the rotation vectors to cartesian velocities
  and B rotates those into the lon/lat system


*/
void assemble_block_f(COMP_PRECISION **f,struct bflt *fault,
		      int nflt,int nrb, struct bck *block,
		      COMP_PRECISION *cfac,
		      my_boolean invert_for_cfac)
{
  int i,j,k,l,icol,irow,n,m,offset,pdim,mf,tsize;
  COMP_PRECISION fac,loc_cfac;
#ifdef BLOCK_SPHERICAL
  COMP_PRECISION *b, *s;
#endif
#ifdef SUPER_DEBUG
  static my_boolean written=FALSE;
#endif
#ifdef DEBUG
  static my_boolean warned=FALSE;
  if(invert_for_cfac)
    if(!warned){
      for(i=0;i<nflt;i++){
	if(fabs(cfac[i]-1.0) > EPS_COMP_PREC){
	  fprintf(stderr,"assemble_block_f: WARNING: at least one coupling factor != 1.0\n");
	  warned = TRUE;
	  break;
	}
      }
    }
#endif
  if((BLOCK_NBASE != 3)||(BLOCK_DIM!=2)){
    fprintf(stderr,"assemble_block_f: dim: %i or nbase: %i not implemented\n",
	    BLOCK_DIM,BLOCK_NBASE);
    exit(-1);
  }
  /* 
     nr of columns 
  */
  n = BLOCK_NBASE * nrb;
#ifdef BLOCK_SPHERICAL
  pdim = 3;  
  /* 
     nr of rows of S 
  */
  m = pdim * nflt;
  /* nr of rows of F  */
  mf = BLOCK_DIM * nflt;
  /* 
     allocate space for S
  */
  my_vecalloc(&s,n*m,"assemble_block_f: s");
#else
  pdim = BLOCK_DIM;
  m = BLOCK_DIM * nflt;
  mf = m;
  tsize = mf * n;
  /* allocate space for F */
  my_vecrealloc(f,tsize,"assemble_block_f: f");
#endif
  for(k=0;k < nrb;k++){		/* loop through blocks */
    for(j=0;j < BLOCK_NBASE;j++){	/* loop through base functions */
      icol = k * BLOCK_NBASE + j;	/* column counter */
      offset = icol * m;
      for(i=0;i < nflt;i++){		/* loop through faults */
	loc_cfac = (invert_for_cfac)?(cfac[i]):(1.0);
	if(block[k].fixed)		/* block is fixed */
	  fac = 0.0;
	else{			/* block has rigid body motion */
	  if(fault[i].block[0] == k) /* left border of fault 
					in block*/
	    fac =  loc_cfac;	/* normally unity */
	  else if(fault[i].block[1] == k) /* right border */
	    fac = -loc_cfac;	/* normally -unity */
	  else
	    fac = 0.0;		/* zero else */
	}
	for(l=0;l < pdim;l++){ /* loop through global directions */
	  irow = pdim * i + l;
	  if(fac != 0.0){
#ifdef BLOCK_SPHERICAL
	    /* 

	       assemble S 

	       |    0  r_z  -r_y |
	       | -r_z    0   r_x |  w = v
	       |  r_y -r_x     0 |
	       
	       NOTE THAT THIS ASSIGMENT WORKS COLUMN BY COLUMN
	       
	    */
	    if(j == 0){	/* first sub-col */
	      if(l == 0)
		*(s+offset+irow)=   0.0;                  /*      0 */
	      else if(l == 1) 
		*(s+offset+irow)=  -fault[i].xc[INT_Z] * fac; /*   -r_z */
	      else
		*(s+offset+irow)=   fault[i].xc[INT_Y] * fac; /*    r_y */
	    }else if(j == 1){	/* second sub-col */
	      if(l == 0)
		*(s+offset+irow)=   fault[i].xc[INT_Z] * fac; /* r_z */
	      else if(l == 1)
		*(s+offset+irow)=   0.0;                  /*   0 */
	      else
		*(s+offset+irow)=  -fault[i].xc[INT_X]* fac; /* -r_x */
	    }else{  /* third sub-col  */
	      if(l == 0)
		*(s+offset+irow)= -fault[i].xc[INT_Y] * fac; /* -r_y */
	      else if(l == 1)
		*(s+offset+irow)=  fault[i].xc[INT_X] * fac; /*  r_x */
	      else
		*(s+offset+irow)=  0.0;                  /*    0 */
	    }
#else
	    /* rotation and translation part */
	    if(j == 0)	/* first sub-col: omega */
	      if(l==0)		/* those are, like in assemble_block_a
				   the projected coordinates of the
				   fault midpoint (with respect to 
				   the genral, not local projection)
				*/
		*(*f+offset+irow)= -fault[i].px[INT_Y] * fac;
	      else 
		*(*f+offset+irow)=  fault[i].px[INT_X] * fac;
	    else if(j == 1)	/* second sub-col: vx_0 */
	      if(l==0)
		*(*f+offset+irow) = fac; /* 1 or -1 */
	      else 
		*(*f+offset+irow) = 0.0;
	    else  /* third sub-col: vy_0  */
	      if(l==0)
		*(*f+offset+irow) = 0.0;
	      else 
		*(*f+offset+irow) = fac; /* 1 or -1 */
#endif
	  }else{
#ifdef BLOCK_SPHERICAL
	    *(s+offset+irow) = 0.0;
#else
	    *(*f+offset+irow) = 0.0;
#endif
	  }
	}
      }
    }
  }
#ifdef BLOCK_SPHERICAL
  /* 

  assemble B which is mf rows by m columns and links the cartesian
  system relative velocities to polar (lon/lat) displacements

  */
  tsize = mf * m;
  /* allocate space for B */
  my_vecalloc(&b,tsize,"assemble_block_f: b");
  /*  have to initialize with zeros  */
  for(i=0;i < tsize;i++)
    b[i] = 0.0;
  /* PS: don't use  calloc to allow memory debugging */
  for(k=0;k < nflt;k++){		/* column flt loop */
    for(j=0;j < pdim;j++){        /* column dimension loop */
      icol = k * pdim + j;
      offset = icol * mf;
      /* assign only one row, i=k, the rest of B is zeroes */
      for(l=0;l < BLOCK_DIM;l++){ /* row dimension loop */
	irow = BLOCK_DIM * k + l;
	if(l==0)// lon or phi direction
	  *(b+offset+irow) =  fault[k].pbase[INT_PHI*3  +j];
	else// lat or -theta direction
	  *(b+offset+irow) = -fault[k].pbase[INT_THETA*3+j];
      }
    }
  }
#ifdef SUPER_DEBUG
  if(!written){
    print_matrix_ftrn_file(s,m,n,"s.dat",FALSE);
    print_matrix_ftrn_file(b,mf,m,"b.dat",FALSE);
    /* `written' gets set to true later */
  }
#endif /* end debug if */
  /* 
     obtain F from B.S 
  */
  /* allocate space for F */
  tsize = n * mf;
  my_vecrealloc(f,tsize,"assemble_block_f: f");  
  calc_AB_ftn(b,mf,m,s,n,*f);
#ifdef MEM_ALLOC_DEBUG
  fprintf(stderr,"assemble_block_f: freeing b and s\n");
#endif /* end mem alloc if */
  free(b);free(s);
#endif /* end block_spherical part if */
#ifdef SUPER_DEBUG
  if(!written){
    print_matrix_ftrn_file(*f,mf,n,"f.dat",FALSE);
    written=TRUE;
  }
#endif
}
/*

  assemble the complete K(m=nrgp+nrsp+nfdamp+nxdamp,
                          n=nbase*nbl) matrix

  -----n-----
  ( A - E  )   E = D.G.F and K2=I.G.F            (m1 = nrgp*2/3   rows)
  (   - K2 )                                     (m2 = nrsp*6     rows) 
  (  g G.F )   if damping of slip   is activated (m3 = nfdamp     rows) 
  (     aI )   if damping of x sol. is activated (m4 = nxdamp(0/n), I being the 
                                                  unity matrix) 

  m1 = nrgp * BLOCK_DIM or nrgp * 3 (for spherical);
  m2 = nrsp * 6;
  m3 = nsnf or nflt rows. if nsnf, full G.F matrix, if nflt only
                          normal motion parts
  m4 = 0 or n depending on nxdamp			  

  on input, can be m1, m1+m2, or m1+m2+m3 

*/
void assemble_block_k(COMP_PRECISION **kmat, int m,int m1,
		      int m2,int n, int nrgp, int nrsp,	int nflt,
		      COMP_PRECISION *a, COMP_PRECISION *e,
		      COMP_PRECISION *k2mat,my_boolean rigid,
		      my_boolean damp_nslip,COMP_PRECISION *gf, 
		      struct bflt *fault, int nsnf,int nslip,
		      int nfdamp, int nxdamp,
		      COMP_PRECISION xdamp)
{
  int i,j,k,nm,o1,offset,offset2,m1m2,icount,kcount,tsize;
  static int k_size=0;
  nm = n * m;
  m1m2 = m1 + m2;

  if(m1m2 + nfdamp + nxdamp != m){
    fprintf(stderr,"assemble_block_k: logic error: m: %i m1: %i m2: %i nfdamp: %i nxdamp: %i\n",
	    m,m1,m2,nfdamp,nxdamp);
    exit(-1);
  }
  tsize = nm;
  if(k_size != tsize){
    // reallocate K matrix
    my_vecrealloc(kmat,tsize,"assemble_block_k: kmat");
    k_size = tsize;
  }
  //
  for(i=0;i < tsize;i++)             /* initialize K with zeroes */
    *(*kmat+i) = 0.0;		
  // 
  // check the damping assignment numbers
  //
  if(damp_nslip){		/* damping part */
    if(rigid){
      fprintf(stderr,"assemble_block_k: logic error: both normal damping and rigid flag set\n");
      exit(-1);
    }
    if(nsnf <= 0){
      fprintf(stderr,"assemble_block_k: error with m1: %i m2: %i m: %i nsnf: %i\n",
	      m1,m2,m,nsnf);
      exit(-1);
    }
    /* prepare j-loop in case of damping:
       if nfdamp is nflt and nslip=2, use every other (normal) 
       component. 
    */
    switch(nslip){
    case 1:
      if(nfdamp != nflt){
	fprintf(stderr,"assemble_block_k: error: nslip: %i nflt: %i nfdamp: %i\n",
		nslip,nflt,nfdamp);
	exit(-1);
      }
      break;
    case 2:
      if((nfdamp != nflt)&&(nfdamp != nsnf)){
	fprintf(stderr,"assemble_block_k: error: nslip: %i nflt: %i nfdamp: %i\n",
		nslip,nflt,nfdamp);
	exit(-1);
      }
      break;
    default:
      fprintf(stderr,"assemble_block_k: error: nslip: %i\n",
	      nslip);
      exit(-1);
      break;
    }
  }else{
    /* 
       no slip damping 
    */
    if(nfdamp != 0){
      fprintf(stderr,"assemble_block_k: logic error: no damping set, but nfdamp %i\n",
	      nfdamp);
      exit(-1);
    }
  }
#ifdef BLOCK_SPHERICAL
  if(m1 != 3 * nrgp){
    fprintf(stderr,"assemble_block_k: error: m1: %i nrgp: %i pdim: %i (spherical)\n",
	    m1,nrgp,3);
    exit(-1);
  }
#else
  if(m1 != BLOCK_DIM * nrgp){
    fprintf(stderr,"assemble_block_k: error: m1: %i nrgp: %i pdim: %i\n",
	    m1,nrgp,BLOCK_DIM);
    exit(-1);
  }
#endif
  if((!rigid)&&(nflt)){
    for(j=o1=0;j < n;j++,o1+=m){ /* j-loop through solution
				    parameters, columns of K */
      //
      // actually subtract interseismic deformation part 
      /* A - E part */
      //
      offset2 = o1;		/* start off at first row */
      /* this is the increment */
      offset  = j * m1;		/* A and E are nrgp*pdim by n */
      for(i=0;i < m1;i++){	/* i loop through GPS data rows  */
	*(*kmat + offset2 + i) =  a[offset + i] - e[offset + i];
      }
      offset2 += m1;		/* shift offset to end of GPS row */
      /* 
	 if there are stress observations, add the stress part, K2
      */
      if(m2){			
	offset   = j * m2;	/* K2 is nrsp * 6 by n */
	for(i=0;i < m2;i++){	/* -K2 = - I.G.F lower part,
				   loop through stress rows
				*/
	  *(*kmat + offset2 + i) = -k2mat[offset + i];
	}
	/* shift to end of stress rows */
	offset2 += m2;
      }
      if(damp_nslip){		
	/* damping of slip motion part */
	offset = j * nsnf;	/* nfdamp may be nflt or
				   nflt*nslip=nsnf however, G.F is
				   nsnf by n */
	for(i=icount=kcount=0;i < nflt;i++){ /* fault loop */
	  for(k=0;k < nslip;k++){ /* slip direction loop */
	    /* only those fault slip motions with damping 
	       have entries other than zero */
	    if(fault[i].use_damp[k]){
	      *(*kmat+ offset2 + icount) = 
		fault[i].sdamp[k]  * gf[offset + kcount];
	      icount++;		/* increases only with damped
				   directions */
	    }
	    kcount++;		/* increases with all possibble
				   slip directions 0...nfns-1*/
	  }
	}
	if((icount != nfdamp)||(kcount != nsnf)){
	  fprintf(stderr,"assemble_block_k: assembly error: icount: %i m3: %i kcount: %i nsnf: %i\n",
		  icount,nfdamp,kcount,nsnf);
	  exit(-1);
	}
	/* shift offset to end of slip damping */
	offset2 += nfdamp;
      }	/* end slip damping loop */
      if(nxdamp){
	/* 
	   solution vector damping toward rigid solution,
	   only the diagonal elements are unity * damping factor,
	   else leave at zero
	*/
	*(*kmat + offset2 + j) = xdamp ;
      }
    } /* end j-loop through the columns of K  */
    /* end non-rigid computation */
  }else{
    /* 
       only rigid computation 
    */
    if(nfdamp || m2){
      fprintf(stderr,"assemble_block_k: logic error: rigid, but m2: %i and nfdamp: %i\n",
	      m2,nfdamp);
      exit(-1);
    }
    for(j=o1=0;j < n;j++,o1 += m){
      /* the upper part of K = A */
      offset2 = o1;
      offset = j * m1;
      for(i=0;i < m1;i++)
	*(*kmat + offset2 + i) = a[offset + i];
      offset2 += m1;		/* shift to end of GPS */
      /* 
	 damping of solution vector toward rigid solution,
	 diagonal elements are 1 * xdamp
      */
      if(nxdamp){
	*(*kmat + offset2 + j) = xdamp;
      }
    }
    fprintf(stderr,"assemble_block_k: WARNING: rigid calculation: no interseismic deformation\n");
  } 
}

/*

  given a fault structure with initialized interaction coefficients,
  and a fault slip vector slip[nslip*nflt] obtained by multiplying
  G.F.x (calc_Ax_ftn(gf,nsnf,n,xsol,fslip)) calculate the stress
  matrix in short (6) vector storage format at point ipnt

  the effect of the stress should be taken negative like
  the velocity field which is block motion - effect of locking
  

*/
void assemble_stress_matrix(COMP_PRECISION *sloc,int ipnt,
			    struct bflt *fault,
			    COMP_PRECISION *fslip,int nflt,
			    int nslip)
{
  int i,j,k,offset,idir;
  for(i=0;i<6;i++)
    sloc[i] = 0.0; /* only upper right half */
  for(j=0;j < nflt;j++)	/* loop through faults */
    for(i = 0;i < nslip;i++){	/* add all slip modes */
      offset = j * nslip + i;	/* stress contributions  */
      idir = block_slip_direction(i,(fault+j));
      for(k=0;k < 6;k++)
	sloc[k] -= (fault+j)->s[ipnt].sc[idir][k] * 
	  fslip[offset]; 
    }
}
/*

  calculate the fit vector component jrow where jrow runs from 
  0 .. m, the number of velocity components (DIM * # GPS obs.)

  if jrow == -1, will produce the whole solution (fit) vector

  if assemble_k is false, the routine simply performs the
  multiplication of row jrow of K with xsol[n] to return 
  the y[jrow] prediction, assuming that kmat has been assembled
  already
  
  if assemble_k is true, then K will be reassembled. for this to 
  work, a, d, and g will have to be precomputed and allocated,
  kmat only has to be allocated

  further input is nrb, the block[] array, nsnf
  

*/
void block_assemble_fit_vector(COMP_PRECISION **kmat,int m, 
			       int ms,int m1,int m2,
			       int n,int nrgp, 
			       int nrsp, int nflt, int nslip,
			       COMP_PRECISION *xsol, int jrow,
			       COMP_PRECISION *y, 
			       COMP_PRECISION *yc,
			       my_boolean assemble_k,
			       COMP_PRECISION *a,
			       COMP_PRECISION **d,
			       COMP_PRECISION *g,
			       COMP_PRECISION **imat,
			       int nrb,struct bck *block,
			       int nsnf, struct bflt **fault,
			       my_boolean rigid,
			       COMP_PRECISION **gf,
			       my_boolean invert_for_ld,
			       my_boolean invert_for_cfac,
			       my_boolean ld_changed,
			       COMP_PRECISION *gx, 
			       COMP_PRECISION *sx,
			       COMP_PRECISION *stress_depths,
			       struct prj projection,
			       my_boolean damp_nslip,
			       int nfdamp, int nxdamp,
			       COMP_PRECISION xdamp,
			       struct bmd *mod)
{
  int pdim;
  /* 

  WARNING: F, E, and K2 will only live here!

  */
  COMP_PRECISION *f,*e,*k2mat;
#ifdef BLOCK_SPHERICAL
  pdim = 3;
#else
  pdim = BLOCK_DIM;
#endif
  f = e = k2mat = NULL;
  if((invert_for_ld) && (ld_changed)){	
    /* if the locking depths change, have to
       recalculate the fault geometry and the
       interaction coefficients
    */
    change_locking_depths(fault,nflt,nrgp,nrsp,nslip,
			  gx, sx,stress_depths,rigid,
			  projection,d,imat,(xsol+n),mod);
    assemble_k = TRUE;
  }
  if(assemble_k){	
    assemble_block_fltdep_matrices(&f,gf,&e,&k2mat,g,*d,*imat,
				   m,n,nflt,nrb,nrgp,nrsp,nsnf,
				   xsol,*fault,block,
				   invert_for_ld,invert_for_cfac);
    assemble_block_k(kmat,ms,m1,m2,n,nrgp,nrsp,nflt,a,e,
		     k2mat,rigid,damp_nslip,*gf,*fault,nsnf,
		     nslip,nfdamp,nxdamp,xdamp);
    if(0)
      fprintf(stderr,"2: |a|: %g |e|: %g |gf|: %g |d|: %g |k2|: %g |k|: %g\n",
	      norm(a,BLOCK_DIM*nrgp*BLOCK_NBASE*nrb),
	      norm(e,pdim*nrgp*BLOCK_NBASE*nrb),
	      norm(*gf,n*nslip*nflt), 
	      norm(*d,pdim * nrgp *nslip* nflt),
	      norm(k2mat,nrsp*6*BLOCK_NBASE*nrb),
	      norm(*kmat,n*ms));
  }
  if(jrow > -1){
    evaluate_block_solution(*kmat,ms,n,xsol,y,yc,mod);
    *(y)=y[jrow];
    my_vecrealloc(&y,1,"block_assemble_fit_vector: y");
#ifdef BLOCK_SPHERICAL
    *(yc)=yc[jrow];
    my_vecrealloc(&yc,1,"block_assemble_fit_vector: yc");
#endif    
  }else{		
    /* 
       produce whole solution vector  with m elements 
    */
    evaluate_block_solution(*kmat,ms,n,xsol,y,yc,mod);
  }
#ifdef MEM_ALLOC_DEBUG
  fprintf(stderr,"block_assemble_fit_vector: freeing f, e, and k2\n");
#endif
  free(f);free(e);free(k2mat);
}
/*
  
  calculate the derivatives with respect to the 
  xsol[1..n...n+nflt] vector at xsol
  
*/
void block_assemble_dyda_matrix(COMP_PRECISION **dyda,
				COMP_PRECISION **kmat,
				int m, int ms,
				int m1,int m2,int n,
				int nrgp,int nrsp,
				int nflt, int nslip,
				COMP_PRECISION *xsol, 
				COMP_PRECISION *a,
				COMP_PRECISION **d,
				COMP_PRECISION *g,
				COMP_PRECISION **gf,
				COMP_PRECISION **imat,
				int nrb,struct bck *block,
				int nsnf, struct bflt **fault,
				my_boolean rigid,
				my_boolean invert_for_ld,
				my_boolean invert_for_cfac,
				COMP_PRECISION *gx, 
				COMP_PRECISION *sx,
				COMP_PRECISION *stress_depths,
				struct prj projection,
				my_boolean damp_nslip,
				int nfdamp,int nxdamp,
				COMP_PRECISION xdamp,
				struct bmd *mod)
{
  int j,os1,na,md,ild[2],icf[2];
  COMP_PRECISION *f1,*f2,*f3,dx;
  static my_boolean init=FALSE;
  my_boolean ld_changed, cfac_changed;
  /* 
     nr of parameters for block solution
  */
  na = n;
  /* 
     add additional degrees of freedom 
     and determine bounds for determining if locking depths
     or coupling factors changed
  */
  if(invert_for_cfac && invert_for_ld){
    /* inversion for both ld and cfac */
    na += 2*nflt;
    ild[0] = n;     ild[1] = n+nflt;
    icf[0] = n+nflt;icf[1] = na;
  }else if(invert_for_cfac){
    na += nflt;
    icf[0] = n;     icf[1] = n+nflt;
  }else if(invert_for_ld){
    ild[0] = n;     ild[1] = n+nflt;
    na += nflt;
  }
  /* 
     number of real data  points + damping 
  */
  md = mod->mgd + mod->m2 + mod->nfdamp + mod->nxdamp;
  /*  

  matrix to hold the derivatives of each data point with respect 
  to changes in the parameters
  
  */
  my_vecrealloc(dyda,md * na,"block_assemble_dyda_matrix: dyda");
  my_vecalloc(&f1,md,"block_assemble_dyda_matrix: f1");
  my_vecalloc(&f2,md,"block_assemble_dyda_matrix: f2");
  my_vecalloc(&f3,md,"block_assemble_dyda_matrix: f3");
  if(!init)
    fprintf(stderr,"calculate_dyda: calculating derivatives of y with resp. to xsol\n");
  //
  // original solution (fit vector) at x = xsol
  //
  block_assemble_fit_vector(kmat,m,ms,m1,m2,n,nrgp,nrsp,nflt,nslip,
			    xsol,-1,f1,mod->vmodc,
			    (invert_for_ld || invert_for_cfac),
			    a,d,g,imat,nrb,block,nsnf,fault,
			    rigid,gf,invert_for_ld,invert_for_cfac,
			    TRUE,gx,sx,stress_depths,projection,
			    damp_nslip,nfdamp,nxdamp,xdamp,mod);
  //
  // loop though all na (either n or n + nflt) parameters
  for(j=os1=0;j < na;j++,os1+=md){	
    /* 
       will the locking depth change? 
    */
    if((invert_for_ld)&&(j>=ild[0])&&(j<ild[1]))
      ld_changed = TRUE;
    else
      ld_changed = FALSE;
    if((invert_for_cfac)&&(j>=icf[0])&&(j<icf[1]))
      cfac_changed = TRUE;
    else
      cfac_changed = FALSE;
    /*
      set finite difference dx
    */
    dx = fabs(xsol[j]) * 0.01;
    if(dx < 0.01)
      dx = 0.01;
    /* 
       forward difference test vector 
    */
    xsol[j] += dx;		
#ifdef SUPER_DEBUG
    /* 
       for debugging
    */
    if(ld_changed)
      fprintf(stderr,"dyda: flt %4i (%4i): old ld: %11g new ld: %11g\n",
	      j-n+1,nflt,xsol[j]-dx,xsol[j]);
#endif
    //
    // obtain the solution for the changed parameters 
    //
    block_assemble_fit_vector(kmat,m,ms,m1,m2,n,nrgp,nrsp,nflt,
			      nslip,xsol,-1,f2,mod->vmodc,
			      (cfac_changed || ld_changed),
			      a,d,g,imat,nrb,block,nsnf,fault,
			      rigid,gf,invert_for_ld,
			      invert_for_cfac,ld_changed,gx,sx,
			      stress_depths,projection,damp_nslip,
			      nfdamp,nxdamp,xdamp,mod);
    if(ld_changed || cfac_changed){
      /* 
	 central differences, non-linear 
      */
      /* take backward step, change sign of dx*/
      xsol[j] -= 2.0*dx;
      block_assemble_fit_vector(kmat,m,ms,m1,m2,n,nrgp,nrsp,nflt,
				nslip,xsol,-1,f3,mod->vmodc,
				(cfac_changed || ld_changed),
				a,d,g,imat,nrb,block,nsnf,fault,
				rigid,gf,invert_for_ld,
				invert_for_cfac,ld_changed,gx,sx,
				stress_depths,projection,
				damp_nslip,nfdamp,nxdamp,xdamp,
				mod);
      /* calculate central differences */
      c_eq_a_minus_b((*dyda+os1),f2,f3,md); 
      scale_vector((*dyda+os1),1.0/(2.0*dx),md);
      /* change sign so that dx gets reset properly */
      dx = -dx;
    }else{
      /* 
	 forward differences
      */
      /* 
	 assign scaled diff as derivative 
      */
      c_eq_a_minus_b((*dyda+os1),f2,f1,md); 
      scale_vector((*dyda+os1),1.0/dx,md);
    }
#ifdef SUPER_DEBUG
    fprintf(stderr,"dyda: i: %3i n: %3i na: %3i |dyda_m|: %11g\n", 
	    j,n,na,norm((*dyda+os1),md)); 
#endif
    /* 
       reset test vector 
    */
    xsol[j] -= dx;		
  }
  free(f1);free(f2);free(f3);
#ifdef SUPER_DEBUG
  if(!init)
    print_matrix_ftrn_file(*dyda,md,na,"dyda.dat",FALSE);
#endif
  init = TRUE;
}


/*

  calculate the weighted (chi2) type deviation between two 
  vectors x and y of dimension n and supplied weights which 
  are given as sigmas and applied as 1/sigma
  
  chi2 = sum_i ((x_i - y_i)/sig_i)^2

  the weight beta are applied for elements

  0 ... n1-1 and n1 ... n-1 respectively

  typically, will be called similarly to:

  chi2[0] = block_chi_square(mod.vmod,mod.v,mod.sigv,mod.mgd,
                             mod.m2,beta,(chi2+1),(chi2+2));

*/
COMP_PRECISION block_chi_square(COMP_PRECISION *x,
				COMP_PRECISION *y, 
				COMP_PRECISION *sigma,
				int n1, int n2,
				COMP_PRECISION beta,
				COMP_PRECISION *vchi2,
				COMP_PRECISION *schi2)
{
  int i,n;
  n = n1 + n2;
  for(*vchi2 = 0.0,i=0;i < n1;i++) /* velocities  */
    *vchi2 += square((x[i]-y[i])/sigma[i]);
  for(*schi2 = 0.0,i=n1;i < n;i++) /* stresses */
    *schi2 += square((x[i]-y[i])/sigma[i]);
  if(n2 == 0){
    return *vchi2;
  }else{
    return ((*vchi2) + beta * (*schi2))/(1.0 + beta);
  }
}
/* 

   calculate the norm of the data vector applying different weights
   for velocities and stress

   is the reference frame was modified, this will add 
   the correction velocities

*/
COMP_PRECISION block_data_norm(COMP_PRECISION *x,
			       COMP_PRECISION *sigma,
			       int n1, int n2,
			       COMP_PRECISION beta,
			       struct bmd *mod)
{
  COMP_PRECISION suma,sumb;
  int i,n;
  n = n1 + n2;
  if(!mod->changed_reference_frame){ 
    for(suma=0.0,i=0;i<n1;i++)
      suma += square(x[i]/sigma[i]);
  }else{/* input vel data was  modified */
    for(suma=0.0,i=0;i<n1;i++)
      suma += square((x[i] + mod->vcorp[i])/sigma[i]);
  }
  for(sumb=0.0,i=n1;i<n;i++)
    sumb += square(x[i]/sigma[i]);
  if(n2 == 0)
    return suma;
  else
    return (suma + beta * sumb)/(1.0 + beta);
}
/* determine the right index for the slip directions
   given the dip of a fault. vertical faults are strike and normal,
   non-vertical strike and dip 
*/
int block_slip_direction(int j,struct bflt *fault)
{
  if(j == 0){		/* first slip direction */
    return STRIKE;
  }else if(j == 1){		/* second direction */
    if(fault->vertical)
      return NORMAL;	/* dip=90, normal direction opening */
    else
      return DIP;		/* dip!=90, dip direction motion
				   with normal component */
  }else{
    fprintf(stderr,"block_slip_direction: nslip error, j: %i\n",
	    j);
    exit(-1);
  }
  return 0;
}
/* 

   assign additional, non-block motion solution given
   an inversion for slip locking factors and/or locking depth


   this routine gets used extensively by the random 
   inversion scheme

*/

/* 
   if defined, uniform distribution, of locking depths 
   else, will use Gaussian
*/
#define UNIFORM_RANDOM_LD

void assign_additional_sol_values(COMP_PRECISION *x,
				  int n, int na,int nflt,
				  struct bflt *fault,
				  my_boolean invert_for_ld,
				  my_boolean invert_for_cfac,
				  long int *seed,
				  int mode)
{
  int i,j;
  static COMP_PRECISION mean_ld=-1.0,ld_range;
  if((!invert_for_ld)&&(!invert_for_cfac)){ /* no additional parameters */
#ifdef DEBUG
    if(na != n){
      fprintf(stderr,"assign_additional_sol_values: na != n but no inversion for cfac or ld\n");
      exit(-1);
    }
#endif
    return;
  }
  /* assign the solution vector such that locking depths come 
     first, then locking factors if both are to be inverted for */
  i = n;
  if(invert_for_ld){		/* assign initial locking depths  */
    if(mode == INIT_ADD_SOL){
      mean_ld = 0.0;
      for(j=0;j < nflt;j++,i++){
	x[i] = fault[j].ld;
	mean_ld += fault[j].ld;
      }
      mean_ld /= (COMP_PRECISION)nflt;
#ifdef UNIFORM_RANDOM_LD
      /* uniform range */
      ld_range = mean_ld * 2.0;
      fprintf(stderr,"assign_additional_sol_values: random lds will range uniformly from %g to %g\n",
	      mean_ld-ld_range/2,mean_ld+ld_range/2);
#else
      /* gaussian */
      ld_range = mean_ld/2.0;
      fprintf(stderr,"assign_additional_sol_values: random lds will have %g Gaussian stddev around %g\n",
	      ld_range, mean_ld);
#endif
    }else if(mode == RANDOM_ADD_SOL){
      if(mean_ld == -1 ){
	fprintf(stderr,"assign_additional_sol_values: error, initalize first\n");
	exit(-1);
      }
      for(j=0;j<nflt;j++,i++){
	/* 
	   vary random locking depth around the mean locking depth
	   with a range ~two times the mean locking depth but not
	   shallower than 0.25 km and not larger than 40 km
	*/
#ifdef UNIFORM_RANDOM_LD
	x[i] = mean_ld + (-.5 + myrand(seed))*ld_range;
#else
	/* Gaussian */
	x[i] = mean_ld + mygauss_randnr(ld_range,seed);
#endif	
      }
    }else{
      fprintf(stderr,"assign_additional_sol_values: mode error\n");
      exit(-1);
    }
  }
  if(invert_for_cfac){		/* assign initial locking factors */
    if(mode == INIT_ADD_SOL){
      for(j=0;j < nflt;j++,i++)
	x[i] = 1.0;
    }else if(mode == RANDOM_ADD_SOL){
      for(j=0;j<nflt;j++,i++)
	x[i] = myrand(seed);	/* have the random factors range between 
				   0 and 1 */
    }else{
      fprintf(stderr,"assign_additional_sol_values: mode error\n");
      exit(-1);
    }
  }
  /* 

  make sure solution is physically reasonable
  
  */
  check_solution_vector(x,n,nflt,invert_for_cfac,invert_for_ld);
#ifdef DEBUG
  if(i != na){
    fprintf(stderr,"assign_additional_sol_values: na != i, parameter mismatch\n");
    exit(-1);
  }
#endif
}
/* 

the fault lockings depths have been changed,
hence we have to recalculate  the fault geometry, 
interaction coefficients, D, and I, if there are stresses


*/
void change_locking_depths(struct bflt **fault, int nflt,
			   int nrgp, int nrsp, int nslip,
			   COMP_PRECISION *gx, 
			   COMP_PRECISION *sx,
			   COMP_PRECISION *stress_depths,
			   my_boolean rigid,
			   struct prj projection,
			   COMP_PRECISION **d,
			   COMP_PRECISION **imat,
			   COMP_PRECISION *new_ld,
			   struct bmd *mod)
{
  int i;
  static unsigned int ncalled=0;
  ncalled++;
  for(i=0;i < nflt;i++){ /* reassign geometry and 
			    interaction coefficients */
    if(new_ld[i] <= 0.0){
      fprintf(stderr,"change_locking_depths: error: call %i: fault %i: locking depth: %g\n",
	      ncalled,i+1,new_ld[i]);
      exit(-1);
    }
    assign_fault_locking_depth_parameters((*fault+i),new_ld[i],
					  projection,FALSE,i);
    get_bflt_intcoeff(fault,i,gx,nrgp,sx,nrsp,stress_depths,
		      rigid);
  }
  assemble_block_d(d,*fault,nflt,nrgp,nslip,mod);
  if(nrsp)
    assemble_block_i(imat,*fault,nflt,nrsp,nslip);
}
/* 
   calculate the RMS of strike and normal fault slip motion
*/
void calc_fault_sn_rms(COMP_PRECISION *vslip,int nflt,
		       int nsdir, COMP_PRECISION *fsrms,
		       COMP_PRECISION *fnrms)
{
  int i,j;
  *fsrms = *fnrms = 0.0;	/* fault strike and normal RMS  */
  for(i=0;i < nflt;i++){
   for(j=0;j < nsdir;j++){	 
     if(j == 0)		/* add to direction type RMSs */
       *fsrms += square(vslip[i*nsdir+j]);
     else
       *fnrms += square(vslip[i*nsdir+j]);
   }
  }
  *fsrms = sqrt(*fsrms/((COMP_PRECISION)nflt));
  *fnrms = sqrt(*fnrms/((COMP_PRECISION)nflt));
}
/*

initialize the block model structure

*/
void init_block_mods(struct bmd **mod)
{
  *mod = (struct bmd *)calloc(1,sizeof(struct bmd));
  if(!(*mod))
    MEMERROR("init_block_mods");
  (*mod)->pbase_init = (*mod)->cvel_init = (*mod)->block_init = FALSE;
  // set pointers to zero so that we can use realloc safely
  (*mod)->kmat = (*mod)->a = (*mod)->d = (*mod)->g = (*mod)->f = (*mod)->gf = 
    (*mod)->e = (*mod)->xsol  = (*mod)->cov = (*mod)->vmod   = (*mod)->sx =
    (*mod)->imat = (*mod)->k2mat = (*mod)->vc = (*mod)->sigvc =  (*mod)->gx =
    (*mod)->vmodc = (*mod)->pbase = (*mod)->gcx = (*mod)->gcx = (*mod)->v =
    (*mod)->sigv = (*mod)->rho = (*mod)->gpx = NULL;
  (*mod)->bcode = NULL;
  (*mod)->vcorp = NULL;
  (*mod)->fault = NULL;
  (*mod)->sigma = NULL;
  (*mod)->fblock_sites = NULL;
}
/* 

   initialize a pointer of n blocks

*/
void init_blocks(struct bmd *mod, int n)
{
  int i,ilim;
  struct bck *ix;
  static int oldn;
  if(mod->block_init){
    /* 
       resize only 
    */
    mod->block = (struct bck *)realloc(mod->block,n*sizeof(struct bck));
    if((!mod->block)&&(n!=0))
      MEMERROR("init_blocks: 1");
    if(n > oldn){
      /* 
	 actually increas the number of blocks 
      */
      /* zero out remainder
	 get zero blocks
      */
      ix = (struct bck *)calloc((n-oldn),sizeof(struct bck));
      if(!ix)
	MEMERROR("init_blocks: 2");
      /* copy the zero blocks to the end of the new blocks */
      ilim = n-oldn;
      for(i=0;i<ilim;i++)
	copy_block((ix+i),(mod->block+oldn+i));
      free(ix);
      for(i=oldn;i < n;i++)
	mod->block[i].fixed = mod->block[i].rot_c = FALSE;
    }
    oldn = n;
  }else{
    /* 
       first call 
    */
    mod->block = (struct bck *)calloc(n,sizeof(struct bck));
    if(!mod->block)
      MEMERROR("init_blocks: orig, n==0?");
    for(i=0;i < n;i++)
      mod->block[i].fixed = mod->block[i].rot_c = FALSE;
    oldn = n;
    mod->block_init = TRUE;
  }
}

/* 
   convert a solution for cartesian velocities to 
   a regular solution that has lon/lat velocity components

*/
void convert_cart_sol(COMP_PRECISION *vmodc,COMP_PRECISION *vmod,
		      struct bmd *mod)
{
  int i,j;
  COMP_PRECISION pv[3];
  /* convert velocities */
  for(i=j=0;i < mod->nrgp;i++,j+=BLOCK_DIM){
    /* go from cartesian to polar (lon/lat)  */
    cv2pv((vmodc+i*3),pv,(mod->pbase+i*9));
    vmod[j+INT_X] =  pv[INT_PHI];
    vmod[j+INT_Y] = -pv[INT_THETA];
  }
  a_equals_b_vector((vmod+mod->mgd),(vmodc+mod->m1),
		    (mod->m2+mod->nfdamp+mod->nxdamp));
}
/* 

check if the solution vector is physically reasonable, 
if not, adjust

*/
void check_solution_vector(COMP_PRECISION *x, 
			   int n, int nflt,
			   my_boolean invert_for_cfac,
			   my_boolean invert_for_ld)
{
  int ild[2]={0,0},icf[2]={0,0},i;
  /* 
     add additional degrees of freedom 
     and determine bounds for determining if locking depths
     or coupling factors changed
  */
  if(invert_for_cfac && invert_for_ld){
    /* inversion for both ld and cfac */
    ild[0] = n;     ild[1] = n+nflt;
    icf[0] = n+nflt;icf[1] = n + 2*nflt;
  }else if(invert_for_cfac){
    icf[0] = n;     icf[1] = n+nflt;
  }else if(invert_for_ld){
    ild[0] = n;     ild[1] = n+nflt;
  }
#ifdef DEBUG
  fprintf(stderr,"check_solution_vector: ld: %i (%i - %i) cfac: %i (%i - %i) n: %i nflt: %i\n",
	  invert_for_ld,ild[0],ild[1],invert_for_cfac,icf[0],icf[1],n,nflt);
#endif
  for(i=ild[0];i < ild[1];i++){	/* locking depths */
    if(x[i] < 0.5)
      x[i] = 0.5;	/* minimum */
    if(x[i] > 50)
      x[i] = 50.0;	/* maximum */
  }
  for(i=icf[0];i < icf[1];i++){	/* coupling factors */
    if(x[i] < 0.0)
      x[i]=0.0;	/* minimum */
    if(x[1] > 2.0)
      x[i]=2.0;	/* maximum */
  }
}
/* 
   b = a for blocks

*/
void copy_block(struct bck *a,struct bck *b)
{
  a_equals_b_vector_3d(b->center,a->center);
  a_equals_b_vector_3d(b->xrigid,a->xrigid);
  b->fixed = a->fixed;
  b->rot_c = a->rot_c;
  b->pcnt = a->pcnt;

}
