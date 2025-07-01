/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu


  this file is an addition to geometry.c and holds the low level
  geometry routines

  also does some linear algebra stuff

  
*/
#include "interact.h"
COMP_PRECISION tensor3d_norm(COMP_PRECISION [3][3]);

/*

  rotates a 3D vector (such as a displacement vector) into a new
  coordinate system whose x-axes is rotated around the 
  z-axis by alpha angle in clockwise direction from east

  WARNING: this routine doesn't take azimuth but direction!

  alternatively: rotates a vector in the old x-y-z coordinate system
  into a new one whose x-axis is rotated by an angle alpha from
  west-east (old x-axis) in a counterclockwise direction

*/
void rotate_vec(COMP_PRECISION *xin, COMP_PRECISION *xout,
		double cos_alpha, double sin_alpha)
{
  xout[INT_X]= ( cos_alpha * xin[INT_X]) + (sin_alpha * xin[INT_Y]);
  xout[INT_Y]= (-sin_alpha * xin[INT_X]) + (cos_alpha * xin[INT_Y]);
  xout[INT_Z]=  xin[INT_Z];
}
// similar for 2-D in X-Y plane
//   WARNING: this routine doesn't take azimuth but direction!
void rotate_vec2d(COMP_PRECISION *xin, COMP_PRECISION *xout,
		  double cos_alpha, double sin_alpha)
{
  xout[INT_X]= ( cos_alpha * xin[INT_X]) + (sin_alpha * xin[INT_Y]);
  xout[INT_Y]= (-sin_alpha * xin[INT_X]) + (cos_alpha * xin[INT_Y]);
}

/* 
   rotate a 3x3 tensor using a general rotation matrix r 
   xout = r . xin . r^T
*/
void rotate_mat(COMP_PRECISION xin[3][3],COMP_PRECISION xout[3][3],
		COMP_PRECISION r[3][3])
{
  COMP_PRECISION tmp[3][3];
  int i,j,k;
  // calculate xin . r^T
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      tmp[i][j] = 0.0;
      for(k=0;k<3;k++)
	tmp[i][j] += xin[i][k] * r[j][k];
    }
  // calculate r . tmp
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      xout[i][j] = 0.0;
      for(k=0;k<3;k++)
	xout[i][j] += r[i][k] * tmp[k][j];
    }
}
/* 
   rotate a SYMMTERIC 3x3 matrix 

   by angle alpha around z axis, cos(alpha) and sin(alpha) are given
   as input
   
   WARNING: this routine doesn't take azimuth but direction!

*/
void rotate_mat_z(COMP_PRECISION xin[3][3],
		  COMP_PRECISION xout[3][3],
		  double cos_alpha, 
		  double sin_alpha)
{
  double cos_square,sin_square,sin_cos;
  /* COMP_PRECISION r[3][3]; */
  /* r[0][0]=1;r[0][1]=0;r[0][2]=0; */
  /* r[1][0]=0;r[1][1]=cos_alpha;r[1][2]=-sin_alpha; */
  /* r[2][0]=0;r[2][1]=sin_alpha;r[2][2]= cos_alpha; */
  /* rotate_mat(xin,xout,r); */
  /* return; */
  //fprintf(stderr,"%20.7e %20.7e\n",cos_alpha,sin_alpha);
  cos_square=cos_alpha*cos_alpha;
  sin_square=sin_alpha*sin_alpha;
  sin_cos=cos_alpha*sin_alpha;

  /* has to be optimized */
  xout[INT_X][INT_X]=(COMP_PRECISION)
    ((double)xin[INT_X][INT_X]*cos_square + (double)xin[INT_Y][INT_Y]*sin_square + 
     ((double)xin[INT_X][INT_Y]*2.0*sin_alpha*cos_alpha));
  xout[INT_X][INT_Y]=xout[INT_Y][INT_X]=(COMP_PRECISION)
    ((double)xin[INT_X][INT_Y]*(cos_square-sin_square) + (-(double)xin[INT_X][INT_X] + (double)xin[INT_Y][INT_Y])*sin_cos);
  xout[INT_X][INT_Z]=xout[INT_Z][INT_X]=(COMP_PRECISION)
    ((double)xin[INT_X][INT_Z]*cos_alpha + (double)xin[INT_Y][INT_Z]*sin_alpha);
  xout[INT_Y][INT_Y]=(COMP_PRECISION)
    ((double)xin[INT_Y][INT_Y]*cos_square -
     (2.0*(double)xin[INT_X][INT_Y]*sin_cos) + 
     (double)xin[INT_X][INT_X]*sin_square);
  xout[INT_Y][INT_Z]=xout[INT_Z][INT_Y]=(COMP_PRECISION)
    ((double)xin[INT_Y][INT_Z]*cos_alpha - (double)xin[INT_X][INT_Z]*sin_alpha);
  xout[INT_Z][INT_Z]=xin[INT_Z][INT_Z];
  //fprintf(stderr,"%20.15e\n",tensor3d_norm(xin)-tensor3d_norm(xout));
}

COMP_PRECISION tensor3d_norm(COMP_PRECISION x[3][3])
{
  int i,j;
  double sum = 0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      sum += x[i][j]*x[i][j];
  return((COMP_PRECISION)sqrt(sum));
    
}
void zero_3x3_matrix(COMP_PRECISION x[3][3])
{
  int i,j;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      x[i][j]=0.0;
}
void scale_3x3_matrix(COMP_PRECISION x[3][3],COMP_PRECISION scale)
{
  int i,j;
  for(i=0;i < 3;i++)
    for(j=0;j < 3;j++)
      x[i][j] *= scale;
}
/* x = x + ay */
void add_ay_to_3x3_matrix(COMP_PRECISION x[3][3],
			  COMP_PRECISION y[3][3],COMP_PRECISION a)
{
  int i,j;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      x[i][j] += a*y[i][j];
}
/* x = x+y */
void add_y_to_3x3_matrix(COMP_PRECISION x[3][3],COMP_PRECISION y[3][3])
{
  int i,j;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      x[i][j] += y[i][j];
}

/* 
   Euclidian norm of vector difference in 3D
*/

COMP_PRECISION distance_3d(COMP_PRECISION *x,COMP_PRECISION *y)
{
  COMP_PRECISION tmp;
  tmp = distance_squared_3d(x,y);
  //  if(tmp <= 0.0)return 0.0;else
  return(sqrt(tmp));
}
/*

  calculate the distance between vectors x[n] and y[n] given 
  in COMP_PRECISIOn precision

*/

COMP_PRECISION distance(COMP_PRECISION *x,
			COMP_PRECISION *y, int n)
{
  COMP_PRECISION tmp;
  tmp = distance_squared(x,y,n);
  //if(tmp <= 0.0)return 0.0;else
  return(sqrt(tmp));
}
/*

  calculate the distance between vectors x[n] and y[n] given 
  in float precision

*/
float distance_float(float *x,float *y, int n)
{
  float tmp;
  tmp = distance_squared_float(x,y,n);
  //if(tmp <= 0.0)return 0.0;else
  return(sqrt(tmp));
}

void assign_to_vector(float *x,float val, int n)
{
  int i;
  for(i=0;i<n;i++)
    x[i]= val;
}
void assign_to_vector_3d(float *x,float val)
{
  x[INT_X] = x[INT_Y] = x[INT_Z] = val;
}

/* 

   Square of Euclidian norm of vector difference in 3D and n
   dimensions

*/
COMP_PRECISION distance_squared_3d(COMP_PRECISION *x,
				   COMP_PRECISION *y)
{
  COMP_PRECISION tmp,tmp1;
  tmp1=x[INT_X]-y[INT_X];tmp  = tmp1*tmp1;
  tmp1=x[INT_Y]-y[INT_Y];tmp += tmp1*tmp1;
  tmp1=x[INT_Z]-y[INT_Z];tmp += tmp1*tmp1;
  return(tmp);
}
COMP_PRECISION distance_squared(COMP_PRECISION *x,
				COMP_PRECISION *y, int n)
{
  COMP_PRECISION tmp,tmp1;
  int i;
  for(tmp=0.0,i=0;i<n;i++){
    tmp1 = x[i]-y[i];
    tmp += tmp1*tmp1;
  }
  return(tmp);
}
float distance_squared_float(float *x,float *y, int n)
{
  float tmp,tmp1;
  int i;
  for(tmp=0.0,i=0;i<n;i++){
    tmp1=x[i]-y[i];
    tmp += tmp1*tmp1;
  }
  return(tmp);
}
/* 
   norm in N dimensions
*/
COMP_PRECISION norm(COMP_PRECISION *x,int n)
{
  COMP_PRECISION tmp;
  tmp = dotp(x,x,n);
  //if(tmp <= 0.0)return 0.0;else
  return(sqrt(tmp));
}
/* 
   norm in 3 dimensions 
*/
COMP_PRECISION norm_3d(COMP_PRECISION *x)
{
 COMP_PRECISION tmp;
 tmp = dotp_3d(x,x);
 //if(tmp <= 0.0)return 0.0;else
 return(sqrt(tmp));
}

// dot products
COMP_PRECISION dotp(COMP_PRECISION *x,COMP_PRECISION *y,
		    int n)
{
  double tmp;
  int i;
  for(tmp = 0.0,i=0;i<n;i++)
    tmp += (double)x[i]*(double)y[i];
  return((COMP_PRECISION)tmp);
}
COMP_PRECISION dotp_3d(COMP_PRECISION *x,COMP_PRECISION *y)
{
  double tmp;
  tmp  = (double)x[INT_X]*(double)y[INT_X];
  tmp += (double)x[INT_Y]*(double)y[INT_Y];
  tmp += (double)x[INT_Z]*(double)y[INT_Z];
  return((COMP_PRECISION)tmp);
}
/*
  calculate y = A . x where
  A is m by n, x is [n], and y is [m]
  AND A IS STORED IN THE FORTRAN CONVENTION
*/
void calc_Ax_ftn(COMP_PRECISION *a, int m, int n, 
		 COMP_PRECISION *x,COMP_PRECISION *y)
{
  int i,j,k;
  for(i=0;i<m;i++)// j=0
    y[i] = a[i] * x[0];
  for(j=1,k=m;j<n;j++,k+=m)// j>0
    for(i=0;i<m;i++)
      y[i] += a[k+i] * x[j];
}
/*
  calculate y = A . x where
  A is m by n, x is [n], and y is [m]
  and A is stored in the C  convention
*/
void calc_Ax(COMP_PRECISION *a, int m, int n, 
	     COMP_PRECISION *x,COMP_PRECISION *y)
{
  int i,j,k;
  for(i=k=0;i<m;i++,k+=n){
    y[i] = 0.0;
    for(j=0;j<n;j++)
      y[i] += a[k+j] * x[j];
  }
}

/*
  calculate C = A . B
  where A is m x n, B is n x p, and C is m X p
  ALL ARE ASSUMED TO BE STORED IN THE FORTRAN CONVENTION,e.g.
  A(i,j) = a[j*m+i]
*/
void calc_AB_ftn(COMP_PRECISION *a, int m, int n,
		 COMP_PRECISION *b, int p, COMP_PRECISION *c)
{
  int i,j,k,o1,o2,o3;
  COMP_PRECISION ctmp;
  for(j=o1=o3=0;j <p;j++,o1+=m,o3+=n)  
    /* C matrix column counter */
    for(i=0;i<m;i++){	/* C matrix row counter */
      ctmp = 0.0;
      for(k=o2=0;k<n;k++,o2+=m)
        ctmp += a[i+o2] * b[o3+k];
      c[o1+i] = ctmp;
    }
}



//
// normalize a vector to unity length
void normalize(COMP_PRECISION *x,int n)
{
  COMP_PRECISION length;
  int i;
  length = norm(x,n);
  if(length != 0.0){
    for(i=0;i<n;i++)
      x[i]/=length;
  }
}
void normalize_3d(COMP_PRECISION *x)
{
  COMP_PRECISION length;
  length = norm_3d(x);
  if(length != 0.0){
    x[INT_X] /= length;x[INT_Y] /= length;x[INT_Z] /= length;
  }
}
//  \vec{a} = \vec{b}
void a_equals_b_vector(COMP_PRECISION *a,COMP_PRECISION *b,
		       int n)
{
  memcpy(a,b,n*sizeof(COMP_PRECISION));
}
void a_equals_b_vector_3d(COMP_PRECISION *a,COMP_PRECISION *b)
{
  //a[INT_X] = b[INT_X];a[INT_Y] = b[INT_Y];a[INT_Z] = b[INT_Z];
  memcpy(a,b,3*sizeof(COMP_PRECISION));
}
void swap_ab_vector_3d(COMP_PRECISION *a,COMP_PRECISION *b)
{
  COMP_PRECISION tmp[3];
  a_equals_b_vector_3d(tmp,a);
  a_equals_b_vector_3d(a,b);
  a_equals_b_vector_3d(b,tmp);
}
//
// do \vec{a} = f * \vec{a}
// where f is a scalar
void scale_vector_3d(COMP_PRECISION *a,COMP_PRECISION f)
{
  a[INT_X] *= f;a[INT_Y] *= f;a[INT_Z] *= f;
}
// do \vec{a} = f * \vec{a}
void scale_vector(COMP_PRECISION *a,COMP_PRECISION f, int n)
{
  int i;
  for(i=0;i<n;i++)
    a[i] *= f;
}
// do \vec{a} = \vec{a} + \vec{b}
void add_b_to_a_vector(COMP_PRECISION *a,COMP_PRECISION *b, int n)
{
  int i;
  for(i=0;i<n;i++)
    a[i] += b[i];
}
// do \vec{a} = \vec{a} + \vec{b} for 3-D
void add_b_to_a_vector_3d(COMP_PRECISION *a,COMP_PRECISION *b)
{
  a[INT_X] += b[INT_X];a[INT_Y] += b[INT_Y];a[INT_Z] += b[INT_Z];
}
// do \vec{a} = \vec{a} - \vec{b}
void sub_b_from_a_vector(COMP_PRECISION *a,COMP_PRECISION *b,int n)
{
  int i;
  for(i=0;i<n;i++)
    a[i] -= b[i];
}
// do \vec{a} = \vec{a} - \vec{b} for 3-D
void sub_b_from_a_vector_3d(COMP_PRECISION *a,COMP_PRECISION *b)
{
  a[INT_X] -= b[INT_X];a[INT_Y] -= b[INT_Y];a[INT_Z] -= b[INT_Z];
}
// \vec{c} = \vec{a} + \vec{b}
void c_eq_a_plus_b_3d(COMP_PRECISION *c,COMP_PRECISION *a,
		      COMP_PRECISION *b)
{
  c[INT_X] = a[INT_X] + b[INT_X];c[INT_Y] = a[INT_Y] + b[INT_Y];c[INT_Z] = a[INT_Z] + b[INT_Z];
}
// \vec{c} = \vec{a} - \vec{b}
void c_eq_a_minus_b_3d(COMP_PRECISION *c,COMP_PRECISION *a,COMP_PRECISION *b)
{
  c[INT_X] = a[INT_X] - b[INT_X];c[INT_Y] = a[INT_Y] - b[INT_Y];c[INT_Z] = a[INT_Z] - b[INT_Z];
}
// \vec{c} = \vec{a} - \vec{b}
void c_eq_a_minus_b(COMP_PRECISION *c,COMP_PRECISION *a,COMP_PRECISION *b, int n)
{
  int i;
  for(i=0;i<n;i++)
    c[i] = a[i] - b[i];
}
/*

  calculate cross product c = a \cross b

*/
void cross_product(COMP_PRECISION *a,COMP_PRECISION *b,
		   COMP_PRECISION *c)
{
  c[INT_X] = b[INT_Z] * a[INT_Y] - b[INT_Y] * a[INT_Z];
  c[INT_Y] = b[INT_X] * a[INT_Z] - b[INT_Z] * a[INT_X];
  c[INT_Z] = b[INT_Y] * a[INT_X] - b[INT_X] * a[INT_Y];
}
/*

  find the maximum value of a vector

 */
COMP_PRECISION find_max_vec(COMP_PRECISION *x,int n)
{
  COMP_PRECISION max = -FLT_MAX;
  int i;
  for(i=0;i<n;i++)
    if(x[i] > max)
      max = x[i];
  return max;
}
/*

  find the maximum absolute value of a vector

 */
COMP_PRECISION find_max_abs_vec(COMP_PRECISION *x,int n)
{
  COMP_PRECISION max = 0.0,tmp;
  int i;
  for(i=0;i<n;i++){
    tmp = fabs(x[i]);
    if(tmp > max)
      max = tmp;
  }
  return max;
}
float find_max_abs_vec_float(float *x,int n)
{
  float max = 0.0,tmp;
  int i;
  for(i=0;i<n;i++){
    tmp = fabs(x[i]);
    if(tmp > max)
      max = tmp;
  }
  return max;
}


/*
  
  find the largest offset from the diagonal of a mxn matrix a
  whose absolute value is larger than acutoff

*/ 

int find_lde_max(A_MATRIX_PREC *a, int m, int n, 
		 A_MATRIX_PREC acutoff)
{
  int i,j,k,tmpi,imaxbw;
  imaxbw = 0;
  for(k=i=0;i<m;i++,k+=n)
    for(j=0;j<n;j++)
      if(fabs(*(a+k+j)) > acutoff){// entry is larger than cutoff
	if((tmpi = abs(j-i)) > imaxbw)// distance from diag > max
	  imaxbw = tmpi;
      }
  return imaxbw;
}

COMP_PRECISION square(COMP_PRECISION x)
{
  return x * x;
}

/*

  trace of matrix stored in 3 x 3 format, fortran style

*/
COMP_PRECISION tracemat9(COMP_PRECISION *s)
{
  return s[INT_XX] + s[INT_YY] + s[INT_ZZ];

}
/*
  
  compare two floating point numbers
  but invert results
*/
int inv_compare_flt(const void *va,const void *vb)
{
  COMP_PRECISION *a = (COMP_PRECISION *)va,
    *b = (COMP_PRECISION *)vb;
  if(*a > *b)
    return -1;
  else if(*a == *b)
    return 0;
  else 
    return 1;
}
/*
  
  compare two numbers

*/
int compare_flt(const void *va,const void *vb)
{
  COMP_PRECISION *a = (COMP_PRECISION *)va,
    *b = (COMP_PRECISION *)vb;
  if(*a < *b)
    return -1;
  else if(*a == *b)
    return 0;
  else 
    return 1;
}

/*

calculate the RMS of a vector
sqrt(sum_i x_i ^2 / n )

*/
COMP_PRECISION rms(COMP_PRECISION *x, int n)
{
  COMP_PRECISION dp;
  dp = dotp(x,x,n);
  //if(dp <= 0.0)
  //return dp;
  //else
  return sqrt(dp/(COMP_PRECISION)n);
}


/*
  
allocate floating point vector

*/
void my_vecalloc(COMP_PRECISION **x,int n,char *message)
{
  *x = (COMP_PRECISION *)malloc((size_t)sizeof(COMP_PRECISION)*(size_t)n);
  if(! (*x))
    MEMERROR(message);
#ifdef MEM_ALLOC_DEBUG
  fprintf(stderr,"my_vecalloc: newly allocating %.4f MB of memory (calling routine: %s)\n",
	  (float)(sizeof(COMP_PRECISION)*(size_t)n)/(float)ONE_MEGABYTE,message);
#endif
}
/* allocate integer vector */
void my_ivecalloc(int **iv,int n,char *message)
{
  *iv = (int *)malloc(sizeof(int)*(size_t)n);
  if(!(*iv))
    MEMERROR(message);
#ifdef MEM_ALLOC_DEBUG
  fprintf(stderr,"my_ivecalloc: newly allocating %.4f MB of memory (calling routine: %s)\n",
	  (float)(sizeof(int)*(size_t)n)/(float)ONE_MEGABYTE,message);
#endif
}
/* reallocate floating point vector */
void my_vecrealloc(COMP_PRECISION **x,int n,char *message)
{
  *x = (COMP_PRECISION *)realloc(*x,
				 sizeof(COMP_PRECISION)*(size_t)n);
  if(!(*x))
    MEMERROR(message);
#ifdef MEM_ALLOC_DEBUG
  fprintf(stderr,"my_vecrealloc: re-allocating %.4f MB of memory (calling routine: %s)\n",
	  (float)(sizeof(COMP_PRECISION)*(size_t)n)/(float)ONE_MEGABYTE,
	  message);
#endif
}
/* re-allocate integer vector */
void my_ivecrealloc(int **iv,int n,char *message)
{
  *iv = (int *)realloc(*iv,sizeof(int)*(size_t)n);
  if(!(*iv))
    MEMERROR(message);
#ifdef MEM_ALLOC_DEBUG
  fprintf(stderr,"my_ivecrealloc: re-allocating %.4f MB of memory (calling routine: %s)\n",
	  (float)(sizeof(int)*(size_t)n)/(float)ONE_MEGABYTE,message);
#endif
}
/* calculate the mean of the entries in a vector  */
COMP_PRECISION mean(COMP_PRECISION *x,int ixoff, int n)
{
  COMP_PRECISION s=0.0;
  int i,os;
  for(i=os=0;i<n;i++,os+=ixoff)
    s += x[os];
  return s/(COMP_PRECISION)n;
}
/* calculate the mean of the absolute entries in a vector  */
COMP_PRECISION mean_abs(COMP_PRECISION *x,int ixoff, int n)
{
  COMP_PRECISION s=0.0;
  int i,os;
  for(i=os=0;i<n;i++,os+=ixoff)
    s += fabs(x[os]);
  return s/(COMP_PRECISION)n;
}

/* 

   calculate the weighted mean given vector x with 
   n elements that are spaced x[i*ixoff] where 0 <= i < n
   and weights in w[] with elements w[i]

*/
COMP_PRECISION wmean(COMP_PRECISION *x, int ixoff, int n,
		     COMP_PRECISION *w)
{
  COMP_PRECISION ws=0.0,s=0.0;
  int i,os;
  for(i=os=0;i<n;i++,os+=ixoff){
    s += w[i] * x[os];
    ws += w[i];
  }
  return s/ws;			/* return NaN if no weights
				   assembled */
}
/* same for wmeand(abs(x)) */
COMP_PRECISION wmean_abs(COMP_PRECISION *x, int ixoff, int n,
			 COMP_PRECISION *w)
{
  COMP_PRECISION ws=0.0,s=0.0;
  int i,os;
  for(i=os=0;i<n;i++,os+=ixoff){
    s += w[i] * fabs(x[os]);
    ws += w[i];
  }
  return s/ws;			/* return NaN if no weights
				   assembled */
}
		     
/*
  
  calculate the mean and standard deviation of a vector
  stat has to be [4]: mean stddev min max
  fast and inaccurate method for stddev

*/
void calc_vec_stat(COMP_PRECISION *x, int n, COMP_PRECISION *stat)
{
  int i;
  stat[0]=0.0;
  stat[1]=0.0;
  stat[2]=FLT_MAX;stat[3]=-FLT_MAX;// mean, STD, min, max
#ifdef DEBUG
  if(n <= 1){
    fprintf(stderr,"calc_vec_stat: n (%i) out of range\n",n);
    exit(-1);
  }
#endif
  for(i=0;i<n;i++){
    stat[0] += x[i];
    stat[1] += SQUARE(x[i]);
    // min
    if(x[i] < stat[2])
      stat[2]=x[i];
    // max
    if(x[i] > stat[3])
      stat[3]=x[i];
  }
  //
  // standard deviation (fast but inaccurate)
  stat[1] = stddev(stat[0],stat[1],n);
  //
  // mean
  //
  stat[0] /= (COMP_PRECISION)n;
}
/* 
   
fast, inaccurate standard deviation formula

sigma = (n * \sum (x_i^2 ) - (\sum (x_i))^2)(n*(n-1))

STD = sqrt(sigma)


*/
COMP_PRECISION stddev(COMP_PRECISION sum, COMP_PRECISION sum_sqrd, int n)
{
  double xn,sigma;
  xn = (double)n;
  sigma =  (xn * (double)sum_sqrd - (double)sum * (double)sum) /(xn * (xn - 1.0));

  //if(sigma <= 0.0)
  //  return 0.0;
  //else
  return sqrt(sigma);
}

/* count the number of zero entries in a vector */
int countzero_vec(COMP_PRECISION *x,int n)
{
  int i,nz=0;
  for(i=0;i<n;i++)
    if(fabs(x[i]) < EPS_COMP_PREC)
      nz ++;
  return nz;

}
/* for output purposes */
COMP_PRECISION reformat_small(COMP_PRECISION x)
{
  if(fabs(x) < EPS_COMP_PREC)
    return 0.0;
  else
    return x;
}
