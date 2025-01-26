
void free_vector(COMP_PRECISION *,long ,long );
COMP_PRECISION *vector(long ,long );
void free_a_vector(A_MATRIX_PREC *,long ,long );
A_MATRIX_PREC *a_vector(long ,long );

COMP_PRECISION **matrix(long ,long ,long ,long );
void free_matrix(COMP_PRECISION **,long ,long ,long ,long );

A_MATRIX_PREC **a_matrix(long ,long ,long ,long );
void free_a_matrix(A_MATRIX_PREC **,long ,long ,long ,long );

int *ivector(long,long);
void free_ivector(int *,long ,long );
COMP_PRECISION pythag(COMP_PRECISION ,COMP_PRECISION );
void svbksb(A_MATRIX_PREC **,A_MATRIX_PREC *,
	    A_MATRIX_PREC **,int ,int ,A_MATRIX_PREC *,
	    A_MATRIX_PREC *);
void svdcmp(A_MATRIX_PREC **,int ,int ,A_MATRIX_PREC *,
	    A_MATRIX_PREC **);

/* numerical recipes routines follow */
void svbksb(A_MATRIX_PREC **u,A_MATRIX_PREC *w,
	    A_MATRIX_PREC **v,int m,int n,
	    A_MATRIX_PREC *b,A_MATRIX_PREC *x)
{
  int jj,j,i;
  COMP_PRECISION s,*tmp;
  
  tmp=vector(1,n);
  for (j=1;j<=n;j++) {
    s=0.0;
    if (w[j]) {
      for (i=1;i<=m;i++) s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=1;j<=n;j++) {
    s=0.0;
    for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
    x[j]=s;
  }
  free_vector(tmp,1,n);
}


void svdcmp(A_MATRIX_PREC **a,int m,int n,
	    A_MATRIX_PREC *w,A_MATRIX_PREC **v)
{
  int flag,i,its,j,jj,k,l,nm;
  COMP_PRECISION anorm,c,f,g,h,s,scale,x,y,z,*rv1,
    tmpdbl;
  
  rv1=vector(1,n);
  g=scale=anorm=0.0;
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<=m;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    tmpdbl=fabs(w[i])+fabs(rv1[i]);
    anorm=MAX(anorm,tmpdbl);
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=MIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
	nm=l-1;
	if ((COMP_PRECISION)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((COMP_PRECISION)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((COMP_PRECISION)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=1;j<=m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 30){
	fprintf(stderr,"no convergence in 30 svdcmp iterations\n");
	exit(-1);
      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free_vector(rv1,1,n);
}



COMP_PRECISION pythag(COMP_PRECISION a,COMP_PRECISION b)
{
  COMP_PRECISION absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQUARE(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQUARE(absa/absb)));
}

void free_vector(COMP_PRECISION *v,
		 long nl,long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}
void free_a_vector(A_MATRIX_PREC *v,
		   long nl,long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}


A_MATRIX_PREC *a_vector(long nl,long nh)
{
  A_MATRIX_PREC *v;
  
  v=(A_MATRIX_PREC *)
    calloc(sizeof(A_MATRIX_PREC),
	   (unsigned int) ((nh-nl+1+NR_END)));
  if (!v){
    fprintf(stderr,"allocation failure in vector()\n");
    exit(-1);
  }
  return v-nl+NR_END;
}

COMP_PRECISION *vector(long nl,long nh)
{
  COMP_PRECISION *v;
  
  v=(COMP_PRECISION *)
    calloc(sizeof(COMP_PRECISION),
	   (unsigned int) ((nh-nl+1+NR_END)));
  if (!v){
    fprintf(stderr,"allocation failure in vector()\n");
    exit(-1);
  }
  return v-nl+NR_END;
}
COMP_PRECISION **matrix(long nrl,long nrh,long ncl,long nch)
/* allocate a COMP_PRECISION matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  COMP_PRECISION **m;
  
  /* allocate pointers to rows */
  m=(COMP_PRECISION **) 
    malloc((unsigned int)((nrow+NR_END)*
			  sizeof(COMP_PRECISION*)));
  if (!m) {
    fprintf(stderr,"allocation failure 1 in matrix()\n");
    exit(-1);}
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(COMP_PRECISION *) 
    calloc(sizeof(COMP_PRECISION),
	   (unsigned int)((nrow*ncol+NR_END)));
  if (!m[nrl]) {
    fprintf(stderr,"allocation failure 2 in matrix()\n");
    exit(-1);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

void free_matrix(COMP_PRECISION **m,
		 long nrl,long nrh,long ncl,long nch)
/* free a COMP_PRECISION matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
A_MATRIX_PREC **a_matrix(long nrl,long nrh,long ncl,long nch)
/* allocate a COMP_PRECISION matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  A_MATRIX_PREC **m;
  
  /* allocate pointers to rows */
  m=(A_MATRIX_PREC **) 
    malloc((unsigned int)((nrow+NR_END)*
			  sizeof(A_MATRIX_PREC*)));
  if (!m) {
    fprintf(stderr,"allocation failure 1 in matrix()\n");
    exit(-1);}
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(A_MATRIX_PREC *) 
    calloc(sizeof(A_MATRIX_PREC),
	   (unsigned int)((nrow*ncol+NR_END)));
  if (!m[nrl]) {
    fprintf(stderr,"allocation failure 2 in matrix()\n");
    exit(-1);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

void free_a_matrix(A_MATRIX_PREC **m,
		 long nrl,long nrh,long ncl,long nch)
/* free a COMP_PRECISION matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

int *ivector(long nl,long nh)
     /* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;
  
  v=(int *)calloc(sizeof(int),
		  (unsigned int) ((nh-nl+1+NR_END)));
  if (!v){
    fprintf(stderr,"allocation failure in ivector()\n");
    exit(-1);
  }
  return v-nl+NR_END;
}
void free_ivector(int *v,long nl,long nh)
     /* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}
