#include <math.h>
#include <stdio.h>
#include "interact.h"
#include "period.h"




#define TWOPID 6.2831853071795865


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

static SPECA_PREC _tmp_sqrarg;
#define SQR(a) ((_tmp_sqrarg=(a)) == 0.0 ? 0.0 : _tmp_sqrarg*_tmp_sqrarg)

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static int _tmp_lmaxarg1,_tmp_lmaxarg2;
#define LMAX(a,b) (_tmp_lmaxarg1=(a),_tmp_lmaxarg2=(b),(_tmp_lmaxarg1) > (_tmp_lmaxarg2) ?\
        (_tmp_lmaxarg1) : (_tmp_lmaxarg2))

static int _tmp_lminarg1,_tmp_lminarg2;
#define LMIN(a,b) (_tmp_lminarg1=(a),_tmp_lminarg2=(b),(_tmp_lminarg1) < (_tmp_lminarg2) ?\
        (_tmp_lminarg1) : (_tmp_lminarg2))

#define NR_MOD(a,b) 	while(a >= b) a -= b;

void period(SPECA_PREC *x,SPECA_PREC *y,int n,SPECA_PREC ofac,SPECA_PREC hifac,SPECA_PREC *px,
	    SPECA_PREC *py,int np,int *nout,int *jmax,SPECA_PREC *prob)
{
  int i,j,n1;
  SPECA_PREC ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
    sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy;
  double arg,wtemp,*wi,*wpi,*wpr,*wr;

  n1=n+1;
  wi=(double *)calloc(n1,sizeof(double));
  wpi=(double *)calloc(n1,sizeof(double));
  wpr=(double *)calloc(n1,sizeof(double));
  wr=(double *)calloc(n1,sizeof(double));
  if(!wi || !wpi || !wpr || !wr)MEMERROR("period");

  *nout=0.5*ofac*hifac*n;

  if (*nout > np){fprintf(stderr,"output arrays too short in period\n");exit(-1);}
  avevar(y,n,&ave,&var);

  xmax=xmin=x[1];
  for (j=1;j<=n;j++) {
    if (x[j] > xmax) xmax=x[j];
    if (x[j] < xmin) xmin=x[j];
  }
  xdif=xmax-xmin;
  xave=0.5*(xmax+xmin);
  pymax=0.0;
  pnow=1.0/(xdif*ofac);
  for (j=1;j<=n;j++) {
    arg=TWOPID*((x[j]-xave)*pnow);
    wpr[j] = -2.0*SQR(sin(0.5*arg));
    wpi[j]=sin(arg);
    wr[j]=cos(arg);
    wi[j]=wpi[j];
  }
  for (i=1;i<=(*nout);i++) {
    px[i]=pnow;
    sumsh=sumc=0.0;
    for (j=1;j<=n;j++) {
      c=wr[j];
      s=wi[j];
      sumsh += s*c;
      sumc += (c-s)*(c+s);
    }
    wtau=0.5*atan2(2.0*sumsh,sumc);
    swtau=sin(wtau);
    cwtau=cos(wtau);
    sums=sumc=sumsy=sumcy=0.0;
    for (j=1;j<=n;j++) {
      s=wi[j];
      c=wr[j];
      ss=s*cwtau-c*swtau;
      cc=c*cwtau+s*swtau;
      sums += ss*ss;
      sumc += cc*cc;
      yy=y[j]-ave;
      sumsy += yy*ss;
      sumcy += yy*cc;
      wr[j]=((wtemp=wr[j])*wpr[j]-wi[j]*wpi[j])+wr[j];
      wi[j]=(wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
    }
    py[i]=0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;
    if (py[i] >= pymax) pymax=py[(*jmax=i)];
    pnow += 1.0/(ofac*xdif);
  }
  expy=exp(-pymax);
  effm=2.0*(*nout)/ofac;
  *prob=effm*expy;
  if (*prob > 0.01) *prob=1.0-pow(1.0-expy,effm);
  free(wr);
  free(wpr);
  free(wpi);
  free(wi);
}
#undef TWOPID

void avevar(SPECA_PREC *data,int n,SPECA_PREC *ave,SPECA_PREC *var)
{
  int j;
  SPECA_PREC s,ep;
  
  for (*ave=0.0,j=1;j<=n;j++) *ave += data[j];
  *ave /= n;
  *var=ep=0.0;
  for (j=1;j<=n;j++) {
    s=data[j]-(*ave);
    ep += s;
    *var += s*s;
  }
  *var=(*var-ep*ep/n)/(n-1);
}


#define MACC 4

void fasper(SPECA_PREC *x,SPECA_PREC *y,int n,SPECA_PREC ofac,SPECA_PREC hifac,
	    SPECA_PREC *wk1,SPECA_PREC *wk2,int nwk,int *nout, int *jmax,
	    SPECA_PREC *prob)
{
  int j,k,ndim,nfreq,nfreqt;
  SPECA_PREC ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt;
  SPECA_PREC hs2wt,hypo,pmax,sterm,swt,var,xdif,xmax,xmin;
  
  *nout=0.5*ofac*hifac*n;
  nfreqt=ofac*hifac*n*MACC;
  nfreq=64;
  while (nfreq < nfreqt) nfreq <<= 1;
  ndim=nfreq << 1;
  if (ndim > nwk) {fprintf(stderr,"workspaces too small in fasper, %i versus %i\n",
			   ndim,nwk);
  exit(-1);}
  avevar(y,n,&ave,&var);
  xmin=x[1];
  xmax=xmin;
  for (j=2;j<=n;j++) {
    if (x[j] < xmin) xmin=x[j];
    if (x[j] > xmax) xmax=x[j];
  }
  xdif=xmax-xmin;
  for (j=1;j<=ndim;j++) wk1[j]=wk2[j]=0.0;
  fac=ndim/(xdif*ofac);
  fndim=ndim;
  for (j=1;j<=n;j++) {
    ck=(x[j]-xmin)*fac;
    NR_MOD(ck,fndim)
      ckk=2.0*(ck++);
    NR_MOD(ckk,fndim)
      ++ckk;
    spread(y[j]-ave,wk1,ndim,ck,MACC);
    spread(1.0,wk2,ndim,ckk,MACC);
  }
  realft(wk1,ndim,1);
  realft(wk2,ndim,1);
  df=1.0/(xdif*ofac);
  pmax = -1.0;
  for (k=3,j=1;j<=(*nout);j++,k+=2) {
    hypo=sqrt(wk2[k]*wk2[k]+wk2[k+1]*wk2[k+1]);
    hc2wt=0.5*wk2[k]/hypo;
    hs2wt=0.5*wk2[k+1]/hypo;
    cwt=sqrt(0.5+hc2wt);
    swt=SIGN(sqrt(0.5-hc2wt),hs2wt);
    den=0.5*n+hc2wt*wk2[k]+hs2wt*wk2[k+1];
    cterm=SQR(cwt*wk1[k]+swt*wk1[k+1])/den;
    sterm=SQR(cwt*wk1[k+1]-swt*wk1[k])/(n-den);
    wk1[j]=j*df;
    wk2[j]=(cterm+sterm)/(2.0*var);
    if (wk2[j] > pmax) pmax=wk2[(*jmax=j)];
  }
  expy=exp(-pmax);
  effm=2.0*(*nout)/ofac;
  *prob=effm*expy;
  if (*prob > 0.01) *prob=1.0-pow(1.0-expy,effm);
}


void spread(SPECA_PREC y,SPECA_PREC *yy,int n,SPECA_PREC x, int m)
{
  int ihi,ilo,ix,j,nden;
  static int nfac[11]={0,1,1,2,6,24,120,720,5040,40320,362880};
  SPECA_PREC fac;
  
  if (m > 10){fprintf(stderr,"factorial table too small in spread\n");exit(-1);}
  ix=(int)x;
  if (x == (SPECA_PREC)ix) yy[ix] += y;
  else {
    ilo=LMIN(LMAX((long)(x-0.5*m+1.0),1),n-m+1);
    ihi=ilo+m-1;
    nden=nfac[m];
    fac=x-ilo;
    for (j=ilo+1;j<=ihi;j++) fac *= (x-j);
    yy[ihi] += y*fac/(nden*(x-ihi));
    for (j=ihi-1;j>=ilo;j--) {
      nden=(nden/(j+1-ilo))*(j-ihi);
      yy[j] += y*fac/(nden*(x-j));
    }
  }
}

void realft(SPECA_PREC *data,int n,int isign)
{
  int i,i1,i2,i3,i4,np3;
  SPECA_PREC c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  
  theta=3.141592653589793/(double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
}


void four1(SPECA_PREC *data,int nn,int isign)
{
  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  SPECA_PREC tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

#undef SWAP
#undef LMIN
#undef LMAX
#undef NR_MOD
#undef MACC
#undef SIGN
#undef SQR
