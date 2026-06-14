#
# linear regression routine
#
# $Id: linreg.awk,v 1.1 2003/01/17 01:12:42 becker Exp becker $
#
 BEGIN{

  n=0;

  if(xlimit == 0)xlimit=9e99;
#  printf("Linear regression, xlimit=%g\n\n",xlimit);
}
{
  if(($1 != "")&&((substr($1,1,1)!="#"))&&($2 != "")&&($3 != -1)&&($1 <= xlimit))
    {
	n++;
	x[n]=$1;
	y[n]=$2;

    };
}
END{
  if(n == 0){
      printf("No matching data points found.\n\n");
  }else{
#      printf("Used %g data pairs.\n\n",n);
      linreg(x,y,n,a);

# outpuf of a b and dev, where y = a + b *x
      printf("%20.10e %20.10e %20.10e\n",a[1],a[2],a[3]);
    }
}

#
# given n x[] and y[] values, compute linear regression and return 
# offset a[1] slope a[2] and misfit a[3]
# 
function linreg(x,y,n,a)
{

    xs=0.0;ys=0.0;
    x2s=0.0;xys=0.0;
    for(i=1;i<=n;i++){
	xs += x[i];
	ys += y[i];
	xys += (x[i]*y[i]);
	x2s += (x[i]*x[i]);
    }
    xm = xs/n;
    ym = ys/n;
    tmp = x2s - xm * xs;
    if(tmp != 0){
	a[2] = (xys - ym * xs) / (x2s - xm * xs);
    }else{
	a[2] = "nan";
    }
    a[1] = ym - a[2] * xm;
# printf("y = x * %g + %g\n\n",a[2],a[1]);
    
    # misfit 
    a[3] = 0.;
    for(i=1;i<=n;i++){
	tp = x[n]*a[2]+a[1];
	a[3] += (y[n]-tp)**2;
    }
    a[3]  = sqrt(a[3]/ n);
}  