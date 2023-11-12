#
# calculate Gutenberg-Richter type statistics 
#
# $(col) is the moment
#
# if cum=1, will do cumulative statistics
#
# there is also a C version
#
#
BEGIN{
    determine_range1 = 0;
    determine_range2 = 0;
    if(xmin == ""){
	xmin= 10;
	determine_range1 = 1;
    }
    if(xmax == ""){
	xmax = -10;
	determine_range2 = 1;
    }
    
# number of boxes
  if(nbox==0)
    nbox=20;
# column
  if(col==0)
    col=1;
  n=0;
# cumulative statistics?
# set cum to unity
  
}
{
# read in data and determine bounds on moments
    if(substr($1,1,1)!="#"){
	n++;
	x[n] = m02mag($col);		# convert to mag
	if(determine_range1)
	    if(x[n] < xmin)
		xmin = x[n];
	if(determine_range2)
	    if(x[n] > xmax)
		xmax = x[n];
    }
}
END{
    xrange=xmax - xmin;
    print("# magnitudes min:",xmin," max:",xmax,
	  " range:",xrange," n:",n);
    
    dx = xrange/(nbox);
# set up boxes
    for(i=1;i <= nbox;i++){
	bb[i] = xmin + dx*(i-1);		# mid bin location
	bm[i] = bb[i] + dx/2;	# left boundary
	be[i] = bb[i] + dx;	# right boundary
	
	nb[i]=0;
    }
# sort the log(m_0) array
    asort(x);
    if(cum == 1){
# cumulative statistics
	nb[1] = nsum = n;
	i=1;j=2;
	while((i <= n)&&(j<=nbox)){
	    while((x[i] < be[j])&&(i<=n)){
		i++;
	    }
	    nb[j] = n-i-1;
	    if(nb[j] < 0)
		nb[j] = 0;
	    j++;
	}
# print them
	print("# cumulative stats: M n(M)/N n(M)");
	for(i=1;i<=nbox;i++)
	    printf("%e %e %i\n",mag2m0(bm[i]),nb[i]/n,nb[i]);
    }else{
# non-cumulative statistics 
# count
	j=1;k=0;
	for(i=1;i<=n;i++){
	    if(x[i] >= be[j]){
		nb[j]=k;
		j++;k=0;nlast=j;
	    }
	    k++;
	}
	nb[nlast] = k;

# print 
	print("# non-cumulative stats: M n(M)/N n(M)");
	k=0;
	for(i=1;i<=nbox;i++){
	    printf("%e %e %i\n",mag2m0(bm[i]),nb[i]/n,nb[i]);
	    k+=nb[i];
	}
    }
}

function m02mag(m0){
    mag = 2./3. * (0.4342944819032518*log(m0)-9.1);
    return mag;
}
function mag2m0(mag){
    m0 = 10**(3./2.*mag + 9.1);
    return m0;
}