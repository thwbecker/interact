BEGIN{
    pi = 3.1415926535897932384626433832795;
    pif=180/pi;
    if(nx=="")
	nx=15;
    if(r=="")
	r=50;
    if(athetaf=="")
	athetaf=10;
    if(nzfac=="")
	nzfac=1;
    print(nx,r,athetaf) > "/dev/stderr"

    aphi=pi/2;
    atheta=pi/2/athetaf;
    
    dazi=pi/2/(nx*2);
    dtheta=atheta/nx/nzfac;
    for(azi=dazi/2;azi < aphi+1e-6;azi+=dazi){
	for(theta=pi/2+dtheta;theta < pi/2+atheta;theta+=dtheta){
	    
	    x0=r*sin(theta)*cos(azi);
	    y0=r*sin(theta)*sin(azi);
	    z0=r*cos(theta);
	    i=1;
	    for(a=azi-dazi/2;a<azi+dazi/2+1e-6;a+=dazi){
		for(t=theta-dtheta/2;t<theta+dtheta/2+1e-6;t+=dtheta){
		    x[i]=r*sin(t)*cos(a);
		    y[i]=r*sin(t)*sin(a);
		    z[i]=r*cos(t);
#		    print(i,a,t,x[i],y[i],z[i])
		    i++;
		}
	    }
	    l[1]=sqrt((x[3]-x[1])**2+(y[3]-y[1])**2);
	    l[2]=sqrt((x[4]-x[2])**2+(y[4]-y[2])**2);
	    lm=(l[1]+l[2])/2;
	    w[1]=sqrt((z[2]-z[1])**2);
	    w[2]=sqrt((z[4]-z[3])**2);
	    wm=(w[1]+w[2])/2;
	    
	    print(x0,y0,z0,(-azi)*pif,theta*pif,lm/2,wm/2,0)
	}
    }

}

