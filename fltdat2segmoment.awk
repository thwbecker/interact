#
BEGIN{
    grpmin=1e10;
    grpmax=-1e5;
}
{
    if(substr($1,1,1)!="#"){
	grp=$13+1;
	if(grp<grpmin)
	    grpmin=grp;
	if(grp>grpmax)
	    grpmax=grp;
	if(grp<1)print("error with group") > "/dev/stderr"
	mom[grp] += $3 * sqrt($6**2+$7**2+$8**2);
    }
}
END{
    for(i=grpmin;i<=grpmax;i++){
	if(mom[i]!=0)
	    printf("%5i %20.7e\n",i-1,mom[i]);
    }
}
