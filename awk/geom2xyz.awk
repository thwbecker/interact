# convert geom.in file into projected fault coordinates
# for plotting with GMT 
# parameters i, j pick the coordinates for projection
#
# $Id: geom2xyz.awk,v 1.1 2003/01/17 01:12:17 becker Exp $
# 
BEGIN{
  f=0.017453293;
  fl=0.3;
  if(i=="")
    i=1;
  if(j=="")
    j=2;
}
{
# read in patch geometry
  x[1]=$1;x[2]=$2;x[3]=$3;

  alpha=90.0-$4;
  alpha *= f;
    
  sin_alpha=sin(alpha);
  cos_alpha=cos(alpha);
  dip=$5*f;
  sin_dip=sin(dip);
  cos_dip=cos(dip);
  len=$6;
  width=$7;

# base vectors
  
  t_strike[1]=  cos_alpha;
  t_strike[2]=  sin_alpha;
  t_strike[3]=  0.0;
  
  t_dip[1]=  -sin_alpha * cos_dip;
  t_dip[2]=   cos_alpha * cos_dip;
  t_dip[3]=   sin_dip;
  

#  calculate corners
  for(l=1;l<=3;l++){
    sx[l]=t_strike[l]*len;
    dx[l]=t_dip[l]   *width;
  }
  for(fac=1.0;fac<=1.0;fac+=0.1){
    for(l=1;l<=4;l++){
      corner1[l]=x[l]-fac*sx[l]-fac*dx[l];
      corner2[l]=x[l]+fac*sx[l]-fac*dx[l];
      corner3[l]=x[l]+fac*sx[l]+fac*dx[l];
      corner4[l]=x[l]-fac*sx[l]+fac*dx[l];
    }
    printf(">\n");
    printf("%g %g\n",corner1[i],corner1[j]);
    printf("%g %g\n",corner2[i],corner2[j]);
    printf("%g %g\n",corner3[i],corner3[j]);
    printf("%g %g\n",corner4[i],corner4[j]);
    printf("%g %g\n",corner1[i],corner1[j]);
    printf(">\n");
  }
}
