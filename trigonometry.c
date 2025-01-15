/*
  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, thwbecker@post.harvard.edu
  
  
  trigonometric formulae for computations of distance and
  such on the surface of a spherical earth

  here, lat refers to latitude (not colatitude)
  
*/
#include "interact.h"
#include <math.h>
/*

  compute the distance on a sphere in radians given location 
  lon1,lat and lon2,lat2 (in radians)

*/
COMP_PRECISION dist_on_sphere(COMP_PRECISION lon1, 
			      COMP_PRECISION lat1,
			      COMP_PRECISION lon2, 
			      COMP_PRECISION lat2)
{
  COMP_PRECISION tmp1,tmp2,tmp3;
  tmp1 = sin((lat1-lat2)/2.0);
  tmp1 *= tmp1;
  tmp2 = sin((lon2-lon1)/2.0);
  tmp2 *= tmp2;
  
  return 2.0*asin(sqrt(tmp1 + cos(lat1) * cos(lat2) * tmp2));
}
/* same for input in degrees */
COMP_PRECISION dist_on_sphere_deg(COMP_PRECISION dlon1, 
				  COMP_PRECISION dlat1,
				  COMP_PRECISION dlon2, 
				  COMP_PRECISION dlat2)
{
  return dist_on_sphere(dlon1*DEG2RAD,dlat1*DEG2RAD,
			dlon2*DEG2RAD,dlat2*DEG2RAD);
}
/*

  calculate the azimuth from point 1 to point 2 at point 1
  given the lon and lats in radians

  the returned azimuth is in radians between 0 .. 2pi

*/
COMP_PRECISION azimuth(COMP_PRECISION lon1, COMP_PRECISION lat1,
		       COMP_PRECISION lon2, COMP_PRECISION lat2)

{
  COMP_PRECISION slat1,clat1,slat2,clat2,azi;

  my_sincos(&slat1,&clat1,lat1);
  my_sincos(&slat2,&clat2,lat2);
  azi = atan2(sin(lon2-lon1)*clat2,
	      clat1*slat2-slat1*clat2*cos(lon1-lon2));
  if(azi < 0)
    azi += TWOPI;
  return azi;
}
/*


  calculate the pole of a great circle as defined by two points
  lon1,lat1 and lon2,lat2

  intput has to be in radians

*/
void get_gc_pole(COMP_PRECISION lon1, COMP_PRECISION lat1,
		 COMP_PRECISION lon2, COMP_PRECISION lat2, 
		 COMP_PRECISION *plon,COMP_PRECISION *plat)
{
  COMP_PRECISION xc1[3],xc2[3],xcp[3];
  // get pole by cross product
  lonlat2xyz(lon1,lat1,xc1);
  lonlat2xyz(lon2,lat2,xc2);
  // do cross product
  cross_product(xc1,xc2,xcp);
  // normalize
  normalize_3d(xcp);
  // go back to lon lat
  xyz2lonlat(xcp,plon,plat);
}
/*

calculate a point at distance fraction f between p1 and p2 along
a great circle. p1 and p2 are given as lon,lat pairs in radians

also returns the total distance between p1 and p2, d
*/
void get_point_on_gc(COMP_PRECISION lon1, COMP_PRECISION lat1,
		     COMP_PRECISION lon2, COMP_PRECISION lat2,
		     COMP_PRECISION f, 
		     COMP_PRECISION *lon,COMP_PRECISION *lat,
		     COMP_PRECISION *d)
{
  COMP_PRECISION slat1,clat1, slat2,clat2, 
    slon1,clon1, slon2,clon2,a,b,x[3],sd;
  if((fabs(lat1+lat2) < EPS_COMP_PREC) && 
     (fabs(fabs(lon1-lon2)-PI) < EPS_COMP_PREC)){
    fprintf(stderr,"get_point_on_gc: error: points antipodal: %g, %g and %g, %g\n",
	    RAD2DEGF(lon1),RAD2DEGF(lat1),
	    RAD2DEGF(lon2),RAD2DEGF(lat2));
    exit(-1);
  }
  *d = dist_on_sphere(lon1,lat1, lon2,lat2);
  sd = sin(*d);
  my_sincos(&slat1,&clat1,lat1);
  my_sincos(&slon1,&clon1,lon1);
  my_sincos(&slat2,&clat2,lat2);
  my_sincos(&slon2,&clon2,lon2);

  a = sin((1.0 - f) * (*d))/sd;
  b = sin(f * (*d))/sd;
  x[INT_X] =  a*clat1*clon1 +  b*clat2*clon2;
  x[INT_Y] =  a*clat1*slon1 +  b*clat2*slon2;
  x[INT_Z] =  a*slat1       +  b*slat2;
  xyz2lonlat(x,lon,lat);
}
/* 
  
calculate lon and lat of a point with distance d from lon1,lat1
and azimuth azi

all in radians

 */
void get_point_on_course(COMP_PRECISION lon1, COMP_PRECISION lat1,
			 COMP_PRECISION d, COMP_PRECISION azi,
			 COMP_PRECISION *lon,COMP_PRECISION *lat)
{
  COMP_PRECISION dlon,sazi,cazi,sd,cd,slat1,clat1;
  my_sincos(&sazi,&cazi,azi);
  my_sincos(&sd,&cd,d);
  my_sincos(&slat1,&clat1,lat1);
  
  *lat = asin(slat1 * cd + clat1 * sd * cazi);
  dlon = atan2(sazi * sd * clat1, cd - slat1 * sin(*lat));
  *lon = -(fmod(-lon1 - dlon + PI, TWOPI) -PI);
}
/* 

   convert from lon - lat to x y z assuming radius = 1
   and lon lat are in radians

*/
void lonlat2xyz(COMP_PRECISION lon, COMP_PRECISION lat,
		COMP_PRECISION *xc)
{
  COMP_PRECISION clon,slon,clat,slat;
  my_sincos(&slon,&clon,lon);
  my_sincos(&slat,&clat,lat);
  xc[INT_X]=clat * clon;
  xc[INT_Y]=clat * slon; 
  xc[INT_Z]=slat;
}
void lonlat2xyz_deg(COMP_PRECISION lon, COMP_PRECISION lat,
		    COMP_PRECISION *xc)
{
  COMP_PRECISION clon,slon,clat,slat;
  my_sincos(&slon,&clon,DEG2RADF(lon));
  my_sincos(&slat,&clat,DEG2RADF(lat));
  xc[INT_X]=clat * clon;
  xc[INT_Y]=clat * slon; 
  xc[INT_Z]=slat;
}
/*
  
  go from xyz to lon, lat (radians)

*/
void xyz2lonlat(COMP_PRECISION *xc, COMP_PRECISION *lon,
		COMP_PRECISION *lat)
{
  *lat = atan2(xc[INT_Z],hypot(xc[INT_X],xc[INT_Y]));
  *lon = atan2(xc[INT_Y],xc[INT_X]);
}
// same in degrees
void xyz2lonlat_deg(COMP_PRECISION *xc, COMP_PRECISION *lon,
		    COMP_PRECISION *lat)
{
  *lat = atan2(xc[INT_Z],hypot(xc[INT_X],xc[INT_Y]));
  *lon = atan2(xc[INT_Y],xc[INT_X]);
  *lat = RAD2DEGF(*lat);
  *lon = RAD2DEGF(*lon);
}

/* convert a spherical system vector to cartesian given 
   existing basis vectors [9] */
void pv2cv(COMP_PRECISION *xp,COMP_PRECISION *xc,
	   COMP_PRECISION *polar_base)
{
  int i;
  // convert vector
  for(i=0;i<3;i++){
    xc[i]  = polar_base[INT_R*3+i]     * xp[INT_R]; /* r contribution */
    xc[i] += polar_base[INT_THETA*3+i] * xp[INT_THETA]; /* theta contribution */
    xc[i] += polar_base[INT_PHI*3+i]   * xp[INT_PHI]; /* phi  contribution */
  }
}

/* convert a cartesian system vector to spherical given 
   existing basis vectors [9] */
void cv2pv(COMP_PRECISION *xc,COMP_PRECISION *xp,
	   COMP_PRECISION *polar_base)
{
  int i,j;
  // convert vector
  for(i=j=0;i<3;i++,j+=3){
    xp[i]  = polar_base[j+INT_X] * xc[INT_X]; /* x contribution */
    xp[i] += polar_base[j+INT_Y] * xc[INT_Y]; /* y contribution */
    xp[i] += polar_base[j+INT_Z] * xc[INT_Z]; /* z contribution */
  }
}


//
// given a location specified as lon lat (in degrees), calculate
// the three spherical basis vectors pr, ptheta, and phi
// polar_base is [9]
//
void calculate_polar_base(COMP_PRECISION lon, COMP_PRECISION lat,
			  COMP_PRECISION *polar_base)
{
  COMP_PRECISION theta, phi,st,ct,sp,cp;
  theta = DEG2RADF(90.0-lat);
  phi =   DEG2RADF(lon);
  // calculate sin and cos of theta and phi
  my_sincos(&st,&ct,theta);
  my_sincos(&sp,&cp,phi);
  // r base vector, R*3+i
  polar_base[INT_R*3+INT_X]= st * cp;
  polar_base[INT_R*3+INT_Y]= st * sp;
  polar_base[INT_R*3+INT_Z]= ct;
  // theta base vector, THETA*3+i
  polar_base[INT_THETA*3+INT_X]= ct * cp;
  polar_base[INT_THETA*3+INT_Y]= ct * sp;
  polar_base[INT_THETA*3+INT_Z]= -st;
  // phi base vector, PHI*3+i
  polar_base[INT_PHI*3+INT_X]= -sp;
  polar_base[INT_PHI*3+INT_Y]=  cp;
  polar_base[INT_PHI*3+INT_Z]= 0.0;
}
