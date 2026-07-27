/*--------------------------------------------------------------------

modified from GMT's project:

routine works like project -Cclon/clat -Aazimuth -Fpq -Q


original comments follow

 *	project.c,v 1.4.4.6 2002/09/27 19:02:08 pwessel Exp
 *
 *	Copyright (c) 1991-2002 by P. Wessel and W. H. F. Smith
 *	See COPYING file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 of the License.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * project.c
 * reads (x,y,[z]) data and writes some combination of (x,y,z,p,q,u,v),
 * where p,q is the distance along,across the track of the projection of (x,y),
 * and u,v are the un-transformed (x,y) coordinates of the projected position.
 * Can also create (x,y) along track.
  
   Author: 	Walter H. F. Smith
   Date:	19 April, 1988.
   Modified:	4 December 1988, to be more flexible.
   Complete rebuild 22 June, 1989 to use vector products and do more things.
   version 2.0
   		23-FEB-1998	PW: Added support for multiple files, multi-segment formats
				and binary i/o.  Old -M renamed -Q.
		03-NOV-1998	PW: Can read any number of data columns; z in -Fz refers to
				all these columns in the output.
   Version:	3.4		PW: Fixed problem with small circle distances
   Version:	3.4.2
*/
#ifdef USE_GEOPROJECT
#include "gmt.h"
#include "myprojectsimple.h"
void my_sincos(COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION );



void myprojectsimple(COMP_PRECISION *xin, COMP_PRECISION *xout,COMP_PRECISION clon,
		     COMP_PRECISION clat, COMP_PRECISION azimuth, int init)
{	
  static COMP_PRECISION	a[3], b[3], xt[3], pole[3], center[3],x_b, y_b;
  static COMP_PRECISION d_to_km;

  //GMT_geographic_in = GMT_geographic_out = TRUE;
  GMT_io.in_col_type[GMT_X] = GMT_io.out_col_type[GMT_X] = GMT_IS_LON;
  GMT_io.in_col_type[GMT_Y] = GMT_io.out_col_type[GMT_Y] = GMT_IS_LAT;
 

  if(init){
    //d_to_km = 0.001 * 2.0 * M_PI * gmtdefs.ellipse[N_ELLIPSOIDS-1].eq_radius / 360.0;
    d_to_km = 0.001 * 2.0 * M_PI * gmtdefs.ref_ellipsoid[gmtdefs.ellipsoid].eq_radius / 360.;
    sphere_project_setup(clat, clon, a, y_b, x_b, b, azimuth, pole, center, FALSE);
  }
  oblique_transform(xin[1], xin[0], &xt[1], &xt[0], pole, center);

  /* At this stage, all values are still in degrees.  */
  xout[0] = xt[0] * d_to_km;
  xout[1] = xt[1] * d_to_km;
  xout[2] = xin[2];
}

COMP_PRECISION oblique_setup (COMP_PRECISION plat, COMP_PRECISION plon,
			      COMP_PRECISION *p, COMP_PRECISION clat,
			      COMP_PRECISION clon, COMP_PRECISION *c, GMT_LONG c_given)
{
	/* routine sets up a unit 3-vector p, the pole of an 
	   oblique projection, given plat, plon, the position 
	   of this pole in the usual coordinate frame.
	   c_given = TRUE means that clat, clon are to be used
	   as the usual coordinates of a point through which the
	   user wants the central meridian of the oblique
	   projection to go.  If such a point is not given, then
	   the central meridian will go through p and the usual
	   N pole.  In either case, a unit 3-vector c is created
	   which is the directed normal to the plane of the central
	   meridian, pointing in the positive normal (east) sense.
	   Latitudes and longitudes are in degrees. */

	COMP_PRECISION cp, sin_lat_to_pole;
	double pd[3],platd,plond,cd[3],sd[3];
	int i;
	platd = (double)plat;
	plond = (double)plon;
	for(i=0;i<3;i++){
	  pd[i] = p[i];
	  cd[i] = c[i];
	}
	sd[0] = sd[1] = 0.0;
	sd[2] = -1.0;
	
	GMT_geo_to_cart(platd, plond, pd, TRUE);

	if (c_given) {	/* s points to user's clat, clon  */
	  GMT_geo_to_cart(clat, clon, sd, TRUE);
	}
	GMT_cross3v(pd, sd, cd);
	GMT_normalize3v(cd);
	cp = GMT_dot3v (pd, sd);
	sin_lat_to_pole = d_sqrt (1.0 - cp * cp);
	return (sin_lat_to_pole);
}

void make_euler_matrix (double *p, double *e, double theta)
{
	/* Routine to fill an euler matrix e with the elements
	   needed to rotate a 3-vector about the pole p through
	   an angle theta (in degrees).  p is a unit 3-vector.
	   Latitudes and longitudes are in degrees. */

	double	cos_theta, sin_theta, one_minus_cos_theta;
	double pxsin, pysin, pzsin, temp;

	sincosd(theta,&sin_theta,&cos_theta);
	one_minus_cos_theta = 1.0 - cos_theta;

	pxsin = p[0] * sin_theta;
	pysin = p[1] * sin_theta;
	pzsin = p[2] * sin_theta;

	temp = p[0] * one_minus_cos_theta;
	e[0] = temp * p[0] + cos_theta;
	e[1] = temp * p[1] - pzsin;
	e[2] = temp * p[2] + pysin;

	temp = p[1] * one_minus_cos_theta;
	e[3] = temp * p[0] + pzsin;
	e[4] = temp * p[1] + cos_theta;
	e[5] = temp * p[2] - pxsin;

	temp = p[2] * one_minus_cos_theta;
	e[6] = temp * p[0] - pysin;
	e[7] = temp * p[1] + pxsin;
	e[8] = temp * p[2] + cos_theta;
}


void	matrix_3v(double *a, double *x, double *b)
{
	/* routine to find b, where Ax = b, A is a 3 by 3 square matrix,
	   and x and b are 3-vectors.  A is stored row wise, that is:
	   
	   A = { a11, a12, a13, a21, a22, a23, a31, a32, a33 }  */
	
	b[0] = x[0]*a[0] + x[1]*a[1] + x[2]*a[2];
	b[1] = x[0]*a[3] + x[1]*a[4] + x[2]*a[5];
	b[2] = x[0]*a[6] + x[1]*a[7] + x[2]*a[8];
}

void	matrix_2v(double *a, double *x, double *b)
{
	/* routine to find b, where Ax = b, A is a 2 by 2 square matrix,
	   and x and b are 2-vectors.  A is stored row wise, that is:
	   
	   A = { a11, a12, a21, a22 }  */
	
	b[0] = x[0]*a[0] + x[1]*a[1];
	b[1] = x[0]*a[2] + x[1]*a[3];
}


void sphere_project_setup (COMP_PRECISION alat, COMP_PRECISION alon,
			   COMP_PRECISION *a, COMP_PRECISION blat,
			   COMP_PRECISION blon, COMP_PRECISION *b,
			   COMP_PRECISION azim, COMP_PRECISION *p, 
			   COMP_PRECISION *c, GMT_LONG two_pts)
{
	/* routine to initialize a pole vector, p, and a central meridian 
	   normal vector, c, for use in projecting points onto a great circle.
	   
	   The great circle is specified in either one of two ways:
	   if (two_pts), then the user has given two points, a and b,
	   which specify the great circle (directed from a to b);
	   if !(two_pts), then the user has given one point, a, and an azimuth,
	   azim, clockwise from north, which defines the projection.

	   The strategy is to use the oblique_transform operations above,
	   in such a way that the great circle of the projection is the
	   equator of an oblique transform, and the central meridian goes
	   through a.  Then the transformed longitude gives the distance
	   along the projection circle, and the transformed latitude gives
	   the distance normal to the projection circle.

	   If (two_pts), then p = normalized(a X b).  If not, we temporarily
	   create p_temp = normalized(a X n), where n is the north pole.
	   p_temp is then rotated about a through the angle azim to give p.
	   After p is found, then c = normalized(p X a).

	   Latitudes and longitudes are in degrees.
	*/

        double 	e[9];	/* Euler rotation matrix, if needed  */
	double ad[3],bd[3],pd[3],cd[3],alatd,alond,blatd,blond;
	int i;
	for(i=0;i<3;i++){
	  ad[i] = (double)a[i];
	  bd[i] = (double)b[i];
	  pd[i] = (double)p[i];
	}
	alond = (double)alon;
	alatd = (double)alat;
	
	blond = (double)blon;
	blatd = (double)blat;
	/* First find p vector  */

	if (two_pts) {
	  GMT_geo_to_cart(alatd, alond, ad, TRUE);
	  GMT_geo_to_cart(blatd, blond, bd, TRUE);
	  GMT_cross3v(ad, bd, pd);
	  GMT_normalize3v(pd);
	}
	else {
	  GMT_geo_to_cart(alatd, alond, ad, TRUE);
	  bd[0] = bd[1] = 0.0;	/* set b to north pole  */
	  bd[2] = 1.0;
	  GMT_cross3v(ad, bd, cd);	/* use c for p_temp  */
	  GMT_normalize3v(cd);
	  make_euler_matrix(ad, e, -azim);
	  matrix_3v(e, cd, pd);	/* c (p_temp) rotates to p  */
	}

	/* Now set c vector  */

	GMT_cross3v(pd, ad, cd);
	GMT_normalize3v(cd);
	for(i=0;i<3;i++)
	  c[i] = cd[i];
}

void oblique_transform (COMP_PRECISION xlat, COMP_PRECISION xlon,
			COMP_PRECISION *x_t_lat, COMP_PRECISION *x_t_lon,
			COMP_PRECISION *p, COMP_PRECISION *c)
{
	/* routine takes the point x at conventional (xlat, xlon) and
	   computes the transformed coordinates (x_t_lat, x_t_lon) in
	   an oblique reference frame specified by the unit 3-vectors
	   p (the pole) and c (the directed normal to the oblique
	   central meridian).  p and c have been computed earlier by
	   the routine oblique_setup().
	   Latitudes and longitudes are in degrees. */

	COMP_PRECISION	 temp1, temp2;
	double xlatd,xlond,xd[3],pd[3],p_cross_xd[3],cd[3];
	int i;
	xlatd = xlat;
	xlond = xlon;
	for(i=0;i<3;i++){
	  pd[i] = p[i];
	  cd[i] = c[i];
	}
	GMT_geo_to_cart(xlatd, xlond, xd, TRUE);

	temp1 = GMT_dot3v(xd,pd);
	*x_t_lat = (COMP_PRECISION)d_asind(temp1);

	GMT_cross3v(pd,xd,p_cross_xd);
	GMT_normalize3v(p_cross_xd);

	temp1 = GMT_dot3v(p_cross_xd, cd);
	temp2 = GMT_dot3v(xd, cd);
	*x_t_lon = (COMP_PRECISION) copysign(d_acosd(temp1), temp2);
}

#endif
