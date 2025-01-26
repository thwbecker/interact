/*

  interact: model fault interactions using dislocations in a 
            halfspace
  (C) Thorsten Becker, becker@eps.harvard.edu

  $Id: get_projected_fault_parameters.c,v 1.7 2003/03/18 03:00:09 becker Exp $



*/
#include "interact.h"
#include "blockinvert.h"
/*

  given fault end points in degrees in lon,lat format in fx,
  and the total depth tdepth (positive in km),
  
  will setup an oblique Mercator projection and 

  return the center[2] and azimuth of the projection (and fault trace)
  in degrees, as well as the projected fault half-length and
  half-width and depth, 

  the center of the fault in the projected system is (0,0,-z) if dip
  is not equal to 90, then the center of the patch has to be moved (as
  in read_bflt.c) OUTSIDE THIS ROUTINE
  


*/
void get_projected_fault_parameters(COMP_PRECISION fx[2][2],
				    COMP_PRECISION lock_depth,
				    COMP_PRECISION *center,
				    COMP_PRECISION *azi,
				    COMP_PRECISION *dip,
				    COMP_PRECISION *l, 
				    COMP_PRECISION *w,
				    COMP_PRECISION *depth)

{
  int i,j;
  COMP_PRECISION fxr[2][2],pfx[2][2],azir,
    centerr[2],gdistance,dummy=0;
  // convert fault endpoint coordinates to radians
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      fxr[i][j] = DEG2RADF(fx[i][j]);
  // get the center point along a great circle and distance
  // in radians
  get_point_on_gc(fxr[0][0],fxr[0][1],fxr[1][0],fxr[1][1],
		  0.5, centerr,(centerr+1),&gdistance);
  // convert center from radians to degrees
  for(i=0;i<2;i++)
    center[i] = RAD2DEGF(centerr[i]);
  /*

    get forward azimuth 

  */
  azir = azimuth(fxr[0][0],fxr[0][1],fxr[1][0],fxr[1][1]); 
  // convert to degrees
  *azi = RAD2DEGF(azir);
  fix_azimuth(azi); // limit to 0 .. 360
  if(*azi >= 180.0)		/* limit further to 0 <= a < 180 */
    *azi = *azi - 180.0;
  //
  // convert fault end points to projected system using fault
  // center and azimuth
  //
  for(i=0;i<2;i++)
    geoproject(&fx[i][0], &pfx[i][0],FLT_ROT_PROJECTION,
	       center[0],center[1],*azi,
	       dummy,dummy,dummy, dummy, (int)FALSE);
  // determine projected half fault length
  *l = 0.5 * distance(&pfx[0][0],&pfx[1][0],2);
  //
  // determine half width and mid-fault depth
  if(fabs(*dip) < EPS_COMP_PREC){
    fprintf(stderr,"get_projected_fault_parameters: error: dip: %g\n",
	    *dip);
    exit(-1);
  }
  *w = lock_depth/(2.0*sin((*dip)*DEG2RAD));
  // 
  // depth of fault mid point (>0)
  *depth = lock_depth/2.0;
}
