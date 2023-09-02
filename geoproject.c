#include "interact.h"
#include "geoproject.h"

/*

WARNING: this routines assumes there is a z-value!

routine to convert geographic input coordinates (lon lat z in degrees
and AU, in[3]) to cartesian, projected coordinates (x y z, out[3], z
just gets copied), or vice-versa. functionality is similar to GMTs
(3.4.2) mapproject but only a limited set of projections are supported

for inverse=false:

go from geographic coordinates to x-y cartesian 

or for inverse=true:

the reverse, from cartesian coordinate to lon lat
the lon-lat output will be in the -180 .. 180 convention


the cartesian units are set to kilometers, and the mapprojection
center is the origin


input: 

PORJECT_AZI: simple projecion given a center and azimuth
             this is probably the only one that really works!!!!

OMERC_AZI: oblique Mercator projection

lon, lat and azimuth (all in degrees) of projection center and azimuth for
oblique Mercator

OMERC_POLE: oblique Mercator projection

lon, lat and of projection center and plon plat of projection pole
(all in degrees)

LCONFORM: Lambert conic conformal 

lon,lat of center and lat1 and lat2 of two standard parallels


this code is basically a clone of the relevant parts of GMTs (3.4.2)
mapproject.c sources, see there for appropriate copyright, comments,
etc


$Id: geoproject.c,v 1.2 2003/02/13 22:45:12 becker Exp $


*/

// this switch is in here for using this routine with the 
// interact package
#ifdef USE_GEOPROJECT

					 /* we cannot include gmt.h because of defined X there.... */


void geoproject(double *in, double *out, int projection,
		double lon, double lat, double azi,
		double plon, double plat, double lat1,
		double lat2, int inverse)
{
  static my_boolean init=FALSE;
  static int oldprojection,unit;
  static double fwd_scale, inv_scale, inch_to_unit, 
    unit_to_inch,oldazi,oldlon,oldlat,
    west = 0.0, east = 1.0, south = 0.0, north = 1.0;
  static char unit_name[80],**dummy;  
  // changing (variables)
  char projection_string[200];  
  double u_scale;
  if(!init){	  
    /* 
       called for first time, initialize GMT stuff 
    */
    dummy=(char **)malloc(sizeof(char*));
    dummy[0]=(char *)malloc(sizeof(char)*200);
    strcpy(dummy[0],"dummy\0");
    GMT_begin (1, dummy);
    /* 
       scaling 
    */
    unit = 1.0; strcpy(unit_name,"km");// kilometers
    init = TRUE;
    // make sure we set up the projection at least once
    oldprojection = -1;oldlat = 100;oldlon=400;oldazi=500;
  }
  if((azi != oldazi) || (lat != oldlat) || (lon != oldlon) ||
     (projection != oldprojection)){
    /* 
       projection has changed: 
       initialize the projection parameters by
       assigning a string in GMT -J style 
    */
    switch(projection){
    case PROJECT_AZI:
      if(inverse){
	fprintf(stderr,"geoproject: error: inverse not implemented for simple project\n");
	exit(-1);
      }
      if(azi < 0)
	azi += 360.0;
      if(azi >= 360)
	azi -= 360.0;
      myprojectsimple(in,out,lon,lat,azi,1); /* initialize */
      break;
    case OMERC_AZI:
      /* oblique Mercator with center 
	 at lon,lat and azimuth azi */
      if(lat == 0){
	fprintf(stderr,"geoproject: OMERC_AZI: lat=0 no good\n");
	exit(-1);
      }
      // go to (0, 180]
      if(azi > 180){
	fprintf(stderr,"geoproject: OMERC_AZI: azi > 180 no good\n");
	exit(-1);
      }

      sprintf(projection_string,"Oa%024.15e/%024.15e/%024.15e/1",
	      lon,lat,azi);
      break;
    case OMERC_POLE:
      /* oblique Mercator with center 
	 at lon,lat and projection pole 
	 at plon,plat: ption -Joc or -JOc */
      sprintf(projection_string,"Oc%024.16e/%024.16e/%024.16e/%024.16e/1",
	      lon,lat,plon,plat);
      break;
    case LCONFORM:		/* Lamber conic conformal */
      sprintf(projection_string,"L%024.16e/%024.16e/%024.16e/%024.16e/1",
	     lon,lat,lat1,lat2);
      break;
    default:
      fprintf(stderr,"geoproject: error: projection %i undefined\n",
	      projection);
      exit(-1);
      break;
    }
    if(projection != PROJECT_AZI){
      /* 
	 
      evaluate the projection string for other projections

      */
      //fprintf(stderr,"geoproject: setting up: -J%s\n",projection_string);
      /* 
	 obtain the projection parameters
      */
      fprintf(stderr,"not implemented yet GMT_map_getproject in geoproject\n");
      exit(-1);

    /*   if(GMT_map_getproject (projection_string)){ */
/* 	fprintf(stderr,"geoproject: GMT projection error: -J%s\n", */
/* 		projection_string); */
/* 	exit(-1); */
/*       } */

      //fprintf(stderr,"geoproject: setting up projection -J%s\n",
      //projection_string);
      /* we move origin to the projection center as in mapproject -C  */
      /* we switch on one-to-one scaling as in mapproject -F */
      /* set up scaling */
      GMT_init_scales(unit, &fwd_scale, &inv_scale, &inch_to_unit, 
		      &unit_to_inch, unit_name);
      /* set up mapping */
      GMT_map_setup(west, east, south, north);
      //fprintf(stderr,"%g %g %g %g\n",fwd_scale,inv_scale,inch_to_unit,unit_to_inch);

    }
    // set for next time
    oldazi = azi;oldlon = lon;oldlat = lat;
    oldprojection = projection;
  }
  if(projection != PROJECT_AZI){
    /* 
       settings that depend on "inverse" switch 
    */
    u_scale = (inverse) ? inv_scale : fwd_scale;  
    //GMT_geographic_in = (inverse) ? FALSE : TRUE;
    //GMT_geographic_out = !GMT_geographic_in;
    if(inverse){	
      /* 
	 DO INVERSE TRANSFORMATION 
      */
      if (unit) {
	in[INT_X] *= u_scale;
	in[INT_Y] *= u_scale;
      }
      in[INT_X] *= project_info.x_scale;
      in[INT_Y] *= project_info.y_scale;
      in[INT_X] += project_info.x0;
      in[INT_Y] += project_info.y0;
      GMT_xy_to_geo (&out[INT_X], &out[INT_Y], in[INT_X], in[INT_Y]);
      // go to -180 ... 180 system
      if(out[INT_X]>180)
	out[INT_X] -= 360;
    } else {		
      /* 
	 DO FORWARD TRANSFORMATION 
      */
      GMT_geo_to_xy (in[INT_X], in[INT_Y], &out[INT_X], &out[INT_Y]);
      out[INT_X] -= project_info.x0;
      out[INT_Y] -= project_info.y0;
      out[INT_X] /= project_info.x_scale;
      out[INT_Y] /= project_info.y_scale;
      if (unit) {
	out[INT_X] *= u_scale;
	out[INT_Y] *= u_scale;
      }
    }
    out[INT_Z] = in[INT_Z];
  }else{
    myprojectsimple(in,out,lon,lat,azi,0); 
  }

}
#endif
