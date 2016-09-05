/*


header file for geoproject, a clone of GMT's (3.4.2) mapproject
for a limited set of projections


$Id: geoproject.h,v 1.1 2003/07/10 02:24:46 becker Exp $

*/

#include <math.h>
#include <stdio.h>
void geoproject(double *, double *, int ,double ,double ,
		double , double, double , double, double,int );
#include "gmt.h" 
#include "myprojectsimple.h"

/* list of projections

don't expect anything to work....
*/

#define PROJECT_AZI 0		/* 
				   simple projection given 
				   clon/clat and azimuth
				   
				 */
#define OMERC_AZI 1		/* oblique Mercator lon/lat/azi */
#define OMERC_POLE 2		/* oblique Mercator lon/lat plon/plat */
#define LCONFORM 3		/* Lambert conic conformal 
				 lon/lat/latp1/latp2 */


#define X 0
#define Y 1
#define Z 2

#define my_boolean unsigned short

