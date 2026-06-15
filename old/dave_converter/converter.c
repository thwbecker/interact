#include <stdio.h>
#include <string.h>
#include "converter.h"

/* file from dave hoang

changes by twb


 */
main( int argc, char * argv[] )
{
   FILE * writeFile;

   char * path2 = "bc.in";   
   char * path3 = "geom.in";

   float xmin, xmax;
   int n; 
   float ymin, ymax; 
   int m;
   float zmin, zmax;
   int o, numObPts;
   
   int fpmode, displacement;
   float dx, dy;  
 
   if( argc != 2 )
   {
     fprintf( stderr, "Usage: %s outFile\n" ,
	      argv[0]);
     exit(1);
   }

   writeFile = fopen( argv[1], "w" );

   fprintf( writeFile, "START\n" );
   fprintf( writeFile, "test\n" );
   fprintf( writeFile, "CPARM\n" );
   fclose( writeFile );
   
   openBC( path2, &xmin, &xmax, &n, &ymin, &ymax, &m, &zmin, &zmax, &o, 
           &fpmode, &displacement );

   numObPts = n*m*o;

   displacement = 0 - displacement;

   openGeom( path3, argv[1], &fpmode, &numObPts, &displacement );

   writeFile = fopen( argv[1], "a" );
   
   fprintf( writeFile, "PRINT\n" );

   if( (n != 1) && (m == 1) && (o == 1) )
   { 
     fprintf( writeFile, "CO2\n" );
     fprintf( writeFile, "%d %d %d %d %d %d\n", numObPts, xmin, ymin, xmax,
              ymax, zmin );
   }
   else if( (m != 1) && (n == 1) && (o == 1) )
   {
     fprintf( writeFile, "CO2\n" );
     fprintf( writeFile, "%d %d %d %d %d %d\n", numObPts, xmin, ymin, xmax,
              ymax, zmin );
   }
   else if( (o != 1) && (m == 1) && (n == 1) )
   {
     fprintf( writeFile, "CO2\n" );
     fprintf( writeFile, "%d %d %d %d %d %d\n", numObPts, xmin, ymin, xmax,
              ymax, zmin );
   }
   else
   {
     
     dx = (xmax-xmin) / (n-1);
     dy = (ymax-ymin) / (n-1); 
     
     
     fprintf( writeFile, "CO3\n" );
     fprintf( writeFile, "LL\n" );
     fprintf( writeFile, "%d %d %d %d %d %d %d\n", 
	      xmin, dx, n, ymin, dy, m, zmin );
   }



   fprintf( writeFile, "DISP\n" );
   fprintf( writeFile, "PRINT\n" );
   fprintf( writeFile, "STOP\n" );
   fclose( writeFile );
      
}
