#include <stdio.h>
#include <string.h>
#include <math.h>

void openGeom( char * path3, char * path2, int * mode,
               int * numPts, int * displace )
{
  FILE * inFile;
  FILE * writeFile;

  int strCounter = 0;
  int counter = 0;

  double strikes, upper, lower, z, width, dips;  
 
  char fileStream[200];
  char * xCoord;
  char * yCoord;
  char * zCoord;
  char * strike;
  char * dip;
  char * half_length;
  char * half_width;
  char * groupNum;

  char * space = " ";
  char * period = ".";

  inFile = fopen( path3, "r" );

  if( inFile == NULL )
  {
    fprintf( stderr, "Error opening inFile %s!  Program Terminating!!!\n", 
	     path3 );
    exit(1);
  }

  writeFile = fopen( path2, "a" );

  if( writeFile == NULL )
  {
    fprintf( stderr, "Error opening outFile %s!  Program Terminating!!!\n", 
	     path2 );
    exit(1);
  }

  while( (fgets( fileStream, 200, inFile ) != NULL) &&
         (strlen(fileStream) > 1)  )
  {
    strCounter++;
  }

  fclose( inFile );

  if( *mode == 0 )
    fprintf( writeFile, "%d %d 0 %d .25 30000 0\n", *numPts, strCounter,
	     strCounter );
  else
    fprintf( writeFile, "%d 0 %d %d .25 30000 0\n", *numPts, strCounter,
             strCounter  );


  fprintf( writeFile, "ITER\n" );
  fprintf( writeFile, "1000 0.00001 1.0 1\n" );
  fprintf( writeFile, "REMOTE\n" );
  fprintf( writeFile, "1\n" );
  fprintf( writeFile, "0. 0. 0. 0. 0. 0.\n" );
  fprintf( writeFile, "0. 0. 0. 0. 0. 0.\n" );
  fprintf( writeFile, "0. 0. 0. 0. 0. 0.\n" );
  fprintf( writeFile, "0. 0. 0. 0. 0. 0. 0.\n" );
  fprintf( writeFile, "PRINT\n" );

  if( *mode == 0 )
    fprintf( writeFile, "FP1\n" );
  else
    fprintf( writeFile, "FP3\n" );

  inFile = fopen( path3, "r" );
  
  while( (fgets( fileStream, 200, inFile ) != NULL) &&
	 (strlen(fileStream) > 1)  )
  { 
    xCoord = strtok( fileStream, space );
    yCoord = strtok( NULL, space );
    zCoord = strtok( NULL, space );
    strike = strtok( NULL, space );
    dip = strtok( NULL, space );
    half_length = strtok( NULL, space );
    half_width = strtok( NULL, space );
    groupNum = strtok( NULL, space );
   
    z = atof( zCoord ); 
    z = -z;

    width = atof( half_width );
    upper = z - width;
    lower = z + width;
        
    strikes = atof( strike );
    strikes = 360 - strikes - 90;
   
    if( *mode == 0 )
    {
      fprintf( writeFile, "1 %s %e %e %s %e %s %s %d 0 0\n", half_length,
	       upper, lower, dip, strikes, xCoord, yCoord, *displace );
    }
    else
    {
      fprintf( writeFile, "1 %s %e %e %s %e %s %s 1 1 0\n", half_length,
               upper, lower, dip, strikes, xCoord, yCoord );
      fprintf( writeFile, "0 0 0 %d 0 0 0 0 0\n", *displace );

    }      
  }

  fclose( inFile );
  fclose( writeFile );
}
